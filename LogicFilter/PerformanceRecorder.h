#ifndef NESTDAQ_LOGICFILTER_PERFORMANCE_RECORDER_H
#define NESTDAQ_LOGICFILTER_PERFORMANCE_RECORDER_H

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>
#include <chrono>
#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include "SubTimeFrameHeader.h"
#include "HeartbeatFrameHeader.h"
#include "ExecProfileData.h"

namespace nestdaq {

// CPU/GPU use exactly the same schema. All elapsed times are host wall time in us.
class LogicFilterPerformanceRecorder {
public:
    using Clock = std::chrono::high_resolution_clock;
    struct Sample {
        ULong64_t iteration{}, input_parts{}, input_bytes{}, hbf_count{}, segment_count{};
        ULong64_t trigger_input_bytes{}, trigger_input_words{}, payload_words{};
        ULong64_t leading_tdc_words{}, trailing_tdc_words{}, tdc_words{}, triggers{};
        ULong64_t output_parts{}, output_bytes{};
        UInt_t tf_id{};
        Int_t receive_result{};
        Bool_t received{}, sent{}, cuda_enabled{}, gpu_exec_enabled{};
        ExecProfileData gpu_exec;
        Double_t conditional_us{}, receive_us{}, cleanup_host_us{}, exec_us{}, trigger_us{};
        Double_t send_us{}, accounting_us{}, processing_us{}, other_us{};
    } sample;

    ~LogicFilterPerformanceRecorder() {
        try { Close(false); } catch (const std::exception& e) { std::cerr << e.what() << '\n'; }
    }
    static double Microseconds(Clock::time_point start, Clock::time_point end) {
        return std::chrono::duration<double, std::micro>(end - start).count();
    }
    bool Enabled() const { return bool(fFile); }

    void Open(const std::string& path, bool cuda, const std::string& signals,
              const std::string& expression, int width, bool execDetails = false, bool execCudaEvents = false) {
        Close(false);
        if (path.empty()) return;
        TDirectory::TContext context;
        fFile.reset(TFile::Open(path.c_str(), "RECREATE"));
        if (!fFile || fFile->IsZombie()) {
            fFile.reset();
            throw std::runtime_error("Cannot open LogicFilter performance output: " + path);
        }
        fCuda = cuda;
        fExecDetails = execDetails;
        fIteration = 0;
        fWriteFailed = false;
        fTree = new TTree("logicfilter", "One ConditionalRun invocation; durations in microseconds");
        fTree->SetDirectory(fFile.get());
        // Periodically save recoverable baskets without writing every event to disk.
        fTree->SetAutoSave(-10000000);
#define LF_BRANCH(name, type) fTree->Branch(#name, &sample.name, #name "/" type)
        LF_BRANCH(iteration, "l"); LF_BRANCH(tf_id, "i");
        LF_BRANCH(received, "O"); LF_BRANCH(sent, "O"); LF_BRANCH(cuda_enabled, "O");
        LF_BRANCH(receive_result, "I");
        LF_BRANCH(gpu_exec_enabled, "O");
#define LF_EXEC_DOUBLE_BRANCH(field) fTree->Branch("gpu_exec_" #field, &sample.gpu_exec.field, "gpu_exec_" #field "/D");
        LF_EXEC_HOST_FIELDS(LF_EXEC_DOUBLE_BRANCH)
        LF_EXEC_DEVICE_FIELDS(LF_EXEC_DOUBLE_BRANCH)
#undef LF_EXEC_DOUBLE_BRANCH
#define LF_EXEC_COUNT_BRANCH(field) fTree->Branch("gpu_exec_" #field, &sample.gpu_exec.field, "gpu_exec_" #field "/l");
        LF_EXEC_COUNT_FIELDS(LF_EXEC_COUNT_BRANCH)
#undef LF_EXEC_COUNT_BRANCH
        LF_BRANCH(input_parts, "l"); LF_BRANCH(input_bytes, "l"); LF_BRANCH(hbf_count, "l");
        LF_BRANCH(segment_count, "l"); LF_BRANCH(trigger_input_bytes, "l");
        LF_BRANCH(trigger_input_words, "l"); LF_BRANCH(payload_words, "l");
        LF_BRANCH(leading_tdc_words, "l"); LF_BRANCH(trailing_tdc_words, "l");
        LF_BRANCH(tdc_words, "l"); LF_BRANCH(triggers, "l");
        LF_BRANCH(output_parts, "l"); LF_BRANCH(output_bytes, "l");
        LF_BRANCH(conditional_us, "D"); LF_BRANCH(receive_us, "D");
        LF_BRANCH(cleanup_host_us, "D"); LF_BRANCH(exec_us, "D"); LF_BRANCH(trigger_us, "D");
        LF_BRANCH(send_us, "D"); LF_BRANCH(accounting_us, "D");
        LF_BRANCH(processing_us, "D"); LF_BRANCH(other_us, "D");
#undef LF_BRANCH
        TNamed("implementation", cuda ? "GPU" : "CPU").Write();
        TNamed("trigger_signals", signals.c_str()).Write();
        TNamed("trigger_expression", expression.c_str()).Write();
        TNamed("clock", "std::chrono::high_resolution_clock").Write();
        TParameter<int>("clock_is_steady", Clock::is_steady).Write();
        TParameter<int>("trigger_width", width).Write();
        TParameter<int>("schema_version", 1).Write();
        TParameter<int>("exec_detail_enabled", execDetails).Write();
        TParameter<int>("exec_cuda_events_enabled", execCudaEvents).Write();
        TParameter<int>("completed", 0).Write();
    }

    // Declare this guard before the FairMQ messages, so they are destroyed before
    // the final timestamp. ROOT Fill itself is deliberately outside conditional_us.
    class Iteration {
    public:
        Iteration(LogicFilterPerformanceRecorder& recorder, Clock::time_point start)
            : fRecorder(recorder), fStart(start) {
            if (!fRecorder.Enabled()) return;
            if (fRecorder.fWriteFailed) throw std::runtime_error("LogicFilter TTree Fill failed");
            fRecorder.sample = {};
            fRecorder.sample.iteration = fRecorder.fIteration++;
            fRecorder.sample.cuda_enabled = fRecorder.fCuda;
            fRecorder.sample.gpu_exec_enabled = fRecorder.fExecDetails;
        }
        ~Iteration() {
            if (!fRecorder.Enabled()) return;
            auto& s = fRecorder.sample;
            s.conditional_us = Microseconds(fStart, Clock::now());
            s.processing_us = s.conditional_us - s.receive_us - s.send_us - s.accounting_us;
            s.other_us = s.processing_us - s.trigger_us;
            if (fRecorder.fTree->Fill() < 0) fRecorder.fWriteFailed = true;
        }
        Iteration(const Iteration&) = delete;
        Iteration& operator=(const Iteration&) = delete;
    private:
        LogicFilterPerformanceRecorder& fRecorder;
        Clock::time_point fStart;
    };

    template<class Segments> void CountSegments(const Segments& segments) {
        for (const auto& seg : segments) {
            ++sample.segment_count;
            sample.trigger_input_bytes += seg.size;
            sample.trigger_input_words += seg.size / sizeof(uint64_t);
            size_t begin = 0;
            if (size_t(seg.size) >= sizeof(HeartbeatFrame::Header)) {
                HeartbeatFrame::Header header;
                std::memcpy(&header, seg.data, sizeof(header));
                if (header.magic == HeartbeatFrame::MAGIC && header.hLength <= seg.size)
                    begin = header.hLength;
            }
            sample.payload_words += (seg.size - begin) / sizeof(uint64_t);
            if (seg.femType != SubTimeFrame::TDC64H && seg.femType != SubTimeFrame::TDC64L &&
                seg.femType != SubTimeFrame::TDC64H_V3 && seg.femType != SubTimeFrame::TDC64L_V3)
                continue;
            for (size_t offset = begin; offset + sizeof(uint64_t) <= size_t(seg.size); offset += sizeof(uint64_t)) {
                uint64_t word;
                std::memcpy(&word, seg.data + offset, sizeof(word));
                const auto type = word >> 58;
                if (type == (0x2c >> 2)) ++sample.leading_tdc_words;
                if (type == (0x34 >> 2)) ++sample.trailing_tdc_words;
            }
        }
        sample.tdc_words = sample.leading_tdc_words + sample.trailing_tdc_words;
    }

    void Close(bool completed) {
        if (!fFile) return;
        TDirectory::TContext context;
        fFile->cd();
        const bool failed = fWriteFailed || fTree->Write("", TObject::kOverwrite) < 0;
        TParameter<int>("completed", completed && !failed).Write("", TObject::kOverwrite);
        fFile->Write();
        const bool ioError = fFile->TestBit(TFile::kWriteError);
        fTree = nullptr;
        fFile->Close();
        fFile.reset();
        if (failed || ioError) throw std::runtime_error("Failed to write LogicFilter performance ROOT file");
    }
private:
    std::unique_ptr<TFile> fFile;
    TTree* fTree = nullptr; // Owned by fFile.
    ULong64_t fIteration = 0;
    bool fCuda = false;
    bool fExecDetails = false;
    bool fWriteFailed = false;
};
} // namespace nestdaq
#endif
