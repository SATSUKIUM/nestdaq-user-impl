#ifndef NESTDAQ_LOGICFILTER_GPU_EXEC_PROFILER_H
#define NESTDAQ_LOGICFILTER_GPU_EXEC_PROFILER_H
#include "LogicFilter/ExecProfileData.h"
#include <cuda_runtime.h>
#include <chrono>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nestdaq {
class GpuExecProfiler {
public:
    using Clock = std::chrono::high_resolution_clock;
    using Field = double ExecProfileData::*;
    ~GpuExecProfiler() { for (auto event : fEvents) cudaEventDestroy(event); }
    void Enable(bool enabled, bool cudaEvents = false) { fEnabled = enabled; fCudaEvents = cudaEvents; }
    bool Active() const { return fActive; }
    ExecProfileData& Data() { return fData; }
    const ExecProfileData& Data() const { return fData; }
    void Begin(size_t segments) {
        fActive = fEnabled;
        if (!fActive) return;
        fStart = Clock::now();
        fData = {};
        fNext = 0;
        fSpans.clear();
        const auto setupStart = Clock::now();
        const size_t capacity = 4 * segments + 12;
        while (fCudaEvents && fEvents.size() < capacity) {
            cudaEvent_t event;
            Check(cudaEventCreate(&event));
            fEvents.push_back(event);
        }
        if (fCudaEvents) fSpans.reserve(2 * segments + 6);
        fData.instrumentation_host_us += Us(setupStart, Clock::now());
    }
    template<class Operation> void Measure(Field hostField, Field deviceField, Operation&& operation) {
        if (!fActive) { std::forward<Operation>(operation)(); return; }
        cudaEvent_t start = nullptr;
        if (deviceField && fCudaEvents) start = Record();
        const auto begin = Clock::now();
        std::forward<Operation>(operation)();
        fData.*hostField += Us(begin, Clock::now());
        if (deviceField && fCudaEvents) {
            auto end = Record();
            fSpans.push_back({start, end, deviceField});
        }
    }
    void End() {
        if (!fActive) return;
        const auto queryStart = Clock::now();
        // One profiling-only wait at the end of Exec, never between scan stages.
        // The existing result copies have already waited for the actual pipeline.
        if (!fSpans.empty()) Check(cudaEventSynchronize(fSpans.back().end));
        for (const auto& span : fSpans) {
            float ms;
            Check(cudaEventElapsedTime(&ms, span.start, span.end));
            fData.*span.field += double(ms) * 1000;
        }
        fData.instrumentation_host_us += Us(queryStart, Clock::now());
        fData.total_host_us = Us(fStart, Clock::now());
        const double classified = fData.buffer_host_us + fData.h2d_host_us + fData.mark_launch_host_us
            + fData.mark_sync_host_us + fData.subtr_launch_host_us + fData.lut_launch_host_us
            + fData.scan_host_us + fData.compact_launch_host_us + fData.count_d2h_host_us
            + fData.hits_resize_host_us + fData.hits_d2h_host_us + fData.instrumentation_host_us;
        fData.other_host_us = fData.total_host_us - classified;
        fActive = false;
    }
private:
    static void Check(cudaError_t result) {
        if (result != cudaSuccess) throw std::runtime_error(cudaGetErrorString(result));
    }
    static double Us(Clock::time_point start, Clock::time_point end) {
        return std::chrono::duration<double, std::micro>(end - start).count();
    }
    cudaEvent_t Record() {
        const auto start = Clock::now();
        auto event = fEvents.at(fNext++);
        Check(cudaEventRecord(event));
        fData.instrumentation_host_us += Us(start, Clock::now());
        return event;
    }
    struct Span { cudaEvent_t start, end; Field field; };
    bool fEnabled = false, fActive = false, fCudaEvents = false;
    size_t fNext = 0;
    Clock::time_point fStart;
    ExecProfileData fData;
    std::vector<cudaEvent_t> fEvents;
    std::vector<Span> fSpans;
};
}
#endif
