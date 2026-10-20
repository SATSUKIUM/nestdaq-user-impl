#include "LogicFilter/PerformanceRecorder.h"
#include <vector>

struct Segment { unsigned char* data; int size; uint32_t femType; };
int main()
{
    nestdaq::LogicFilterPerformanceRecorder recorder;
    std::vector<uint64_t> words(6);
    HeartbeatFrame::Header header;
    header.length = words.size() * sizeof(uint64_t);
    // A header word that looks like leading TDC must not count as payload TDC.
    header.type = 0x2c00;
    std::memcpy(words.data(), &header, sizeof(header));
    words[2] = uint64_t(0x0b) << 58;
    words[3] = uint64_t(0x0d) << 58;
    words[4] = uint64_t(0x1c) << 58;
    words[5] = 0;
    uint64_t unsupported = uint64_t(0x0b) << 58;
    std::vector<Segment> segments = {
        {reinterpret_cast<unsigned char*>(words.data()), int(words.size() * 8), SubTimeFrame::TDC64H_V3},
        {reinterpret_cast<unsigned char*>(&unsupported), 8, SubTimeFrame::INPUT_REG}
    };
    recorder.CountSegments(segments);
    const auto& s = recorder.sample;
    if (s.segment_count != 2 || s.trigger_input_bytes != 56 || s.trigger_input_words != 7 ||
        s.payload_words != 5 || s.leading_tdc_words != 1 || s.trailing_tdc_words != 1 || s.tdc_words != 2)
        throw std::runtime_error("TDC accounting regression");
    // Check all supported FEM formats and accumulation over multiple HBFs.
    for (auto type : {SubTimeFrame::TDC64H, SubTimeFrame::TDC64L, SubTimeFrame::TDC64L_V3}) {
        segments[0].femType = type;
        recorder.CountSegments(segments);
    }
    if (s.tdc_words != 8 || s.segment_count != 8) throw std::runtime_error("TDC accumulation regression");
    std::cout << "PASS: header exclusion, TDC types, unsupported FEM, accumulation\n";
}
