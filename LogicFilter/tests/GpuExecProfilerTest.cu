#include "LogicFilter/GPU/Trigger.cxx"
void require(bool value) { if (!value) throw std::runtime_error("GPU Exec profile regression"); }
int main() {
    int devices = 0;
    if (cudaGetDeviceCount(&devices) != cudaSuccess || !devices) return 77;
    Trigger trigger;
    trigger.SetTimeRegion(64);
    std::vector<cudamap::FEEAddr> addresses;
    std::vector<cudamap::TriggerConfig> configs;
    cudamap::stageItem(1, 0, 1, 3, 0, cudamap::NO_SUBGROUP, addresses, configs);
    cudamap::stageItem(2, 0, 1, 3, 2, cudamap::NO_SUBGROUP, addresses, configs);
    trigger.InitializeFEEAddrDecoder(addresses, configs);
    std::string expression = "RPN 0";
    trigger.MakeTable(expression);
    for (uint32_t type : {SubTimeFrame::TDC64H, SubTimeFrame::TDC64L,
                         SubTimeFrame::TDC64H_V3, SubTimeFrame::TDC64L_V3}) {
        unsigned shift = type == SubTimeFrame::TDC64L ? 18 : type == SubTimeFrame::TDC64L_V3 ? 17 : 12;
        uint64_t word = (uint64_t(TDC64H::T_TDC) << 58) | (uint64_t(10) << shift);
        uint64_t secondWord = (uint64_t(TDC64H::T_TDC) << 58) | (uint64_t(20) << shift);
        std::vector<HBFIndex> segments;
        // Reuse the device input buffer for 39 boards and TimeRegion for 4 HBFs.
        for (int board = 0; board < 39; ++board)
            segments.push_back({0, 0, type, 1, reinterpret_cast<unsigned char*>(board % 2 ? &secondWord : &word), 8});
        for (int mode : {0, 1, 2}) {
            const bool enabled = mode > 0;
            trigger.EnableExecDetails(enabled, mode == 2);
            for (int hbf = 0; hbf < 4; ++hbf) {
                trigger.CleanUpTimeRegion();
                require(*trigger.Exec(segments) == std::vector<uint32_t>({9, 19}));
            }
            if (enabled) {
                const auto& p = trigger.GetExecProfile();
                require(p.mark_calls == 39 && p.h2d_bytes == 39 * 8);
                require(p.mark_sync_host_us == 0 && p.mark_launches == 1);
                require(p.metadata_h2d_bytes == 39 * sizeof(MarkSegment));
                require(p.count_d2h_bytes == 4 && p.hits_d2h_bytes == 8);
                require(p.total_host_us > 0 && p.other_host_us >= 0);
                require(mode == 2 ? p.mark_device_us > 0 : p.mark_device_us == 0);
            }
        }
    }
    // Mixed FEM/type, unequal sizes, partial blocks, multiple blocks per segment,
    // skipped empty inputs, and unknown FEM. Compare the original Mark path.
    std::vector<uint64_t> large(1025, 0);
    large.back() = (uint64_t(TDC64H::T_TDC) << 58) | (uint64_t(10) << 12);
    uint64_t lowWord = (uint64_t(TDC64H::T_TDC) << 58) | (uint64_t(20) << 18);
    std::vector<HBFIndex> mixed = {
        {0, 0, SubTimeFrame::TDC64H, 1, reinterpret_cast<unsigned char*>(large.data()), int(large.size() * 8)},
        {0, 0, SubTimeFrame::TDC64L, 2, reinterpret_cast<unsigned char*>(&lowWord), 8},
        {0, 0, SubTimeFrame::TDC64H, 1, nullptr, 0},
        {0, 0, SubTimeFrame::TDC64L, 99, reinterpret_cast<unsigned char*>(&lowWord), 8}};
    trigger.EnableExecDetails(false);
    trigger.CleanUpTimeRegion();
    for (auto& seg : mixed) trigger.Mark(seg.data, seg.size, seg.femId, seg.femType);
    const auto reference = *trigger.Scan();
    require(reference == std::vector<uint32_t>({9, 21}));
    for (int mode : {0, 1, 2}) {
        trigger.EnableExecDetails(mode > 0, mode == 2);
        trigger.CleanUpTimeRegion();
        require(*trigger.Exec(mixed) == reference);
        if (mode > 0) require(trigger.GetExecProfile().mark_launches == 1);
    }
    trigger.CleanUpTimeRegion();
    std::vector<HBFIndex> empty;
    require(trigger.Exec(empty)->empty());
    const auto& p = trigger.GetExecProfile();
    require(p.mark_calls == 0 && p.mark_launches == 0 && p.metadata_h2d_bytes == 0 && p.h2d_bytes == 0 && p.count_d2h_bytes == 4 && p.hits_d2h_bytes == 0);
    std::cout << "PASS: profiling on/off results, four TDC formats, mixed FEM/type and multi-block segments, empty HBF, transfer counters\n";
}
