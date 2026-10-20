#include "LogicFilter/GPU/Trigger.cxx"

void require(bool value) {
    if (!value) throw std::runtime_error("GPU batched Mark regression");
}

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
        uint64_t first = (uint64_t(TDC64H::T_TDC) << 58) | (uint64_t(10) << shift);
        uint64_t second = (uint64_t(TDC64H::T_TDC) << 58) | (uint64_t(20) << shift);
        std::vector<HBFIndex> segments;
        for (int board = 0; board < 39; ++board)
            segments.push_back({type, 1,
                reinterpret_cast<unsigned char*>(board % 2 ? &second : &first), 8});
        for (int hbf = 0; hbf < 4; ++hbf) {
            trigger.CleanUpTimeRegion();
            require(*trigger.Exec(segments) == std::vector<uint32_t>({9, 19}));
        }
    }

    // Different IP/type/length, multiple blocks, a partial block, empty input,
    // and an unregistered IP. The two expected rising edges are bins 9 and 21.
    std::vector<uint64_t> large(1025, 0);
    large.back() = (uint64_t(TDC64H::T_TDC) << 58) | (uint64_t(10) << 12);
    uint64_t low = (uint64_t(TDC64H::T_TDC) << 58) | (uint64_t(20) << 18);
    std::vector<HBFIndex> mixed = {
        {SubTimeFrame::TDC64H, 1, reinterpret_cast<unsigned char*>(large.data()), int(large.size() * 8)},
        {SubTimeFrame::TDC64L, 2, reinterpret_cast<unsigned char*>(&low), 8},
        {SubTimeFrame::TDC64H, 1, nullptr, 0},
        {SubTimeFrame::TDC64L, 99, reinterpret_cast<unsigned char*>(&low), 8}};
    const std::vector<uint32_t> reference = {9, 21};
    trigger.CleanUpTimeRegion();
    require(*trigger.Exec(mixed) == reference);
    trigger.CleanUpTimeRegion();
    std::vector<HBFIndex> empty;
    require(trigger.Exec(empty)->empty());
    // All subgroup bits must clear before its main-region bit is set.
    // Exercise the 32-bit mask boundary and mappings copied during setup.
    Trigger subgroup;
    subgroup.SetTimeRegion(64);
    std::vector<cudamap::FEEAddr> subAddresses;
    std::vector<cudamap::TriggerConfig> subConfigs;
    std::vector<uint64_t> subWords;
    for (unsigned channel = 0; channel < 32; ++channel) {
        cudamap::stageItem(1, channel, uint32_t(1) << channel, 3, 0, 0, subAddresses, subConfigs);
        subWords.push_back((uint64_t(TDC64H::T_TDC) << 58) | (uint64_t(channel) << 51) | (uint64_t(10) << 12));
    }
    subgroup.InitializeFEEAddrDecoder(subAddresses, subConfigs);
    uint32_t entries = 32, mapping = 1;
    subgroup.SetSubTimeRegion(64, &entries, 1);
    subgroup.SetfiSubTimeRegionMainTimeRegion(&mapping, 1);
    subgroup.MakeTable(expression);
    std::vector<HBFIndex> subSegments = {{SubTimeFrame::TDC64H, 1,
        reinterpret_cast<unsigned char*>(subWords.data()), int(subWords.size() * 8)}};
    subgroup.CleanUpTimeRegion();
    subgroup.CleanUpSubTimeRegion(&entries, 1);
    require(*subgroup.Exec(subSegments) == std::vector<uint32_t>({9}));
    subSegments[0].size -= 8;
    subgroup.CleanUpTimeRegion();
    subgroup.CleanUpSubTimeRegion(&entries, 1);
    require(subgroup.Exec(subSegments)->empty());
    std::cout << "PASS: four TDC formats, 39 boards x 4 HBFs, mixed FEM/type, multi-block, empty input and 32-signal subgroup\n";
}
