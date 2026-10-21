#include "LogicFilter/GPU/Trigger.cxx"
#include <random>

void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
uint64_t word(uint32_t type, unsigned time, unsigned channel = 0) {
    unsigned shift = type == SubTimeFrame::TDC64L ? 18 : type == SubTimeFrame::TDC64L_V3 ? 17 : 12;
    unsigned channelShift = type == SubTimeFrame::TDC64L_V3 ? 50 : 51;
    return (uint64_t(TDC64H::T_TDC) << 58) | (uint64_t(channel) << channelShift) | (uint64_t(time) << shift);
}
void clear(Trigger& trigger, unsigned nHbf) {
    trigger.PrepareHbfBatch(nHbf);
    trigger.CleanUpTimeRegion();
}
std::vector<uint32_t> sequential(Trigger& trigger, std::vector<std::vector<HBFIndex>>& frames, unsigned size) {
    std::vector<uint32_t> expected;
    for (size_t frame = 0; frame < frames.size(); ++frame) {
        clear(trigger, 1);
        for (uint32_t bin : *trigger.Exec(frames[frame])) expected.push_back(frame * size + bin);
    }
    return expected;
}
int main() {
    int devices = 0;
    if (cudaGetDeviceCount(&devices) != cudaSuccess || !devices) return 77;
    constexpr unsigned size = 64;
    const std::vector<uint32_t> types = {SubTimeFrame::TDC64H, SubTimeFrame::TDC64L,
        SubTimeFrame::TDC64H_V3, SubTimeFrame::TDC64L_V3};
    Trigger trigger;
    trigger.SetTimeRegion(size, 4);
    trigger.InitializeScanBuffers();
    std::vector<cudamap::FEEAddr> addresses;
    std::vector<cudamap::TriggerConfig> configs;
    cudamap::stageItem(1, 0, 1, 3, 0, cudamap::NO_SUBGROUP, addresses, configs);
    cudamap::stageItem(2, 0, 1, 3, 2, cudamap::NO_SUBGROUP, addresses, configs);
    trigger.InitializeFEEAddrDecoder(addresses, configs);
    std::string expression = "RPN 0";
    trigger.MakeTable(expression);

    std::vector<std::vector<uint64_t>> data = {
        {word(types[0], 62), word(types[0], 63)}, // End window: 61..63; 63 is rejected.
        {word(types[1], 2)}, {},
        {word(types[3], 0), word(types[3], 10)}};
    std::vector<std::vector<HBFIndex>> frames(4);
    for (unsigned frame = 0; frame < 4; ++frame)
        if (!data[frame].empty()) frames[frame].push_back({types[frame], 1,
            reinterpret_cast<unsigned char*>(data[frame].data()), int(data[frame].size() * 8)});
    const std::vector<uint32_t> boundaryExpected = {61, size + 1, 3 * size + 9};
    require(sequential(trigger, frames, size) == boundaryExpected, "Incorrect boundary reference");
    clear(trigger, 4);
    require(*trigger.ExecBatch(frames) == boundaryExpected, "HBF boundary window or rising edge crossed frames");

    // Unequal segment lengths and types, one partial/multiple block segment,
    // 39 segments per frame, empty frames and a deterministic random payload.
    std::mt19937 random(20261021);
    std::vector<std::vector<uint64_t>> randomData(4);
    for (unsigned frame = 0; frame < 4; ++frame)
        for (unsigned i = 0; i < 1025 + frame; ++i)
            randomData[frame].push_back(word(types[frame], random() % (size + 8)));
    for (unsigned frame = 0; frame < 4; ++frame) {
        frames[frame].clear();
        for (unsigned board = 0; board < 39; ++board)
            frames[frame].push_back({types[frame], board % 2 ? 1u : 2u,
                reinterpret_cast<unsigned char*>(randomData[frame].data()),
                int((board ? 1 : randomData[frame].size()) * 8)});
    }
    for (unsigned n : {4, 2, 1, 3, 4}) {
        auto input = std::vector<std::vector<HBFIndex>>(frames.begin(), frames.begin() + n);
        const auto expected = sequential(trigger, input, size);
        clear(trigger, n);
        require(*trigger.ExecBatch(input) == expected, "Batched/serial results or workspace reuse differ");
    }
    clear(trigger, 4);
    std::vector<std::vector<HBFIndex>> empty(4);
    require(trigger.ExecBatch(empty)->empty(), "Old hits leaked into empty batch");
    try { trigger.PrepareHbfBatch(0); require(false, "Zero HBF count accepted"); }
    catch (const std::invalid_argument&) {}
    try { trigger.PrepareHbfBatch(5); require(false, "Oversized batch accepted"); }
    catch (const std::invalid_argument&) {}

    // AND subgroup signals in different HBFs must never combine.
    Trigger subgroup;
    subgroup.SetTimeRegion(size, 4);
    subgroup.InitializeScanBuffers();
    std::vector<cudamap::FEEAddr> subAddresses;
    std::vector<cudamap::TriggerConfig> subConfigs;
    cudamap::stageItem(1, 0, 1, 3, 0, 0, subAddresses, subConfigs);
    cudamap::stageItem(2, 0, 2, 3, 0, 0, subAddresses, subConfigs);
    subgroup.InitializeFEEAddrDecoder(subAddresses, subConfigs);
    uint32_t entries = 2, mapping = 1;
    subgroup.SetSubTimeRegion(size, &entries, 1);
    subgroup.SetfiSubTimeRegionMainTimeRegion(&mapping, 1);
    subgroup.MakeTable(expression);
    uint64_t ten = word(types[0], 10), twenty = word(types[0], 20);
    auto seg = [&](unsigned fem, uint64_t* data) {
        return HBFIndex{types[0], fem, reinterpret_cast<unsigned char*>(data), 8};
    };
    std::vector<std::vector<HBFIndex>> subFrames = {
        {seg(1, &ten)}, {seg(2, &ten)}, {seg(1, &ten), seg(2, &ten)}, {seg(1, &ten), seg(2, &twenty)}};
    for (unsigned n : {4, 2, 4}) {
        clear(subgroup, n);
        subgroup.CleanUpSubTimeRegion();
        auto input = std::vector<std::vector<HBFIndex>>(subFrames.begin(), subFrames.begin() + n);
        require(*subgroup.ExecBatch(input) == (n == 4 ? std::vector<uint32_t>{2 * size + 9} : std::vector<uint32_t>{}),
            "Subgroup state leaked across HBFs");
    }

    // Maximum rising-edge density, including odd region sizes.
    for (unsigned length : {1u, 5u, 64u}) {
        Trigger dense;
        dense.SetTimeRegion(length, 4);
        dense.InitializeScanBuffers();
        std::vector<cudamap::FEEAddr> a;
        std::vector<cudamap::TriggerConfig> c;
        cudamap::stageItem(1, 0, 1, 1, 0, cudamap::NO_SUBGROUP, a, c);
        dense.InitializeFEEAddrDecoder(a, c);
        dense.MakeTable(expression);
        std::vector<uint64_t> words;
        for (unsigned bin = 1; bin < length; bin += 2) words.push_back(word(types[0], bin));
        HBFIndex payload{types[0], 1, reinterpret_cast<unsigned char*>(words.data()), int(words.size() * 8)};
        std::vector<std::vector<HBFIndex>> input(4, std::vector<HBFIndex>{payload});
        std::vector<uint32_t> expected;
        for (unsigned frame = 0; frame < 4; ++frame)
            for (unsigned bin = 1; bin < length; bin += 2) expected.push_back(frame * length + bin);
        clear(dense, 4);
        require(*dense.ExecBatch(input) == expected, "Maximum hit capacity or odd HBF size");
    }
    std::cout << "PASS: four HBFs, boundaries, subgroups, formats, 39 segments, workspace reuse and hit capacity\n";
}
