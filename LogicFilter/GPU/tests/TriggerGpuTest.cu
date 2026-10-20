// nvcc -std=c++17 -arch=sm_86 -rdc=true -I. -ILogicFilter/GPU -I$NESTDAQ/include \
//   LogicFilter/GPU/tests/TriggerGpuTest.cu -L$NESTDAQ/lib -lcudamap -o /tmp/trigger-gpu-test
#include "LogicFilter/GPU/Trigger.cxx"
#include <cstdlib>

void require(bool condition) { if (!condition) throw std::runtime_error("GPU regression failed"); }
// Exercise the complete public Mark -> LUT -> scan -> compact -> host path.
void testCompaction()
{
    Trigger trigger;
    std::vector<cudamap::FEEAddr> addresses;
    std::vector<cudamap::TriggerConfig> configs;
    cudamap::stageItem(1, 0, 1, 1, 0, cudamap::NO_SUBGROUP, addresses, configs);
    cudamap::stageItem(1, 1, 4, 1, 0, cudamap::NO_SUBGROUP, addresses, configs);
    trigger.InitializeFEEAddrDecoder(addresses, configs);
    // Resize exercises scratch-buffer reallocation; sizes straddle block boundaries.
    for (uint32_t size : {1u, 2u, 3u, 255u, 256u, 257u, 1025u}) {
        trigger.SetTimeRegion(size);
        for (std::string formula : {"RPN 0", "RPN 0 !"}) {
            trigger.MakeTable(formula);
            TriggerMap reference;
            reference.MakeTable(formula);
            for (unsigned pattern = 0; pattern < 6; ++pattern) {
                trigger.CleanUpTimeRegion();
                std::vector<uint32_t> bins(size, 0);
                std::vector<uint64_t> words;
                uint32_t random = 12345;
                for (uint32_t bin = 0; bin < size; ++bin) {
                    random = random * 1664525u + 1013904223u;
                    uint32_t value = pattern == 0 ? 0 : pattern == 1 ? 1 :
                        pattern == 2 ? (bin & 1) : pattern == 3 ? !(bin & 1) :
                        pattern == 4 ? uint32_t(bin == size - 1) : (random >> 30);
                    // Values >= 2 deliberately exercise out-of-range LUT lookup.
                    bins[bin] = value >= 2 ? 4 : value;
                    if (bins[bin]) words.push_back((uint64_t(TDC64H::T_TDC) << 58)
                        | (uint64_t(value >= 2) << 51) | (uint64_t(bin) << 12));
                }
                trigger.Mark(reinterpret_cast<unsigned char*>(words.data()),
                    words.size() * sizeof(uint64_t), 1, SubTimeFrame::TDC64H);
                std::vector<uint32_t> expected;
                for (uint32_t bin = 1; bin < size; ++bin)
                    if (!reference.LookUp(bins[bin - 1]) && reference.LookUp(bins[bin]))
                        expected.push_back(bin);
                require(expected.size() <= size / 2);
                require(*trigger.Scan() == expected);
                require(*trigger.Scan() == expected); // repeat without stale output
            }
        }
    }
}
int main()
{
    int devices = 0;
    if (cudaGetDeviceCount(&devices) != cudaSuccess || !devices) {
        std::cerr << "SKIP: CUDA device unavailable\n";
        return 77;
    }
    constexpr uint32_t size = 64;
    Trigger trigger;
    trigger.SetTimeRegion(size);
    uint32_t counts[] = {1, 2, 32};
    uint32_t mapping[] = {1u << 3, 1u << 7, 1u << 9};
    trigger.SetSubTimeRegion(size, counts, 3);
    trigger.SetfiSubTimeRegionMainTimeRegion(mapping, 3);
    uint32_t* sub;
    uint32_t* n;
    CheckCuda(cudaMalloc(&sub, 3 * size * sizeof(uint32_t)));
    CheckCuda(cudaMalloc(&n, sizeof(counts)));
    CheckCuda(cudaMemcpy(n, counts, sizeof(counts), cudaMemcpyHostToDevice));
    devCleanUpSubTimeRegion<<<1, 256>>>(sub, size, n, 3);
    uint32_t initial[3 * size];
    CheckCuda(cudaMemcpy(initial, sub, sizeof(initial), cudaMemcpyDeviceToHost));
    for (int i = 0; i < size; ++i) {
        require(initial[i] == 1);
        require(initial[size + i] == 3);
        require(initial[2 * size + i] == UINT32_MAX);
    }
    cudaFree(sub);
    cudaFree(n);
    std::vector<cudamap::FEEAddr> addresses;
    std::vector<cudamap::TriggerConfig> configs;
    cudamap::stageItem(1, 0, 1, 5, 0, 1, addresses, configs);
    cudamap::stageItem(1, 1, 2, 5, 0, 1, addresses, configs);
    cudamap::stageItem(1, 2, 1u << 2, 5, 0, cudamap::NO_SUBGROUP, addresses, configs);
    trigger.InitializeFEEAddrDecoder(addresses, configs);
    std::string formula = "RPN 7";
    trigger.MakeTable(formula);
    const uint32_t types[] = {SubTimeFrame::TDC64H, SubTimeFrame::TDC64L,
        SubTimeFrame::TDC64H_V3, SubTimeFrame::TDC64L_V3};
    for (auto type : types) {
        trigger.CleanUpTimeRegion();
        trigger.CleanUpSubTimeRegion(counts, 3);
        std::vector<uint64_t> words;
        for (int repeat = 0; repeat < 100; ++repeat) {
            for (uint64_t ch = 0; ch < 3; ++ch) {
                unsigned shift = type == SubTimeFrame::TDC64L_V3 ? 50 : 51;
                unsigned timeShift = type == SubTimeFrame::TDC64L ? 18 :
                    type == SubTimeFrame::TDC64L_V3 ? 17 : 12;
                words.push_back((uint64_t(TDC64H::T_TDC) << 58) | (ch << shift) | (uint64_t(10) << timeShift));
            }
        }
        trigger.Mark(reinterpret_cast<unsigned char*>(words.data()), words.size() * 8, 1, type);
        require(*trigger.Scan() == std::vector<uint32_t>{8});
        for (unsigned i = 0; i < size; ++i)
            require(trigger.GetTimeRegion()[i] == ((i >= 8 && i <= 12) ? ((1u << 7) | (1u << 2)) : 0));
    }
    trigger.CleanUpTimeRegion();
    trigger.CleanUpSubTimeRegion(counts, 3);
    trigger.ScanSubTRandMarkMainTR();
    for (unsigned i = 0; i < size; ++i) require(trigger.GetTimeRegion()[i] == 0);
    testCompaction();
    std::cout << "PASS\n";
}
