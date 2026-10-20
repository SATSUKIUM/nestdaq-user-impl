#include <cstdint>
#include <string>
#include <vector>
#include <limits>
#include <climits>
#include <stdexcept>

#include "UnpackTdc.h"
#include "SubTimeFrameHeader.h"
#include "TriggerMap.cxx"

#include <cuda_runtime.h>
#include <cuda/cmath>
#include <cudamap.h>
#include <cub/device/device_scan.cuh>

struct HBFIndex {
	uint32_t femType;
	uint32_t femId;
	unsigned char* data;
	int size;
};

__global__ void devCleanUpTimeRegion(uint32_t* timeRegion, uint32_t timeRegionSize)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < timeRegionSize) {
		timeRegion[idx] = 0;
	}
}

__global__ void devCleanUpSubTimeRegion(uint32_t* subTimeRegions, const uint32_t subTimeRegionSize, const uint32_t* nEntryInSubTimeRegion, const uint32_t nSubTimeRegions)
{
	size_t workerIndex = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
	size_t totalWorkers = (size_t)nSubTimeRegions * (size_t)subTimeRegionSize;

	if (workerIndex < totalWorkers) {
		uint32_t iSubTimeRegion = workerIndex / subTimeRegionSize;
		uint32_t nEntry = nEntryInSubTimeRegion[iSubTimeRegion];
		subTimeRegions[workerIndex] = nEntry == 32 ? UINT32_MAX : (uint32_t(1u) << nEntry) - uint32_t(1u);
	}
}

__global__ void devScanSubTRandMarkMainTR(const uint32_t timeRegionSize, const uint32_t* subTimeRegions, uint32_t* timeRegion, const uint32_t* iSubTimeRegionMainTimeRegion, const uint32_t nSubTimeRegions)
{
	size_t workerIndex = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
	size_t totalWorkers = (size_t)nSubTimeRegions * (size_t)timeRegionSize;

	if (workerIndex < totalWorkers) {
		uint32_t iSubTimeRegion = workerIndex / timeRegionSize;
		uint32_t iiSubTimeRegion = workerIndex % timeRegionSize;
		if(subTimeRegions[workerIndex] == 0){
			uint32_t mainTimeRegionBit = iSubTimeRegionMainTimeRegion[iSubTimeRegion];
			atomicOr(&timeRegion[iiSubTimeRegion], mainTimeRegionBit);
		}
	}
}

inline void CheckCuda(cudaError_t result)
{
    if (result != cudaSuccess) throw std::runtime_error(cudaGetErrorString(result));
}

// One worker per input word. Formats match UnpackTdc.h (leading TDC only).
__device__ void markWord(const uint64_t* words, size_t count, uint32_t fem,
    uint32_t type, const cudamap::FEEAddrDecoder* decoder,
    uint32_t* mainTR, uint32_t* subTR, uint32_t size, uint32_t nSub, size_t i)
{
    if (i >= count) return;
    uint64_t word = words[i];
    if ((word >> 58) != TDC64H::T_TDC) return;
    uint8_t ch;
    int time;
    if (type == SubTimeFrame::TDC64H || type == SubTimeFrame::TDC64H_V3) {
        ch = (word >> 51) & 0x7f;
        time = (word & 0x1fffffffULL) >> 12;
    } else if (type == SubTimeFrame::TDC64L) {
        ch = (word >> 51) & 0x7f;
        time = (word & 0x00000007fffc0000ULL) >> 18;
    } else if (type == SubTimeFrame::TDC64L_V3) {
        ch = (word >> 50) & 0xff;
        time = (word & 0x00000003fffe0000ULL) >> 17;
    } else return;
    uint32_t key;
    if (!decoder->getKey(fem, ch, key)) return;
    auto config = decoder->getTriggerConfig(key);
    if (!config.isValid()) return;
    int hit = time + config.delay;
    int halfWidth = config.width / 2;
    // Preserve the CPU window admission rule; clip individual bins at zero.
    if (hit < 0 || hit >= int(size) - halfWidth) return;
    for (int bin = hit - halfWidth; bin <= hit + halfWidth; ++bin) {
        if (bin < 0 || bin >= int(size)) continue;
        if (config.iSubTimeRegion == cudamap::NO_SUBGROUP)
            atomicOr(mainTR + bin, config.hitBit);
        else if (config.iSubTimeRegion < nSub)
            atomicAnd(subTR + size_t(config.iSubTimeRegion) * size + bin, ~config.hitBit);
    }
}

struct MarkSegment {
    size_t offset, count; // Offset and length in uint64_t words, not bytes.
    uint32_t fem, type;
    uint32_t firstBlock; // First CUDA block assigned to this segment.
};

// Each block belongs to one segment, so all its threads share FEM/type.
// A large segment receives multiple blocks; small segments still get one.
__global__ void devMarkBatch(const uint64_t* words, const MarkSegment* segments,
    uint32_t nSegments, const cudamap::FEEAddrDecoder* decoder,
    uint32_t* mainTR, uint32_t* subTR, uint32_t size, uint32_t nSub)
{
    uint32_t lo = 0, hi = nSegments;
    while (lo + 1 < hi) {
        uint32_t mid = lo + (hi - lo) / 2;
        if (segments[mid].firstBlock <= blockIdx.x) lo = mid;
        else hi = mid;
    }
    const auto seg = segments[lo];
    markWord(words + seg.offset, seg.count, seg.fem, seg.type, decoder,
        mainTR, subTR, size, nSub,
        size_t(blockIdx.x - seg.firstBlock) * blockDim.x + threadIdx.x);
}

// Keep LookUp's out-of-range=false behavior, including bits outside the LUT.
__device__ bool devLookUp(uint32_t value, const uint8_t* lut, size_t lutSize)
{
    return value < lutSize && lut[value] != 0;
}

__global__ void devFindRisingEdges(const uint32_t* timeRegion, uint32_t size,
    const uint8_t* lut, size_t lutSize, uint32_t* flags)
{
    size_t bin = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (bin >= size) return;
    // The CPU scan does not report bin 0, even if its LUT value is true.
    flags[bin] = bin > 0 && devLookUp(timeRegion[bin], lut, lutSize)
        && !devLookUp(timeRegion[bin - 1], lut, lutSize);
}

__global__ void devCompactRisingEdges(const uint32_t* flags, const uint32_t* offsets,
    uint32_t size, uint32_t* hits, uint32_t* hitCount)
{
    size_t bin = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (bin >= size) return;
    if (flags[bin]) hits[offsets[bin]] = uint32_t(bin);
    if (bin == size - 1) *hitCount = offsets[bin] + flags[bin];
}

class Trigger
{
public:
    Trigger() = default;
    ~Trigger();
    bool SetTimeRegion(uint32_t size);
    bool SetSubTimeRegion(uint32_t subTimeRegionSize, const uint32_t* nEntryInSubTimeRegion, int nSubTimeRegions);
    void SetfiSubTimeRegionMainTimeRegion(const uint32_t* mapping, int nSubTimeRegions);
    void InitializeFEEAddrDecoder(const std::vector<cudamap::FEEAddr>& feeAddrs,
        const std::vector<cudamap::TriggerConfig>& triggerConfigs);
    void MakeTable(std::string& formula);
    void CleanUpTimeRegion();
    void CleanUpSubTimeRegion(const uint32_t* nEntryInSubTimeRegion, int nSubTimeRegions);
    std::vector<uint32_t>* Exec(const std::vector<HBFIndex>& hbf);

private:
    void SetfnEntryInSubTimeRegion(const uint32_t* nEntryInSubTimeRegion, int nSubTimeRegions);
    void ScanSubTRandMarkMainTR();
    std::vector<uint32_t>* Scan();

    uint32_t fTimeRegionSize = 0;
    uint32_t fSubTimeRegionSize = 0;
    uint32_t fNSubTimeRegions = 0;
    uint32_t* devTimeRegion = nullptr;
    uint32_t* devSubTimeRegions = nullptr;
    uint32_t* devNEntryInSubTimeRegion = nullptr;
    uint32_t* deviSubTimeRegionMainTimeRegion = nullptr;
    std::vector<uint32_t> fHits;

    uint8_t* devLut = nullptr;
    size_t fLutSize = 0;
    uint32_t* devEdgeFlags = nullptr;
    uint32_t* devEdgeOffsets = nullptr;
    uint32_t* devHits = nullptr;
    uint32_t* devHitCount = nullptr;
    void* devScanTemp = nullptr;
    size_t fScanTempBytes = 0;

    cudamap::FEEAddrDecoder* devDecoderPtr = nullptr;
    cudamap::TriggerConfig* devTriggerConfigs = nullptr;
    uint64_t* devMarkWords = nullptr;
    size_t fMarkWordCapacity = 0;
    MarkSegment* devMarkSegments = nullptr;
    size_t fMarkSegmentCapacity = 0;
    std::vector<MarkSegment> fMarkSegments;

    static constexpr int numThreadsPerBlock = 256;
};

Trigger::~Trigger()
{
	if (devDecoderPtr != nullptr) {
		cudaFree(devDecoderPtr);
	}
	if (devTriggerConfigs != nullptr) {
		cudaFree(devTriggerConfigs);
	}

    cudaFree(devLut);
    cudaFree(devEdgeFlags);
    cudaFree(devEdgeOffsets);
    cudaFree(devHits);
    cudaFree(devHitCount);
    cudaFree(devScanTemp);
    cudaFree(devMarkWords);
    cudaFree(devMarkSegments);
    cudaFree(devTimeRegion);
    cudaFree(devSubTimeRegions);
    cudaFree(devNEntryInSubTimeRegion);
    cudaFree(deviSubTimeRegionMainTimeRegion);

	return;
}

void Trigger::SetfnEntryInSubTimeRegion(const uint32_t* nEntryInSubTimeRegion, const int nSubTimeRegions){
	if (nSubTimeRegions < 0) throw std::invalid_argument("Negative subgroup count");
    for (int i = 0; i < nSubTimeRegions; ++i)
        if (nEntryInSubTimeRegion[i] == 0 || nEntryInSubTimeRegion[i] > 32)
            throw std::invalid_argument("Subgroup must contain 1..32 signals");
    cudaFree(devNEntryInSubTimeRegion);
    devNEntryInSubTimeRegion = nullptr;
    fNSubTimeRegions = nSubTimeRegions;
    if (!nSubTimeRegions) return;

	CheckCuda(cudaMalloc(&devNEntryInSubTimeRegion, fNSubTimeRegions * sizeof(uint32_t)));
	CheckCuda(cudaMemcpy(devNEntryInSubTimeRegion, nEntryInSubTimeRegion, fNSubTimeRegions * sizeof(uint32_t), cudaMemcpyHostToDevice));
}

void Trigger::SetfiSubTimeRegionMainTimeRegion(const uint32_t* iSubTimeRegionMainTimeRegion, const int nSubTimeRegions){
    if (nSubTimeRegions != int(fNSubTimeRegions)) throw std::invalid_argument("Subgroup count mismatch");
    cudaFree(deviSubTimeRegionMainTimeRegion);
    deviSubTimeRegionMainTimeRegion = nullptr;
    if (!nSubTimeRegions) return;

	CheckCuda(cudaMalloc(&deviSubTimeRegionMainTimeRegion, fNSubTimeRegions * sizeof(uint32_t)));
	CheckCuda(cudaMemcpy(deviSubTimeRegionMainTimeRegion, iSubTimeRegionMainTimeRegion, fNSubTimeRegions * sizeof(uint32_t), cudaMemcpyHostToDevice));
}

void Trigger::MakeTable(std::string & formula)
{
    TriggerMap table;
    table.MakeTable(formula);
    auto bytes = table.GetLutBytes();
    CheckCuda(cudaFree(devLut));
    devLut = nullptr;
    fLutSize = 0;
    CheckCuda(cudaMalloc(&devLut, bytes.size()));
    CheckCuda(cudaMemcpy(devLut, bytes.data(), bytes.size(), cudaMemcpyHostToDevice));
    fLutSize = bytes.size();
	return;
}

void Trigger::InitializeFEEAddrDecoder(
	const std::vector<cudamap::FEEAddr>& feeAddrs,
	const std::vector<cudamap::TriggerConfig>& triggerConfigs)
{
    cudamap::FEEAddrDecoder decoder;
	if (!decoder.initialize(feeAddrs, triggerConfigs))
        throw std::runtime_error("FEEAddrDecoder initialization failed");

	const uint32_t nConfigs = decoder.getSizeSpaceOfKey();
	CheckCuda(cudaMalloc(&devTriggerConfigs, nConfigs * sizeof(cudamap::TriggerConfig)));
	CheckCuda(cudaMemcpy(
		devTriggerConfigs,
		decoder.getTriggerConfigArray(),
		nConfigs * sizeof(cudamap::TriggerConfig),
		cudaMemcpyHostToDevice));

	cudamap::FEEAddrDecoder devDecoder = decoder;
	devDecoder.setTriggerConfigArray(devTriggerConfigs);

	CheckCuda(cudaMalloc(&devDecoderPtr, sizeof(cudamap::FEEAddrDecoder)));
	CheckCuda(cudaMemcpy(
		devDecoderPtr,
		&devDecoder,
		sizeof(cudamap::FEEAddrDecoder),
		cudaMemcpyHostToDevice));

}

bool Trigger::SetTimeRegion(uint32_t size)
{
	if (!size || size > uint32_t(INT_MAX)) throw std::invalid_argument("TimeRegion size must be 1..INT_MAX");
    cudaFree(devTimeRegion);
    devTimeRegion = nullptr;
	fTimeRegionSize = size;

	CheckCuda(cudaMalloc(&devTimeRegion, fTimeRegionSize * sizeof(uint32_t)));

    CheckCuda(cudaFree(devEdgeFlags)); devEdgeFlags = nullptr;
    CheckCuda(cudaFree(devEdgeOffsets)); devEdgeOffsets = nullptr;
    CheckCuda(cudaFree(devHits)); devHits = nullptr;
    CheckCuda(cudaFree(devHitCount)); devHitCount = nullptr;
    CheckCuda(cudaFree(devScanTemp)); devScanTemp = nullptr;
    CheckCuda(cudaMalloc(&devEdgeFlags, size_t(size) * sizeof(uint32_t)));
    CheckCuda(cudaMalloc(&devEdgeOffsets, size_t(size) * sizeof(uint32_t)));
    // Rising edges cannot be adjacent; bin 0 is excluded, so floor(N/2) suffices.
    if (size / 2) CheckCuda(cudaMalloc(&devHits, size_t(size / 2) * sizeof(uint32_t)));
    CheckCuda(cudaMalloc(&devHitCount, sizeof(uint32_t)));
    fScanTempBytes = 0;
    CheckCuda(cub::DeviceScan::ExclusiveSum(nullptr, fScanTempBytes,
        devEdgeFlags, devEdgeOffsets, int(size)));
    if (fScanTempBytes) CheckCuda(cudaMalloc(&devScanTemp, fScanTempBytes));
    fHits.clear();
    fHits.reserve(size / 2);

	const int numBlocks = cuda::ceil_div(fTimeRegionSize, numThreadsPerBlock);
	devCleanUpTimeRegion<<<numBlocks, numThreadsPerBlock>>>(devTimeRegion, fTimeRegionSize);
    CheckCuda(cudaGetLastError());
	CheckCuda(cudaDeviceSynchronize());

	return true;
}

bool Trigger::SetSubTimeRegion(uint32_t subTimeRegionSize, const uint32_t* nEntryInSubTimeRegion, const int nSubTimeRegions){
    if (subTimeRegionSize != fTimeRegionSize) throw std::invalid_argument("TimeRegion sizes must match");
    fSubTimeRegionSize = subTimeRegionSize;
    SetfnEntryInSubTimeRegion(nEntryInSubTimeRegion, nSubTimeRegions);
    cudaFree(devSubTimeRegions);
    devSubTimeRegions = nullptr;
    if (nSubTimeRegions) {
        CheckCuda(cudaMalloc(&devSubTimeRegions, size_t(nSubTimeRegions) * subTimeRegionSize * sizeof(uint32_t)));
        CleanUpSubTimeRegion(nEntryInSubTimeRegion, nSubTimeRegions);
    }

	return true;

}

void Trigger::CleanUpTimeRegion()
{
	int numBlocks = cuda::ceil_div(fTimeRegionSize, numThreadsPerBlock);
	devCleanUpTimeRegion<<<numBlocks, numThreadsPerBlock>>>(devTimeRegion, fTimeRegionSize);
    CheckCuda(cudaGetLastError());
	return;
}

void Trigger::CleanUpSubTimeRegion(const uint32_t* nEntryInSubTimeRegion, const int nSubTimeRegions)
{
	if (nSubTimeRegions != int(fNSubTimeRegions)) throw std::invalid_argument("Subgroup count mismatch");
    if (!nSubTimeRegions) return;
    CheckCuda(cudaMemcpy(devNEntryInSubTimeRegion, nEntryInSubTimeRegion, nSubTimeRegions * sizeof(uint32_t), cudaMemcpyHostToDevice));
    int numBlocks = cuda::ceil_div(size_t(fSubTimeRegionSize) * nSubTimeRegions, size_t(numThreadsPerBlock));
	devCleanUpSubTimeRegion<<<numBlocks, numThreadsPerBlock>>>(devSubTimeRegions, fSubTimeRegionSize, devNEntryInSubTimeRegion, nSubTimeRegions);
    CheckCuda(cudaGetLastError());
	return;
}

std::vector<uint32_t> *Trigger::Scan()
{
    if (!devLut) throw std::runtime_error("Trigger LUT is not initialized");
    if (!fTimeRegionSize) throw std::runtime_error("TimeRegion is not initialized");
    ScanSubTRandMarkMainTR();
    const auto blocks = (size_t(fTimeRegionSize) + numThreadsPerBlock - 1) / numThreadsPerBlock;
    devFindRisingEdges<<<blocks, numThreadsPerBlock>>>(devTimeRegion,
        fTimeRegionSize, devLut, fLutSize, devEdgeFlags);
    CheckCuda(cudaGetLastError());
    CheckCuda(cub::DeviceScan::ExclusiveSum(devScanTemp, fScanTempBytes,
        devEdgeFlags, devEdgeOffsets, int(fTimeRegionSize)));
    devCompactRisingEdges<<<blocks, numThreadsPerBlock>>>(devEdgeFlags,
        devEdgeOffsets, fTimeRegionSize, devHits, devHitCount);
    CheckCuda(cudaGetLastError());
    uint32_t count = 0;
    CheckCuda(cudaMemcpy(&count, devHitCount, sizeof(count), cudaMemcpyDeviceToHost));
    if (count > fTimeRegionSize / 2) throw std::runtime_error("Rising edge count exceeds capacity");
    fHits.resize(count);
    if (count) CheckCuda(cudaMemcpy(fHits.data(), devHits,
        size_t(count) * sizeof(uint32_t), cudaMemcpyDeviceToHost));
    return &fHits;
}

void Trigger::ScanSubTRandMarkMainTR()
{
    if (fNSubTimeRegions) {
        if (!deviSubTimeRegionMainTimeRegion) throw std::runtime_error("Subgroup bit mapping is not initialized");
        size_t count = size_t(fNSubTimeRegions) * fTimeRegionSize;
        devScanSubTRandMarkMainTR<<<(count + numThreadsPerBlock - 1) / numThreadsPerBlock, numThreadsPerBlock>>>(
            fTimeRegionSize, devSubTimeRegions, devTimeRegion,
            deviSubTimeRegionMainTimeRegion, fNSubTimeRegions);
        CheckCuda(cudaGetLastError());
    }
}

std::vector<uint32_t> *Trigger::Exec(const std::vector<HBFIndex>& hbf)
{
    size_t totalWords = 0, totalBlocks = 0;
    fMarkSegments.clear();
    fMarkSegments.reserve(hbf.size());
    for (const auto& seg : hbf) {
        if (seg.size < int(sizeof(uint64_t))) continue;
        const size_t count = size_t(seg.size) / sizeof(uint64_t);
        const size_t blocks = (count + numThreadsPerBlock - 1) / numThreadsPerBlock;
        if (totalBlocks + blocks > size_t(INT_MAX) ||
            count > std::numeric_limits<size_t>::max() / sizeof(uint64_t) - totalWords)
            throw std::overflow_error("Batched Mark input exceeds capacity");
        fMarkSegments.push_back({totalWords, count, seg.femId, seg.femType, uint32_t(totalBlocks)});
        totalWords += count;
        totalBlocks += blocks;
    }
    if (totalWords && !devDecoderPtr) throw std::runtime_error("FEEAddrDecoder is not initialized");
    if (totalWords > fMarkWordCapacity) {
        CheckCuda(cudaFree(devMarkWords));
        devMarkWords = nullptr;
        fMarkWordCapacity = 0;
        CheckCuda(cudaMalloc(&devMarkWords, totalWords * sizeof(uint64_t)));
        fMarkWordCapacity = totalWords;
    }
    if (fMarkSegments.size() > fMarkSegmentCapacity) {
        CheckCuda(cudaFree(devMarkSegments));
        devMarkSegments = nullptr;
        fMarkSegmentCapacity = 0;
        CheckCuda(cudaMalloc(&devMarkSegments, fMarkSegments.size() * sizeof(MarkSegment)));
        fMarkSegmentCapacity = fMarkSegments.size();
    }
    // All copies and the single Mark launch share the default stream.
    // Scan waits for the result before buffers can be reused for the next HBF.
    if (totalWords) {
        size_t index = 0;
        for (const auto& seg : hbf) {
            if (seg.size < int(sizeof(uint64_t))) continue;
            const auto& dst = fMarkSegments[index++];
            CheckCuda(cudaMemcpyAsync(devMarkWords + dst.offset, seg.data,
                dst.count * sizeof(uint64_t), cudaMemcpyHostToDevice, 0));
        }
        CheckCuda(cudaMemcpyAsync(devMarkSegments, fMarkSegments.data(),
            fMarkSegments.size() * sizeof(MarkSegment), cudaMemcpyHostToDevice, 0));
        devMarkBatch<<<uint32_t(totalBlocks), numThreadsPerBlock>>>(devMarkWords,
            this->devMarkSegments, uint32_t(fMarkSegments.size()), devDecoderPtr,
            devTimeRegion, devSubTimeRegions, fTimeRegionSize, fNSubTimeRegions);
        CheckCuda(cudaGetLastError());
    }

    return Trigger::Scan();
}
