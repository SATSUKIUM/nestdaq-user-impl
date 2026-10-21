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

inline void ErrorHandleCuda(cudaError_t result)
{
    if (result != cudaSuccess) throw std::runtime_error(cudaGetErrorString(result));
}

// One worker per input word. Formats match UnpackTdc.h (leading TDC only).
__device__ void markWord(const uint64_t* words, size_t count, uint32_t fem, uint32_t type, const cudamap::FEEAddrDecoder* decoder, uint32_t* mainTR, uint32_t* subTR, uint32_t size, uint32_t nSub, size_t i)
{
    if (i >= count) return;
    uint64_t word = words[i];
    if ((word >> 58) != TDC64H::T_TDC) return;
    uint8_t ch;
    int time;
    if (type == SubTimeFrame::TDC64H || type == SubTimeFrame::TDC64H_V3) {
        // 6bit | 7bit | 22bit | 29bit
        // Data type | Channel | TOT | TDC
        // TDC LSB = 0.9765625 ps
        ch = (word >> 51) & 0x7f;
        time = (word & 0x1fffffffULL) >> 12; // &lower 29 bits, then convert to 4 ns bin index
    } else if (type == SubTimeFrame::TDC64L) {
        // 6bit | 7bit | 16bit | 19bit | 16bit
        // Data type | Channel | TOT | TDC | Zero padding
        ch = (word >> 51) & 0x7f;
        // TDC: bits [34:16], 1 ns ticks; padding: bits [15:0].
        // Keep bits [34:18], discarding the lowest 2 TDC bits; shift by 16 + 2 for 4 ns bins.
        time = (word & 0x00000007fffc0000ULL) >> 18;
    } else if (type == SubTimeFrame::TDC64L_V3) {
        // 6bit | 8bit | 16bit | 19bit | 15bit
        // Data type | Channel | TOT | TDC | Zero padding
        ch = (word >> 50) & 0xff;
        time = (word & 0x00000003ffff8000ULL) >> 17; // &(19bits 1)(15bits 0), then convert to 4 ns bin index
    } else return;
    uint32_t key;
    if (!decoder->getKey(fem, ch, key)) return;
    auto config = decoder->getTriggerConfig(key);
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
} // __device__ markWord(const uint64_t*, size_t, uint32_t, uint32_t, const cudamap::FEEAddrDecoder*, uint32_t*, uint32_t*, uint32_t, uint32_t, size_t)

struct MarkSegment {
    size_t offset, count; // Offset and length in uint64_t words, not bytes.
    uint32_t fem, type;
    uint32_t firstBlock; // First CUDA block assigned to this segment.
};


__global__ void devMarkBatch(const uint64_t* words, const MarkSegment* segments, uint32_t nSegments, const cudamap::FEEAddrDecoder* decoder, uint32_t* mainTR, uint32_t* subTR, uint32_t size, uint32_t nSub)
{
    // Binary search for the segment that contains this block.
    uint32_t lo = 0, hi = nSegments;
    while (lo + 1 < hi) {
        uint32_t mid = lo + (hi - lo) / 2;
        if (segments[mid].firstBlock <= blockIdx.x) lo = mid;
        else hi = mid;
    }
    const auto seg = segments[lo];

    // Mark where (&words + seg.offset + (blockIdx.x - seg.firstBlock) * blockDim.x + threadIdx.x)
    // markWord(first word in this segment, # of words in this segment, femId, tdcType, decoder, mainTR, subTR, sizeTR, # of subTRs, word index in this segment)
    markWord(words + seg.offset, seg.count, seg.fem, seg.type, decoder, mainTR, subTR, size, nSub, size_t(blockIdx.x - seg.firstBlock) * blockDim.x + threadIdx.x);
} // __global__ void devMarkBatch(const uint64_t*, const MarkSegment*, uint32_t, const cudamap::FEEAddrDecoder*, uint32_t*, uint32_t*, uint32_t, uint32_t, uint32_t)


__device__ bool devLookUp(uint32_t value, const uint8_t* lut, size_t lutSize)
{
    return value < lutSize && lut[value] != 0;
}

__global__ void devFindRisingEdges(const uint32_t* timeRegion, uint32_t size, const uint8_t* lut, size_t lutSize, uint32_t* flags)
{
    size_t bin = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (bin >= size) return;
    flags[bin] = (bin > 0) && (devLookUp(timeRegion[bin], lut, lutSize)) && (!devLookUp(timeRegion[bin - 1], lut, lutSize));
    return;
}

__global__ void devCompactRisingEdges(const uint32_t* flags, const uint32_t* offsets, uint32_t size, uint32_t* hits, uint32_t* hitCount)
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
    void InitializeScanBuffers();
    bool SetSubTimeRegion(uint32_t subTimeRegionSize, const uint32_t* nEntryInSubTimeRegion, int nSubTimeRegions);
    void SetfiSubTimeRegionMainTimeRegion(const uint32_t* mapping, int nSubTimeRegions);
    void InitializeFEEAddrDecoder(const std::vector<cudamap::FEEAddr>& feeAddrs,
        const std::vector<cudamap::TriggerConfig>& triggerConfigs);
    void MakeTable(std::string& formula);
    void CleanUpTimeRegion();
    void CleanUpSubTimeRegion();
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
    uint32_t fScanRegionSize = 0;

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
    for (int i = 0; i < nSubTimeRegions; ++i){
        if (nEntryInSubTimeRegion[i] == 0 || nEntryInSubTimeRegion[i] > 32){
            throw std::invalid_argument("Subgroup must contain 1..32 signals");
        }
    }
    cudaFree(devNEntryInSubTimeRegion);
    devNEntryInSubTimeRegion = nullptr;
    fNSubTimeRegions = nSubTimeRegions;
    if (!nSubTimeRegions) return;

	ErrorHandleCuda(cudaMalloc(&devNEntryInSubTimeRegion, fNSubTimeRegions * sizeof(uint32_t)));
	ErrorHandleCuda(cudaMemcpy(devNEntryInSubTimeRegion, nEntryInSubTimeRegion, fNSubTimeRegions * sizeof(uint32_t), cudaMemcpyHostToDevice));
}

void Trigger::SetfiSubTimeRegionMainTimeRegion(const uint32_t* iSubTimeRegionMainTimeRegion, const int nSubTimeRegions){
    if (nSubTimeRegions != int(fNSubTimeRegions)) throw std::invalid_argument("Subgroup count mismatch");
    cudaFree(deviSubTimeRegionMainTimeRegion);
    deviSubTimeRegionMainTimeRegion = nullptr;
    if (!nSubTimeRegions) return;

	ErrorHandleCuda(cudaMalloc(&deviSubTimeRegionMainTimeRegion, fNSubTimeRegions * sizeof(uint32_t)));
	ErrorHandleCuda(cudaMemcpy(deviSubTimeRegionMainTimeRegion, iSubTimeRegionMainTimeRegion, fNSubTimeRegions * sizeof(uint32_t), cudaMemcpyHostToDevice));
}

void Trigger::MakeTable(std::string & formula)
{
    TriggerMap table;
    table.MakeTable(formula);
    auto bytes = table.GetLutBytes();
    ErrorHandleCuda(cudaFree(devLut));
    devLut = nullptr;
    fLutSize = 0;
    ErrorHandleCuda(cudaMalloc(&devLut, bytes.size()));
    ErrorHandleCuda(cudaMemcpy(devLut, bytes.data(), bytes.size(), cudaMemcpyHostToDevice));
    fLutSize = bytes.size();
	return;
}

void Trigger::InitializeFEEAddrDecoder(const std::vector<cudamap::FEEAddr>& feeAddrs, const std::vector<cudamap::TriggerConfig>& triggerConfigs)
{
    cudamap::FEEAddrDecoder decoder;
	if (!decoder.initialize(feeAddrs, triggerConfigs)){
        throw std::runtime_error("FEEAddrDecoder initialization failed");
    }
	const uint32_t nConfigs = decoder.getSizeSpaceOfKey();
	ErrorHandleCuda(cudaMalloc(&devTriggerConfigs, nConfigs * sizeof(cudamap::TriggerConfig)));
	ErrorHandleCuda(cudaMemcpy(devTriggerConfigs, decoder.getTriggerConfigArray(), nConfigs * sizeof(cudamap::TriggerConfig), cudaMemcpyHostToDevice));

	cudamap::FEEAddrDecoder devDecoder = decoder;
	devDecoder.setTriggerConfigArray(devTriggerConfigs);

	ErrorHandleCuda(cudaMalloc(&devDecoderPtr, sizeof(cudamap::FEEAddrDecoder)));
	ErrorHandleCuda(cudaMemcpy(devDecoderPtr, &devDecoder, sizeof(cudamap::FEEAddrDecoder), cudaMemcpyHostToDevice));

    return;
}

bool Trigger::SetTimeRegion(uint32_t size)
{
	if (!size || size > uint32_t(INT_MAX)) throw std::invalid_argument("TimeRegion size must be 1..INT_MAX");
    cudaFree(devTimeRegion);
    devTimeRegion = nullptr;
	fTimeRegionSize = size;

	ErrorHandleCuda(cudaMalloc(&devTimeRegion, fTimeRegionSize * sizeof(uint32_t)));

	const int numBlocks = cuda::ceil_div(fTimeRegionSize, numThreadsPerBlock);
	devCleanUpTimeRegion<<<numBlocks, numThreadsPerBlock>>>(devTimeRegion, fTimeRegionSize);
    ErrorHandleCuda(cudaGetLastError());
	ErrorHandleCuda(cudaDeviceSynchronize());

	return true;
}

// Allocate rising-edge flags, exclusive-scan workspace and compacted-hit storage.
void Trigger::InitializeScanBuffers()
{
    if (!fTimeRegionSize) throw std::runtime_error("SetTimeRegion must precede InitializeScanBuffers");
    fScanRegionSize = 0;
    const uint32_t size = fTimeRegionSize;
    ErrorHandleCuda(cudaFree(devEdgeFlags)); devEdgeFlags = nullptr;
    ErrorHandleCuda(cudaFree(devEdgeOffsets)); devEdgeOffsets = nullptr;
    ErrorHandleCuda(cudaFree(devHits)); devHits = nullptr;
    ErrorHandleCuda(cudaFree(devHitCount)); devHitCount = nullptr;
    ErrorHandleCuda(cudaFree(devScanTemp)); devScanTemp = nullptr;
    ErrorHandleCuda(cudaMalloc(&devEdgeFlags, size_t(size) * sizeof(uint32_t)));
    ErrorHandleCuda(cudaMalloc(&devEdgeOffsets, size_t(size) * sizeof(uint32_t)));
    // Rising edges cannot be adjacent; bin 0 is excluded, so floor(N/2) suffices.
    if (size / 2 > 0) ErrorHandleCuda(cudaMalloc(&devHits, size_t(size / 2) * sizeof(uint32_t)));
    ErrorHandleCuda(cudaMalloc(&devHitCount, sizeof(uint32_t)));

    // Query the CUB::ExcluSiveSum workspace size, allocate it
    fScanTempBytes = 0; // CUB::ExcluSiveSum workspace size in bytes
    ErrorHandleCuda(cub::DeviceScan::ExclusiveSum(nullptr, fScanTempBytes, devEdgeFlags, devEdgeOffsets, int(size)));
    if (fScanTempBytes) ErrorHandleCuda(cudaMalloc(&devScanTemp, fScanTempBytes));

    fHits.clear();
    fHits.reserve(size / 2);
    fScanRegionSize = size;
}

bool Trigger::SetSubTimeRegion(uint32_t subTimeRegionSize, const uint32_t* nEntryInSubTimeRegion, const int nSubTimeRegions){
    if (subTimeRegionSize != fTimeRegionSize) throw std::invalid_argument("TimeRegion sizes must match");
    fSubTimeRegionSize = subTimeRegionSize;
    SetfnEntryInSubTimeRegion(nEntryInSubTimeRegion, nSubTimeRegions);
    cudaFree(devSubTimeRegions);
    devSubTimeRegions = nullptr;
    if (nSubTimeRegions) {
        ErrorHandleCuda(cudaMalloc(&devSubTimeRegions, size_t(nSubTimeRegions) * subTimeRegionSize * sizeof(uint32_t)));
        CleanUpSubTimeRegion();
    }
	return true;
}

void Trigger::CleanUpTimeRegion()
{
	int numBlocks = cuda::ceil_div(fTimeRegionSize, numThreadsPerBlock);
	devCleanUpTimeRegion<<<numBlocks, numThreadsPerBlock>>>(devTimeRegion, fTimeRegionSize);
    ErrorHandleCuda(cudaGetLastError());
	return;
}

void Trigger::CleanUpSubTimeRegion()
{
    if (!fNSubTimeRegions) return;

    const int numBlocks = cuda::ceil_div(size_t(fSubTimeRegionSize) * fNSubTimeRegions, size_t(numThreadsPerBlock));
    devCleanUpSubTimeRegion<<<numBlocks, numThreadsPerBlock>>>(devSubTimeRegions, fSubTimeRegionSize, devNEntryInSubTimeRegion, fNSubTimeRegions);
    ErrorHandleCuda(cudaGetLastError());
}

std::vector<uint32_t> *Trigger::Scan()
{
    if (!devLut) throw std::runtime_error("Trigger LUT is not initialized");
    if (!fTimeRegionSize) throw std::runtime_error("TimeRegion is not initialized");
    if (fScanRegionSize != fTimeRegionSize) throw std::runtime_error("InitializeScanBuffers must be called for the current TimeRegion size");

    // Scan SubTR and mark MainTR
    ScanSubTRandMarkMainTR();

    // Find rising edges in mainTR
    const auto blocks = (size_t(fTimeRegionSize) + numThreadsPerBlock - 1) / numThreadsPerBlock;
    devFindRisingEdges<<<blocks, numThreadsPerBlock>>>(devTimeRegion, fTimeRegionSize, devLut, fLutSize, devEdgeFlags);
    ErrorHandleCuda(cudaGetLastError());

    // Exclusive scan of rising-edge flags
    // cub::DeviceScan::ExclusiveSum(workspace, workspaceSize, input, output, numItems)
    ErrorHandleCuda(cub::DeviceScan::ExclusiveSum(devScanTemp, fScanTempBytes, devEdgeFlags, devEdgeOffsets, int(fTimeRegionSize)));

    // Compact rising edges into hits array
    // (flags, offsets, size, return hits, return hitCount)
    devCompactRisingEdges<<<blocks, numThreadsPerBlock>>>(devEdgeFlags, devEdgeOffsets, fTimeRegionSize, devHits, devHitCount);
    ErrorHandleCuda(cudaGetLastError());

    // Copy hit count and hits to host
    uint32_t count = 0;
    ErrorHandleCuda(cudaMemcpy(&count, devHitCount, sizeof(count), cudaMemcpyDeviceToHost));

    if (count > fTimeRegionSize / 2) throw std::runtime_error("Rising edge count exceeds capacity");
    fHits.resize(count);
    if (count){
        ErrorHandleCuda(cudaMemcpy(fHits.data(), devHits, size_t(count) * sizeof(uint32_t), cudaMemcpyDeviceToHost));
    }
    return &fHits;
}

void Trigger::ScanSubTRandMarkMainTR()
{
    if (fNSubTimeRegions) {
        if (!deviSubTimeRegionMainTimeRegion) throw std::runtime_error("Subgroup bit mapping is not initialized");
        size_t count = size_t(fNSubTimeRegions) * fTimeRegionSize;
        devScanSubTRandMarkMainTR<<<(count + numThreadsPerBlock - 1) / numThreadsPerBlock, numThreadsPerBlock>>>(fTimeRegionSize, devSubTimeRegions, devTimeRegion, deviSubTimeRegionMainTimeRegion, fNSubTimeRegions);
        ErrorHandleCuda(cudaGetLastError());
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

    // check if we need to reallocate the mark buffers
    if (totalWords > fMarkWordCapacity) {
        // Allow 2x headroom to avoid reallocating for small increases in HBF size.
        constexpr size_t safetyFactor = 2;
        if (totalWords > std::numeric_limits<size_t>::max() / sizeof(uint64_t) / safetyFactor)
            throw std::overflow_error("Batched Mark buffer capacity overflow");
        const size_t capacity = totalWords * safetyFactor;
        ErrorHandleCuda(cudaFree(devMarkWords));
        devMarkWords = nullptr;
        fMarkWordCapacity = 0;
        ErrorHandleCuda(cudaMalloc(&devMarkWords, capacity * sizeof(uint64_t)));
        fMarkWordCapacity = capacity;
    }

    // check if we need to reallocate the mark segment buffer
    if (fMarkSegments.size() > fMarkSegmentCapacity) {
        ErrorHandleCuda(cudaFree(devMarkSegments));
        devMarkSegments = nullptr;
        fMarkSegmentCapacity = 0;
        ErrorHandleCuda(cudaMalloc(&devMarkSegments, fMarkSegments.size() * sizeof(MarkSegment)));
        fMarkSegmentCapacity = fMarkSegments.size();
    }

    // All copies and the single Mark launch share the default stream.
    // Scan waits for the result before buffers can be reused for the next HBF.
    if (totalWords) {
        size_t index = 0;
        for (const auto& seg : hbf) {
            if (seg.size < int(sizeof(uint64_t))) continue;
            const auto& dst = fMarkSegments[index++];
            ErrorHandleCuda(cudaMemcpyAsync(devMarkWords + dst.offset, seg.data, dst.count * sizeof(uint64_t), cudaMemcpyHostToDevice, 0));
        }
        ErrorHandleCuda(cudaMemcpyAsync(devMarkSegments, fMarkSegments.data(), fMarkSegments.size() * sizeof(MarkSegment), cudaMemcpyHostToDevice, 0));
        devMarkBatch<<<uint32_t(totalBlocks), numThreadsPerBlock>>>(devMarkWords, this->devMarkSegments, uint32_t(fMarkSegments.size()), devDecoderPtr, devTimeRegion, devSubTimeRegions, fTimeRegionSize, fNSubTimeRegions);
        ErrorHandleCuda(cudaGetLastError());
    }

    return Trigger::Scan();
}
