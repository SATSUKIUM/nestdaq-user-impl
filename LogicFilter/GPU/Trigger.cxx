/*
 *
 *
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <limits>
#include <climits>
#include <stdexcept>

#include <string.h>
#include <assert.h>

#include <bitset> // for std::bitset, debugging print

#include "UnpackTdc.h"
#include "SubTimeFrameHeader.h"
#include "TriggerMap.cxx"
#include "SignalParser.h"

#include <cuda_runtime.h>
#include <cuda/cmath>
#include <cudamap.h>
#include <cub/device/device_scan.cuh>


struct HBFIndex {
	int msg_index;
	uint32_t timeFrameId;
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
		return;
	}
	else{
		return;
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
		return;
	}
	else{
		return;
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
		return;
	}
	else{
		return;
	}
}

inline void CheckCuda(cudaError_t result)
{
    if (result != cudaSuccess) throw std::runtime_error(cudaGetErrorString(result));
}

// One worker per input word. Formats match UnpackTdc.h (leading TDC only).
__global__ void devMark(const uint64_t* words, size_t count, uint32_t fem,
    uint32_t type, const cudamap::FEEAddrDecoder* decoder,
    uint32_t* mainTR, uint32_t* subTR, uint32_t size, uint32_t nSub)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
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
	Trigger();
	virtual ~Trigger();
	void InitParam();
	bool SetTimeRegion(uint32_t size);
	bool SetSubTimeRegion(uint32_t subTimeRegionSize, const uint32_t* nEntryInSubTimeRegion, const int nSubTimeRegions);
	void CleanUpTimeRegion();
	void CleanUpSubTimeRegion(const uint32_t* nEntryInSubTimeRegion, const int nSubTimeRegions);
	uint32_t *GetTimeRegion();
	uint32_t GetTimeRegionSize();
	void Entry(uint32_t, int, int); // fem, ch, offset
	void Entry(uint32_t, int, int, uint32_t, uint32_t); // fem, ch, offset, leftwidth, rightwidth
	bool EntryTo(uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, int); // group_id, subgroup_id, iSubTimeRegion, fem, ch, offset
	bool EntryTo(uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, int, uint32_t, uint32_t); // group_id, subgroup_id, iSubTimeRegion, fem, ch, offset,
	void ClearEntry();
	bool CheckEntryFEM(uint32_t);
	void Mark(unsigned char *, int, int, uint32_t);
	void ScanSubTRandMarkMainTR();
	std::vector<uint32_t> *Scan();
	std::vector<uint32_t> *Exec(std::vector<struct HBFIndex> &);
	void SetMarkLen(int val) {fMarkLen = val;};
	int GetMarkLen() {return fMarkLen;};
	//void SetLogic(int);
	void MakeTable(std::string &);
	void SetfnEntryInSubTimeRegion(const uint32_t* nEntryInSubTimeRegion, const int nSubTimeRegions);
	void SetfiSubTimeRegionMainTimeRegion(const uint32_t* iSubTimeRegionMainTimeRegion, const int nSubTimeRegions);
	void InitializeFEEAddrDecoder(
		const std::vector<cudamap::FEEAddr>& feeAddrs,
		const std::vector<cudamap::TriggerConfig>& triggerConfigs);
protected:
private:
	//std::vector<struct CoinCh> fEntry;
	std::map< uint32_t, std::vector<int> > fEntryCh;
	std::map< uint32_t, std::vector<int> > fEntryChDelay;
	std::map< uint32_t, std::vector<uint32_t> > fEntryChBit;
	std::map< uint32_t, std::vector<uint32_t> > fEntryChLeftWidth; // [T - leftwidth, T + rightwidth]
	std::map< uint32_t, std::vector<uint32_t> > fEntryChRightWidth; // [T - leftwidth, T + rightwidth]
	std::map<std::pair<uint32_t, uint32_t>, int> fEntryCounts; // key: (gp_id, subgp_id), value: count of entry
	std::map< uint32_t, std::vector<uint32_t> > fEntryChiSubTimeRegion; // key: fem, value: list of iSubTimeRegion
	uint32_t* fiSubTimeRegionMainTimeRegion = nullptr; // key: iSubTimeRegion, value: iMainTimeRegion
	uint32_t* deviSubTimeRegionMainTimeRegion = nullptr; // key: iSubTimeRegion, value: iMainTimeRegion

	uint32_t* fNEntryInSubTimeRegion = nullptr; // key: iSubTimeRegion, value: nEntryInSubTimeRegion
	uint32_t* devNEntryInSubTimeRegion = nullptr; // key: iSubTimeRegion, value: nEntryInSubTimeRegion


	int fEntryCounts_old = 0;
	uint32_t fEntryMask = 0;

	int fNentry = 0;
	uint32_t fTimeRegionSize = 0; // size of main TimeRegion
	uint32_t *fTimeRegion = nullptr; // main TimeRegion
	uint32_t fSubTimeRegionSize = 0; // size of each SubTimeRegion
	// std::vector<uint32_t*> fSubTimeRegions; // vector of SubTimeRegions
	uint32_t* fSubTimeRegions = nullptr;
	uint32_t fNSubTimeRegions = 0; // number of SubTimeRegions
	uint32_t* devTimeRegion = nullptr;
	uint32_t* devSubTimeRegions = nullptr;

	//int fMarkCount = 0;
	//uint32_t fMarkMask = 0;
	std::vector<uint32_t> fHits;
	int fMarkLen = 5;

	TriggerMap fTMap;
    uint8_t* devLut = nullptr;
    size_t fLutSize = 0;
    uint32_t* devEdgeFlags = nullptr;
    uint32_t* devEdgeOffsets = nullptr;
    uint32_t* devHits = nullptr;
    uint32_t* devHitCount = nullptr;
    void* devScanTemp = nullptr;
    size_t fScanTempBytes = 0;

	cudamap::FEEAddrDecoder fFEEAddrDecoder;
	cudamap::FEEAddrDecoder* devDecoderPtr = nullptr;
	cudamap::TriggerConfig* devTriggerConfigs = nullptr;
    uint64_t* devMarkWords = nullptr;
    size_t fMarkWordCapacity = 0;

	static constexpr int numThreadsPerBlock = 256;
	int numBlocks = 0;
};

Trigger::Trigger()
{
	return;
}

Trigger::~Trigger()
{
	if (devDecoderPtr != nullptr) {
		cudaFree(devDecoderPtr);
	}
	if (devTriggerConfigs != nullptr) {
		cudaFree(devTriggerConfigs);
	}
	if (fTimeRegion != nullptr) {
		cudaFreeHost(fTimeRegion);
		fTimeRegion = nullptr;
	}

    cudaFree(devLut);
    cudaFree(devEdgeFlags);
    cudaFree(devEdgeOffsets);
    cudaFree(devHits);
    cudaFree(devHitCount);
    cudaFree(devScanTemp);
    cudaFree(devMarkWords);
    cudaFree(devTimeRegion);
    cudaFree(devSubTimeRegions);
    cudaFree(devNEntryInSubTimeRegion);
    cudaFree(deviSubTimeRegionMainTimeRegion);
    cudaFreeHost(fSubTimeRegions);
    cudaFreeHost(fNEntryInSubTimeRegion);
    cudaFreeHost(fiSubTimeRegionMainTimeRegion);

	return;
}

void Trigger::SetfnEntryInSubTimeRegion(const uint32_t* nEntryInSubTimeRegion, const int nSubTimeRegions){
	if (nSubTimeRegions < 0) throw std::invalid_argument("Negative subgroup count");
    for (int i = 0; i < nSubTimeRegions; ++i)
        if (nEntryInSubTimeRegion[i] == 0 || nEntryInSubTimeRegion[i] > 32)
            throw std::invalid_argument("Subgroup must contain 1..32 signals");
    cudaFreeHost(fNEntryInSubTimeRegion);
    cudaFree(devNEntryInSubTimeRegion);
    fNEntryInSubTimeRegion = nullptr;
    devNEntryInSubTimeRegion = nullptr;
    fNSubTimeRegions = nSubTimeRegions;
    if (!nSubTimeRegions) return;
	CheckCuda(cudaMallocHost(&fNEntryInSubTimeRegion, fNSubTimeRegions * sizeof(uint32_t)));
	for(int i = 0; i < nSubTimeRegions; i++){
		fNEntryInSubTimeRegion[i] = nEntryInSubTimeRegion[i];
	}
	CheckCuda(cudaMalloc(&devNEntryInSubTimeRegion, fNSubTimeRegions * sizeof(uint32_t)));
	CheckCuda(cudaMemcpy(devNEntryInSubTimeRegion, fNEntryInSubTimeRegion, fNSubTimeRegions * sizeof(uint32_t), cudaMemcpyHostToDevice));
}

void Trigger::SetfiSubTimeRegionMainTimeRegion(const uint32_t* iSubTimeRegionMainTimeRegion, const int nSubTimeRegions){
    if (nSubTimeRegions != int(fNSubTimeRegions)) throw std::invalid_argument("Subgroup count mismatch");
    cudaFreeHost(fiSubTimeRegionMainTimeRegion);
    cudaFree(deviSubTimeRegionMainTimeRegion);
    fiSubTimeRegionMainTimeRegion = nullptr;
    deviSubTimeRegionMainTimeRegion = nullptr;
    if (!nSubTimeRegions) return;
	CheckCuda(cudaMallocHost(&fiSubTimeRegionMainTimeRegion, fNSubTimeRegions * sizeof(uint32_t)));
	for(int i = 0; i < nSubTimeRegions; i++){
		fiSubTimeRegionMainTimeRegion[i] = iSubTimeRegionMainTimeRegion[i];
	}
	CheckCuda(cudaMalloc(&deviSubTimeRegionMainTimeRegion, fNSubTimeRegions * sizeof(uint32_t)));
	CheckCuda(cudaMemcpy(deviSubTimeRegionMainTimeRegion, fiSubTimeRegionMainTimeRegion, fNSubTimeRegions * sizeof(uint32_t), cudaMemcpyHostToDevice));
}

//void Trigger::SetLogic(int nsignal, std::string formula)
void Trigger::MakeTable(std::string & formula)
{
	fTMap.MakeTable(formula);
    auto bytes = fTMap.GetLutBytes();
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
	if (!fFEEAddrDecoder.initialize(feeAddrs, triggerConfigs))
        throw std::runtime_error("FEEAddrDecoder initialization failed");

	const uint32_t nConfigs = fFEEAddrDecoder.getSizeSpaceOfKey();
	CheckCuda(cudaMalloc(&devTriggerConfigs, nConfigs * sizeof(cudamap::TriggerConfig)));
	CheckCuda(cudaMemcpy(
		devTriggerConfigs,
		fFEEAddrDecoder.getTriggerConfigArray(),
		nConfigs * sizeof(cudamap::TriggerConfig),
		cudaMemcpyHostToDevice));

	cudamap::FEEAddrDecoder devDecoder = fFEEAddrDecoder;
	devDecoder.setTriggerConfigArray(devTriggerConfigs);

	CheckCuda(cudaMalloc(&devDecoderPtr, sizeof(cudamap::FEEAddrDecoder)));
	CheckCuda(cudaMemcpy(
		devDecoderPtr,
		&devDecoder,
		sizeof(cudamap::FEEAddrDecoder),
		cudaMemcpyHostToDevice));

}

void Trigger::InitParam()
{
	//fMarkCount = 0;
	//fMarkMask = 0;
	fHits.clear();
	fHits.resize(0);
	if (fTimeRegion != nullptr) {
		memset(fTimeRegion, 0, fTimeRegionSize * sizeof(uint32_t));
	}

	return;
}

bool Trigger::SetTimeRegion(uint32_t size)
{
	if (!size || size > uint32_t(INT_MAX)) throw std::invalid_argument("TimeRegion size must be 1..INT_MAX");
    cudaFree(devTimeRegion);
    devTimeRegion = nullptr;
	fTimeRegionSize = size;
	if (fTimeRegion != nullptr) {
		cudaFreeHost(fTimeRegion);
		fTimeRegion = nullptr;
	}
	CheckCuda(cudaMallocHost(&fTimeRegion, fTimeRegionSize * sizeof(uint32_t)));
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

	numBlocks = cuda::ceil_div(fTimeRegionSize, numThreadsPerBlock);
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

} // bool Trigger::SetSubTimeRegion(int subTimeRegionSize, const std::map<uint32_t, uint32_t>& nEntryInSubTimeRegion)

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
} // void Trigger::CleanUpSubTimeRegion(const std::map<uint32_t, uint32_t>& nEntryInSubTimeRegion)

uint32_t *Trigger::GetTimeRegion()
{
    // Explicit diagnostic access only; the normal Scan path copies compacted hits.
    CheckCuda(cudaMemcpy(fTimeRegion, devTimeRegion, size_t(fTimeRegionSize) * sizeof(uint32_t), cudaMemcpyDeviceToHost));
	return fTimeRegion;
}

uint32_t Trigger::GetTimeRegionSize()
{
	return fTimeRegionSize;
}


//if (Trig::CheckEntryFEM(lsubtimeframe.femId)) Trig::Mark(pdata);

#if 0
void Trigger::Entry(uint64_t fem, int ch)
{
	bool match = false;
	for (auto ent ; fEntry) {
		if (ent.femId == fem) {
			fEntry.Ch.emplace_back(ch);
			match = true;
			fNentry++;
		}
	} 
	if (match) {
		struct CoinCh newentry;
		newentry.fenId = fem;
		newentry.Ch.emplace_back(ch);
		fEntry.emplace_back(newentry);
		fNentry++;
	}

	return;
}

bool Trigger::CheckEntryFEM(uint64_t fem)
{
	bool rval = false;
	for (auto ent ; fEntry) if (ent.femId == fem) rval = true;
	return rval;
}

#endif

void Trigger::Entry(uint32_t fem, int ch, int offset)
{

	fEntryCh[fem].emplace_back(ch);
	fEntryChDelay[fem].emplace_back(offset);
	fEntryChBit[fem].emplace_back(0x0000001 << fEntryCounts_old);
	fEntryChLeftWidth[fem].emplace_back(fMarkLen / 2); // if mq-param give a signal(fem, ch, offset) with 3 parameters, leftwidth and rightwidth are set to default value (MarkLen / 2)
	fEntryChRightWidth[fem].emplace_back(fMarkLen / 2); // if mq-param give a signal(fem, ch, offset) with 3 parameters, leftwidth and rightwidth are set to default value (MarkLen / 2)
	fEntryMask |= 0x00000001 << fEntryCounts_old;
	fEntryCounts_old++;

	#if 0
	std::cout << "#D Trig Entry : Module: " << fem << " Ch: " << ch << std::endl;
	#endif
	if (static_cast<unsigned int>(fEntryCounts_old) > (sizeof(uint32_t) * 8)) {
		std::cerr << "Entry Ch. exceed " << sizeof(uint32_t) * 8<< std::endl;
	}
	assert(fEntryCounts_old <= static_cast<int>(sizeof(uint32_t) * 8));

	return;
}

bool Trigger::EntryTo(uint32_t group_id, uint32_t subgroup_id, uint32_t iSubTimeRegion, uint32_t fem, uint32_t ch, int offset)
{
	constexpr int bitCount = std::numeric_limits<uint32_t>::digits;
	if(group_id >= bitCount) return false;
	if(subgroup_id != SignalParser::NO_SUBGROUP){
		auto entry = fEntryCounts.find(std::make_pair(group_id, subgroup_id));
		// SetSubTimeRegion/CleanUpSubTimeRegion also shift by the entry count.
		if(entry != fEntryCounts.end() && entry->second >= bitCount - 1) return false;
	}

	fEntryCh[fem].emplace_back(ch);
	fEntryChDelay[fem].emplace_back(offset);
	if(subgroup_id == SignalParser::NO_SUBGROUP){
		fEntryChBit[fem].emplace_back(uint32_t(1u) << group_id);
	}
	else{
		fEntryChBit[fem].emplace_back(uint32_t(1u) << fEntryCounts[std::make_pair(group_id, subgroup_id)]++);
	}

	fEntryChLeftWidth[fem].emplace_back(fMarkLen / 2);
	fEntryChRightWidth[fem].emplace_back(fMarkLen / 2);
	
	if(subgroup_id == SignalParser::NO_SUBGROUP){
		fEntryChiSubTimeRegion[fem].emplace_back(UINT32_MAX); // なにか入れておかないと他のエントリーとずれてしまうので、とりあえず。
	}
	else{
		fEntryChiSubTimeRegion[fem].emplace_back(iSubTimeRegion);
	}

	if(iSubTimeRegion != UINT32_MAX) fiSubTimeRegionMainTimeRegion[iSubTimeRegion] = uint32_t(1u) << group_id;
	if(subgroup_id == SignalParser::NO_SUBGROUP){
		std::cout << "[Trigger::EntryTo] This signal (group_id=" << group_id << ", subgroup_id=NO_SUBGROUP) is assigned to MainTimeRegion bit " << std::bitset<32>(fEntryChBit[fem].back()) << std::endl;
	}
	else{
		std::cout << "[Trigger::EntryTo] This signal (group_id=" << group_id << ", subgroup_id=" << subgroup_id << ") is assigned to SubTimeRegion " << iSubTimeRegion << " and its bit " << std::bitset<32>(fEntryChBit[fem].back()) << ". And to be assigned to MainTimeRegion bit " << std::bitset<32>(fiSubTimeRegionMainTimeRegion[iSubTimeRegion]) << std::endl;
	}

	fEntryMask |= uint32_t(1u) << group_id;
	
	return true;
} // bool Trigger::EntryTo(uint32_t group_id, uint32_t subgroup_id, uint32_t iSubTimeRegion, uint32_t fem, uint32_t ch, int offset)

void Trigger::ClearEntry()
{
	fEntryCh.clear();
	fEntryChDelay.clear();
	fEntryChBit.clear();
	fEntryMask = 0x00000000;
	fEntryCounts.clear();
	fEntryChLeftWidth.clear();
	fEntryChRightWidth.clear();
	fEntryChiSubTimeRegion.clear();

	return;
}

bool Trigger::CheckEntryFEM(uint32_t fem)
{
	bool rval = false;
	if (fEntryCh.count(fem) >= 1) rval = true;
	return rval;
}

void Trigger::Mark(unsigned char *pdata, int len, int fem, uint32_t type)
{
    if (len <= 0) return;
    if (!devDecoderPtr) throw std::runtime_error("FEEAddrDecoder is not initialized");
    size_t count = size_t(len) / sizeof(uint64_t);
    if (!count) return;
    if (count > fMarkWordCapacity) {
        CheckCuda(cudaFree(devMarkWords));
        devMarkWords = nullptr;
        fMarkWordCapacity = 0;
        CheckCuda(cudaMalloc(&devMarkWords, count * sizeof(uint64_t)));
        fMarkWordCapacity = count;
    }
    CheckCuda(cudaMemcpy(devMarkWords, pdata, count * sizeof(uint64_t), cudaMemcpyHostToDevice));
    devMark<<<(count + numThreadsPerBlock - 1) / numThreadsPerBlock, numThreadsPerBlock>>>(
        devMarkWords, count, uint32_t(fem), type, devDecoderPtr,
        devTimeRegion, devSubTimeRegions, fTimeRegionSize, fNSubTimeRegions);
    CheckCuda(cudaGetLastError());
    CheckCuda(cudaDeviceSynchronize());
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

std::vector<uint32_t> *Trigger::Exec(std::vector<struct HBFIndex> &hbf)
{
	for (auto &seg : hbf) {
		//fTrig->Mark(
		//	reinterpret_cast<unsigned char *>(inParts[mindex].GetData()),
		//	inParts[mindex].GetSize(),
		//	vfemid, dbl->Type);
		Trigger::Mark(seg.data, seg.size, seg.femId, seg.femType);
	}

	#if 0
	uint32_t *tr = fTrig->GetTimeRegion();
	std::cout << "####DDDD Hit TimeRegion: ";
	for (uint32_t ii = 0 ; ii < fTrig->GetTimeRegionSize() ; ii++) {
		if (tr[ii] != 0) {
		std::cout << " " << std::dec << i << ":"
			<< std::hex << std::setw(4) << std::setfill('0')
			<< tr[ii];
		}
	}
	std::cout << std::endl;
	#endif

	#if 0
	std::cout << "# HB: " << std::dec << i;
	for (size_t iifem = 0 ; iifem < block_map.size() ; iifem++) {
		struct DataBlock *dbl = &block_map[iifem][i];
		uint64_t vfemid = dbl->femId;
		uint64_t vhbframe = dbl->HBFrame;
		//std::cout << "# HB: " << std::dec << i
		//<< " FEM: " << std::hex << vfemid
		std::cout << " " << std::dec << (vfemid  & 0xff) << ":" << vhbframe;
	}
	std::cout << std::endl;
	#endif

	return Trigger::Scan();
}




#ifdef TEST_MAIN_TRIG
int main(int argc, char* argv[])
{
	Trigger trig;

	static char cbuf[16];
	uint64_t *pdata = reinterpret_cast<uint64_t *>(cbuf);

	uint64_t *buf = new uint64_t[1024*1024*8];


	int time_region = 1024 * 256; //18 bit
	trig.SetTimeRegion(time_region);

	trig.Entry(1, 26, 10);
	trig.Entry(1, 29, 0);


	//struct TDC64H::tdc64 tdc;
	int i = 0;
	while (true) {
		std::cin.read(cbuf, 8);
		//std::cin >> *pdata; 
		if (std::cin.eof()) break;
		buf[i++] = *pdata;

		std::cout << "\r "  << i << ": " << *pdata << "  " << std::flush;
	
	}
	//int len = i * sizeof(uint64_t);
	std::cout << std::endl;

	unsigned char *pcurr = reinterpret_cast<unsigned char *>(buf);
	unsigned char *pend = reinterpret_cast<unsigned char *>(buf + i);
	unsigned char *pnext = nullptr;
	while (true) {
		int len = TDC64H::GetHBFrame(pcurr, pend, &pnext);
		if (len <= 0) break;

		std::cout << "#D HB frame size: " << std::dec << len
			<< " curr: " << std::hex
			<< reinterpret_cast<uintptr_t>(pcurr)
			<< " next: "
			<< reinterpret_cast<uintptr_t>(pnext) << std::endl;

		pdata = reinterpret_cast<uint64_t *>(pcurr);
		#if 0
		for (unsigned int j = 0 ; j < len / sizeof(uint64_t) ; j++) {
			tdc64h_dump(pdata[j]);
		}
		#endif

		trig.InitParam();
		trig.Mark(reinterpret_cast<unsigned char *>(pdata), len, 1, 0);
		std::vector<uint32_t> *nhits = trig.Scan();
		std::cout << "# Hits: ";
		for (auto hit : *nhits) std::cout << " " << hit;
		std::cout << std::endl;
		std::cout << "Nhits: " << std::dec << nhits->size()
			<< " / " << time_region << std::endl;




		pcurr = pnext;
	}


	#if 0
	trig.Mark(reinterpret_cast<unsigned char *>(buf), len, 1);
	std::vector<uint32_t> *nhits = trig.Scan();
	std::cout << "# ";
	for (auto hit : *nhits) std::cout << " " << hit;
	std::cout << std::endl;
	std::cout << "Nhits: " << std::dec << nhits->size()
		<< " / " << time_region << std::endl;
	#endif

	return 0;
}
#endif
