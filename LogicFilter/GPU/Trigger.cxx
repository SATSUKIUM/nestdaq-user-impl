/*
 *
 *
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <limits>

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
		uint32_t iiSubTimeRegion = workerIndex % subTimeRegionSize;
		uint32_t nEntry = nEntryInSubTimeRegion[iSubTimeRegion];
		subTimeRegions[workerIndex] = (uint32_t(1u) << nEntry) - uint32_t(1u);
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

__global__ void testDecoder(const cudamap::FEEAddrDecoder* decoder)
{
	uint32_t key;

	if (decoder->getKey(0xc0a802a9, 8, key)) {
		printf("%d.%d.%d.%d, %d -> key=%u\n", 0xc0, 0xa8, 0x02, 0xa9, 8, key);
		auto config = decoder->getTriggerConfig(key);
		config.print();
	}
	else {
		printf("Failed to get key for 192.168.2.169, ch 8\n");
	}

	if (decoder->getKey(0xc0a802aa, 16, key)) {
		printf("%d.%d.%d.%d, %d -> key=%u\n", 0xc0, 0xa8, 0x02, 0xaa, 16, key);
		auto config = decoder->getTriggerConfig(key);
		config.print();
	}
	else {
		printf("Failed to get key for 192.168.2.170, ch 16\n");
	}
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
	uint32_t* fiSubTimeRegionMainTimeRegion; // key: iSubTimeRegion, value: iMainTimeRegion
	uint32_t* deviSubTimeRegionMainTimeRegion = nullptr; // key: iSubTimeRegion, value: iMainTimeRegion

	uint32_t* fNEntryInSubTimeRegion; // key: iSubTimeRegion, value: nEntryInSubTimeRegion
	uint32_t* devNEntryInSubTimeRegion; // key: iSubTimeRegion, value: nEntryInSubTimeRegion


	int fEntryCounts_old = 0;
	uint32_t fEntryMask = 0;

	int fNentry = 0;
	uint32_t fTimeRegionSize; // size of main TimeRegion
	uint32_t *fTimeRegion = nullptr; // main TimeRegion
	uint32_t fSubTimeRegionSize; // size of each SubTimeRegion
	// std::vector<uint32_t*> fSubTimeRegions; // vector of SubTimeRegions
	uint32_t* fSubTimeRegions;
	uint32_t fNSubTimeRegions = 0; // number of SubTimeRegions
	uint32_t* devTimeRegion;
	uint32_t* devSubTimeRegions;

	//int fMarkCount = 0;
	//uint32_t fMarkMask = 0;
	std::vector<uint32_t> fHits;
	int fMarkLen = 5;

	TriggerMap fTMap;

	cudamap::FEEAddrDecoder fFEEAddrDecoder;
	cudamap::FEEAddrDecoder* devDecoderPtr = nullptr;
	cudamap::TriggerConfig* devTriggerConfigs = nullptr;

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
		delete[] fTimeRegion;
		fTimeRegion = nullptr;
	}

	return;
}

void Trigger::SetfnEntryInSubTimeRegion(const uint32_t* nEntryInSubTimeRegion, const int nSubTimeRegions){
	fNSubTimeRegions = nSubTimeRegions;
	cudaMallocHost(&fNEntryInSubTimeRegion, fNSubTimeRegions * sizeof(uint32_t));
	for(int i = 0; i < nSubTimeRegions; i++){
		fNEntryInSubTimeRegion[i] = nEntryInSubTimeRegion[i];
	}
	cudaMalloc(&devNEntryInSubTimeRegion, fNSubTimeRegions * sizeof(uint32_t));
	cudaMemcpy(devNEntryInSubTimeRegion, fNEntryInSubTimeRegion, fNSubTimeRegions * sizeof(uint32_t), cudaMemcpyHostToDevice);
}

void Trigger::SetfiSubTimeRegionMainTimeRegion(const uint32_t* iSubTimeRegionMainTimeRegion, const int nSubTimeRegions){
	cudaMallocHost(&fiSubTimeRegionMainTimeRegion, fNSubTimeRegions * sizeof(uint32_t));
	for(int i = 0; i < nSubTimeRegions; i++){
		fiSubTimeRegionMainTimeRegion[i] = iSubTimeRegionMainTimeRegion[i];
	}
	cudaMalloc(&deviSubTimeRegionMainTimeRegion, fNSubTimeRegions * sizeof(uint32_t));
	cudaMemcpy(deviSubTimeRegionMainTimeRegion, fiSubTimeRegionMainTimeRegion, fNSubTimeRegions * sizeof(uint32_t), cudaMemcpyHostToDevice);
}

//void Trigger::SetLogic(int nsignal, std::string formula)
void Trigger::MakeTable(std::string & formula)
{
	fTMap.MakeTable(formula);
	return;
}

void Trigger::InitializeFEEAddrDecoder(
	const std::vector<cudamap::FEEAddr>& feeAddrs,
	const std::vector<cudamap::TriggerConfig>& triggerConfigs)
{
	fFEEAddrDecoder.initialize(feeAddrs, triggerConfigs);

	const uint32_t nConfigs = fFEEAddrDecoder.getSizeSpaceOfKey();
	cudaMalloc(&devTriggerConfigs, nConfigs * sizeof(cudamap::TriggerConfig));
	cudaMemcpy(
		devTriggerConfigs,
		fFEEAddrDecoder.getTriggerConfigArray(),
		nConfigs * sizeof(cudamap::TriggerConfig),
		cudaMemcpyHostToDevice);

	cudamap::FEEAddrDecoder devDecoder = fFEEAddrDecoder;
	devDecoder.setTriggerConfigArray(devTriggerConfigs);

	cudaMalloc(&devDecoderPtr, sizeof(cudamap::FEEAddrDecoder));
	cudaMemcpy(
		devDecoderPtr,
		&devDecoder,
		sizeof(cudamap::FEEAddrDecoder),
		cudaMemcpyHostToDevice);

	#if 1
	// Test the decoder on GPU
	testDecoder<<<1, 1>>>(devDecoderPtr);
	cudaDeviceSynchronize();
	#endif
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
	fTimeRegionSize = size;
	if (fTimeRegion != nullptr) {
		delete[] fTimeRegion;
		fTimeRegion = nullptr;
	}
	cudaMallocHost(&fTimeRegion, fTimeRegionSize * sizeof(uint32_t));
	cudaMalloc(&devTimeRegion, fTimeRegionSize * sizeof(uint32_t));

	numBlocks = cuda::ceil_div(fTimeRegionSize, numThreadsPerBlock);
	devCleanUpTimeRegion<<<numBlocks, numThreadsPerBlock>>>(devTimeRegion, fTimeRegionSize);
	cudaDeviceSynchronize();

	return true;
}

bool Trigger::SetSubTimeRegion(uint32_t subTimeRegionSize, const uint32_t* nEntryInSubTimeRegion, const int nSubTimeRegions){
	// set the size of each SubTimeRegion
	fSubTimeRegionSize = subTimeRegionSize;
	// set the number of SubTimeRegions
	fNSubTimeRegions = nSubTimeRegions;
	// allocate memory for fSubTimeRegions and devSubTimeRegions
	cudaMallocHost(&fSubTimeRegions, nSubTimeRegions * fTimeRegionSize * sizeof(uint32_t*));
	cudaMalloc(&devSubTimeRegions, nSubTimeRegions * fTimeRegionSize * sizeof(uint32_t*));
	// allocate memory for devNEntryInSubTimeRegion and copy nEntryInSubTimeRegion to device
	cudaMalloc(&devNEntryInSubTimeRegion, nSubTimeRegions * sizeof(uint32_t));
	cudaMemcpy(devNEntryInSubTimeRegion, nEntryInSubTimeRegion, nSubTimeRegions * sizeof(uint32_t), cudaMemcpyHostToDevice);

	numBlocks = cuda::ceil_div(nSubTimeRegions * fTimeRegionSize, numThreadsPerBlock);
	devCleanUpSubTimeRegion<<<numBlocks, numThreadsPerBlock>>>(devSubTimeRegions, fSubTimeRegionSize, nEntryInSubTimeRegion, nSubTimeRegions);
	cudaDeviceSynchronize();

	// for(size_t iSubTimeRegion = 0; iSubTimeRegion < nSubTimeRegions; ++iSubTimeRegion){
	// 	uint32_t nEntryInSubTimeRegion = nEntryInSubTimeRegion[iSubTimeRegion];
	// 	#if 0
	// 	std::cout << "#D SubTimeRegion " << iSubTimeRegion << " has " << nEntryInSubTimeRegion << " entries." << std::endl;
	// 	#endif
	// 	uint32_t initialValue = (uint32_t(1u) << nEntryInSubTimeRegion) - uint32_t(1u); // nEntryInSubTimeRegion個の下位ビットが立っている値

	// 	uint32_t* subTimeRegion = new uint32_t[subTimeRegionSize];
	// 	for( uint32_t i = 0 ; i < subTimeRegionSize ; i++) {
	// 		subTimeRegion[i] = initialValue;
	// 	}
	// 	fSubTimeRegions[iSubTimeRegion] = subTimeRegion;
	// }

	return true;

} // bool Trigger::SetSubTimeRegion(int subTimeRegionSize, const std::map<uint32_t, uint32_t>& nEntryInSubTimeRegion)

void Trigger::CleanUpTimeRegion()
{
	int numBlocks = cuda::ceil_div(fTimeRegionSize, numThreadsPerBlock);
	devCleanUpTimeRegion<<<numBlocks, numThreadsPerBlock>>>(devTimeRegion, fTimeRegionSize);
	return;
}

void Trigger::CleanUpSubTimeRegion(const uint32_t* nEntryInSubTimeRegion, const int nSubTimeRegions)
{
	int numBlocks = cuda::ceil_div(fSubTimeRegionSize * nSubTimeRegions, numThreadsPerBlock);
	devCleanUpSubTimeRegion<<<numBlocks, numThreadsPerBlock>>>(devSubTimeRegions, fSubTimeRegionSize, nEntryInSubTimeRegion, nSubTimeRegions);
	return;
} // void Trigger::CleanUpSubTimeRegion(const std::map<uint32_t, uint32_t>& nEntryInSubTimeRegion)

uint32_t *Trigger::GetTimeRegion()
{
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
	if (fEntryCh.count(fem) >= 1) {
		uint64_t *tdcval;
		tdcval = reinterpret_cast<uint64_t *>(pdata);

		for (unsigned int i = 0 ; i < fEntryCh[fem].size() ; i++) {
			int ch = fEntryCh[fem][i];
			int delay = fEntryChDelay[fem][i];
			uint32_t leftwidth = fEntryChLeftWidth[fem][i];
			uint32_t rightwidth = fEntryChRightWidth[fem][i];
			uint32_t markbit = fEntryChBit[fem][i];
			uint32_t iSubTimeRegion = fEntryChiSubTimeRegion[fem][i];
			bool isMemberOfSubGroup = (iSubTimeRegion != UINT32_MAX); // true for member of subgroup, false for not of subgroup

			#if 0
			std::cout << "#DD Trigger::Mark " 
				<< " FEM: " << std::hex << fem
				<< " Ch: " << std::dec << ch
				<< " MarkBit: " << markbit << std::endl;
			#endif

			for (unsigned int j = 0 ; j < (len / sizeof(uint64_t)) ; j++) {

				if (type == SubTimeFrame::TDC64H) {
					struct TDC64H::tdc64 tdc;
					if (TDC64H::Unpack(tdcval[j], &tdc) == TDC64H::T_TDC) {
						if (tdc.ch == ch) {
							uint32_t hit = tdc.tdc4n + delay;

							#if 0
							std::cout << "#D Mark"
								<< " FEM: " << std::hex << fem
								<< " Ch: " << std::dec << ch
								<< " Hit: " << hit
								<< std::endl;
							#endif

							if (hit < fTimeRegionSize - leftwidth) {
								for (int k = -static_cast<int>(leftwidth) ; k < static_cast<int>(rightwidth) + 1 ; k++) {
									if ((hit + k) < fTimeRegionSize) {
										if(isMemberOfSubGroup){
											fSubTimeRegions[iSubTimeRegion][hit + k] &= ~markbit; // サブTimeRegionのビットを降ろす
										}
										else{
											fTimeRegion[hit + k] |= markbit;
										}
									} else if ((static_cast<int>(hit) + k) >= 0) {
										std::cout << "#E Over range hit!"
											<< " FEM: " << std::hex << fem
											<< " Ch: " << std::dec << ch
											<< " Hit: " << hit
											<< std::endl;
									}
								}
							}
						}
					}
				} else
				if (type == SubTimeFrame::TDC64L) {
					struct TDC64L::tdc64 tdc;
					if (TDC64L::Unpack(tdcval[j], &tdc) == TDC64L::T_TDC) {
						if (tdc.ch == ch) {
							uint32_t hit = tdc.tdc4n + delay;

							//std::cout << "#D Mark Ch: " << std::dec << ch
							//	<< " Hit: " << hit << std::endl;

							if (hit < fTimeRegionSize - leftwidth) {
								for (int k = -static_cast<int>(leftwidth) ; k < static_cast<int>(rightwidth) + 1 ; k++) {
									if ((hit + k) < fTimeRegionSize) {
										if(isMemberOfSubGroup){
											fSubTimeRegions[iSubTimeRegion][hit + k] &= ~markbit; // サブTimeRegionのビットを降ろす
										}
										else{
											fTimeRegion[hit + k] |= markbit;
										}
									} else if ((static_cast<int>(hit) + k) >= 0) {
										std::cout << "#E Over range hit!"
											<< " FEM: " << std::hex << fem
											<< " Ch: " << std::dec << ch
											<< " Hit: " << hit
											<< std::endl;
									}
								}
							}
						}
					}
				} else

				if (type == SubTimeFrame::TDC64H_V3) {
					struct TDC64H_V3::tdc64 tdc;
					if (TDC64H_V3::Unpack(tdcval[j], &tdc) == TDC64H_V3::T_TDC) {
						if (tdc.ch == ch) {
							uint32_t hit = tdc.tdc4n + delay;

							#if 0
							std::cout << "#D Mark"
								<< " FEM: " << std::hex << fem
								<< " Ch: " << std::dec << ch
								<< " Hit: " << hit
								<< std::endl;
							#endif

							if (hit < fTimeRegionSize - rightwidth) {
								for (int k = -static_cast<int>(leftwidth) ; k < static_cast<int>(rightwidth) + 1 ; k++) {
									if ((hit + k) < fTimeRegionSize) {
										if(isMemberOfSubGroup){
											fSubTimeRegions[iSubTimeRegion][hit + k] &= ~markbit; // サブTimeRegionのビットを降ろす
										}
										else{
											fTimeRegion[hit + k] |= markbit;
										}
									} else if ((static_cast<int>(hit) + k) >= 0) {
										std::cout << "#E Over range hit!"
											<< " FEM: " << std::hex << fem
											<< " Ch: " << std::dec << ch
											<< " Hit: " << hit
											<< " Mark: " << static_cast<int>(hit) + k
											<< std::endl;
									}
								}
							}
						}
					}
				} else
				if (type == SubTimeFrame::TDC64L_V3) {
					struct TDC64L_V3::tdc64 tdc;
					if (TDC64L_V3::Unpack(tdcval[j], &tdc) == TDC64L_V3::T_TDC) {
						if (tdc.ch == ch) {
							uint32_t hit = tdc.tdc4n + delay;

							//std::cout << "#D Mark Ch: " << std::dec << ch
							//	<< " Hit: " << hit << std::endl;

							if (hit < fTimeRegionSize - leftwidth) {
								for (int k = -static_cast<int>(leftwidth) ; k < static_cast<int>(rightwidth) + 1 ; k++) {
									if ((hit + k) < fTimeRegionSize) {
										if(isMemberOfSubGroup){
											fSubTimeRegions[iSubTimeRegion][hit + k] &= ~markbit; // サブTimeRegionのビットを降ろす
										}
										else{
											fTimeRegion[hit + k] |= markbit;
										}
									} else if ((static_cast<int>(hit) + k) >= 0) {
										std::cout << "#E Over range hit!"
											<< " FEM: " << std::hex << fem
											<< " Ch: " << std::dec << ch
											<< " Hit: " << hit
											<< std::endl;
									}
								}
							}
						}
					}
				}

			}

			//fMarkMask |= (0x1 << fMarkCount);
			//fMarkCount++;
			//std::cout << "#D Trig Mark entry : Module: " << fem << " Ch: " << ch << std::endl;
			//if (static_cast<unsigned int>(fMarkCount) > sizeof(uint32_t)) {
			//	std::cerr << "Entry Ch. exceed " << sizeof(uint32_t) << std::endl;
			//}
			//assert(fMarkCount <= static_cast<int>(sizeof(uint32_t)));
		}
	}


	return;
}

std::vector<uint32_t> *Trigger::Scan()
{
	// --------------------------------
	// Process AND logic
	// --------------------------------
	Trigger::ScanSubTRandMarkMainTR();
	// --------------------------------
	// Finish processing AND logic, and marking MainTimeRegion
	// --------------------------------
	//std::cout << "#D Scan fMarkMask: " << std::hex << fMarkMask << std::endl;
	//std::cout << "#D Scan fEntryMask: " << std::hex << fEntryMask << std::endl;
	fHits.clear();
	fHits.resize(0);

	for (unsigned int i = 0 ; i < fTimeRegionSize - 1; i++) {
		#if 0
		if (((fEntryMask & fTimeRegion[i]) != fEntryMask)
		&& ((fEntryMask & fTimeRegion[i + 1]) == fEntryMask)) {
			fHits.emplace_back(i + 1);
		}
		#else
		if ((! fTMap.LookUp(fTimeRegion[i]))
		 && (fTMap.LookUp(fTimeRegion[i + 1]))) {
			fHits.emplace_back(i + 1);
		}
		#endif

		#if 0
		if (fTimeRegion[i] != 0) {
			std::cout << "#D Scan Time: " << std::dec << i
				<< " Bits: " << std::hex << fTimeRegion[i]
				//<< " Mask: " << fMarkMask << std::endl;
				<< " Mask: " << fEntryMask << std::endl;
		}
		#endif
	}

	return &fHits;
}

void Trigger::ScanSubTRandMarkMainTR()
{
	for(const auto &p : fiSubTimeRegionMainTimeRegion){
		uint32_t iSubTimeRegion = p.first; // index of SubTimeRegion
		uint32_t iMainTimeRegion = p.second; // a bit in MainTimeRegion reserved by this SubTimeRegion

		for(uint32_t i = 0; i < fSubTimeRegionSize; ++i){
			if(fSubTimeRegions[iSubTimeRegion][i] == uint32_t(0u)){
				fTimeRegion[i] |= iMainTimeRegion;
			}
		}
	}
	return;
} // void Trigger::ScanSubTRandMarkMainTR()

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
