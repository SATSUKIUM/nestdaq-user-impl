/*
 *
 *
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

#include <string.h>
#include <assert.h>

#include <chrono>
#include <cmath>

#include "UnpackTdc.h"
#include "SubTimeFrameHeader.h"
#include "TriggerMap.h"
#include "LogiCalc.h"

#include <bitset>


struct HBFIndex {
	int msg_index;
	uint32_t timeFrameId;
	uint32_t femType;
	uint32_t femId;
	unsigned char* data;
	int size;
};

class Trigger
{
public:
	Trigger();
	virtual ~Trigger();
	virtual void InitParam() = 0;
	virtual bool SetTimeRegion(int) = 0;
	virtual void CleanUpTimeRegion() = 0;
	// uint32_t *GetTimeRegion();
	uint32_t GetTimeRegionSize();
	virtual void Entry(uint32_t, int, int) = 0; // fem, ch, offset
	virtual void Entry(uint32_t, int, int, uint32_t, uint32_t) = 0; // fem, ch, offset, leftwidth, rightwidth
	virtual void ClearEntry();
	bool CheckEntryFEM(uint32_t);
	virtual void Mark(unsigned char *, int, int, uint32_t) = 0;
	virtual std::vector<uint32_t> *Scan() = 0;
	std::vector<uint32_t> *Exec(std::vector<struct HBFIndex> &);
	void SetMarkLen(int val) {fMarkLen = val;};
	int GetMarkLen() {return fMarkLen;};
	//void SetLogic(int);
	virtual void SetTExpression(std::string &) = 0;
protected:
	//std::vector<struct CoinCh> fEntry;
	std::map< uint32_t, std::vector<int> > fEntryCh;
	std::map< uint32_t, std::vector<int> > fEntryChDelay;
	// std::map< uint32_t, std::vector<uint32_t> > fEntryChBit;
	std::map< uint32_t, std::vector<uint32_t> > fEntryChLeftWidth; // [T - leftwidth, T + rightwidth]
	std::map< uint32_t, std::vector<uint32_t> > fEntryChRightWidth; // [T - leftwidth, T + rightwidth]
	int fEntryCounts = 0;

	uint32_t fTimeRegionSize;
	int fMarkLen = 5; // default set value ( to be changed by mq-param ), common for all channels

	int fNentry = 0;
	
	std::vector<uint32_t> fHits;
	TriggerMap fTMap; // これ持つのTrigger32だけでええかもな...
}; // class Trigger

class Trigger32 : public Trigger, public TriggerMap{
public:
	Trigger32() : Trigger() {};
	virtual ~Trigger32() {};
	void InitParam() override;
	void SetTExpression(std::string &tx) override; // call TriggerMap::MakeTable()
	bool SetTimeRegion(int size) override;
	void CleanUpTimeRegion() override;
	void Mark(unsigned char *, int, int, uint32_t) override;
	std::vector<uint32_t> *Scan() override;
	void Entry(uint32_t fem, int ch, int offset) override;
	void Entry(uint32_t fem, int ch, int offset, uint32_t leftwidth, uint32_t rightwidth) override;
	void ClearEntry() override;

private:
	uint32_t* fTimeRegion = nullptr;
	uint32_t fEntryMask = 0;

	std::map< uint32_t, std::vector<uint32_t> > fEntryChBit;
}; // class Trigger32 : public Trigger, public TriggerMap


class TriggerBitSet : public Trigger, public LogiCalc {
public:
	TriggerBitSet() : Trigger(), LogiCalc() {};
	virtual ~TriggerBitSet() {};
	void InitParam() override; // もともと誰からも呼ばれてなかったが、一応。
	constexpr static int defaultSizeFTimeRegion = 512;
	bool SetTimeRegion(int size) override;
	void CleanUpTimeRegion() override;
	void Mark(unsigned char *, int, int, uint32_t) override;
	std::vector<uint32_t> *Scan() override;
	bool DetCoin(std::bitset<defaultSizeFTimeRegion> &flagcollection);
	void Entry(uint32_t fem, int ch, int offset) override;
	void Entry(uint32_t fem, int ch, int offset, uint32_t leftwidth, uint32_t rightwidth) override;
	void ClearEntry() override;
	bool isUseUserDefinedCoin = false; // false: default, true: user defined coin (to be implemented)
	void SetTExpression(std::string &tx) override;
	
	using LogiCalc::Calc;
	bool Calc(std::bitset<defaultSizeFTimeRegion> &flagcollection);

private:
	// std::string tx;
	void SetWorkFlowCalc();
	using OpFunc = void(*)(TriggerBitSet*, const std::bitset<defaultSizeFTimeRegion> &, int);
	struct OpWithArg{
		OpFunc exec;
		int sig;
	};
	std::vector<OpWithArg> fWorkFlowCalc;
	using LogiCalc::ExtractBit;
	bool ExtractBit(const std::bitset<defaultSizeFTimeRegion> &flagcollection, int digit);

	std::bitset<defaultSizeFTimeRegion>* fTimeRegion = nullptr;
	std::bitset<defaultSizeFTimeRegion> fEntryMask = ~std::bitset<defaultSizeFTimeRegion>(0x00000000);

	std::map< uint32_t, std::vector<std::bitset<defaultSizeFTimeRegion>> > fEntryChBit;

	static void OpPush(TriggerBitSet* self, const std::bitset<defaultSizeFTimeRegion>& flagcollection, int sig){
		self->fStack.push(self->ExtractBit(flagcollection, sig));
	};
	static void OpAnd(TriggerBitSet* self, const std::bitset<defaultSizeFTimeRegion>& flagcollection, int){
		auto v1 = self->fStack.top(); self->fStack.pop();
		auto v2 = self->fStack.top(); self->fStack.pop();
		self->fStack.push(v1 & v2);
	};
	static void OpOr(TriggerBitSet* self, const std::bitset<defaultSizeFTimeRegion>& flagcollection, int){
		auto v1 = self->fStack.top(); self->fStack.pop();
		auto v2 = self->fStack.top(); self->fStack.pop();
		self->fStack.push(v1 | v2);
	};
	static void OpSwap(TriggerBitSet* self, const std::bitset<defaultSizeFTimeRegion>& flagcollection, int){
		auto vTop = self->fStack.top(); self->fStack.pop();
		auto vBelowTheTop = self->fStack.top(); self->fStack.pop();
		self->fStack.push(vTop);
		self->fStack.push(vBelowTheTop);
	};
	static void OpNot(TriggerBitSet* self, const std::bitset<defaultSizeFTimeRegion>& flagcollection, int){
		auto vTop = self->fStack.top(); self->fStack.pop();
		self->fStack.push(!vTop);
	};
	static void OpPop(TriggerBitSet* self, const std::bitset<defaultSizeFTimeRegion>& flagcollection, int){
		self->fStack.pop();
	};



}; // class TriggerBitSet : public Trigger, public LogiCalc




Trigger::Trigger()
{
	return;
}

Trigger::~Trigger()
{
	return;
}

void Trigger::InitParam() {};
bool Trigger::SetTimeRegion(int) {};
void Trigger::CleanUpTimeRegion() {};
void Trigger::Entry(uint32_t, int, int) {};
void Trigger::Entry(uint32_t, int, int, uint32_t, uint32_t) {};
void Trigger::Mark(unsigned char *, int, int, uint32_t) {};
std::vector<uint32_t> *Trigger::Scan() {};


void Trigger32::SetTExpression(std::string &tx)
{
	#if 1
	std::cout << "\n[Trigger32::SetTExpression] start makeing LUT" << std::endl;
	#endif
	auto t0 = std::chrono::high_resolution_clock::now();
	MakeTable(tx); // convert mq-param to RPN, and make LUT
	auto t1 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
	#if 1
	std::cout << "[Trigger32::SetTExpression] finished makeing LUT in " << duration << " ms" << std::endl;
	double lut_size = std::pow(2.0, fEntryCounts) / (8 * 1024.0); // LUT size in KB
	std::cout << "[Trigger32::SetTExpression] LUT size: " << lut_size << " KB" << std::endl;
	#endif

	return;
}

void TriggerBitSet::SetTExpression(std::string &tx)
{
	fCommands = SetFormula(tx); // class LogiCalc member variable std::vector<std::string> fCommands, LogiCalc::SetFormula()
	#if 1
	std::cout << "\n[TriggerBitSet::SetTExpression] Commands:";
	for (auto & com : fCommands) std::cout << " " << com;
	std::cout << std::endl;
	#endif
	SetWorkFlowCalc(); // make job list for evaluation of a set of flags of hitmap, stored in TriggerBitSet::fWorkFlowCalc
	return;
}


void Trigger32::InitParam()
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

void TriggerBitSet::InitParam()
{
	//fMarkCount = 0;
	//fMarkMask = 0;
	fHits.clear();
	fHits.resize(0);
	if (fTimeRegion != nullptr) {
		for (uint32_t i = 0 ; i < fTimeRegionSize ; i++) fTimeRegion[i].reset();
	}

	// SetFormula
	// fCommands = SetFormula(tx);

	return;
}

bool Trigger32::SetTimeRegion(int size)
{
	fTimeRegionSize = size;
	if (fTimeRegion != nullptr) {
		delete[] fTimeRegion;
		fTimeRegion = nullptr;
	}
	fTimeRegion = new uint32_t[size];

	return true;
} // bool Trigger32::SetTimeRegion(int size)

bool TriggerBitSet::SetTimeRegion(int size)
{
	fTimeRegionSize = size;
	if (fTimeRegion != nullptr) {
		delete[] fTimeRegion;
		fTimeRegion = nullptr;
	}
	fTimeRegion = new std::bitset<defaultSizeFTimeRegion>[size];

	return true;
} // bool TriggerBitSet::SetTimeRegion(int size)

void Trigger32::CleanUpTimeRegion()
{
	for (uint32_t i = 0 ; i < fTimeRegionSize ; i++) fTimeRegion[i] = 0;
	return;
}

void TriggerBitSet::CleanUpTimeRegion()
{
	for (uint32_t i = 0 ; i < fTimeRegionSize ; i++) fTimeRegion[i].reset();
	return;
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

void Trigger32::Entry(uint32_t fem, int ch, int offset)
{

	fEntryCh[fem].emplace_back(ch);
	fEntryChDelay[fem].emplace_back(offset);
	fEntryChBit[fem].emplace_back(0x0000001 << fEntryCounts);
	fEntryChLeftWidth[fem].emplace_back(fMarkLen / 2); // if mq-param give a signal(fem, ch, offset) with 3 parameters, leftwidth and rightwidth are set to default value (MarkLen / 2)
	fEntryChRightWidth[fem].emplace_back(fMarkLen / 2); // if mq-param give a signal(fem, ch, offset) with 3 parameters, leftwidth and rightwidth are set to default value (MarkLen / 2)
	fEntryMask |= 0x00000001 << fEntryCounts;
	fEntryCounts++;

	#if 0
	std::cout << "#D Trig Entry : Module: " << fem << " Ch: " << ch << std::endl;
	#endif
	if (static_cast<unsigned int>(fEntryCounts) > (sizeof(uint32_t) * 8)) {
		std::cerr << "Entry Ch. exceed " << sizeof(uint32_t) * 8<< std::endl;
	}
	assert(fEntryCounts <= static_cast<int>(sizeof(uint32_t) * 8));

	return;
} // void Trigger32::Entry(uint32_t fem, int ch, int offset)

void Trigger32::Entry(uint32_t fem, int ch, int offset, uint32_t leftwidth, uint32_t rightwidth)
{

	fEntryCh[fem].emplace_back(ch);
	fEntryChDelay[fem].emplace_back(offset);
	fEntryChBit[fem].emplace_back(0x0000001 << fEntryCounts);
	fEntryChLeftWidth[fem].emplace_back(leftwidth);
	fEntryChRightWidth[fem].emplace_back(rightwidth);
	fEntryMask |= 0x00000001 << fEntryCounts;
	fEntryCounts++;

	#if 0
	std::cout << "#D Trig Entry : Module: " << fem << " Ch: " << ch << std::endl;
	#endif
	if (static_cast<unsigned int>(fEntryCounts) > (sizeof(uint32_t) * 8)) {
		std::cerr << "Entry Ch. exceed " << sizeof(uint32_t) * 8<< std::endl;
	}
	assert(fEntryCounts <= static_cast<int>(sizeof(uint32_t) * 8));

	return;
} // void Trigger32::Entry(uint32_t fem, int ch, int offset, uint32_t leftwidth, uint32_t rightwidth)

void TriggerBitSet::Entry(uint32_t fem, int ch, int offset)
{
	fEntryCh[fem].emplace_back(ch);
	fEntryChDelay[fem].emplace_back(offset);
	fEntryChBit[fem].emplace_back(std::bitset<defaultSizeFTimeRegion>(0x1 << fEntryCounts));
	fEntryChLeftWidth[fem].emplace_back(fMarkLen / 2); // if mq-param give a signal(fem, ch, offset) with 3 parameters, leftwidth and rightwidth are set to default value (MarkLen / 2)
	fEntryChRightWidth[fem].emplace_back(fMarkLen / 2); // if mq-param give a signal(fem, ch, offset) with 3 parameters, leftwidth and rightwidth are set to default value (MarkLen / 2)
	fEntryMask |= std::bitset<defaultSizeFTimeRegion>(0x1 << fEntryCounts);
	fEntryCounts++;

	#if 0
	std::cout << "#D Trig Entry : Module: " << fem << " Ch: " << ch << std::endl;
	#endif

	return;
} // void TriggerBitSet::Entry(uint32_t fem, int ch, int offset)

void TriggerBitSet::Entry(uint32_t fem, int ch, int offset, uint32_t leftwidth, uint32_t rightwidth)
{
	fEntryCh[fem].emplace_back(ch);
	fEntryChDelay[fem].emplace_back(offset);
	fEntryChBit[fem].emplace_back(std::bitset<defaultSizeFTimeRegion>(0x1 << fEntryCounts));
	fEntryChLeftWidth[fem].emplace_back(leftwidth);
	fEntryChRightWidth[fem].emplace_back(rightwidth);
	fEntryMask |= std::bitset<defaultSizeFTimeRegion>(0x1 << fEntryCounts);
	fEntryCounts++;

	#if 0
	std::cout << "#D Trig Entry : Module: " << fem << " Ch: " << ch << std::endl;
	#endif
	if (static_cast<unsigned int>(fEntryCounts) > (sizeof(uint32_t) * 8)) {
		std::cerr << "Entry Ch. exceed " << sizeof(uint32_t) * 8<< std::endl;
	}
	assert(fEntryCounts <= static_cast<int>(sizeof(uint32_t) * 8));

	return;
} // void TriggerBitSet::Entry(uint32_t fem, int ch, int offset, uint32_t leftwidth, uint32_t rightwidth)

void Trigger::ClearEntry() {};
void Trigger32::ClearEntry()
{
	fEntryCh.clear();
	fEntryChDelay.clear();
	fEntryChBit.clear();
	fEntryMask = 0x00000000;
	fEntryCounts = 0;
	fEntryChLeftWidth.clear();
	fEntryChRightWidth.clear();

	return;
} // void Trigger32::ClearEntry()

void TriggerBitSet::ClearEntry()
{
	fEntryChBit.clear();
	fEntryChDelay.clear();
	fEntryMask = std::bitset<defaultSizeFTimeRegion>(0x00000000);
	fEntryCounts = 0;
	fEntryChLeftWidth.clear();
	fEntryChRightWidth.clear();

	return;
} // void TriggerBitSet::ClearEntry()

bool Trigger::CheckEntryFEM(uint32_t fem)
{
	bool rval = false;
	if (fEntryCh.count(fem) >= 1) rval = true;
	return rval;
}

void Trigger32::Mark(unsigned char *pdata, int len, int fem, uint32_t type)
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
								for (int k = -1 * leftwidth ; k < (rightwidth + 1) ; k++) {
									if ((hit + k) < fTimeRegionSize) {
										fTimeRegion[hit + k] |= markbit;
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
								for (int k = -1 * leftwidth ; k < (rightwidth + 1) ; k++) {
									if ((hit + k) < fTimeRegionSize) {
										fTimeRegion[hit + k] |= markbit;
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

							if (hit < fTimeRegionSize - leftwidth) {
								for (int k = -1 * leftwidth ; k < (rightwidth + 1) ; k++) {
									if ((hit + k) < fTimeRegionSize) {
										fTimeRegion[hit + k] |= markbit;
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
								for (int k = -1 * leftwidth ; k < (rightwidth + 1) ; k++) {
									if ((hit + k) < fTimeRegionSize) {
										fTimeRegion[hit + k] |= markbit;
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
} // void Trigger32::Mark(unsigned char *pdata, int len, int fem, uint32_t type)

void TriggerBitSet::Mark(unsigned char *pdata, int len, int fem, uint32_t type)
{
	if (fEntryCh.count(fem) >= 1) {
		uint64_t *tdcval;
		tdcval = reinterpret_cast<uint64_t *>(pdata);

		for(unsigned int i = 0; i < fEntryCh[fem].size(); i++){
			int ch = fEntryCh[fem][i];
			int delay = fEntryChDelay[fem][i];
			uint32_t leftwidth = fEntryChLeftWidth[fem][i];
			uint32_t rightwidth = fEntryChRightWidth[fem][i];
			std::bitset<defaultSizeFTimeRegion> markbit = fEntryChBit[fem][i];

			for(unsigned int j = 0; j < (len / sizeof(uint64_t)); j++){
				if(type == SubTimeFrame::TDC64H){
					struct TDC64H::tdc64 tdc;
					if(TDC64H::Unpack(tdcval[j], &tdc) == TDC64H::T_TDC){
						if(tdc.ch == ch){
							uint32_t hit = tdc.tdc4n + delay;

							if(hit < fTimeRegionSize - leftwidth){
								for(int k = -1 * leftwidth; k < (rightwidth + 1); k++){
									if((hit + k) < fTimeRegionSize){
										fTimeRegion[hit + k] |= markbit;
									}else if((static_cast<int>(hit) + k) >= 0){
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
				} // endif(type == SubTimeFrame::TDC64H)
				else if(type == SubTimeFrame::TDC64L){
					struct TDC64L::tdc64 tdc;
					if(TDC64L::Unpack(tdcval[j], &tdc) == TDC64L::T_TDC){
						if(tdc.ch == ch){
							uint32_t hit = tdc.tdc4n + delay;

							if(hit < fTimeRegionSize - leftwidth){
								for(int k = -1 * leftwidth; k < (rightwidth + 1); k++){
									if((hit + k) < fTimeRegionSize){
										fTimeRegion[hit + k] |= markbit;
									}else if((static_cast<int>(hit) + k) >= 0){
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
				} // endif(type == SubTimeFrame::TDC64L)
				else if(type == SubTimeFrame::TDC64H_V3){
					struct TDC64H_V3::tdc64 tdc;
					if(TDC64H_V3::Unpack(tdcval[j], &tdc) == TDC64H_V3::T_TDC){
						if(tdc.ch == ch){
							uint32_t hit = tdc.tdc4n + delay;

							if(hit < fTimeRegionSize - leftwidth){
								for(int k = -1 * leftwidth; k < (rightwidth + 1); k++){
									if((hit + k) < fTimeRegionSize){
										fTimeRegion[hit + k] |= markbit;
									}else if((static_cast<int>(hit) + k) >= 0){
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
				} // endif(type == SubTimeFrame::TDC64H_V3)
				else if(type == SubTimeFrame::TDC64L_V3){
					struct TDC64L_V3::tdc64 tdc;
					if(TDC64L_V3::Unpack(tdcval[j], &tdc) == TDC64L_V3::T_TDC){
						if(tdc.ch == ch){
							uint32_t hit = tdc.tdc4n + delay;

							if(hit < fTimeRegionSize - leftwidth){
								for(int k = -1 * leftwidth; k < (rightwidth + 1); k++){
									if((hit + k) < fTimeRegionSize){
										fTimeRegion[hit + k] |= markbit;
									}else if((static_cast<int>(hit) + k) >= 0){
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
		}
	}

	return;
} // void TriggerBitSet::Mark(unsigned char *pdata, int len, int fem, uint32_t type)

std::vector<uint32_t> *Trigger32::Scan()
{
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
} // std::vector<uint32_t> *Trigger32::Scan()

std::vector<uint32_t> *TriggerBitSet::Scan()
{
	fHits.clear();
	fHits.resize(0);
	
	if(isUseUserDefinedCoin){ // user defined coin (to be implemented)
		// to be implemented
	}
	else{ // default coin
		for(unsigned int i = 0; i < fTimeRegionSize; i++){
			if ((! DetCoin(fTimeRegion[i]))
			&& (DetCoin(fTimeRegion[i + 1]))) {
				fHits.emplace_back(i + 1);
			}
		}
	}

	return &fHits;

} // std::vector<uint32_t> *TriggerBitSet::Scan()

bool TriggerBitSet::DetCoin(std::bitset<defaultSizeFTimeRegion> &flagcollection){
	for(auto& step : fWorkFlowCalc){
		step.exec(this, flagcollection, step.sig);
	}

	bool ret;
	if(fStack.size() > 0){
		ret = fStack.top(); fStack.pop();
	}
	else{
		std::cout << "#E empty stack result" << std::endl;
		ret = false;
	}
	return ret;
} // bool TriggerBitSet::DetCoin(std::bitset<defaultSizeFTimeRegion> &flagcollection)

std::vector<uint32_t> *Trigger::Exec(std::vector<struct HBFIndex> &hbf)
{
	Trigger::CleanUpTimeRegion();

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
} // std::vector<uint32_t> *Trigger::Exec(std::vector<struct HBFIndex> &hbf)

bool TriggerBitSet::Calc(std::bitset<defaultSizeFTimeRegion> &flagcollection){
	
	return false;
} // bool TriggerBitSet::Calc(std::bitset<defaultSizeFTimeRegion> &flagcollection)

bool TriggerBitSet::ExtractBit(const std::bitset<defaultSizeFTimeRegion> &flagcollection, int digit){
	if (digit < defaultSizeFTimeRegion){
		bool bit = flagcollection[digit];
		return bit;
	}
	else{
		return false;
	}
} // bool TriggerBitSet::ExtractBit(std::bitset<defaultSizeFTimeRegion> &flagcollection, int digit)

void TriggerBitSet::SetWorkFlowCalc(){
	while(!fStack.empty()) fStack.pop();

	int nStack = 0;
	bool handled = false;
	for(auto& com : fCommands){
		handled = false;
		if(std::all_of(com.cbegin(), com.cend(), isdigit)){
			int sig = atoi(com.c_str());
			fWorkFlowCalc.push_back({OpPush, sig});
			nStack++;
			handled = true;
		}
		else{
			if(nStack > 1){
				if((com == "&") || (com == "*")){
					fWorkFlowCalc.push_back({OpAnd, 0});
					nStack--;
				}
				else if((com == "|") || (com == "+")){
					fWorkFlowCalc.push_back({OpOr, 0});
					nStack--;
				}
				else if(com == "x"){
					fWorkFlowCalc.push_back({OpSwap, 0});
				}
				else{
					std::cout << "\n#E Invalid command: " << com << " not enough operands in stack" << std::endl;
				}
				handled = true;
			}
			if((com == "!") || (com == "^")){
				fWorkFlowCalc.push_back({OpNot, 0});
				handled = true;
			}
			else if((com == "d") || (com == "p")){
				if(nStack > 0){
					fWorkFlowCalc.push_back({OpPop, 0});
					nStack--;
				}
				else{
					std::cout << "\n#E empty stack" << std::endl;
				}
				handled = true;
			}
			else if(!handled){
				std::cout << "\n#E Unknown command: " << com << std::endl;
			}
		}
	}

	return;
} // void TriggerBitSet::SetWorkFlowCalc





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
