#ifndef TRIGGER_H
#define TRIGGER_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

#include <string.h>
#include <assert.h>

#include <chrono>
#include <cmath>


#include "SubTimeFrameHeader.h"
#include "TriggerMap.h"
#include "LogiCalc.h"

#include <bitset>

#define BENCHMARK_MORE32 1
#include <chrono>


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
    virtual void SetIsUseUserDefinedCoin(bool val) = 0;
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
	
	bool DetCoin_as_you_like(std::bitset<defaultSizeFTimeRegion> &flagcollection);
	void SetTExpression(std::string &tx) override;
	
	using LogiCalc::Calc;
	bool Calc(std::bitset<defaultSizeFTimeRegion> &flagcollection);
    void SetIsUseUserDefinedCoin(bool val){isUseUserDefinedCoin = val;};

private:
    bool isUseUserDefinedCoin = false; // false: default, true: user defined coin (to be implemented)
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

#endif // TRIGGER_H