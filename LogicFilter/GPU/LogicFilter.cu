/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *										*
 *	      This software is distributed under the terms of the		*
 *	      GNU Lesser General Public Licence (LGPL) version 3,		*
 *		  copied verbatim in the file "LICENSE"				*
 ********************************************************************************/

#include <fairmq/Device.h>
#include <fairmq/runDevice.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <stdexcept>
#include <unordered_set>
#include <cassert>

#include <cstring>

#include "utility/MessageUtil.h"

#include "SubTimeFrameHeader.h"
#include "TimeFrameHeader.h"
#include "HeartbeatFrameHeader.h"
#include "FilterHeader.h"

#include "SignalParser.h"
#include "Trigger.cxx"
#include "KTimer.cxx"

#include <cudamap.h> // dependency: https://github.com/SATSUKIUM/light-map-in-cuda.git

namespace bpo = boost::program_options;

struct LogicFilter : fair::mq::Device
{
	struct OptionKey {
		static constexpr std::string_view InputChannelName   {"in-chan-name"};
		static constexpr std::string_view OutputChannelName  {"out-chan-name"};
		static constexpr std::string_view DQMChannelName     {"dqm-chan-name"};
		static constexpr std::string_view DataSuppress       {"data-suppress"};
		static constexpr std::string_view RemoveHB           {"remove-hb"};
		static constexpr std::string_view PollTimeout        {"poll-timeout"};
		static constexpr std::string_view SplitMethod        {"split"};

		static constexpr std::string_view TriggerSignals     {"trigger-signals"};
		static constexpr std::string_view TriggerExpression  {"trigger-expression"};
		static constexpr std::string_view TriggerWidth       {"trigger-width"};
	};

	struct DataBlock {
		uint32_t femId;
		uint32_t Type;
		bool is_HB;
		int msg_index;
	};

	LogicFilter()
	{
		// register a handler for data arriving on "data" channel

		fTrig = new Trigger();
		fKt1 = new KTimer(1000);
		fKt3 = new KTimer(1000);
		fKt4 = new KTimer(1000);
	}

    ~LogicFilter() override {
        delete fTrig;
        delete fKt1;
        delete fKt3;
        delete fKt4;
        cudaFreeHost(fNEntryInSubTimeRegion);
        cudaFreeHost(fiSubTimeRegionMainTimeRegion);
    }
	void InitTask() override;

	
	bool ConditionalRun() override;
	void PostRun() override;

private:
	int IsHartBeat(uint64_t, uint32_t);
	void CheckMultiPart(FairMQParts &);
	int MakeBlockMap(FairMQParts &,
		std::vector< std::vector<struct DataBlock> > &);
	int BuildHBF(
		std::vector<struct HBFIndex> &,
		std::vector< std::vector<struct DataBlock> > &,
		FairMQParts &,
		int);
	int AddFilterMessage(
		FairMQParts &,
		std::vector< std::vector<uint32_t> > &,
		uint32_t,
		uint32_t);
	
	int MarkFlagSending(
		std::vector< std::vector<struct DataBlock> > &block_map,
		int i,
		int nhits,
		std::vector<bool> &flag_sending);

	std::string fInputChannelName;
	std::string fOutputChannelName;
	std::string fDQMChannelName;

	int fNumDestination {0};
	uint32_t fDirection {0};
	int fPollTimeoutMS  {0};
	int fSplitMethod    {2};

	uint32_t fId {0};
	Trigger *fTrig;
	bool fIsDataSuppress = true;
	bool fIsRemoveHB = false;

	int fHBflag = 0;

	KTimer *fKt1;
	KTimer *fKt3;
	KTimer *fKt4;

	uint32_t* fNEntryInSubTimeRegion = nullptr; // key: iSubTR, value: # of signals in iSubTR
	int fNSubTimeRegions {0};
	uint32_t* fiSubTimeRegionMainTimeRegion = nullptr; // key: iSubTR, value: iMainTR

};

void LogicFilter::InitTask()
{
	using opt = OptionKey;

	fInputChannelName  = fConfig->GetValue<std::string>(opt::InputChannelName.data());
	fOutputChannelName = fConfig->GetValue<std::string>(opt::OutputChannelName.data());
	fDQMChannelName    = fConfig->GetValue<std::string>(opt::DQMChannelName.data());

	LOG(info) << "InitTask: Input Channel : " << fInputChannelName
		<< " Output Channel : " << fOutputChannelName;

	fNumDestination = GetNumSubChannels(fOutputChannelName);
	fPollTimeoutMS  = std::stoi(fConfig->GetProperty<std::string>(opt::PollTimeout.data()));

	const auto name = fConfig->GetProperty<std::string>("id");
	std::istringstream ss(name.substr(name.rfind("-") + 1));
	ss >> fId;

	fSplitMethod    = std::stoi(
		fConfig->GetProperty<std::string>(opt::SplitMethod.data()));
	LOG(info) << "InitTask: SplitMethod : " << fSplitMethod;

	std::string sIsDataSuppress = fConfig->GetValue<std::string>(opt::DataSuppress.data());
	if (sIsDataSuppress == "true") {
		fIsDataSuppress = true;
	} else {
		fIsDataSuppress = false;
	}
	LOG(info) << "InitTask: DataSuppress : " << fIsDataSuppress;

	std::string sIsRemoveHB = fConfig->GetValue<std::string>(opt::RemoveHB.data());
	if (sIsRemoveHB == "true") {
		fIsRemoveHB = true;
	} else {
		fIsRemoveHB = false;
	}
	LOG(info) << "InitTask: RemoveHB : " << fIsRemoveHB;

	constexpr int timeRegionSize = 1024 * 128;
	fTrig->SetTimeRegion(timeRegionSize);

	std::string str_signals = fConfig->GetProperty<std::string>(opt::TriggerSignals.data());
	std::string expression = fConfig->GetProperty<std::string>(opt::TriggerExpression.data());
	int window_width = std::stoi(fConfig->GetProperty<std::string>(opt::TriggerWidth.data()));

	LOG(info) << "Trigger windows width: " << window_width;
	if (window_width <= 0 || window_width > UINT8_MAX)
        throw std::invalid_argument("Trigger width must be 1..255");

	std::vector<struct SignalParser::psig> signals = SignalParser::Parsing(str_signals);

	// --------------------------------
	// distribute signals into groups, and scan max group id
	// --------------------------------
	std::cout << "\n[LogicFilter::InitTask] Distributing signals into groups... " << std::endl;
	std::vector<std::vector<struct SignalParser::psig>> groups;
	uint32_t max_group_id = 0;
    std::unordered_set<uint64_t> registeredAddresses;
	for(auto &sig : signals){
		if (!registeredAddresses.insert((uint64_t(sig.femId) << 32) | sig.channel).second)
            throw std::invalid_argument("FEEAddrDecoder supports only one configuration per FEM/channel");
        if (sig.group_id >= 32 || sig.channel > UINT8_MAX || sig.offset < INT16_MIN || sig.offset > INT16_MAX)
            throw std::invalid_argument("Signal exceeds decoder field limits");
        if(sig.group_id > max_group_id) max_group_id = sig.group_id;
	}
	groups.resize(max_group_id + 1);
	for(auto &sig : signals){
		groups[sig.group_id].emplace_back(sig);
	}
	std::cout << "\tMax group id: " << max_group_id << std::endl;

	// --------------------------------
	// sort signals in each group by subgroup_id for making SubTimeRegion in order of subgroup_id
	// --------------------------------
	std::cout << "\n[LogicFilter::InitTask] Sorting signals in each group by subgroup_id... " << std::endl;
	for(auto &group : groups){
		if(group.size() == 0) continue; // skip empty group
		std::sort(group.begin(), group.end(),
			[](const struct SignalParser::psig &left, const struct SignalParser::psig &right)
			{
				if(left.subgroup_id == right.subgroup_id){
					return left.femId < right.femId;
				}
				else{
					return left.subgroup_id < right.subgroup_id;
				}
			});
	}
	std::cout << "\tDone." << std::endl;

	// --------------------------------
	// check and print group information
	// --------------------------------
	std::cout << "\n[LogicFilter::InitTask] Check group information: " << std::endl;
	for(size_t i = 0; i < max_group_id + 1; i++){
		if(groups[i].size() > 0){
			std::cout << "\tgroup " << i << " has " << groups[i].size() << std::endl;
			std::cout << "\t\tSignals: ";
			for(const auto &sig : groups[i]){
				if(sig.subgroup_id != SignalParser::NO_SUBGROUP){
					std::cout << sig.group_id << "-" << sig.subgroup_id << " ";
				}
				else{
					std::cout << sig.group_id << " ";
				}
			}
			std::cout << std::endl;
		} // endif(groups[i].size() > 0)
		else{
			std::cout << "\tgroup " << i << " is empty" << std::endl;
		}
	} // for(size_t i = 0; i < max_group_id + 1; i++)
	std::cout << "\tDone." << std::endl;
	std::cout << "\n[LogicFilter::InitTask] Trigger expression: " << expression << std::endl;

	// --------------------------------
	// initialize cudamap
	// --------------------------------
	std::cout << "\n[LogicFilter::InitTask] Initializing cudamap... " << std::endl;
    std::vector<cudamap::FEEAddr> feeAddrs;
    std::vector<cudamap::TriggerConfig> triggerConfigs;
	int iSubTimeRegion_ = -1; // subgroupが見つかれば、0から順番に割り当てる
	std::pair<uint32_t, uint32_t> last_subgroup_ = std::make_pair(UINT32_MAX, UINT32_MAX); // 最初にsubgroupが見つかれば必ず、このペアとは異なる
	std::vector<uint32_t> nEntryInSubTimeRegion; // key: iSubTR, value: # of signals in iSubTR
	std::vector<uint32_t> iSubTimeRegionMainTimeRegion; // key: iSubTR, value: iMainTR
	for(size_t iGroupId = 0; iGroupId < max_group_id + 1; ++iGroupId){
		auto &group = groups[iGroupId];
		if(group.size() == 0){
			std::cout << "\tgroup " << iGroupId << " is empty, skip." << std::endl;
			continue; // skip empty group
		}
		for(const auto &sig : group){
			if(sig.subgroup_id != SignalParser::NO_SUBGROUP){ // subgroupに属していれば
				if(last_subgroup_ != std::make_pair(sig.group_id, sig.subgroup_id)){ // 前のgroup-subgroupのペアと異なれば
					iSubTimeRegion_++;
					nEntryInSubTimeRegion.push_back(1);
					iSubTimeRegionMainTimeRegion.push_back(uint32_t(1u) << iGroupId);
					last_subgroup_ = std::make_pair(sig.group_id, sig.subgroup_id);
					std::cout << "\t\tgroup " << sig.group_id << " has subgroup " << sig.subgroup_id << " assigned to SubTimeRegion " << iSubTimeRegion_ << std::endl;
				}
				else{ // 前のgroup-subgroupのペアと同じであれば
					nEntryInSubTimeRegion[iSubTimeRegion_]++;
					std::cout << "\t\tgroup " << sig.group_id << " has subgroup " << sig.subgroup_id << " assigned to SubTimeRegion " << iSubTimeRegion_ << " (entry count in this SubTimeRegion: " << nEntryInSubTimeRegion[iSubTimeRegion_] << "th)" << std::endl;
				}
				std::cout << "\t\tRegistering signal: group_id=" << sig.group_id
					<< " subgroup_id=" << sig.subgroup_id
					<< " iSubTimeRegion=" << (sig.subgroup_id != SignalParser::NO_SUBGROUP ? iSubTimeRegion_ : SignalParser::NO_SUBGROUP)
					<< " femId=" << std::hex << sig.femId
					<< " channel=" << std::dec <<sig.channel
					<< " offset=" << sig.offset
					<< " window_width=" << window_width
					<< std::endl;
				if (nEntryInSubTimeRegion[iSubTimeRegion_] > 32 || iSubTimeRegion_ >= cudamap::NO_SUBGROUP)
                    throw std::invalid_argument("Too many subgroup signals or subgroups");
                cudamap::stageItem(sig.femId, static_cast<uint8_t>(sig.channel), static_cast<uint32_t>(1u << (nEntryInSubTimeRegion[iSubTimeRegion_] - 1)), static_cast<uint8_t>(window_width), static_cast<int16_t>(sig.offset), static_cast<uint8_t>(iSubTimeRegion_), feeAddrs, triggerConfigs);
			} // if(sig.subgroup_id != SignalParser::NO_SUBGROUP)
			else{ // subgroupに属していなければ
				std::cout << "\t\tRegistering signal: group_id=" << sig.group_id
					<< " subgroup_id=NO_SUBGROUP"
					<< " iSubTimeRegion=NO_SUBGROUP"
					<< " femId=" << std::hex << sig.femId
					<< " channel=" << std::dec <<sig.channel
					<< " offset=" << sig.offset
					<< " window_width=" << window_width
					<< std::endl;
				cudamap::stageItem(sig.femId, static_cast<uint8_t>(sig.channel), static_cast<uint32_t>(1u << sig.group_id), static_cast<uint8_t>(window_width), static_cast<int16_t>(sig.offset), static_cast<uint8_t>(cudamap::NO_SUBGROUP), feeAddrs, triggerConfigs);
			} // else(sig.subgroup_id != SignalParser::NO_SUBGROUP)
		} // for(const auto &sig : group)
	} // for(size_t iGroupId = 0; iGroupId < max_group_id + 1; ++iGroupId)

	fNSubTimeRegions = iSubTimeRegion_ + 1;
	// check for development
	if(fNSubTimeRegions != nEntryInSubTimeRegion.size()){
		std::cerr << "#E nSubTimeRegions != nEntryInSubTimeRegion.size()" << std::endl;
		std::cerr << "#E nSubTimeRegions: " << fNSubTimeRegions << std::endl;
		std::cerr << "#E nEntryInSubTimeRegion.size(): " << nEntryInSubTimeRegion.size() << std::endl;
		throw std::runtime_error("fNSubTimeRegions != nEntryInSubTimeRegion.size()");
	}
	if(fNSubTimeRegions != iSubTimeRegionMainTimeRegion.size()){
		std::cerr << "#E nSubTimeRegions != iSubTimeRegionMainTimeRegion.size()" << std::endl;
		std::cerr << "#E nSubTimeRegions: " << fNSubTimeRegions << std::endl;
		std::cerr << "#E iSubTimeRegionMainTimeRegion.size(): " << iSubTimeRegionMainTimeRegion.size() << std::endl;
		throw std::runtime_error("fNSubTimeRegions != iSubTimeRegionMainTimeRegion.size()");
	}
	cudaMallocHost(&fNEntryInSubTimeRegion, sizeof(uint32_t) * fNSubTimeRegions);
	for(size_t i=0; i<nEntryInSubTimeRegion.size(); ++i){
		fNEntryInSubTimeRegion[i] = nEntryInSubTimeRegion[i];
	}
	cudaMallocHost(&fiSubTimeRegionMainTimeRegion, sizeof(uint32_t) * fNSubTimeRegions);
	for(size_t i=0; i<iSubTimeRegionMainTimeRegion.size(); ++i){
		fiSubTimeRegionMainTimeRegion[i] = iSubTimeRegionMainTimeRegion[i];
	}
	std::cout << "\tNumber of SubTimeRegions: " << fNSubTimeRegions << std::endl;

	fTrig->InitializeFEEAddrDecoder(feeAddrs, triggerConfigs);

	std::cout << "\tDone." << std::endl;

	// --------------------------------
	// Initialize SubTimeRegions
	// --------------------------------
	std::cout << "\n[LogicFilter::InitTask] Setting SubTimeRegions... " << std::endl;
	fTrig->SetSubTimeRegion(timeRegionSize, fNEntryInSubTimeRegion, fNSubTimeRegions);
    fTrig->SetfiSubTimeRegionMainTimeRegion(fiSubTimeRegionMainTimeRegion, fNSubTimeRegions);
	std::cout << "\tDone." << std::endl;

	// --------------------------------
	// Make LUT
	// --------------------------------
	LOG(info) << "Trigger Expression: " << expression;
	fTrig->MakeTable(expression);
	std::cout << "\tMaking table Done." << std::endl;

}

int LogicFilter::IsHartBeat(uint64_t val, uint32_t type)
{
	int hbflag = -1;
	int hbframe = -1;
	if (type == SubTimeFrame::TDC64H) {
		struct TDC64H::tdc64 tdc;
		if (TDC64H::Unpack(val, &tdc) == TDC64H::T_HB) {
			hbframe = tdc.hartbeat;
			hbflag = tdc.flag;
		}
	} else if (type == SubTimeFrame::TDC64L) {
		struct TDC64L::tdc64 tdc;
		if (TDC64L::Unpack(val, &tdc) == TDC64L::T_HB) {
			hbframe = tdc.hartbeat;
			hbflag = tdc.flag;
		}
	} else if (type == SubTimeFrame::TDC64L_V1) {
		struct TDC64L::tdc64 tdc;
		if (TDC64L::v1::Unpack(val, &tdc) == TDC64L::v1::T_HB) {
			hbframe = tdc.hartbeat;
			hbflag = tdc.flag;
		}
	} else if (type == SubTimeFrame::TDC64H_V3) {
		struct TDC64H_V3::tdc64 tdc;
		if (TDC64H_V3::Unpack(val, &tdc) == TDC64H_V3::T_HB1) {
			hbframe = tdc.hartbeat;
			hbflag = tdc.flag;
		}
	} else if (type == SubTimeFrame::TDC64L_V3) {
		struct TDC64L_V3::tdc64 tdc;
		if (TDC64L_V3::Unpack(val, &tdc) == TDC64L_V3::T_HB1) {
			hbframe = tdc.hartbeat;
			hbflag = tdc.flag;
		}
	} else {
		std::cout << "Unknown device : " << std::hex << type << std::endl;
		const uint64_t* data = reinterpret_cast<const uint64_t*>(&val);
		for (size_t i = 0; i < sizeof(val) / sizeof(uint64_t); ++i) {
			std::cout << "Data[" << i << "]: " << std::hex << data[i] << std::endl;
		}
	}

	if (hbflag > 0) fHBflag |= hbflag;
	if (fKt4->Check()) {
		if (fHBflag > 0) {
			if ((fHBflag & 0x200) == 0x200) LOG(warn) << "HB Data lost";
			if ((fHBflag & 0x100) == 0x100) LOG(warn) << "HB Data confiliction";
			if ((fHBflag & 0x080) == 0x080) LOG(warn) << "HB LFN mismatch";
			if ((fHBflag & 0x040) == 0x040) LOG(warn) << "HB GFN mismatch";
		}
		fHBflag = 0;
	}

	return hbframe;
}

void LogicFilter::CheckMultiPart(FairMQParts &inParts)
{
	std::cout << "#Nmsg: " << std::dec << inParts.Size() << std::endl;
	std::cout << "# " << std::dec;
	for(int i = 0 ; i < inParts.Size() ; i++) {
		uint64_t *top = reinterpret_cast<uint64_t *>(inParts[i].GetData());
		struct TimeFrame::Header *tbh
			= reinterpret_cast<TimeFrame::Header *>(inParts[i].GetData());
		struct SubTimeFrame::Header *stbh
			= reinterpret_cast<SubTimeFrame::Header *>(inParts[i].GetData());
		struct HeartbeatFrame::Header *hbh
			= reinterpret_cast<HeartbeatFrame::Header *>(inParts[i].GetData());
		if (tbh->magic == TimeFrame::MAGIC) {
			std::cout << "TF Id: "
				<< tbh->timeFrameId  << " Nsrc: " << tbh->numSource;
		} else
		if (stbh->magic == SubTimeFrame::MAGIC) {
			std::cout << std::endl << "* STF: TFid: 0x" << std::hex
				<< std::setw(8) << stbh->timeFrameId << " :"; 
		} else
		if (hbh->magic == HeartbeatFrame::MAGIC) {
			std::cout << " HB:" << std::dec
				<< std::setw(8) << inParts[i].GetSize() << " :"; 
		} else {
			std::cout << "UK " << std::setw(8) << std::setfill('0')
				<< IsHartBeat(top[0], SubTimeFrame::TDC64H);
		}
	}
	std::cout << std::dec << std::endl;

	return;
}

int LogicFilter::MakeBlockMap(
	FairMQParts &inParts,
	std::vector< std::vector<struct DataBlock> > &block_map)
{
	uint64_t femid = 0;
	uint64_t devtype = 0;
	std::vector<struct DataBlock> blocks;

	for(int i = 0 ; i < inParts.Size() ; i++) {

		auto tfHeader = reinterpret_cast<struct TimeFrame::Header *>(inParts[i].GetData());
		auto stfHeader = reinterpret_cast<struct SubTimeFrame::Header *>(inParts[i].GetData());
		struct DataBlock dblock;

		if (tfHeader->magic == TimeFrame::MAGIC) {
			dblock.is_HB = false;
			dblock.msg_index = 0;
		} else
		if (stfHeader->magic == SubTimeFrame::MAGIC) {
			femid = stfHeader->femId;
			devtype = stfHeader->femType;
			if (blocks.size() > 0) block_map.push_back(blocks);
			dblock.is_HB = false;
			dblock.msg_index = -1;
			blocks.clear();
			blocks.resize(0);

		} else {
			// make block map;

			uint64_t *data = reinterpret_cast<uint64_t *>(inParts[i].GetData());
			int hbframe = IsHartBeat(data[0], devtype);

			//if (hbframe < 0) {
			if ((data[0] == HeartbeatFrame::MAGIC)
				|| (hbframe < 0)) {
				dblock.femId = femid;
				dblock.Type = devtype;
				dblock.is_HB = false;
				dblock.msg_index = i;
			} else {
				if (fSplitMethod == 1) {
					//data ga nakattatokimo push_back dummy
					if (   (blocks.size() == 0)
						|| ((blocks.size() > 0)
						&& (blocks.back().is_HB == true)) ) {

						dblock.femId = femid;
						dblock.Type = SubTimeFrame::NULDEV;
						dblock.is_HB = false;
						dblock.msg_index = -2;
						blocks.push_back(dblock);
					}
				}

				dblock.femId = femid;
				dblock.Type = devtype;
				dblock.is_HB = true;
				dblock.msg_index = i;
			}
			blocks.push_back(dblock);
		}
	} /// end of the for loop
	block_map.push_back(blocks);

	size_t bsize_min = block_map[0].size();
	for (auto& blk : block_map) {
		if (blk.size() < bsize_min) {
			LOG(warn) << "Unmatched number of stf in TF "
				<< bsize_min << " " << blk.size();
			bsize_min = blk.size();
		}
	}

	return 0;
}

int LogicFilter::BuildHBF(
	std::vector<struct HBFIndex> &hbf_list,
	std::vector< std::vector<struct DataBlock> > &block_map,
	FairMQParts &inParts,
	int imsg
)
{

	for (size_t iifem = 0 ; iifem < block_map.size() ; iifem++) {
		struct DataBlock *dbl = &block_map[iifem][imsg];
		uint64_t vfemid = dbl->femId;
		if (!(dbl->is_HB)) {
			int mindex = dbl->msg_index;

			if (mindex > 0) {
				struct HBFIndex hbf_seg;
				hbf_seg.femType = dbl->Type;
				hbf_seg.femId = vfemid;
				hbf_seg.data = reinterpret_cast<unsigned char *>(
					inParts[mindex].GetData());
				hbf_seg.size = inParts[mindex].GetSize();

				hbf_list.emplace_back(hbf_seg);
			}
		}
	}

	return 0;
}

int LogicFilter::AddFilterMessage(
	FairMQParts &outParts,
	std::vector< std::vector<uint32_t> > &fltdata,
	uint32_t elapse,
	uint32_t tf_id)
{
	int flt_data_len = 0;

	auto tp = std::chrono::system_clock::now() ;
	auto d = tp.time_since_epoch();
	uint64_t sec = std::chrono::duration_cast
		<std::chrono::seconds>(d).count();
	uint64_t usec = std::chrono::duration_cast
		<std::chrono::microseconds>(d).count();

	uint32_t flt_datasize = 0;
	uint32_t totalhits = 0;
	for (auto &v : fltdata) {
		flt_datasize += v.size() * sizeof(uint32_t)
			+ sizeof(struct Filter::TrgTimeHeader);
		totalhits += v.size();
	}
	uint64_t flt_len = flt_datasize + sizeof(struct Filter::Header);

	//add FLT Header
	auto fltHeader = std::make_unique<struct Filter::Header>();
	fltHeader->magic = Filter::MAGIC;
	fltHeader->length = flt_len;
	fltHeader->timeFrameId = tf_id;
	fltHeader->numTrigs = totalhits;
	fltHeader->workerId = fId;
	fltHeader->numMessages = 1 + fltdata.size();
	fltHeader->elapseTime = elapse;
	fltHeader->timeSec = sec;
	fltHeader->timeUSec = usec;
	outParts.AddPart(MessageUtil::NewMessage(*this, std::move(fltHeader)));

	flt_data_len += sizeof(struct Filter::Header);

	//add FLT data
	for (auto &v : fltdata) {
		//outParts.AddPart(MessageUtil::NewMessage(*this, std::move(trgtdc)));

		struct Filter::TrgTimeHeader trg_time_header;
		trg_time_header.length = sizeof(Filter::TrgTimeHeader)
			+ v.size() * sizeof(Filter::TrgTime);

		std::vector<uint32_t> vv;
		for (auto &h : trg_time_header.u32data) vv.emplace_back(h);
		for (auto &tdc : v) {
			uint32_t trg_type = 0xaa000000;
			vv.emplace_back(tdc);
			vv.emplace_back(trg_type);
		}
	
		outParts.AddPart(MessageUtil::NewMessage
			(*this, std::make_unique<std::vector<uint32_t>>(std::move(vv))));

		flt_data_len += outParts[outParts.Size() - 1].GetSize();
	}

	return flt_data_len;
}

int LogicFilter::MarkFlagSending(
	std::vector< std::vector<struct DataBlock> > &block_map,
	int i,
	int nhits,
	std::vector<bool> &flag_sending)
{
	for (size_t iifem = 0 ; iifem < block_map.size() ; iifem++) {
		if (nhits == 0) {
			int mindex = block_map[iifem][i].msg_index;
			bool is_HB = block_map[iifem][i].is_HB;
			uint32_t dtype = block_map[iifem][i].Type;

			if ((mindex > 0) && (! is_HB)) {
				if (dtype != SubTimeFrame::NULDEV) {
					flag_sending[mindex] = false;
				}

				if ((fSplitMethod == 0) || (fSplitMethod == 1)) {
					//// kokoha mondai gaaru ////
					if (fIsRemoveHB) {
						flag_sending[mindex + 1] = false;
					}
				}
			}

		} else {
		}
	}

	return 0;
}

bool LogicFilter::ConditionalRun()
{
	//Receive
	FairMQParts inParts;

	int i_tf_msg_index = 0;
	uint32_t tf_tf_id = 0;
	std::chrono::system_clock::time_point sw_start, sw_end;

	if (Receive(inParts, fInputChannelName, 0, 1000) > 0) {
		assert(inParts.Size() >= 2);
		sw_start = std::chrono::system_clock::now();

		if (fKt1->Check()) CheckMultiPart(inParts);

		std::vector< std::vector<struct DataBlock> > block_map;
		std::vector<bool> flag_sending;
		for(int i = 0 ; i < inParts.Size() ; i++) flag_sending.push_back(true);

		for(int i = 0 ; i < inParts.Size() ; i++) {
			auto tfHeader = reinterpret_cast<struct TimeFrame::Header *>(inParts[i].GetData());
			if (tfHeader->magic == TimeFrame::MAGIC) {
				i_tf_msg_index = i;
				tf_tf_id = tfHeader->timeFrameId;
			}
		}
		MakeBlockMap(inParts, block_map);

		size_t bsize_min = block_map[0].size();
		for (auto& blk : block_map) {
			if (blk.size() < bsize_min) {
				LOG(warn) << "Unmatched number of stf in TF "
					<< bsize_min << " " << blk.size();
				bsize_min = blk.size();
			}
		}

		// Trigger processing in unit of HBF
		std::vector< std::vector<uint32_t> > fltdata;
		int totalhits = 0;
		for (size_t i = 0 ; i < bsize_min ; i++) {

			//HBF building
			std::vector<struct HBFIndex> hbf_list;
			BuildHBF(hbf_list, block_map, inParts, i);

			//Trigger process
			fTrig->CleanUpTimeRegion();
			fTrig->CleanUpSubTimeRegion(fNEntryInSubTimeRegion, fNSubTimeRegions);
			std::vector<uint32_t> *hits = fTrig->Exec(hbf_list);
			fltdata.emplace_back(*hits);
			int nhits = hits->size();

			// mark sending flag for STF
			MarkFlagSending(block_map, i, nhits, flag_sending);

			totalhits += nhits;
		}

		sw_end = std::chrono::system_clock::now();
		uint32_t elapse = std::chrono::duration_cast<std::chrono::microseconds>(
			sw_end - sw_start).count();
		static uint64_t int_hits = 0;
		static uint64_t int_processed_hbf = 0;
		double trig_ratio;
		int_hits += totalhits;
		int_processed_hbf += bsize_min / 2;
		if (fKt3->Check()) {
			if (int_processed_hbf > 0) {
				trig_ratio = static_cast<double>(int_hits)
					/ static_cast<double>(int_processed_hbf);
			}

			std::cout << "#Elapse: " << std::dec << elapse << " us"
				<< " Hits: " << totalhits
				<< " T.Ratio: " << trig_ratio
				<< " Hits(inte): " << int_hits
				<< " HBF(inte): " << int_processed_hbf
				<< std::endl;

			//int_hits = 0;
			//int_processed_hbf = 0;
		}

		FairMQParts outParts;
  
		//Copy TimeFrameHeader
		if (i_tf_msg_index != 0) {
			LOG(error) << "Wrong TimeFrame Header locataion : " << i_tf_msg_index;
		}
		FairMQMessagePtr msgCopyTF(fTransportFactory->CreateMessage());
		msgCopyTF->Copy(inParts.AtRef(0));
		outParts.AddPart(std::move(msgCopyTF));

		//FilterHeader
		AddFilterMessage(outParts, fltdata, elapse, tf_tf_id);

		//Copy SubTimeFrame
		unsigned int msg_size = inParts.Size();
		for (unsigned int ii = 1 ; ii < msg_size ; ii++) {
			if (flag_sending[ii] || (! fIsDataSuppress)) {
				FairMQMessagePtr msgCopy(fTransportFactory->CreateMessage());
				msgCopy->Copy(inParts.AtRef(ii));
				outParts.AddPart(std::move(msgCopy));
			} else {
				// add blank SubTime Frame
				struct HeartbeatFrame::Header *hbf
					= reinterpret_cast<struct HeartbeatFrame::Header *>(inParts[ii].GetData());
				//	= reinterpret_cast<uint64_t *>(inParts[ii].GetData());
				if (hbf->magic == HeartbeatFrame::MAGIC) {
				//if (hbf_u64[0] == HeartbeatFrame::MAGIC) {
					#pragma pack(4)
					struct BlankHBF{
						struct HeartbeatFrame::Header Header;
						uint64_t Delimiter[2];
					};
					#pragma pack()
					auto blank = std::make_unique<struct BlankHBF>();
					uint64_t *d = reinterpret_cast<uint64_t *>(
						reinterpret_cast<char *>(inParts[ii].GetData())
						+ inParts[ii].GetSize() - (sizeof(uint64_t) * 2));
					std::memcpy(&(blank->Header), inParts[ii].GetData(),
						sizeof(struct HeartbeatFrame::Header));
					blank->Header.length = sizeof(HeartbeatFrame::Header)
						+ (sizeof(uint64_t) * 2);
					std::memcpy(&(blank->Delimiter), d, sizeof(uint64_t) * 2);
					outParts.AddPart(MessageUtil::NewMessage(*this, std::move(blank)));
					//flag_sending[ii] = true;

				}
			}
		}

		//Modify SubTimeFrameHeader
		if (fIsDataSuppress) {
			for (int ii = 0 ; ii < outParts.Size() ; ii++) {
				auto stfh = reinterpret_cast<struct SubTimeFrame::Header *>
					(outParts[ii].GetData());
				if (stfh->magic == SubTimeFrame::MAGIC) {
					uint32_t len_stf = 0;
					uint32_t nmsg_stf = 0;
					int kk = ii + 1;
					for (int jj = ii + 1 ; jj < outParts.Size() ; jj++) {
						auto sstf =
							reinterpret_cast<struct SubTimeFrame::Header *>
							(outParts[jj].GetData());
						if (sstf->magic == SubTimeFrame::MAGIC) {
							kk = jj;
							break;
						} else
						if (flag_sending[jj]) {
							len_stf += outParts[jj].GetSize();
							nmsg_stf++;
						}
					}
					// if (len_stf == 0) { // STF header kara nankunaru? OK?
					// 	flag_sending[ii] = false;
					// }
					stfh->length
						= len_stf + sizeof(struct SubTimeFrame::Header);
					stfh->numMessages = nmsg_stf;
					ii = kk - 1;
				}
			}
		}

		// rewrite TimeFrameHeader length
		auto tfHeader = reinterpret_cast<struct TimeFrame::Header *>(
			outParts[0].GetData());
		// Calcurate TimeFrame.length
		uint32_t tf_len = 0;
		for (auto &msg : outParts) tf_len += msg->GetSize();
		tfHeader->length = tf_len;

		FairMQParts dqmParts;
		bool dqmSocketExists = fChannels.count(fDQMChannelName);
		if (dqmSocketExists) {
			for (auto & m : outParts) {
				FairMQMessagePtr msgCopy(fTransportFactory->CreateMessage());
				msgCopy->Copy(*m);
				dqmParts.AddPart(std::move(msgCopy));
			}

			if (Send(dqmParts, fDQMChannelName) < 0) {
				if (NewStatePending()) {
					LOG(info) << "Device is not RUNNING";
				} else {
					LOG(error) << "Failed to enqueue dqm-channel";
				}
			} else {
			}
		} else {
		}

		//Send
		auto poller = NewPoller(fOutputChannelName);
		while (!NewStatePending()) {
			auto direction = (fDirection++) % fNumDestination;
			poller->Poll(fPollTimeoutMS);
			if (poller->CheckOutput(fOutputChannelName, direction)) {
				if (Send(outParts, fOutputChannelName, direction) > 0) {
					// successfully sent
					break;
				} else {
					LOG(error) << "Failed to queue output-channel";
				}
			}
			if (fNumDestination==1) {
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
			}
		}

	}

	return true;
}

void LogicFilter::PostRun()
{
	LOG(info) << "Post Run";
	return;
}

void addCustomOptions(bpo::options_description& options)
{
	using opt = LogicFilter::OptionKey;

	options.add_options()
		//("max-iterations", bpo::value<uint64_t>()->default_value(0),
		//"Maximum number of iterations of Run/ConditionalRun/OnData (0 - infinite)")
		(opt::InputChannelName.data(),
			bpo::value<std::string>()->default_value("in"),
			"Name of the input channel")
		(opt::OutputChannelName.data(),
			bpo::value<std::string>()->default_value("out"),
			"Name of the output channel")
		(opt::DQMChannelName.data(),
			bpo::value<std::string>()->default_value("dqm"),
			"Name of the data quality monitoring channel")
		(opt::DataSuppress.data(),
			bpo::value<std::string>()->default_value("false"),
			"Data suppression enable")
		(opt::RemoveHB.data(),
			bpo::value<std::string>()->default_value("false"),
			"Remove HB without hit")
		(opt::PollTimeout.data(), 
			bpo::value<std::string>()->default_value("1"),
			"Timeout of polling (in msec)")
		(opt::SplitMethod.data(),
			bpo::value<std::string>()->default_value("2"),
			"STF split method")

		(opt::TriggerSignals.data(),
			bpo::value<std::string>()->default_value(
			"(0 0xc0a802a9 0 0) (1 0xc0a802a9 1 0)"),
			"Triger signals (module_IP Channel_number Offset)")
		(opt::TriggerExpression.data(),
			bpo::value<std::string>()->default_value("RPN 0 1 &"),
			"Trigger expression")
		(opt::TriggerWidth.data(),
			bpo::value<std::string>()->default_value("10"),
			"Trigger window width (4 ns unit)")

    		;
}

std::unique_ptr<fair::mq::Device> getDevice(fair::mq::ProgOptions& /*config*/)
{
	return std::make_unique<LogicFilter>();
}
