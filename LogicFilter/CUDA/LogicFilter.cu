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
#include <vector>

#include <unordered_map>
#include <unordered_set>

#include <cstring>

#include <fstream>
#include <sstream>
#include <ctime>

#include "utility/MessageUtil.h"

//#include "HulStrTdcData.h"
#include "SubTimeFrameHeader.h"
#include "TimeFrameHeader.h"
#include "HeartbeatFrameHeader.h"
#include "FilterHeader.h"

// #include "SignalParser.cxx"
// #include "KTimer.cxx"
// #include "Trigger.cxx"



//std::atomic<int> gQdepth = 0;

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
		static constexpr std::string_view TriggerExpression     {"trigger-expression"};
		static constexpr std::string_view TriggerWidth       {"trigger-width"};
	};

	struct DataBlock {
		uint32_t femId;
		uint32_t Type;
		uint32_t timeFrameId;
		int HBFrame;
		bool is_HB;
		int msg_index;
		int nTrig;
	};

	LogicFilter()
	{

	}

	void InitTask() override;
	bool ConditionalRun() override;
	void PostRun() override;

private:
	std::string fInputChannelName;
	std::string fOutputChannelName;
	std::string fDQMChannelName;
	std::string fName;

	int fNumDestination {0};
	uint32_t fDirection {0};
	int fPollTimeoutMS  {0};
	int fSplitMethod    {2};

};

__global__ void testKernel()
{
    int workerIndex = blockIdx.x * blockDim.x + threadIdx.x;
    printf("Hello from worker %d\n", workerIndex);
}


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

	fName = fConfig->GetProperty<std::string>("id");
	std::istringstream ss(fName.substr(fName.rfind("-") + 1));
	ss >> fId;

	fSplitMethod    = std::stoi(
		fConfig->GetProperty<std::string>(opt::SplitMethod.data()));
	LOG(info) << "InitTask: SplitMethod : " << fSplitMethod;

} // void LogicFilter::InitTask()



bool LogicFilter::ConditionalRun()
{
    const std::string_view funcName = "[LogicFilter::ConditionalRun] ";
	//Receive
	FairMQParts inParts;
    Receive(inParts, fInputChannelName, 0, 1000);

    std::cout << std::endl;
    std::cout << funcName << "Called" << std::endl;
    testKernel<<<1, 16>>>();
    cudaDeviceSynchronize();
    std::cout << funcName << "Kernel finished" << std::endl;


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