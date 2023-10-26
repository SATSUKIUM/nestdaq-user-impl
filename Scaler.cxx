#include <chrono>
#include <thread>
#include <future>
#include <map>
#include <iostream>
#include <boost/format.hpp>
#include <cmath>

#include <fairmq/runDevice.h>

#include "TimeFrameHeader.h"
#include "SubTimeFrameHeader.h"
#include "STFBuilder.h"
#include "AmQStrTdcData.h"
#include "utility/Reporter.h"
#include "utility/HexDump.h"
#include "utility/MessageUtil.h"

#include "Scaler.h"


namespace bpo = boost::program_options;

//______________________________________________________________________________

void addCustomOptions(bpo::options_description& options)
{
  {
    using opt = Scaler::OptionKey;
  options.add_options()
    (opt::NumSource.data(),          bpo::value<std::string>()->default_value("1"), "Number of source endpoint")
    (opt::BufferTimeoutInMs.data(),  bpo::value<std::string>()->default_value("100000"), "Buffer timeout in milliseconds")
    (opt::RunNumber.data(),          bpo::value<std::string>()->default_value("1"), "Run number (integer) given by DAQ plugin")
    (opt::InputChannelName.data(),   bpo::value<std::string>()->default_value("in"), "Name of the input channel")
    (opt::OutputChannelName.data(),  bpo::value<std::string>()->default_value("out"), "Name of the output channel")    
    (opt::UpdateInterval.data(),     bpo::value<std::string>()->default_value("1000"), "Canvas update rate in milliseconds")
    (opt::ServerUri.data(),          bpo::value<std::string>()->default_value("tcp://127.0.0.1:6379/3"), "Redis server URI (if empty, the same URI of the service registry is used.)")
    ;
  }
   {   // FileUtil's options ------------------------------
     using fu = nestdaq::FileUtil;
     //using fopt = fu::OptionKey;
     //using sopt = fu::SplitOption;
     //using oopt = fu::OpenmodeOption;
         auto &desc = fu::GetDescriptions();	
         fu::AddOptions(options);
   }
}


//______________________________________________________________________________
std::unique_ptr<fair::mq::Device> getDevice(fair::mq::ProgOptions& /*config*/)
{
   return std::make_unique<Scaler>();
}

//______________________________________________________________________________
Scaler::Scaler()
    : FairMQDevice()
{
    
}

//______________________________________________________________________________
bool Scaler::HandleData(FairMQParts& parts, int index)
{
  namespace STF = SubTimeFrame;
  namespace Data = AmQStrTdc::Data;
  using Bits = Data::Bits;

  auto stfHeader = reinterpret_cast<STF::Header*>(parts.At(0)->GetData());
  int nCh = 0;
  if(stfHeader->FEMType == 1 || stfHeader->FEMType == 3){
    nCh = 128;
  }else if(stfHeader->FEMType == 2){
    nCh = 64;
  }
  if (nCh != hScaler->GetNBins()){
    hScaler->Reset();
    hScaler->SetNBins(nCh);
    hScaler->SetMinimum(0.);
    hScaler->SetMaximum((double)nCh);
    hScalerPrev->Reset();
    hScalerPrev->SetNBins(nCh);
    hScalerPrev->SetMinimum(0.);
    hScalerPrev->SetMaximum((double)nCh);
  }
  
  (void)index;
  assert(parts.Size()>=2);
  Reporter::AddInputMessageSize(parts);
  {
    auto dt = std::chrono::steady_clock::now() - fPrevUpdate;
    if (std::chrono::duration_cast<std::chrono::milliseconds>(dt).count() > fUpdateIntervalInMs) {
      try {
        fTsHeartbeatCounterKey = join({"scaler", fdevId, "heartbeatCounter"}, fSeparator);
	data_store->ts_add(fTsHeartbeatCounterKey, "*", std::to_string(tsHeartbeatCounter));
	fTsHeartbeatFlagKey = join({"scaler", fdevId, "heartbeatFlag"}, fSeparator);
	data_store->ts_add(fTsHeartbeatFlagKey, "*", std::to_string(tsHeartbeatFlag));
	
	std::string fCountsVsCh = join({"scaler", fdevId, "counts-vs-ch"}, fSeparator);
	data_store->write(fCountsVsCh,Slowdashify(*hScaler));
	
	uint64_t deltaHeartbeatCounter = tsHeartbeatCounter - fPrevHeartbeatCounter;
	double   dt_in_sec = deltaHeartbeatCounter * 0.000524; /* 1 heart beat = 512 us */
	if (dt_in_sec == 0.) {
	  dt_in_sec = 1.;
	}
	UH1Book hCountRate(*hScaler);
	hCountRate.Subtract(*hScalerPrev);
	for (int i=1; i <= hScaler->GetNBins(); i++) {
	  hCountRate.SetBinContent(i, hCountRate.GetBinContent(i)/dt_in_sec);
	}
	
	std::string fRateVsCh = join({"scaler", fdevId, "rate-vs-ch"}, fSeparator);
	data_store->write(fRateVsCh,Slowdashify(hCountRate));
	
	// for checking the flag bit for AmQ
	std::string fFlagBit = join({"flag", fdevId, "bit"}, fSeparator);
	data_store->write(fFlagBit,Slowdashify(*hFlag));
	
	LOG(info) << "tsHeartbeatCounter: " << tsHeartbeatCounter;
	LOG(info) << "tsHeartbeatFlag   : " << tsHeartbeatFlag;
	
      } catch (const std::exception& e) {
      	LOG(error) << "Scaler " << __FUNCTION__ << " exception : what() " << e.what();
      } catch (...) {
      	LOG(error) << "Scaler " << __FUNCTION__ << " exception : unknown ";
      }

      fPrevUpdate = std::chrono::steady_clock::now();
      fPrevHeartbeatCounter = tsHeartbeatCounter;
      *hScalerPrev = *hScaler;
      hFlag->Print();
      hFlag->Draw();
      hScaler->Print();
      hScaler->Draw();
    }
  }
  
  //  Check(stfs);
  
  auto nmsg      = parts.Size();
  for(int imsg = 1; imsg < nmsg; ++imsg){

    const auto& msg = parts.At(imsg);
    
    auto wb =  reinterpret_cast<Bits*>(msg->GetData());
    if(fDebug){
      LOG(debug) << " word =" << std::hex << wb->raw << std::dec;
      LOG(debug) << " head =" << std::hex << wb->head << std::dec;
    }
      
    switch (wb->head) {
    case Data::Heartbeat:
      {
	//      LOG(info) << "HBF comes:  " << std::hex << wb->raw;
	if(fDebug){
	  LOG(debug) << "== Data::Heartbeat --> ";
	  LOG(debug) << "hbframe: " << std::hex << wb->hbframe << std::dec;
	  LOG(debug) << "hbspill#: " << std::hex << wb->hbspilln << std::dec;
	  LOG(debug) << "hbfalg: " << std::hex << wb->hbflag << std::dec;
	  LOG(debug) << "header: " << std::hex << wb->htype << std::dec;
	  LOG(debug) << "head =" << std::hex << wb->head << std::dec;
	  LOG(debug) << "== scaler --> ";
	}
	if(fDebug){
	  LOG(debug) << "============================";
	  LOG(debug) << "HB frame : " << wb->hbframe;		
	  LOG(debug) << "timeFrameId : " << stfHeader->timeFrameId;
	  LOG(debug) << "# of HB: " << wb->hbframe - stfHeader->timeFrameId;
	  LOG(debug) << "hbflag: "  << wb->hbflag;
	}
	
	tsHeartbeatCounter++;
	
	for(int i=0; i<10; i++){
	  auto bit_sum = (wb->hbflag >> i) & 0x01;
	  FlagSum[i] = fpreFlagSum[i] + bit_sum;
	  if (bit_sum) {
	    hFlag->Fill(i);
	  }
	}      
	tsHeartbeatFlag = wb->hbflag;      

	for(int i=0; i<10; i++)
	  fpreFlagSum[i] = FlagSum[i];
            
	//      LOG(info) << "HBF Count:  " << tsHeartbeatCounter;
	//      LOG(info) << "HBF Flag:  " << tsHeartbeatFlag;
	
	break;
      }
    case Data::SpillOn: 
      LOG(info) << "SpillOn:  timeframeId->" << stfHeader->timeFrameId << " " << std::hex << stfHeader->FEMId ;
      break;

    case Data::SpillEnd: 
      LOG(info) << "SpillEnd:  timeframeId->" << stfHeader->timeFrameId << " " << std::hex << stfHeader->FEMId ;
      break;
	
    case Data::Data:
      {
	auto msgBegin = reinterpret_cast<Data::Word*>(msg->GetData());	
	auto msgSize  = msg->GetSize();
	auto nWord    = msgSize / sizeof(uint64_t);
	
	for(long unsigned int i = 0; i < nWord; ++i){	  
	  auto iwb = reinterpret_cast<Bits*>(msgBegin+i);	
	  
	  if(fDebug){
	    if( (stfHeader->FEMType == 1) || (stfHeader->FEMType==3) )
	      LOG(debug) << "LRtdc: "   << std::hex << iwb->tdc << std::dec;
	    if(stfHeader->FEMType == 2)	  
	      LOG(debug) << "HRtdc: " << std::hex << iwb->hrtdc<< std::dec;
	  }
	
	  if( (stfHeader->FEMType == 1) || (stfHeader->FEMType==3) ) {
	    hScaler->Fill(static_cast<int>(iwb->ch)+1);
	  }else if(stfHeader->FEMType==2) {
	    hScaler->Fill(static_cast<int>(iwb->ch)+1);
	  }
	    
	}//
      }
	
      break;	

    default:
      LOG(error) << " unknown Head : " << std::hex << wb->head << std::dec
		 << " FEMId: " << std::hex << stfHeader->FEMId;
                   
      break;
    }
  }


  /* make scaler data */

  auto outdata  = std::make_unique<std::vector<uint64_t>>();
  auto scHeader = std::make_unique<STF::Header>();
  
  for ( int ich = 0; ich < nCh; ++ich){
    outdata->push_back(hScaler->GetBinContent(ich+1));	      
  }
      
  /* send data to Filesink */
  FairMQParts outParts;
  fFEMId = stfHeader->FEMId;
  
  //  auto femid = stfHeader->FEMId & 0xff;
  scHeader->FEMType     = 200;  
  scHeader->FEMId       = stfHeader->FEMId;
  scHeader->length      = outdata->size()*sizeof(uint64_t) + sizeof(STF::Header);  
  scHeader->timeFrameId = stfHeader->timeFrameId;
  scHeader->numMessages = stfHeader->numMessages;
  scHeader->time_sec    = stfHeader->time_sec;
  scHeader->time_usec   = stfHeader->time_usec;     
    
    
  FairMQMessagePtr tmsg = MessageUtil::NewMessage(*this, std::move(scHeader));

#if 0
  auto n = tmsg->GetSize()/sizeof(uint64_t);      
  LOG(debug) << "tmsg: " << n << "  " << tmsg->GetSize();      
  std::for_each(reinterpret_cast<uint64_t*>(tmsg->GetData()),
		reinterpret_cast<uint64_t*>(tmsg->GetData()) + n,
		::HexDump{4});
#endif
      
  outParts.AddPart(std::move(tmsg));
  //  
  FairMQMessagePtr smsg;       
  if(outdata->size() > 0){
    smsg = MessageUtil::NewMessage(*this, std::move(outdata));

#if 0
    auto m = smsg->GetSize()/sizeof(uint64_t);	
    LOG(debug) << "smsg: " << m << "  " << smsg->GetSize();	
    std::for_each(reinterpret_cast<uint64_t*>(smsg->GetData()),
		  reinterpret_cast<uint64_t*>(smsg->GetData()) + m,
		  ::HexDump{4});
#endif
	
    outParts.AddPart(std::move(smsg));	
  }	

      
  while(Send(outParts, fOutputChannelName) < 0){
    if (NewStatePending()) {
      LOG(info) << "Device is not RUNNING";
      return false;
    }
    LOG(error) << "Failed to enqueue OutputChannel"; 
  }
	
  return true;
}

//______________________________________________________________________________
void Scaler::Init()
{
    Reporter::Instance(fConfig);
}

//______________________________________________________________________________
void Scaler::InitTask()
{
    using opt = OptionKey;
    
    fBufferTimeoutInMs  = std::stoi(fConfig->GetProperty<std::string>(opt::BufferTimeoutInMs.data()));
    LOG(debug) << "fBufferTimeoutInMs: "<< fBufferTimeoutInMs ;    

    fInputChannelName    = fConfig->GetProperty<std::string>(opt::InputChannelName.data());
    fOutputChannelName   = fConfig->GetProperty<std::string>(opt::OutputChannelName.data());    
    
    auto numSubChannels = GetNumSubChannels(fInputChannelName);
    fNumSource = 0;
    for (auto i=0u; i<numSubChannels; ++i) {
      fNumSource += GetNumberOfConnectedPeers(fInputChannelName, i);
    }
    LOG(debug) << "fNumSource: "<< fNumSource;

    LOG(debug) << " input channel : name = " << fInputChannelName
               << " num = " << GetNumSubChannels(fInputChannelName)
               << " num peer = " << GetNumberOfConnectedPeers(fInputChannelName,0);

    LOG(debug) << " output channel : name = " << fOutputChannelName;
    
    LOG(debug) << " number of source = " << fNumSource;

    fUpdateIntervalInMs = std::stoi(fConfig->GetProperty<std::string>(opt::UpdateInterval.data()));
    LOG(debug) << "fUpdateIntervalInMs: "<< fUpdateIntervalInMs ;      

    fHbc.resize(fNumSource);

    using opt = Scaler::OptionKey;
    std::string serverUri = fConfig->GetProperty<std::string>(opt::ServerUri.data());
    LOG(debug) << "serverUri: " << serverUri;
    data_store = new RedisDataStore(serverUri);
    
    fdevId        = fConfig->GetProperty<std::string>("id");
    LOG(debug) << "fdevId: " << fdevId;
    fTopPrefix   = "scaler";
    fSeparator   = fConfig->GetProperty<std::string>("separator");
    LOG(debug) << "fSeparator: " << fSeparator;
    const auto fCreatedTimeKey = join({fTopPrefix, opt::CreatedTimePrefix.data()},   fSeparator);
    const auto fHostNameKey    = join({fTopPrefix, opt::HostnamePrefix.data()},      fSeparator);
    const auto fIpAddressKey   = join({fTopPrefix, opt::HostIpAddressPrefix.data()}, fSeparator);
    data_store->hset(fCreatedTimeKey, fdevId, "n/a");
    data_store->hset(fHostNameKey,    fdevId, fConfig->GetProperty<std::string>("hostname"));
    data_store->hset(fIpAddressKey,   fdevId, fConfig->GetProperty<std::string>("host-ip"));

    Reporter::Reset();

    fFile = std::make_unique<nestdaq::FileUtil>();
    fFile->Init(fConfig->GetVarMap());
    fFile->Print();
    fFileExtension = fFile->GetExtension();

    //
    OnData(fInputChannelName, &Scaler::HandleData);
    
}

//______________________________________________________________________________
void Scaler::PreRun()
{
    using opt = Scaler::OptionKey;
    if (fConfig->Count(opt::RunNumber.data())) {
        fRunNumber = std::stoll(fConfig->GetProperty<std::string>(opt::RunNumber.data()));
        LOG(debug) << " RUN number = " << fRunNumber;
    }
    fFile->SetRunNumber(fRunNumber);
    fFile->ClearBranch();
    fFile->Open();

    for(int i=0; i<10; ++i){
      FlagSum[i] = 0;
      fpreFlagSum[i] = 0;
    }
    hScaler     = new UH1Book("ScalerHisto",128,0.,128.);
    hScalerPrev = new UH1Book("ScalerHistoPrev",128,0.,128.);
    hFlag       = new UH1Book("FlagHisto",10,0.,10.);
    hFlagPrev   = new UH1Book("FlagHistoPrev",10,0.,10.);
}
//______________________________________________________________________________
void Scaler::PostRun()
{

    int nrecv=0;
    if (fChannels.count(fInputChannelName) > 0) {
        auto n = fChannels.count(fInputChannelName);

        for (auto i = 0u; i < n; ++i) {
	    if (fDebug) {
	        LOG(debug) << " #i : "<< i;
	    }
            while(true) {

                FairMQParts part;
                if (Receive(part, fInputChannelName, i, 1000) <= 0) {
		  //                    LOG(debug) << __func__ << " no data received " << nrecv;
                    ++nrecv;
                    if (nrecv > 10) {
                        break;
                    }
                    std::this_thread::sleep_for(std::chrono::milliseconds(200));
                } else {
		  //                    LOG(debug) << __func__ << " data comes..";
                }
            }
        }
    }// for clean up

    std::stringstream ss;
    ss << "Scaler data at PostRun()" << std::endl;
    ss << "== scaler --> module:  "  << std::hex << fFEMId << std::endl;
    for (int i = 0; i < hScaler->GetNBins(); i++){
      ss << "ch: " << i << "     sc: " << hScaler->GetBinContent(i+1) << std::endl;
    }
    ss << " --> scaler ==" << std::endl;
    
    fFile->Write(reinterpret_cast<const char *>(ss.str().c_str()), ss.str().length());
    ss.str("");
    ss.clear(std::stringstream::goodbit);
    fFile->Close();
    fFile->ClearBranch();

    delete hScaler;
    delete hScalerPrev;
    delete hFlag;
    delete hFlagPrev;  
    LOG(debug) << __func__ << " done";
}
