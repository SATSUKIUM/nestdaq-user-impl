/**
 * @file FilterTimeFrameSliceABC.h
 * @brief Abstract base class for filtering time frame slice
 * @date Created : 2024-05-04 12:27:57 JST
 *       Last Modified : 2024-05-23 07:19:25 JST
 *
 * @author Shinsuke OTA <ota@rcnp.osaka-u.ac.jp>
 *
 */


#ifndef NESTDAQ_FILTERTIMEFRAMESLICEABC_H
#define NESTDAQ_FILTERTIMEFRAMESLICEABC_H

#include "fairmq/Device.h"
#include "KTimer.cxx"
#include "SubTimeFrameHeader.h"
#include "TimeFrameHeader.h"
#include "FilterHeader.h"
#include "HeartbeatFrameHeader.h"
#include "FrameContainer.h"

#include <TFile.h>
#include <TTree.h>

#include <chrono>

#include <chmap/item.hpp>

namespace nestdaq {
   class FilterTimeFrameSliceABC;
}


class nestdaq::FilterTimeFrameSliceABC : public fair::mq::Device {
public:
   FilterTimeFrameSliceABC();
   virtual ~FilterTimeFrameSliceABC() override = default;

   virtual void PreRun() override;
   virtual void InitTask() override;
   virtual bool ConditionalRun() override;
   virtual void PostRun() override;
   virtual void ResetTask() override;


   struct OptionKey {
      static constexpr std::string_view InputChannelName {"in-chan-name"};
      static constexpr std::string_view OutputChannelName {"out-chan-name"};
      static constexpr std::string_view DQMChannelName {"out-chan-name"};
      static constexpr std::string_view PollTimeout        {"poll-timeout"};
      static constexpr std::string_view SplitMethod        {"split"};
   };

protected:

   virtual bool ParseMessages(FairMQParts& inParts);

   virtual bool ProcessSlice(TTF& ) { return true; }
   
   std::string fInputChannelName;
   std::string fOutputChannelName;
   std::string fName;
   uint32_t fID;


   // control
   uint32_t fNextIdx;

   std::vector<KTimer> fKTimer;
   bool fDoCheck;

   std::vector<TTF> fTFs; // time frame

   // ================================
   // for catching the Coincidence TDC from LogicFilter block
   // ================================
   TLF fLF;
   std::vector<TLF> fLFs;
   uint32_t fLFTDC4n; // 4ns-unit TDC of Coincidence

   // output
   int fNumDestination {0};
   uint32_t fDirection {0};
   int fPollTimeoutMS  {0};
   int fSplitMethod    {0};
   
   // ================================
   // Throughput measurement
   // ================================
   std::chrono::high_resolution_clock::time_point fStartTime, fEndTime;
   std::chrono::high_resolution_clock::time_point fBeforeReceive, fAfterReceive;
   double fElapsedTime{0.0};
   double fElapsedTimeReceive{0.0};
   TFile* fRootFile{nullptr};
   TTree* fRootTree1{nullptr};
   const std::string fRootFileName{"./fileout/tracking/ConditionalRun_throughput.root"};

   double fTree_elapsed_time{0.0};
   int fTree_nslice{0};
   int fTree_npassslice{0};


};


#endif  // NESTDAQ_FILTERTIMEFRAMESLICEABC_H
