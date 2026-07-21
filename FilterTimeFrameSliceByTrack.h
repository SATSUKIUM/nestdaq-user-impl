/**
 * @file FilterTimeFrameSliceByTrack.h
 * @brief Abstract base class for filtering time frame slice
 * @date Created : 2026-07-14 13:41:20 JST
 *       Last Modified : 2026-07-14 13:41:20 JST
 *
 * @author Shinsuke OTA <ota@rcnp.osaka-u.ac.jp>
 * @author Shunichi KASHIMA <kashima@rcnp.osaka-u.ac.jp>
 *
 */


#ifndef NESTDAQ_TIMEFRAMESLICERBYTRACK_H
#define NESTDAQ_TIMEFRAMESLICERBYTRACK_H

#include "fairmq/Device.h"
#include "KTimer.cxx"
#include "SubTimeFrameHeader.h"
#include "TimeFrameHeader.h"
#include "FilterHeader.h"
#include "HeartbeatFrameHeader.h"
#include "FrameContainer.h"
#include "FilterTimeFrameSliceABC.h"

#include <chmap/item.hpp>

namespace nestdaq {
   class FilterTimeFrameSliceByTrack;
   struct DCRawHit;
}


class nestdaq::FilterTimeFrameSliceByTrack : public nestdaq::FilterTimeFrameSliceABC {
public:
   FilterTimeFrameSliceByTrack();
   virtual ~FilterTimeFrameSliceByTrack() override = default;

   virtual bool ProcessSlice(TTF& ) override;

};

struct nestdaq::DCRawHit {
   const chmap::DETIdItem* detid;
   uint32_t tdc;
}; // struct nestdaq::DCRawHit


// channel-map provides (uint8_t DetectorNameIdx, uint8_t PlaneNameIdx, uint8_t SegmentNumber, uint16_t ChannelNumber, uint8_t ChannelNameIdx) or (std::string DetectorName, std::string PlaneName, uint8_t SegmentNumber, uint16_t ChannelNumber, std::string ChannelName)
// example: ("KLDC", " U", 1, 50, " 0")

// DCRawHit: TDC, DriftLength
// DCHitContainer: PlaneName, ChannelName, std::vector<std::vector<DCRawHit>>(<- 1st index is WireID, 2nd index is the hit number)
// DCHitContainerContainer: SegmentNumber, DetectorName, SegmentNumber, std::vector<DCHitContainer>(<- index is PlaneName)
// DCHitContainerContainerContainer: std::vector<DCHitContainerContainer>(<- index is SegmentNumber)

// class nestdaq::DCRawHit {
//    public:
//       DCRawHit() : TDC(0), DriftLength(0.0) {}
//       DCRawHit(uint32_t tdc) : TDC(tdc), DriftLength(0.0) {}
//       DCRawHit(uint32_t tdc, double driftLength) : TDC(tdc), DriftLength(driftLength) {}
//       ~DCRawHit() = default;

//       void SetDriftLength(double driftLength) {
//          this->DriftLength = driftLength;
//       }
//       double GetDriftLength() const {
//          return this->DriftLength;
//       }

//    private:
//       uint32_t TDC;
//       double DriftLength;
// }; // class nestdaq::DCRawHit

// class nestdaq::DCHitContainer{
//    public:   
//       DCHitContainer() : PlaneName(""), ChannelName(""), ChannelNumber(0) {}
//       DCHitContainer(const std::string& planeName, const std::string& channelName, uint16_t channelNumber) : PlaneName(planeName), ChannelName(channelName), ChannelNumber(channelNumber) {}
         
//       void Clear() {
//          this->Hits.clear();
//       }
//       void AlocHC(size_t n) {
//          this->Hits.resize(n);
//       }

//       void AddHitTo(uint16_t WireID, uint32_t tdc, double driftLength) {
//          Hits[WireID].emplace_back(DCRawHit(tdc, driftLength));
//       }
//       void AddHitTo(uint16_t WireID, uint32_t tdc) {
//          Hits[WireID].emplace_back(DCRawHit(tdc));
//       }
   
//    private:
//       uint8_t SegmentNumber;
//       std::vector<std::vector<DCRawHit>> Hits;
//       std::string PlaneName;
//       std::string ChannelName;
//       uint16_t ChannelNumber;

// };

// class nestdaq::DCHitContainerContainer {
//    public:
//       void Clear() {
//          std::for_each(this->HitContainers.begin(), this->HitContainers.end(),
//                        [] (DCHitContainer& c) {
//                           c.Clear();
//                        });
//          this->HitContainers.clear();
//       }

//       void AlocHC(size_t n) {
//          for (size_t i = 0; i < n; ++i) {
//             this->HitContainers.emplace_back(DCHitContainer());
//          }
//       }
//       void AddHitContainer(const std::string& planeName, const std::string& channelName, uint16_t channelNumber) {
//          DCHitContainer container(planeName, channelName, channelNumber);
//          this->HitContainers.emplace_back(std::move(container));
//       }
//       void SetDetectorName(const std::string& detectorName) {
//          this->DetectorName = detectorName;
//       }

//    private:
//       std::vector<DCHitContainer> HitContainers;
//       std::string DetectorName;
// };




#endif  // NESTDAQ_TIMEFRAMESLICERBYTRACK_H
