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

#include <chmap/channel_map_dopeness.hpp>
#include <chmap/item.hpp>

namespace nestdaq {
   class FilterTimeFrameSliceByTrack;
   struct DCRawHit;
   struct temporary_geometry;
   struct temporary_dctdccalib;
   struct temporary_dcdriftparam;
   class DCHit;
}

struct nestdaq::DCRawHit {
   DCRawHit(chmap::DETIdItem* detid, uint32_t tdc) : detid(detid), tdc(tdc) {};
   const chmap::DETIdItem* detid;
   uint32_t tdc;
}; // struct nestdaq::DCRawHit

class nestdaq::DCHit {
public:
   DCHit();
   DCHit(double wirePos, double wireAngle, const chmap::DETIdItem* detid) {
      this->wirePos = wirePos;
      this->wireAngle = wireAngle;
      this->detid = detid;
   };
   ~DCHit() = default;

   void AddHit(double tdc){
      TDCs.push_back(tdc);
      return;
   }

   const chmap::DETIdItem* GetDETIdItem() const { return detid; };
   int Clear(){
      int n = TDCs.size();
      TDCs.clear();
      DriftTimes.clear();
      DriftLengths.clear();
      return n;
   };

   void CalcDriftTimes();
   void CalcDriftLengths();

   double GetWirePos() const { return wirePos; };
   double GetWireAngle() const { return wireAngle; };
   double GetDriftLength(int nth) const { return DriftLengths[nth]; };

private:
   double wirePos;
   double wireAngle;
   
   std::vector<double> TDCs;
   std::vector<double> DriftTimes;
   std::vector<double> DriftLengths;
   const chmap::DETIdItem* detid;
}; // class nestdaq::DCHit

struct nestdaq::temporary_geometry {
   int detectoridentifier;
   std::string detectorname;
   double x, y, z;
   double tiltangle, rotationangle1, rotationangle2;
   double length, resolution, wirecenternumber, wirepitch, offset;
}; // struct nestdaq::temporary_geometry

struct nestdaq::temporary_dctdccalib {
   int detectoridentifier;
   int wireidentifier;
   double offset, scale;
}; // struct nestdaq::temporary_dctdccalib

struct nestdaq::temporary_dcdriftparam {
   int detectoridentifier;
   int approxOrder;
   std::vector<double> coefficients;
}; // struct nestdaq::temporary_dcdriftparam


class nestdaq::FilterTimeFrameSliceByTrack : public nestdaq::FilterTimeFrameSliceABC {
public:
   struct OptionKey : FilterTimeFrameSliceABC::OptionKey {
      static constexpr std::string_view ChannelMapDataFile {"chmap-data-file"};
      static constexpr std::string_view GeometryConfigFile {"geometry-file"};
      static constexpr std::string_view DCTdcCalibConfigFile {"dctdc-calib-file"};
      static constexpr std::string_view DCDriftParamConfigFile {"dc-drift-param-file"};
   };

   FilterTimeFrameSliceByTrack();
   virtual ~FilterTimeFrameSliceByTrack() override = default;

   void InitTask() override;
   virtual bool ProcessSlice(TTF& ) override;

protected:
   // ================================
   // channel map
   // ================================
   bool isCreateInvMap{true};
   chmap::ChannelMapDopeness* fChMap{nullptr};
   std::string fChMapDataFile;

   // ================================
   // detector configuration
   // ================================
   int LoadDetectorConfig_Geometry(std::string_view filename);
   int LoadDetectorConfig_DCTdcCalib(std::string_view filename);
   int LoadDetectorConfig_DCDriftParam(std::string_view filename);

   std::vector<temporary_geometry> fTemporaryGeometries;
   std::vector<temporary_dctdccalib> fTemporaryDCTdcCalibs;
   std::vector<temporary_dcdriftparam> fTemporaryDCDriftParams;

   std::map<int, std::string> detectorNameMap;
   std::map<int, std::string> detectorPlaneMap;
   std::map<int, uint8_t> detectorSegmentMap;
   void DefineDetectorIdMap(); // for reading detector configuration files from AnalyzerT103

   bool RegisterDetectorConfig_Geometry();
   bool RegisterDetectorConfig_DCTdcCalib();
   bool RegisterDetectorConfig_DCDriftParam(); 

   // ================================
   // TDC value calculation
   // ================================
   const double tdc64h_to_ns = 0.9765625*0.001; // HighResolution TDC LSB unit: 0.9762265 ps, convert to ns

   // ================================
   // File I/O for debugging
   // ================================
   std::ofstream fDebugFile;
   std::string fDebugFileName;


}; // class nestdaq::FilterTimeFrameSliceByTrack




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
