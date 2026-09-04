/**
 * @file FilterTimeFrameSliceByTrack.cxx
 * @brief Slice Timeframe by Logic timing for NestDAQ
 * @date Created : 2026-07-14 13:41:20 JST
 *       Last Modified : 2026-07-14 13:41:20 JST
 *
 * @author Shinsuke OTA <ota@rcnp.osaka-u.ac.jp>
 * @author Shunichi KASHIMA <kashima@rcnp.osaka-u.ac.jp>
 *
 */
#include "FilterTimeFrameSliceByTrack.h"
#include "FilterTimeFrameSliceABC.icxx"
#include "fairmq/runDevice.h"

#include "utility/MessageUtil.h"
#include "UnpackTdc.h"

#include "SubTimeFrameHeader.h"
#include "TimeFrameHeader.h"
#include "FrameContainer.h" // TTF, TSTF, THBF, TLF(parse raw data)

// for the Parser of detector configurations
#include <fstream>
#include <sstream>

// for registering the detector configurations
#include <map>

// for cout
#include <iomanip>
#include <iostream>

// for debugging
#include "FilterTimeFrameSliceByTrackDebugger.h"

#define DEBUG 0

using nestdaq::FilterTimeFrameSliceByTrack;
using nestdaq::DCRawHit;
using nestdaq::DCHit;
using nestdaq::KLDCHitContainer;
namespace bpo = boost::program_options;




FilterTimeFrameSliceByTrack::FilterTimeFrameSliceByTrack()
{
}

void FilterTimeFrameSliceByTrack::InitTask()
{
   const std::string_view funcname = "[FilterTimeFrameSliceByTrack::InitTask] ";
   FilterTimeFrameSliceABC::InitTask();
   using opt = OptionKey;

   // ================================
   // channel map initialization
   // ================================
   fChMapDataFile = fConfig->GetProperty<std::string>(opt::ChannelMapDataFile.data());
   chmap::ChannelMapDopeness& chmap = chmap::ChannelMapDopeness::get_instance();
   chmap.initialize(fChMapDataFile, isCreateInvMap);

   fChMap = &chmap;
   std::cout << "[FilterTimeFrameSliceABC::InitTask] ChannelMapDopeness initialized with " << fChMapDataFile << std::endl;
   std::cout << "\t# of channels: " << std::dec << fChMap->getNumberOfChannels() << std::endl;

   #if CHECK_COUT_CHMAP
   uint32_t dopeKey_FEtoDET;
   uint32_t dopeKey_DETtoFE;
   chmap::DETIdItem detiditem;
   chmap::FEAddrItem feaddritem;

   // test t1 right channel
   uint8_t test_ip3rd_T1right = 0x02;
   uint8_t test_ip4th_T1right = 0xAA;
   uint8_t test_ch_T1right = 12;
   std::cout << "\n\tchecking T1 right DETIdItem search..." << std::endl;
   bool _FOUND_FEtoDET = fChMap->getDopeKey_FEtoDET(test_ip3rd_T1right, test_ip4th_T1right, test_ch_T1right, dopeKey_FEtoDET);
   if(_FOUND_FEtoDET == true){
      std::cout << "\t\t-> found." << std::endl;
      detiditem = fChMap->getDETIdItem(dopeKey_FEtoDET);
      detiditem.decode();
   }
   else{
      std::cout << "\t\t-> not found." << std::endl;
   }

   // test kldc 2 U 95
   uint8_t test_ip3rd_kldc2u95 = 0x02;
   uint8_t test_ip4th_kldc2u95 = 0xB2;
   uint8_t test_ch_kldc2u95 = 96;

   // check existance of kldc 2 U 95 DETIdItem
   std::cout << "\n\t" << "checking kldc 2 U 95 DETIdItem search" << std::endl;
   _FOUND_FEtoDET = fChMap->getDopeKey_FEtoDET(test_ip3rd_kldc2u95, test_ip4th_kldc2u95, test_ch_kldc2u95, dopeKey_FEtoDET);
   if(_FOUND_FEtoDET == true){
      std::cout << "\t" << funcname << "-> found." << std::endl;
      detiditem = fChMap->getDETIdItem(dopeKey_FEtoDET);
      detiditem.decode();
   }
   else{
      std::cout << "\t" << funcname << "-> not found." << std::endl;
   }

   // check the counterpart of kldc 2 U 95 DETidItem
   std::cout << "\t" << "checking kldc 2 U 95 FEAddrItem search..." << std::endl;
   bool _FOUND_DETtoFE = fChMap->getDopeKey_DETtoFE(detiditem, dopeKey_DETtoFE);
   if(_FOUND_DETtoFE == true){
      std::cout << "\t" << funcname << "-> found." << std::endl;
      feaddritem = fChMap->getFEAddrItem(dopeKey_DETtoFE);
      feaddritem.decode();
   }
   else{
      std::cout << "\t" << funcname << "-> not found when searching the counterpart of\n";
      detiditem.decode();
   }

   // check existance of utof right DETIdItem
   const std::string_view test_detectorname_utofright = "utof";
   const std::string_view test_planename_utofright = "0";
   const uint8_t test_segment_utofright = 0;
   const uint16_t test_channelnumber_utofright = 0;
   const std::string_view test_channelname_utofright = "right";
   std::cout << "\n\t" << "checking utof right FEAddrItem search..." << std::endl;
   _FOUND_DETtoFE = fChMap->getDopeKey_DETtoFE(test_detectorname_utofright, test_planename_utofright, test_segment_utofright, test_channelname_utofright, test_channelnumber_utofright, dopeKey_DETtoFE);
   if(_FOUND_DETtoFE == true){
      std::cout << "\t" << "utof right FEAddrItem -> found." << std::endl;
      feaddritem = fChMap->getFEAddrItem(dopeKey_DETtoFE);
      feaddritem.decode();
      _FOUND_FEtoDET = fChMap->getDopeKey_FEtoDET(feaddritem, dopeKey_FEtoDET);
      if(_FOUND_FEtoDET == true){
         std::cout << "\t" << "utof right DETIdItem -> found." << std::endl;
         detiditem = fChMap->getDETIdItem(dopeKey_FEtoDET);
         detiditem.decode();
      }
   }

   // check existance of utof left DETIdItem
   const std::string_view test_detectorname_utofleft = "utof";
   const std::string_view test_planename_utofleft = "0";
   const uint8_t test_segment_utofleft = 0;
   const uint16_t test_channelnumber_utofleft = 0;
   const std::string_view test_channelname_utofleft = "left";
   std::cout << "\n\t" << "checking utof left FEAddrItem search..." << std::endl;
   _FOUND_DETtoFE = fChMap->getDopeKey_DETtoFE(test_detectorname_utofleft, test_planename_utofleft, test_segment_utofleft, test_channelname_utofleft, test_channelnumber_utofleft, dopeKey_DETtoFE);
   if(_FOUND_DETtoFE == true){
      std::cout << "\t" << "utof left FEAddrItem -> found." << std::endl;
      feaddritem = fChMap->getFEAddrItem(dopeKey_DETtoFE);
      feaddritem.decode();
      _FOUND_FEtoDET = fChMap->getDopeKey_FEtoDET(feaddritem, dopeKey_FEtoDET);
      if(_FOUND_FEtoDET == true){
         std::cout << "\t" << "utof left DETIdItem -> found." << std::endl;
         detiditem = fChMap->getDETIdItem(dopeKey_FEtoDET);
         detiditem.decode();
      }
   }
   #endif

   // ================================
   // detector configuration initialization
   // ================================

   const auto geometryFile = fConfig->GetProperty<std::string>(opt::GeometryConfigFile.data());
   std::cout << "[FilterTimeFrameSliceByTrack::InitTask] Geometry configuration file: " << geometryFile << std::endl;
   const auto dctdcCalibFile = fConfig->GetProperty<std::string>(opt::DCTdcCalibConfigFile.data());
   std::cout << "[FilterTimeFrameSliceByTrack::InitTask] DC TDC calibration configuration file: " << dctdcCalibFile << std::endl;
   const auto dcDriftParamFile = fConfig->GetProperty<std::string>(opt::DCDriftParamConfigFile.data());
   std::cout << "[FilterTimeFrameSliceByTrack::InitTask] DC drift parameter configuration file: " << dcDriftParamFile << std::endl;

   if (!geometryFile.empty()) {
      if (!LoadDetectorConfig_Geometry(geometryFile)) {
         std::cerr << "[FilterTimeFrameSliceByTrack::InitTask] Failed to parse geometry configuration file: " << geometryFile << std::endl;
      }
   }
   if (!dctdcCalibFile.empty()) {
      if (!LoadDetectorConfig_DCTdcCalib(dctdcCalibFile)) {
         std::cerr << "[FilterTimeFrameSliceByTrack::InitTask] Failed to parse DC TDC calibration configuration file: " << dctdcCalibFile << std::endl;
      }
   }
   if (!dcDriftParamFile.empty()) {
      if (!LoadDetectorConfig_DCDriftParam(dcDriftParamFile)) {
         std::cerr << "[FilterTimeFrameSliceByTrack::InitTask] Failed to parse DC drift parameter configuration file: " << dcDriftParamFile << std::endl;
      }
   }

   DefineDetectorIdMap(); // for reading detector configuration files from AnalyzerT103
   RegisterDetectorConfig_Geometry();
   RegisterDetectorConfig_DCTdcCalib();
   RegisterDetectorConfig_DCDriftParam();
   #if CHECK_COUT_DETCONF
   // kldc 1 U 63
   std::cout << "\n\t" << "checking kldc 1 U 63 DETIdItem search..." << std::endl;
   const std::string_view test_detectorname_kldc1u63 = "kldc";
   const std::string_view test_planename_kldc1u63 = "U";
   const uint8_t test_segment_kldc1u63 = 1;
   const uint16_t test_channelnumber_kldc1u63 = 63;
   const std::string_view test_channelname_kldc1u63 = "0";
   _FOUND_DETtoFE = fChMap->getDopeKey_DETtoFE(test_detectorname_kldc1u63, test_planename_kldc1u63, test_segment_kldc1u63, test_channelname_kldc1u63, test_channelnumber_kldc1u63, dopeKey_DETtoFE);
   if(_FOUND_DETtoFE == true){
      std::cout << "\t" << "-> found." << std::endl;
      feaddritem = fChMap->getFEAddrItem(dopeKey_DETtoFE);
      feaddritem.decode();
      _FOUND_FEtoDET = fChMap->getDopeKey_FEtoDET(feaddritem, dopeKey_FEtoDET);
      if(_FOUND_FEtoDET == true){
         std::cout << "\t" << "-> found." << std::endl;
         detiditem = fChMap->getDETIdItem(dopeKey_FEtoDET);
         detiditem.decode();
         // check detector configuration for kldc 1 U 63
         const chmap::GeomItemDC* retrieved_geomitemdc = dynamic_cast<const chmap::GeomItemDC*>(detiditem.detconf->membername_geom.get());
         if(retrieved_geomitemdc != nullptr){
            std::cout << "\t" << "-> found detector configuration for kldc 1 U 63." << std::endl;
         }
         else{
            std::cout << "\t" << "-> not found detector configuration for kldc 1 U 63." << std::endl;
         }
      } // if(_FOUND_FEtoDET == true)
   } // if(_FOUND_DETtoFE == true)
   else{
      std::cout << "\t" << "-> not found." << std::endl;
   }
   // plane name index check
   const std::string test_planename_up = "Up";
   uint8_t test_planename_up_index = 0;
   bool _FOUND_Index = fChMap->plane_dictionary.StringToIndex(test_planename_up, test_planename_up_index);
   if(_FOUND_Index == true){
      std::cout << "\t" << "-> found plane name index for Up: " << std::hex << static_cast<int>(test_planename_up_index) << std::dec << std::endl;
   }
   else{
      std::cout << "\t" << "-> not found plane name index for Up." << std::endl;
   }

   const std::string test_planename_v = "V";
   uint8_t test_planename_v_index = 0;
   _FOUND_Index = fChMap->plane_dictionary.StringToIndex(test_planename_v, test_planename_v_index);
   if(_FOUND_Index == true){
      std::cout << "\t" << "-> found plane name index for V: " << std::hex << static_cast<int>(test_planename_v_index) << std::dec << std::endl;
   }
   else{
      std::cout << "\t" << "-> not found plane name index for V." << std::endl;
   }

   const std::string test_planename_u = "U";
   uint8_t test_planename_u_index = 0;
   _FOUND_Index = fChMap->plane_dictionary.StringToIndex(test_planename_u, test_planename_u_index);
   if(_FOUND_Index == true){
      std::cout << "\t" << "-> found plane name index for U: " << std::hex << static_cast<int>(test_planename_u_index) << std::dec << std::endl;
   }
   else{
      std::cout << "\t" << "-> not found plane name index for U." << std::endl;
   }

   const std::string test_planename_vp = "Vp";
   uint8_t test_planename_vp_index = 0;
   _FOUND_Index = fChMap->plane_dictionary.StringToIndex(test_planename_vp, test_planename_vp_index);
   if(_FOUND_Index == true){
      std::cout << "\t" << "-> found plane name index for Vp: " << std::hex << static_cast<int>(test_planename_vp_index) << std::dec << std::endl;
   }
   else{
      std::cout << "\t" << "-> not found plane name index for Vp." << std::endl;
   }
   #endif

   // ================================
   // File I/O for debugging
   // ================================
   #if FILEOUT_LFTDC
   fDebugFileName = "./fileout/tracking/FilterTimeFrameSliceByTrack_debug.txt";
   fDebugFile.open(fDebugFileName, std::ios::out);
   if (!fDebugFile.is_open()) {
      std::cerr << "[FilterTimeFrameSliceByTrack::InitTask] Failed to open debug file: " << fDebugFileName << std::endl;
   }
   #endif
} // void FilterTimeFrameSliceByTrack::InitTask()

bool FilterTimeFrameSliceByTrack::ProcessSlice(TTF& tf)
{
   const std::string_view funcname = "[FilterTimeFrameSliceByTrack::ProcessSlice] ";
   #if DEBUG_LFTDC
   std::cout << funcname << "Function called" << std::endl;
   std::cout << "\tchecking TLF TDC 4ns unit: " << std::dec << std::setw(10) << fLFTDC4n << " -> " << std::setw(10) << fLFTDC4n * 4 << " [ns]" << std::endl;
   #endif

   auto tfHeader = tf.GetHeader();
   auto numSTF = tfHeader->numSource;

   uint32_t femId = 0;
   uint16_t ch = 0; // 8でも十分なんだけど、255よりでかい謎のエントリーを除外するために16で受ける
   double ftdc = 0; // unit: ns
   uint32_t ftdc_int = 0; // unit: ns
   uint32_t ftot_int = 0; // unit: ns
   const uint32_t lftdc = fLFTDC4n * 4; // FilterTimeFrameSliceABCの持ってるfield lftdc4n
   std::vector<double> utof_left_times;

   chmap::DETIdItem* detiditem = nullptr;

   // ================================
   // Scan, searching for UTOF
   // ================================
   const uint32_t femId_utof = 0xc0a802a9;
   const uint16_t ch_utof_right = 10;
   const uint16_t ch_utof_left = 8;
   const uint32_t tdc_min = lftdc - 100; // unit: ns
   const uint32_t tdc_max = lftdc + 100; // unit: ns
   int nTDC_utof_right = 0;
   int nTDC_utof_left = 0;
   uint32_t time_utof_right = 0;
   uint32_t time_utof_left = 0;

   #if FILEOUT_LFTDC
   double min_diff_utof_right = 1e6;
   double min_diff_utof_left = 1e6;
   #endif

   for(auto& stf : tf){
      auto stfHeader = stf->GetHeader();
      femId = stfHeader->femId;
      if(femId != femId_utof){
         continue;
      }
      // auto numHBF = stfHeader->numMessages; // tabun, 1

      TDC64H_V3::tdc64 tdc64_h;
      // TDC64L_V3::tdc64 tdc64_l;
      uint32_t nTDC = 0;
      if(stfHeader->femType == SubTimeFrame::TDC64H_V3){
         auto& hbf = stf->at(0);
         nTDC = hbf->GetNumData();
         if(nTDC == 0){
            continue; // no TDC data
         }
         #if DEBUG_LFTDC
         std::cout << funcname << "HR TDC FEE found. femId = " << std::hex << femId << std::dec << std::endl;
         std::cout << "\tnTDC: " << nTDC << std::endl;
         #endif
         for(uint32_t iTDC=0; iTDC<nTDC; ++iTDC){
            TDC64H_V3::Unpack(hbf->UncheckedAt(iTDC), &tdc64_h);
            ch = tdc64_h.ch;
            // if(ch == ch_utof_right){
            //    ftdc = tdc64_h.tdc * 0.9765625 * 0.001; // HR TDCのLSBが0.9765625 ps = 1/2^10 nsなので、(0.9765625 * 0.001)を掛ける
            //    // ftdc_int = tdc64_h.tdc>>10; // HR TDCのLSBが0.9765625 ps = 1/2^10 nsなので、(0.9765625 * 0.001)を掛ける代わりに2^10を掛ける
            //    #if DEBUG_LFTDC
            //    std::cout << funcname << "(double) utof right: ch = " << ch << ", tdc     - lftdc = " << ftdc - lftdc << std::endl;
            //    std::cout << funcname << "(int)    utof right: ch = " << ch << ", tdc_int - lftdc = " << static_cast<int>(ftdc_int) - static_cast<int>(lftdc) << std::endl;
            //    #endif

            //    #if FILEOUT_LFTDC
            //    if(fabs(ftdc - lftdc) < min_diff_utof_right){
            //       min_diff_utof_right = ftdc - lftdc;
            //    }
            //    #endif

            //    if(ftdc >= tdc_min && ftdc <= tdc_max){
            //       nTDC_utof_right++;
            //       time_utof_right = ftdc;
            //    }
            // } // if(ch == ch_utof_right)
            if(ch == ch_utof_left){
               ftdc = tdc64_h.tdc * tdc64h_to_ns; // HR TDCのLSBが0.9765625 ps = 1/2^10 nsなので、(0.9765625 * 0.001)を掛ける
               // ftdc_int = tdc64_h.tdc>>10; // HR TDCのLSBが0.9765625 ps = 1/2^10 nsなので、(0.9765625 * 0.001)を掛ける代わりに2^10を掛ける
               #if DEBUG_LFTDC
               std::cout << funcname << "(double) utof left: ch = " << ch << ", tdc     - lftdc = " << ftdc - lftdc << std::endl;
               std::cout << funcname << "(int)    utof left: ch = " << ch << ", tdc_int - lftdc = " << static_cast<int>(ftdc_int) - static_cast<int>(lftdc) << std::endl;
               #endif
               
               #if FILEOUT_LFTDC
               if(fabs(ftdc - lftdc) < min_diff_utof_left){
                  min_diff_utof_left = ftdc - lftdc;
               }
               #endif

               utof_left_times.push_back(ftdc);
               // if(ftdc >= tdc_min && ftdc <= tdc_max){
               //    nTDC_utof_left++;
               //    time_utof_left = ftdc;
               // }
            } // if(ch == ch_utof_left)
         } // for(uint32_t iTDC=0; iTDC<nTDC; ++iTDC)
      } // if(stfHeader->femType == SubTimeFrame::TDC64H_V3)
   } // for(auto& stf : tf)

   #if FILEOUT_LFTDC
   if(min_diff_utof_right < 1e6){
      fDebugFile << 10 << " " << min_diff_utof_right << std::endl;
   }
   if(min_diff_utof_left < 1e6){
      fDebugFile << 8 << " " << min_diff_utof_left << std::endl;
   }
   #endif


   #if DEBUG_LFTDC
   if(nTDC_utof_right > 0 || nTDC_utof_left > 0){
      std::cout << funcname << "UTOF hit in range [lftdc-2, lftdc+2] = [" << lftdc << " - 2, " << lftdc <<  " + 2]" << std::endl;
      std::cout << "\tutof right: nTDC = " << nTDC_utof_right << ", time = " << time_utof_right << std::endl;
      std::cout << "\tutof left: nTDC = " << nTDC_utof_left << ", time = " << time_utof_left << std::endl;
   }
   #endif

   // ================================
   // Scan, searching for KLDC with hit time after the coincidence TDC from LogicFilter block
   // ================================
   // define the standard time for the start time of drift time calculation
   int nStandardTime = utof_left_times.size();
   int standardTime = 0;
   if(nStandardTime > 0){
      standardTime = static_cast<int>(*std::min_element(utof_left_times.begin(), utof_left_times.end()));
   }
   else{
      return false;
   }

   // scan
   uint8_t femId_ip3rd = 0;
   uint8_t femId_ip4th = 0;
   uint32_t keyFEtoDET;
   const uint8_t kldc_detname_index = 0x06;
   std::vector<DCRawHit> kldcRawHits;
   for(auto& stf : tf){
      auto stfHeader = stf->GetHeader();
      if(stfHeader->femType == SubTimeFrame::TDC64L_V3){
         // std::cout << "TDC64L_V3 FEM found. femId = " << std::hex << femId << std::dec << std::endl;
         auto& hbf = stf->at(0);
         // auto hbfHeader = hbf->GetHeader();
         uint32_t nData = hbf->GetNumData();
         femId = stfHeader->femId;
         femId_ip3rd = (femId >> 8) & 0xff;
         femId_ip4th = femId & 0xff;

         TDC64L_V3::tdc64 tdc;
         for(uint32_t i=0; i<nData; ++i){
            TDC64L_V3::Unpack(hbf->UncheckedAt(i), &tdc);
            ch = tdc.ch;
            ftdc_int = tdc.tdc;
            ftot_int = tdc.tot;
            bool isFound = fChMap->getDopeKey_FEtoDET(femId_ip3rd, femId_ip4th, ch, keyFEtoDET);
            if(isFound){
               detiditem = &fChMap->getDETIdItem(keyFEtoDET);
               if(detiditem->name == kldc_detname_index){
                  if(ftdc_int >= standardTime){
                     DCRawHit hit(detiditem, ftdc_int, ftot_int);
                     kldcRawHits.push_back(hit);
                  } // if(ftdc_int >= standardTime)
               } // if(detiditem->name == kldc_detname_index)
            } // if(isFound)
         } // for(uint32_t i=0; i<nData; ++i)
      } // else if(stfHeader->femType == SubTimeFrame::TDC64L_V3)
      else{
         continue; // skip other FEM types
      } // if(stfHeader->femType == SubTimeFrame::TDC64L_V3)
   } // for(auto& stf : tf)

   #if 0
   std::cout << funcname << "Number of KLDC hits after standard time: " << kldcRawHits.size() << std::endl;
   #endif

   // sort the raw hits, first by segment, then by plane, and finally by channel number
   std::sort(kldcRawHits.begin(), kldcRawHits.end(), [](const DCRawHit& left, const DCRawHit& right) {
      if(left.detid->segment != right.detid->segment) {
         return left.detid->segment < right.detid->segment;
      } else if(left.detid->plane != right.detid->plane) {
         return left.detid->plane < right.detid->plane;
      } else {
         return left.detid->channel_number < right.detid->channel_number;
      }
   });

   #if 0
   // print the sorted KLDC hits for debugging
   std::cout << funcname << "\tKLDC hit:" << std::endl;
   for(const auto& hit : kldcRawHits){
      std::cout << "\t\tsegment: " << std::dec << std::setfill(' ') << std::setw(2) << static_cast<int>(hit.detid->segment)
                << ", plane: " << std::setw(2) << static_cast<int>(hit.detid->plane)
                << ", channel: " << std::setw(3) << hit.detid->channel_number
                << ", tdc: " << std::setw(10) << hit.tdc
                << std::endl;
   } // for(const auto& hit : kldcRawHits)
   #endif

   // distribute the KLDC hits to its container
   const uint8_t u_plane_index = 0x01; // (std::string)"U"
   const uint8_t v_plane_index = 0x03; // (std::string)"V"
   const uint8_t up_plane_index = 0x02; // (std::string)"Up"
   const uint8_t vp_plane_index = 0x04; // (std::string)"Vp"

   // 
   for(auto& hit : kldcRawHits){
      const chmap::DETIdItem* detiditem = hit.detid;
      uint8_t segment = detiditem->segment;
      uint16_t wireNumber = detiditem->channel_number;
      uint32_t tdc = hit.tdc;
      uint32_t tot = hit.tot;
      double wirePos = 0.0;
      double wireAngle = 0.0;

      if(detiditem->segment == 1){ // KLDC1
         if(detiditem->plane == u_plane_index || detiditem->plane == up_plane_index){
            if(!(fKLDCHitContainer[0].empty()) && fKLDCHitContainer[0].back().GetDETIdItem()->channel_number == wireNumber){
               fKLDCHitContainer[0].back().AddHit(tdc, tot);
            }
            else{
               if(detiditem->detconf != nullptr){
                  const chmap::GeomItemDC* geomitemdc = dynamic_cast<const chmap::GeomItemDC*>(detiditem->detconf->membername_geom.get());
                  if(geomitemdc != nullptr){
                     wirePos = geomitemdc->GetWirePosition();
                     wireAngle = geomitemdc->GetTiltAngle();
                  }
                  DCHit h(wirePos, wireAngle, detiditem, tdc, tot);
                  fKLDCHitContainer[0].push_back(h);
               } // if(detiditem->detconf != nullptr)
            }
         } // if(detiditem->plane == u_plane_index || detiditem->plane == up_plane_index)
         else if(detiditem->plane == v_plane_index || detiditem->plane == vp_plane_index){
            if(!(fKLDCHitContainer[1].empty()) && fKLDCHitContainer[1].back().GetDETIdItem()->channel_number == wireNumber){
               fKLDCHitContainer[1].back().AddHit(tdc, tot);
            }
            else{
               if(detiditem->detconf != nullptr){
                  const chmap::GeomItemDC* geomitemdc = dynamic_cast<const chmap::GeomItemDC*>(detiditem->detconf->membername_geom.get());
                  if(geomitemdc != nullptr){
                     wirePos = geomitemdc->GetWirePosition();
                     wireAngle = geomitemdc->GetTiltAngle();
                  }
                  DCHit h(wirePos, wireAngle, detiditem, tdc, tot);
                  fKLDCHitContainer[1].push_back(h);
               } // if(detiditem->detconf != nullptr)
            }
         } // if(detiditem->plane == v_plane_index || detiditem->plane == vp_plane_index)
      } // if(detiditem->segment == 1)
      else if(detiditem->segment == 2){ // KLDC2
         if(detiditem->plane == u_plane_index || detiditem->plane == up_plane_index){
            if(!(fKLDCHitContainer[2].empty()) && fKLDCHitContainer[2].back().GetDETIdItem()->channel_number == wireNumber){
               fKLDCHitContainer[2].back().AddHit(tdc, tot);
            }
            else{
               if(detiditem->detconf != nullptr){
                  const chmap::GeomItemDC* geomitemdc = dynamic_cast<const chmap::GeomItemDC*>(detiditem->detconf->membername_geom.get());
                  if(geomitemdc != nullptr){
                     wirePos = geomitemdc->GetWirePosition();
                     wireAngle = geomitemdc->GetTiltAngle();
                  }
                  DCHit h(wirePos, wireAngle, detiditem, tdc, tot);
                  fKLDCHitContainer[2].push_back(h);
               } // if(detiditem->detconf != nullptr)
            }
         } // if(detiditem->plane == u_plane_index || detiditem->plane == up_plane_index)
         else if(detiditem->plane == v_plane_index || detiditem->plane == vp_plane_index){
            if(!(fKLDCHitContainer[3].empty()) && fKLDCHitContainer[3].back().GetDETIdItem()->channel_number == wireNumber){
               fKLDCHitContainer[3].back().AddHit(tdc, tot);
            }
            else{
               if(detiditem->detconf != nullptr){
                  const chmap::GeomItemDC* geomitemdc = dynamic_cast<const chmap::GeomItemDC*>(detiditem->detconf->membername_geom.get());
                  if(geomitemdc != nullptr){
                     wirePos = geomitemdc->GetWirePosition();
                     wireAngle = geomitemdc->GetTiltAngle();
                  }
                  DCHit h(wirePos, wireAngle, detiditem, tdc, tot);
                  fKLDCHitContainer[3].push_back(h);
               } // if(detiditem->detconf != nullptr)
            }
         } // if(detiditem->plane == v_plane_index || detiditem->plane == vp_plane_index)
      } // if(detiditem->segment == 2)
   } // for(auto& hit : kldcRawHits)

   // ================================
   // tracking for each UTOF hit
   // ================================
   std::cout << funcname << "Number of UTOF left hits: " << nStandardTime << std::endl;
   for(int i=0; i<nStandardTime; ++i){
      int standardTime = static_cast<int>(utof_left_times[i]);
      fKLDCHitContainer.SetStandardTime(standardTime, fDCTimeRange);
   }

#if 0
   int doKeep = false;

   if (tf[0]->GetHeader()->femid == /* some module */) {
      auto& hbf = tf[0]->at(0)); // 
      uint64_t nData = hbf->GetNumData();
      for (int i = 0; i < nData; ++i) {
         // unpack data
         // Unpack(hbf->UncheckedAt(i),tdc);
      //
      }
      // judge
      doKeep = true;
   }

   if (!doKeep) {
      return false;
   }
#endif   
   return false;
}

int FilterTimeFrameSliceByTrack::LoadDetectorConfig_Geometry(std::string_view filename)
{
   /*
   input file format
   - ignore lines starting with '#'
   - each line contains the following fields:
     detectorID detectorname x y z tiltangle rotationangle1 rotationangle2 length resolution wirecenternumber wirepitch offset
      - fields are separated by whitespace
   - wirecenternumber is the number of wires where the center of Drift Chamber is located
      - 4.5 means the center is between wire 4 and wire 5
   */
   const std::string_view funcName = "[FilterTimeFrameSliceByTrack::LoadDetectorConfig_Geometry] ";
   std::ifstream ifs(filename.data());
   if(!ifs.is_open()){
      std::cerr << funcName << "Failed to open file: " << filename << std::endl;
      return -1;
   }

   std::string line;
   while(std::getline(ifs, line)){
      if(line.empty() || line.front() == '#'){
         continue;
      }

      std::istringstream iss(line);

      int detectoridentifier{0};
      std::string detectorname;
      double x{0.0}, y{0.0}, z{0.0};
      double tiltangle{0.0}, rotationangle1{0.0}, rotationangle2{0.0};
      double length{0.0}, resolution{0.0}, wirecenternumber{0.0}, wirepitch{0.0}, offset{0.0};

      if(!(iss >> detectoridentifier >> detectorname >> x >> y >> z >> tiltangle >> rotationangle1 >> rotationangle2 >> length >> resolution >> wirecenternumber >> wirepitch >> offset)){
         std::cerr << funcName << "Failed to parse line: " << line << std::endl;
         continue;
      }
      fTemporaryGeometries.push_back({detectoridentifier, detectorname, x, y, z, tiltangle, rotationangle1, rotationangle2, length, resolution, wirecenternumber, wirepitch, offset});
   } // while(std::getline(ifs, line))

   std::cout << funcName << "Parsed " << fTemporaryGeometries.size() << " geometry entries from file: " << filename << std::endl;

   return fTemporaryGeometries.size();
} // int FilterTimeFrameSliceByTrack::LoadDetectorConfig_Geometry(std::string_view filename)

int FilterTimeFrameSliceByTrack::LoadDetectorConfig_DCTdcCalib(std::string_view filename)
{
   const std::string_view funcName = "[FilterTimeFrameSliceByTrack::LoadDetectorConfig_DCTdcCalib] ";
   std::ifstream ifs(filename.data());
   if(!ifs.is_open()){
      std::cerr << funcName << "Failed to open file: " << filename << std::endl;
      return -1;
   }


   std::string line;
   while(std::getline(ifs, line)){
      if(line.empty() || line.front() == '#'){
         continue;
      }

      std::istringstream iss(line);

      int detectoridentifier{0};
      int wireidentifier{0};
      double offset{0.0}, scale{0.0};

      if(!(iss >> detectoridentifier >> wireidentifier >> offset >> scale)){
         std::cerr << funcName << "Failed to parse line: " << line << std::endl;
         continue;
      }
      fTemporaryDCTdcCalibs.push_back({detectoridentifier, wireidentifier, offset, scale});
   } // while(std::getline(ifs, line))

   std::cout << funcName << "Parsed " << fTemporaryDCTdcCalibs.size() << " TDC calibration entries from file: " << filename << std::endl;

   return fTemporaryDCTdcCalibs.size();
} // int FilterTimeFrameSliceByTrack::LoadDetectorConfig_DCTdcCalib(std::string_view filename)

int FilterTimeFrameSliceByTrack::LoadDetectorConfig_DCDriftParam(std::string_view filename)
{
   const std::string_view funcName = "[FilterTimeFrameSliceByTrack::LoadDetectorConfig_DCDriftParam] ";
   std::ifstream ifs(filename.data());
   if(!ifs.is_open()){
      std::cerr << funcName << "Failed to open file: " << filename << std::endl;
      return -1;
   }

   std::string line;
   while(std::getline(ifs, line)){
      if(line.empty() || line.front() == '#'){
         continue;
      }

      std::istringstream iss(line);

      int detectoridentifier{0};
      int wireidentifier{0};
      int approxOrder{0};
      std::vector<double> coefficients;

      int buf;

      if(!(iss >> detectoridentifier >> buf >> approxOrder >> buf)){
         std::cerr << funcName << "Failed to parse line: " << line << std::endl;
         continue;
      }
      coefficients.resize(approxOrder);
      for(int i = 0; i < approxOrder; ++i){
         if(!(iss >> coefficients[i])){
            std::cerr << funcName << "Failed to parse coefficient " << i << " in line: " << line << std::endl;
            continue;
         }
      }
      fTemporaryDCDriftParams.push_back({detectoridentifier, approxOrder, coefficients});
      #if 0
      std::cout << "{detectoridentifier, approxOrder, coefficients} = {" << detectoridentifier << ", " << approxOrder << ", {";
      for(int i = 0; i < approxOrder; ++i){
         std::cout << coefficients[i];
         if(i < approxOrder - 1){
            std::cout << ", ";
         }
      }
      std::cout << "}}" << std::endl;
      #endif
   } // while(std::getline(ifs, line))

   std::cout << funcName << "Parsed " << fTemporaryDCDriftParams.size() << " drift parameter entries from file: " << filename << std::endl;

   return fTemporaryDCDriftParams.size();
} // int FilterTimeFrameSliceByTrack::LoadDetectorConfig_DCDriftParam(std::string_view filename)

void FilterTimeFrameSliceByTrack::DefineDetectorIdMap()
{
   {
      detectorNameMap[101] = "bdc";
      detectorNameMap[102] = "bdc";
      detectorNameMap[103] = "bdc";
      detectorNameMap[104] = "bdc";
      detectorNameMap[105] = "bdc";
      detectorNameMap[106] = "bdc";
      detectorNameMap[107] = "bdc";
      detectorNameMap[108] = "bdc";

      detectorNameMap[201] = "kldc";
      detectorNameMap[202] = "kldc";
      detectorNameMap[203] = "kldc";
      detectorNameMap[204] = "kldc";
      detectorNameMap[205] = "kldc";
      detectorNameMap[206] = "kldc";
      detectorNameMap[207] = "kldc";
      detectorNameMap[208] = "kldc";

      detectorNameMap[301] = "bft";
      detectorNameMap[302] = "bft";
      detectorNameMap[303] = "bft";
      detectorNameMap[304] = "bft";
      detectorNameMap[305] = "bft";
      detectorNameMap[306] = "bft";
   } // detectorNameMap

   {
      detectorPlaneMap[101] = "X";
      detectorPlaneMap[102] = "Xp";
      detectorPlaneMap[103] = "U";
      detectorPlaneMap[104] = "V";
      detectorPlaneMap[105] = "X";
      detectorPlaneMap[106] = "Xp";
      detectorPlaneMap[107] = "U";
      detectorPlaneMap[108] = "V";

      detectorPlaneMap[201] = "V";
      detectorPlaneMap[202] = "Vp";
      detectorPlaneMap[203] = "U";
      detectorPlaneMap[204] = "Up";
      detectorPlaneMap[205] = "Up";
      detectorPlaneMap[206] = "U";
      detectorPlaneMap[207] = "Vp";
      detectorPlaneMap[208] = "V";

      detectorPlaneMap[301] = "X";
      detectorPlaneMap[302] = "U";
      detectorPlaneMap[303] = "V";
      detectorPlaneMap[304] = "X";
      detectorPlaneMap[305] = "U";
      detectorPlaneMap[306] = "V";
   } // detectorPlaneMap

   {
      detectorSegmentMap[101] = 1;
      detectorSegmentMap[102] = 1;
      detectorSegmentMap[103] = 1;
      detectorSegmentMap[104] = 1;
      detectorSegmentMap[105] = 2;
      detectorSegmentMap[106] = 2;
      detectorSegmentMap[107] = 2;
      detectorSegmentMap[108] = 2;

      detectorSegmentMap[201] = 1;
      detectorSegmentMap[202] = 1;
      detectorSegmentMap[203] = 1;
      detectorSegmentMap[204] = 1;
      detectorSegmentMap[205] = 2;
      detectorSegmentMap[206] = 2;
      detectorSegmentMap[207] = 2;
      detectorSegmentMap[208] = 2;

      detectorSegmentMap[301] = 1;
      detectorSegmentMap[302] = 1;
      detectorSegmentMap[303] = 1;
      detectorSegmentMap[304] = 2;
      detectorSegmentMap[305] = 2;
      detectorSegmentMap[306] = 2;
   } // detectorSegmentMap
} // void FilterTimeFrameSliceByTrack::DefineDetectorIdMap()

bool FilterTimeFrameSliceByTrack::RegisterDetectorConfig_Geometry()
{
   const std::string_view funcname = "[FilterTimeFrameSliceByTrack::RegisterDetectorConfig_Geometry] ";
   std::cout << "\n" << funcname << "Registering geometry configurations...(# of items: " << fTemporaryGeometries.size() << ")" << std::endl;

   int registered_count_geomitemdc_kldc = 0;
   // Loop over loaded temporary geometries
   for(const auto& geom : fTemporaryGeometries){
      #if CHECK_COUT_DETCONF_REGISTERING
      std::cout << funcname << "Registering geometry for detector ID: " << geom.detectoridentifier << std::endl;
      #endif

      // geom -> local variables
      int detectorId = geom.detectoridentifier;
      double x = geom.x;
      double y = geom.y;
      double z = geom.z;
      double tiltAngle = geom.tiltangle;
      double rotationAngle1 = geom.rotationangle1;
      double rotationAngle2 = geom.rotationangle2;
      double length = geom.length;
      double resolution = geom.resolution;
      double wireCenterNumber = geom.wirecenternumber;
      double wirePitch = geom.wirepitch;
      double offset = geom.offset;

      // prepare channel map information
      std::string DetectorName = detectorNameMap[detectorId];
      std::string PlaneName = detectorPlaneMap[detectorId];
      uint8_t SegmentNumber = detectorSegmentMap[detectorId];
      std::string ChannelName = "0";

      int missing_count_FEAddrItem = 0;
      int missing_count_DETIdItem = 0;
      // Register KLDC
      if(DetectorName == "kldc"){
         for(int i=0+32; i<128-32; ++i){
            int ChannelNumber = i;
            std::unique_ptr<chmap::GeomItemDC> geomitemdc = std::make_unique<chmap::GeomItemDC>();
            geomitemdc->SetGlobalPosition(x, y, z);
            geomitemdc->SetResolution(resolution, resolution, resolution);
            geomitemdc->SetRotationAngles(tiltAngle, rotationAngle1, rotationAngle2);
            geomitemdc->SetWireGeometry(wireCenterNumber, wirePitch, offset);
            geomitemdc->CalcWirePosition(ChannelNumber);

            uint32_t dopeKey_DETtoFE;
            bool found_DETtoFE = fChMap->getDopeKey_DETtoFE(DetectorName, PlaneName, SegmentNumber, std::string("0"), static_cast<uint16_t>(ChannelNumber), dopeKey_DETtoFE);
            if(found_DETtoFE){
               chmap::FEAddrItem feaddritem = fChMap->getFEAddrItem(dopeKey_DETtoFE);
               uint32_t dopeKeyFEtoDET;
               bool found_FEtoDET = fChMap->getDopeKey_FEtoDET(feaddritem.ip3rd, feaddritem.ip4th, feaddritem.ch, dopeKeyFEtoDET);
               if(found_FEtoDET){
                  fChMap->registerDETConfSubItem<chmap::GeomItem, chmap::GeomItemDC>(dopeKeyFEtoDET, std::move(geomitemdc), &chmap::DETConfItem::membername_geom);

                  chmap::DETIdItem detiditem = fChMap->getDETIdItem(dopeKeyFEtoDET);
                  const chmap::GeomItemDC* retrieved_geomitemdc = dynamic_cast<const chmap::GeomItemDC*>(detiditem.detconf->membername_geom.get());
                  #if 0
                  std::cout << "pointer address of geomitemdc after move: " << geomitemdc.get() << std::endl;
                  std::cout << "pointer address of detiditem.detconf: " << detiditem.detconf << std::endl;
                  std::cout << "pointer address of detiditem.detconf->membername_geom.get(): " << detiditem.detconf->membername_geom.get() << std::endl;
                  std::cout << "pointer address of retrieved_geomitemdc: " << retrieved_geomitemdc << std::endl;
                  #endif
                  if(retrieved_geomitemdc != nullptr){
                     registered_count_geomitemdc_kldc++;
                  }
                  else{
                     std::cerr << funcname << "Failed to register GeomItemDC for DetectorName: " << DetectorName << ", PlaneName: " << PlaneName << ", SegmentNumber: " << static_cast<int>(SegmentNumber) << ", ChannelNumber: " << ChannelNumber << std::endl;
                  }
               } // if(found_FEtoDET)
               else{
                  ++missing_count_DETIdItem;
               }
            } // if(found_DETtoFE)
            else{
               ++missing_count_FEAddrItem;
            }
         } // for(int i=0+32; i<128-32; ++i)
         if(missing_count_DETIdItem > 0){
            std::cout << funcname << "Missing FEAddrItem count: " << missing_count_FEAddrItem << std::endl;
         }
         if(missing_count_DETIdItem > 0){
            std::cout << funcname << "Missing DETIdItem count: " << missing_count_DETIdItem << std::endl;
         }
      } // if(DetectorName == "kldc")
   } // for(const auto& geom : fTemporaryGeometries)

   std::cout << funcname << "Registered " << registered_count_geomitemdc_kldc << " geometry entries of KLDC" << std::endl;
   return true;
} // bool FilterTimeFrameSliceByTrack::RegisterDetectorConfig_Geometry()

bool FilterTimeFrameSliceByTrack::RegisterDetectorConfig_DCTdcCalib()
{
   const std::string_view funcname = "[FilterTimeFrameSliceByTrack::RegisterDetectorConfig_DCTdcCalib] ";
   std::cout << "\n" << funcname << "Registering TDC calibration configurations...(# of items: " << fTemporaryDCTdcCalibs.size() << ")" << std::endl;

   int missing_count_FEAddrItem = 0;
   int missing_count_DETIdItem = 0;
   int registered_count_dctdccalib_kldc = 0;
   bool isFirstEntry = true;

   // Loop over loaded temporary DC TDC calibration entries
   for(const auto& calib : fTemporaryDCTdcCalibs){
      // calib -> local variables
      int detectorId = calib.detectoridentifier;
      int wireId = calib.wireidentifier; // 1-indexed
      double offset = calib.offset;
      double scale = calib.scale;

      // prepare channel map information
      std::string DetectorName = detectorNameMap[detectorId];
      std::string PlaneName = detectorPlaneMap[detectorId];
      uint8_t SegmentNumber = detectorSegmentMap[detectorId];
      std::string ChannelName = std::string("0");
      if(isFirstEntry){
         std::cout << funcname << "Registering first tdccalib for KLDC detector ID: " << detectorId << std::endl;
         std::cout << "\t" << DetectorName << ", " << PlaneName << ", Segment: " << static_cast<int>(SegmentNumber) << ",wire: " << wireId << std::endl;
         isFirstEntry = false;
      }

      // Create CalibrationItem_DCTdcCalib and set its properties
      std::unique_ptr<chmap::CalibrationItem_DCTdcCalib> calibitem_dctdccalib = std::make_unique<chmap::CalibrationItem_DCTdcCalib>();
/* CLASS DEFINITION
    class CalibrationItem_DCTdcCalib : public CalibrationItem {
        public:
            virtual ~CalibrationItem_DCTdcCalib() = default;
            void SetTdcCalibration(double offset_, double scale_) {
                offset = offset_;
                scale = scale_;
            }
        private:
            double offset, scale; // relative time = (absolute time * scale) + offset
    };
*/
      calibitem_dctdccalib->SetTdcCalibration(offset, scale);


      // Register KLDC
      if(DetectorName == "kldc"){
         int ChannelNumber = wireId - 1; // wireId is 1-indexed, ChannelNumber is 0-indexed
         uint32_t dopeKey_DETtoFE;
         bool found_DETtoFE = fChMap->getDopeKey_DETtoFE(DetectorName, PlaneName, SegmentNumber, ChannelName, static_cast<uint16_t>(ChannelNumber), dopeKey_DETtoFE);
         if(found_DETtoFE){
            chmap::FEAddrItem feaddritem = fChMap->getFEAddrItem(dopeKey_DETtoFE);
            uint32_t dopeKeyFEtoDET;
            bool found_FEtoDET = fChMap->getDopeKey_FEtoDET(feaddritem.ip3rd, feaddritem.ip4th, feaddritem.ch, dopeKeyFEtoDET);
            if(found_FEtoDET){
               fChMap->registerDETConfSubItem<chmap::CalibrationItem, chmap::CalibrationItem_DCTdcCalib>(dopeKeyFEtoDET, std::move(calibitem_dctdccalib), &chmap::DETConfItem::membername_calib_dctdccalib);
               registered_count_dctdccalib_kldc++;
            } // if(found_FEtoDET)
            else{
               ++missing_count_DETIdItem;
            }
         } // if(found_DETtoFE)
         else{
            ++missing_count_FEAddrItem;
         }
      } // if(DetectorName == "kldc")
   } // for(const auto& calib : fTemporaryDCTdcCalibs)

   if(missing_count_FEAddrItem > 0){
      std::cout << funcname << "Missing FEAddrItem count: " << missing_count_FEAddrItem << std::endl;
   }
   if(missing_count_DETIdItem > 0){
      std::cout << funcname << "Missing DETIdItem count: " << missing_count_DETIdItem << std::endl;
   }

   std::cout << funcname << "Registered " << registered_count_dctdccalib_kldc << " TDC calibration entries of KLDC." << std::endl;
   return true;
} // bool FilterTimeFrameSliceByTrack::RegisterDetectorConfig_DCTdcCalib()

bool FilterTimeFrameSliceByTrack::RegisterDetectorConfig_DCDriftParam()
{
   const std::string_view funcname = "[FilterTimeFrameSliceByTrack::RegisterDetectorConfig_DCDriftParam] ";
   std::cout << "\n" << funcname << "Registering drift parameter configurations...(# of items: " << fTemporaryDCDriftParams.size() << ")" << std::endl;

   int registered_count_geomitemdc_kldc = 0;
   // Loop over loaded temporary DC drift length parameter entries
   for(const auto& drift : fTemporaryDCDriftParams){
      #if 0
      std::cout << funcname << "detectorID: " << drift.detectoridentifier << ", approxOrder: " << drift.approxOrder << ", coefficients: {";
      for(int i = 0; i < drift.approxOrder; ++i){
         std::cout << drift.coefficients[i];
         if(i < drift.approxOrder - 1){
            std::cout << ", ";
         }
      }
      std::cout << "}" << std::endl;
      #endif
      // drift -> local variables
      int detectorId = drift.detectoridentifier;
      int approxOrder = drift.approxOrder;
      const std::vector<double>& coefficients = drift.coefficients;

      std::string DetectorName = detectorNameMap[detectorId];
      std::string PlaneName = detectorPlaneMap[detectorId];
      uint8_t SegmentNumber = detectorSegmentMap[detectorId];

      // prepare channel map information
      std::string ChannelName = std::string("0");

      // Create CalibrationItem_DCDriftLength and set its properties
      std::unique_ptr<chmap::CalibrationItem_DCDriftLength> calibitem_dcdriftlength = std::make_unique<chmap::CalibrationItem_DCDriftLength>();
/* CLASS DEFINITION
    class CalibrationItem_DCDriftLength : public CalibrationItem {
        public:
            virtual ~CalibrationItem_DCDriftLength() = default;
            void SetApproximation(int approxOrder_, const std::vector<double>& coeffs_) {
                approxOrder = approxOrder_;
                coeffs = coeffs_;
            }
        private:
            int approxOrder; // the number of coefficients for polynomial approximation of drift length
            std::vector<double> coeffs; // coefficients for polynomial approximation of drift length
            // dLen(t) = coeffs[0] + coeffs[1]*t + coeffs[2]*t^2 + ... + coeffs[n]*t^n
    };
*/
      calibitem_dcdriftlength->SetApproximation(approxOrder, coefficients);

      int missing_count_FEAddrItem = 0;
      int missing_count_DETIdItem = 0;
      // Register KLDC
      if(DetectorName == "kldc"){
         for(int i=0+32; i<128-32; ++i){
            int ChannelNumber = i;
            uint32_t dopeKey_DETtoFE;
            bool found_DETtoFE = fChMap->getDopeKey_DETtoFE(DetectorName, PlaneName, static_cast<uint8_t>(SegmentNumber), ChannelName, static_cast<uint8_t>(ChannelNumber), dopeKey_DETtoFE);
            if(found_DETtoFE){
               chmap::FEAddrItem feaddritem = fChMap->getFEAddrItem(dopeKey_DETtoFE);
               uint32_t dopeKeyFEtoDET;
               bool found_FEtoDET = fChMap->getDopeKey_FEtoDET(feaddritem.ip3rd, feaddritem.ip4th, feaddritem.ch, dopeKeyFEtoDET);
               if(found_FEtoDET){
                  fChMap->registerDETConfSubItem<chmap::CalibrationItem, chmap::CalibrationItem_DCDriftLength>(dopeKeyFEtoDET, std::move(calibitem_dcdriftlength), &chmap::DETConfItem::membername_calib_dcdriftlen);
                  registered_count_geomitemdc_kldc++;
               } // if(found_FEtoDET)
               else{
                  ++missing_count_DETIdItem;
               }
            } // if(found_DETtoFE)
            else{
               ++missing_count_FEAddrItem;
            }
         } // for(int i=0+32; i<128-32; ++i)
         if(missing_count_FEAddrItem > 0){
            std::cout << funcname << "Missing FEAddrItem count: " << missing_count_FEAddrItem << std::endl;
         }
         if(missing_count_DETIdItem > 0){
            std::cout << funcname << "Missing DETIdItem count: " << missing_count_DETIdItem << std::endl;
         }
      } // if(DetectorName == "kldc")
   } // for(const auto& drift : fTemporaryDCDriftParams)

   std::cout << funcname << "Registered " << registered_count_geomitemdc_kldc << " drift parameter entries." << std::endl;
   return true;
} // bool FilterTimeFrameSliceByTrack::RegisterDetectorConfig_DCDriftParam()

////////////////////////////////////////////////////
// override runDevice
////////////////////////////////////////////////////

void addCustomOptions(bpo::options_description& options)
{
   using opt = FilterTimeFrameSliceByTrack::OptionKey;

   options.add_options()
      (opt::InputChannelName.data(),
       bpo::value<std::string>()->default_value("in"),
       "Name of the input channel")
      (opt::OutputChannelName.data(),
       bpo::value<std::string>()->default_value("out"),
       "Name of the output channel")
      (opt::DQMChannelName.data(),
       bpo::value<std::string>()->default_value("dqm"),
       "Name of the data quality monitoring channel")
      (opt::PollTimeout.data(),
       bpo::value<std::string>()->default_value("1"),
       "Timeout of polling (in msec)")
      (opt::SplitMethod.data(),
       bpo::value<std::string>()->default_value("1"),
       "STF split method")
      (opt::ChannelMapDataFile.data(),
       bpo::value<std::string>()->default_value("channel_map.dat"),
       "Path to the channel map data file")
      (opt::GeometryConfigFile.data(),
       bpo::value<std::string>()->default_value(""),
       "Geometry configuration file")
      (opt::DCTdcCalibConfigFile.data(),
       bpo::value<std::string>()->default_value(""),
       "DC TDC calibration file")
      (opt::DCDriftParamConfigFile.data(),
       bpo::value<std::string>()->default_value(""),
       "DC drift parameter file")
      ;
   
}



std::unique_ptr<fair::mq::Device> getDevice(fair::mq::ProgOptions& /*config*/)
{
    return std::make_unique<FilterTimeFrameSliceByTrack>();
}

void KLDCHitContainer::SetStandardTime(double standardTime, const DCTimeRange& DCTimeRange){
   for(auto& std_vector_dchit : *this){
      for(auto& dchit : std_vector_dchit){
         dchit.CalcDriftTimes(standardTime, DCTimeRange);
      } // for(auto& dchit : std_vector_dchit)
   }
} // void nestdaq::FilterTimeFrameSliceByTrack::KLDCHitContainer::SetStandardTime(int standardTime)

bool DCHit::CalcDriftTimes(double standardTime, const DCTimeRange& DCTimeRange){
   int dctdc_min = DCTimeRange.lower_bound;
   int dctdc_max = DCTimeRange.upper_bound;
   int dctot_min = DCTimeRange.tot_min;

   if(detid == nullptr){
      return false;
   }
   else{
      if(detid->detconf == nullptr){
         return false;
      }
      else{
         const chmap::CalibrationItem_DCTdcCalib* calibitem_dctdccalib = dynamic_cast<const chmap::CalibrationItem_DCTdcCalib*>(detid->detconf->membername_calib_dctdccalib.get());
         if(calibitem_dctdccalib == nullptr){
            return false;
         }
         double offset = calibitem_dctdccalib->GetOffset();
         double scale = calibitem_dctdccalib->GetScale();

         for(size_t i=0; i<TDCs.size(); ++i){ // already ensured that TDCs.size() == TOTs.size()
            uint32_t tdc = TDCs[i];
            uint32_t tot = TOTs[i];
            if((tot > dctot_min) && (tdc >= dctdc_min) && (tdc <= dctdc_max)){
               double driftTime = scale * (tdc - standardTime) + offset;
               DriftTimes.push_back(driftTime);             
            }
         } // for(const auto& tdc : tdcs)
      } // if(detid->detconf == nullptr)
   } // if(detid == nullptr)

   return this->CalcDriftLengths();
} // bool nestdaq::FilterTimeFrameSliceByTrack::DCHit::CalcDriftTimes(double standardTime)

bool DCHit::CalcDriftLengths(){
   if(detid == nullptr){
      return false;
   }
   else{
      if(detid->detconf == nullptr){
         return false;
      }
      else{
         const chmap::CalibrationItem_DCDriftLength* calibitem_dcdriftlen = dynamic_cast<const chmap::CalibrationItem_DCDriftLength*>(detid->detconf->membername_calib_dcdriftlen.get());
         if(calibitem_dcdriftlen == nullptr){
            return false;
         }
         for(const auto& driftTime : DriftTimes){
            double driftLength = calibitem_dcdriftlen->GetDriftLength(driftTime);
            DriftLengths.push_back(driftLength);
         } // for(const auto& driftTime : driftTimes)
      } // if(detid->detconf == nullptr)
   } // if(detid == nullptr)
   return true;
} // bool nestdaq::FilterTimeFrameSliceByTrack::DCHit::CalcDriftLengths()

