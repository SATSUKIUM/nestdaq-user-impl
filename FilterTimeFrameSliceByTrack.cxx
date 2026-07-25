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

// for the Parser of detector configurations
#include <fstream>
#include <sstream>

// for registering the detector configurations
#include <map>

// for debugging
#include "FilterTimeFrameSliceByTrackDebugger.h"

#define DEBUG 0

using nestdaq::FilterTimeFrameSliceByTrack;
namespace bpo = boost::program_options;




FilterTimeFrameSliceByTrack::FilterTimeFrameSliceByTrack()
{
}

void FilterTimeFrameSliceByTrack::InitTask()
{
   FilterTimeFrameSliceABC::InitTask();
   using opt = OptionKey;

   // ================================
   // channel map initialization
   // ================================
   fChMapDataFile = fConfig->GetProperty<std::string>(opt::ChannelMapDataFile.data());
   chmap::ChannelMapDopeness& chmap = chmap::ChannelMapDopeness::get_instance();
   chmap.initialize(fChMapDataFile, isCreateInvMap);

   fChMap = &chmap;

   #if CHECK_COUT_CHMAP
   std::cout << "[FilterTimeFrameSliceABC::InitTask] ChannelMapDopeness initialized with " << fChMapDataFile << std::endl;
   std::cout << "\t# of channels: " << std::dec << fChMap->getNumberOfChannels() << std::endl;

   // test t1 right channel
   uint8_t test_ip3rd_T1right = 0x02;
   uint8_t test_ip4th_T1right = 0xAA;
   uint8_t test_ch_T1right = 12;

   uint32_t dope_key;
   chmap::DETIdItem detitem;
   std::cout << "\tchecking T1 right search" << std::endl;
   bool found = fChMap->getDopeKey_FEtoDET(test_ip3rd_T1right, test_ip4th_T1right, test_ch_T1right, dope_key);
   if(found == true){
      std::cout << "\t\tfound." << std::endl;
      detitem = fChMap->getDETIdItem(dope_key);
      detitem.decode();
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
      if (!ParseDetectorConfig_Geometry(geometryFile)) {
         std::cerr << "[FilterTimeFrameSliceByTrack::InitTask] Failed to parse geometry configuration file: " << geometryFile << std::endl;
      }
   }
   if (!dctdcCalibFile.empty()) {
      if (!ParseDetectorConfig_DCTdcCalib(dctdcCalibFile)) {
         std::cerr << "[FilterTimeFrameSliceByTrack::InitTask] Failed to parse DC TDC calibration configuration file: " << dctdcCalibFile << std::endl;
      }
   }
   if (!dcDriftParamFile.empty()) {
      if (!ParseDetectorConfig_DCDriftParam(dcDriftParamFile)) {
         std::cerr << "[FilterTimeFrameSliceByTrack::InitTask] Failed to parse DC drift parameter configuration file: " << dcDriftParamFile << std::endl;
      }
   }

   DefineDetectorIdMap(); // for reading detector configuration files from AnalyzerT103
   RegisterDetectorConfig_Geometry();
   RegisterDetectorConfig_DCTdcCalib();
   RegisterDetectorConfig_DCDriftParam();
} // void FilterTimeFrameSliceByTrack::InitTask()

bool FilterTimeFrameSliceByTrack::ProcessSlice(TTF& tf)
{
   // std::cout << "[FilterTimeFrameSliceByTrack::ProcessSlice] Function called" << std::endl;
   // std::cout << "\tchecking TLF TDC 4ns unit: " << fLFTDC4n << std::endl;

   

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
   return true;
}

bool FilterTimeFrameSliceByTrack::ParseDetectorConfig_Geometry(std::string_view filename)
{
   std::ifstream ifs(filename.data());
   if(!ifs.is_open()){
      std::cerr << "[FilterTimeFrameSliceByTrack::ParseDetectorConfig_Geometry] Failed to open file: " << filename << std::endl;
      return false;
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
         std::cerr << "[FilterTimeFrameSliceByTrack::ParseDetectorConfig_Geometry] Failed to parse line: " << line << std::endl;
         continue;
      }
      fTemporaryGeometries.push_back({detectoridentifier, detectorname, x, y, z, tiltangle, rotationangle1, rotationangle2, length, resolution, wirecenternumber, wirepitch, offset});
   } // while(std::getline(ifs, line))

   std::cout << "[FilterTimeFrameSliceByTrack::ParseDetectorConfig_Geometry] Parsed " << fTemporaryGeometries.size() << " geometry entries from file: " << filename << std::endl;

   // to be implemented
   return true;
}

bool FilterTimeFrameSliceByTrack::ParseDetectorConfig_DCTdcCalib(std::string_view filename)
{
   std::ifstream ifs(filename.data());
   if(!ifs.is_open()){
      std::cerr << "[FilterTimeFrameSliceByTrack::ParseDetectorConfig_DCTdcCalib] Failed to open file: " << filename << std::endl;
      return false;
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
         std::cerr << "[FilterTimeFrameSliceByTrack::ParseDetectorConfig_DCTdcCalib] Failed to parse line: " << line << std::endl;
         continue;
      }
      fTemporaryDCTdcCalibs.push_back({detectoridentifier, wireidentifier, offset, scale});
   } // while(std::getline(ifs, line))

   std::cout << "[FilterTimeFrameSliceByTrack::ParseDetectorConfig_DCTdcCalib] Parsed " << fTemporaryDCTdcCalibs.size() << " TDC calibration entries from file: " << filename << std::endl;

   // to be implemented
   return true;
}

bool FilterTimeFrameSliceByTrack::ParseDetectorConfig_DCDriftParam(std::string_view filename)
{
   std::ifstream ifs(filename.data());
   if(!ifs.is_open()){
      std::cerr << "[FilterTimeFrameSliceByTrack::ParseDetectorConfig_DCDriftParam] Failed to open file: " << filename << std::endl;
      return false;
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
         std::cerr << "[FilterTimeFrameSliceByTrack::ParseDetectorConfig_DCDriftParam] Failed to parse line: " << line << std::endl;
         continue;
      }
      coefficients.resize(approxOrder);
      for(int i = 0; i < approxOrder; ++i){
         if(!(iss >> coefficients[i])){
            std::cerr << "[FilterTimeFrameSliceByTrack::ParseDetectorConfig_DCDriftParam] Failed to parse coefficient " << i << " in line: " << line << std::endl;
            continue;
         }
      }
      fTemporaryDCDriftParams.push_back({detectoridentifier, approxOrder, coefficients});
   } // while(std::getline(ifs, line))

   std::cout << "[FilterTimeFrameSliceByTrack::ParseDetectorConfig_DCDriftParam] Parsed " << fTemporaryDCDriftParams.size() << " drift parameter entries from file: " << filename << std::endl;

   // to be implemented
   return true;
}

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
   std::cout << "\n[FilterTimeFrameSliceByTrack::RegisterDetectorConfig_Geometry] Registering geometry configurations..." << std::endl;
   for(const auto& geom : fTemporaryGeometries){
      #if CHECK_COUT_DETCONF_REGISTERING
      std::cout << "[FilterTimeFrameSliceByTrack::RegisterDetectorConfig_Geometry] Registering geometry for detector ID: " << geom.detectoridentifier << std::endl;
      #endif

      int detectorId = geom.detectoridentifier;
      std::string detectorName = geom.detectorname;
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

      std::string DetectorName = detectorNameMap[detectorId];
      std::string PlaneName = detectorPlaneMap[detectorId];
      int SegmentNumber = detectorSegmentMap[detectorId];
      std::string ChannelName;
      uint8_t ChannelNameIdx = chmap::dictionary::queryIndex_readout_channel("0");
      std::cout << "[FilterTimeFrameSliceByTrack::RegisterDetectorConfig_Geometry] ChannelNameIdx: " << static_cast<int>(ChannelNameIdx) << std::endl;
      fChMap->detname_dictionary.invIndex(ChannelNameIdx, ChannelName);
      #if CHECK_COUT_DETCONF_REGISTERING
      std::cout << "\tDetectorName: " << DetectorName << std::endl;
      std::cout << "\tPlaneName: " << PlaneName << std::endl;
      std::cout << "\tSegmentNumber: " << SegmentNumber << std::endl;
      std::cout << "\tChannelName: " << ChannelName << std::endl;
      #endif

      // Create GeomItemDC and set its properties
      std::unique_ptr<chmap::GeomItemDC> geomitemdc = std::make_unique<chmap::GeomItemDC>();
      geomitemdc->SetGlobalPosition(x, y, z);
      geomitemdc->SetResolution(resolution, resolution, resolution);
      geomitemdc->SetRotationAngles(tiltAngle, rotationAngle1, rotationAngle2);
      geomitemdc->SetWireGeometry(wireCenterNumber, wirePitch, offset);

      // Register
      for(int i=0; i<128; ++i){
         int ChannelNumber = i;
         uint32_t dopeKey_DETtoFE;
         bool found_DETtoFE = fChMap->getDopeKey_DETtoFE(DetectorName, PlaneName, static_cast<uint8_t>(SegmentNumber), ChannelName, static_cast<uint8_t>(ChannelNumber), dopeKey_DETtoFE);
         if(found_DETtoFE){
            chmap::FEAddrItem feaddritem = fChMap->getFEAddrItem(dopeKey_DETtoFE);
            uint32_t dopeKeyFEtoDET;
            bool found_FEtoDET = fChMap->getDopeKey_FEtoDET(feaddritem.ip3rd, feaddritem.ip4th, feaddritem.ch, dopeKeyFEtoDET);
            if(found_FEtoDET){
               #if 1
               fChMap->getDETIdItem(dopeKeyFEtoDET).decode();
               #endif
               fChMap->registerDETConfSubItem<chmap::GeomItem, chmap::GeomItemDC>(dopeKeyFEtoDET, std::move(geomitemdc), &chmap::DETConfItem::geom);
            } // if(found_FEtoDET)
            else{
               std::cout << "[FilterTimeFrameSliceByTrack::RegisterDetectorConfig_Geometry] no DETIdItem found" << std::endl;
            }
         } // if(found_DETtoFE)
         else{
            std::cout << "[FilterTimeFrameSliceByTrack::RegisterDetectorConfig_Geometry] no FEAddrItem found" << std::endl;
         }
      } // for(int i=0; i<128; ++i)
   } // for(const auto& geom : fTemporaryGeometries)

   std::cout << "[FilterTimeFrameSliceByTrack::RegisterDetectorConfig_Geometry] Registered " << fTemporaryGeometries.size() << " geometry entries." << std::endl;
   return true;
} // bool FilterTimeFrameSliceByTrack::RegisterDetectorConfig_Geometry()

bool FilterTimeFrameSliceByTrack::RegisterDetectorConfig_DCTdcCalib()
{
   std::cout << "\n[FilterTimeFrameSliceByTrack::RegisterDetectorConfig_DCTdcCalib] Registering TDC calibration configurations..." << std::endl;
   for(const auto& calib : fTemporaryDCTdcCalibs){
      #if CHECK_COUT_DETCONF_REGISTERING
      std::cout << "[FilterTimeFrameSliceByTrack::RegisterDetectorConfig_DCTdcCalib] Registering TDC calibration for detector ID: " << calib.detectoridentifier << ", wire ID: " << calib.wireidentifier << std::endl;
      #endif
      int detectorId = calib.detectoridentifier;
      int wireId = calib.wireidentifier;
      double offset = calib.offset;
      double scale = calib.scale;

      std::string DetectorName = detectorNameMap[detectorId];
      std::string PlaneName = detectorPlaneMap[detectorId];
      int SegmentNumber = detectorSegmentMap[detectorId];
   } // for(const auto& calib : fTemporaryDCTdcCalibs)

   std::cout << "[FilterTimeFrameSliceByTrack::RegisterDetectorConfig_DCTdcCalib] Registered " << fTemporaryDCTdcCalibs.size() << " TDC calibration entries." << std::endl;
   return true;
} // bool FilterTimeFrameSliceByTrack::RegisterDetectorConfig_DCTdcCalib()

bool FilterTimeFrameSliceByTrack::RegisterDetectorConfig_DCDriftParam()
{
   std::cout << "\n[FilterTimeFrameSliceByTrack::RegisterDetectorConfig_DCDriftParam] Registering drift parameter configurations..." << std::endl;
   for(const auto& drift : fTemporaryDCDriftParams){
      #if CHECK_COUT_DETCONF_REGISTERING
      std::cout << "[FilterTimeFrameSliceByTrack::RegisterDetectorConfig_DCDriftParam] Registering drift parameter for detector ID: " << drift.detectoridentifier << ", approx order: " << drift.approxOrder << std::endl;
      #endif
      int detectorId = drift.detectoridentifier;
      int approxOrder = drift.approxOrder;
      const std::vector<double>& coefficients = drift.coefficients;

      std::string DetectorName = detectorNameMap[detectorId];
      std::string PlaneName = detectorPlaneMap[detectorId];
      int SegmentNumber = detectorSegmentMap[detectorId];
   } // for(const auto& drift : fTemporaryDCDriftParams)

   std::cout << "[FilterTimeFrameSliceByTrack::RegisterDetectorConfig_DCDriftParam] Registered " << fTemporaryDCDriftParams.size() << " drift parameter entries." << std::endl;
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

