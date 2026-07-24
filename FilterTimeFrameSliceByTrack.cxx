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
#include <stdexcept>

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
   chmap::ChannelMapSimpleItem_DET detitem;
   std::cout << "\tchecking T1 right search" << std::endl;
   bool found = fChMap->getDopeKey_FEtoDET(test_ip3rd_T1right, test_ip4th_T1right, test_ch_T1right, dope_key);
   if(found == true){
      std::cout << "\t\tfound." << std::endl;
      detitem = fChMap->getDETItem(dope_key);
      detitem.decode();
   }
   #endif

   // ================================
   // detector configuration initialization
   // ================================

   const auto geometryFile = fConfig->GetProperty<std::string>(opt::GeometryConfigFile.data());
   const auto dctdcCalibFile = fConfig->GetProperty<std::string>(opt::DCTdcCalibConfigFile.data());
   const auto dcDriftParamFile = fConfig->GetProperty<std::string>(opt::DCDriftParamConfigFile.data());

   if (!geometryFile.empty()) {
      if (!ParseDetectorConfig_Geometry(geometryFile)) {
         throw std::runtime_error("Failed to parse geometry configuration");
      }
   }

   if (!dctdcCalibFile.empty()) {
      if (!ParseDetectorConfig_DCTdcCalib(dctdcCalibFile)) {
         throw std::runtime_error("Failed to parse DC TDC calibration configuration");
      }
   }

   if (!dcDriftParamFile.empty()) {
      if (!ParseDetectorConfig_DCDriftParam(dcDriftParamFile)) {
         throw std::runtime_error("Failed to parse DC drift parameter configuration");
      }
   }
}

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

   struct temporary_geometry {
      int detectoridentifier;
      std::string detectorname;
      double x, y, z;
      double tiltangle, rotationangle1, rotationangle2;
      double length, resolution, wirecenternumber, wirepitch, offset;
   };
   std::vector<temporary_geometry> temp_geometries;

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
      temp_geometries.push_back({detectoridentifier, detectorname, x, y, z, tiltangle, rotationangle1, rotationangle2, length, resolution, wirecenternumber, wirepitch, offset});
   } // while(std::getline(ifs, line))

   std::cout << "[FilterTimeFrameSliceByTrack::ParseDetectorConfig_Geometry] Parsed " << temp_geometries.size() << " geometry entries from file: " << filename << std::endl;

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

   struct temporary_dctdccalib {
      int detectoridentifier;
      int wireidentifier;
      double offset, scale;
   };
   std::vector<temporary_dctdccalib> temp_dctdccalibs;

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
      temp_dctdccalibs.push_back({detectoridentifier, wireidentifier, offset, scale});
   } // while(std::getline(ifs, line))

   std::cout << "[FilterTimeFrameSliceByTrack::ParseDetectorConfig_DCTdcCalib] Parsed " << temp_dctdccalibs.size() << " TDC calibration entries from file: " << filename << std::endl;

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

   struct temporary_dcdriftparam {
      int detectoridentifier;
      int approxOrder;
      std::vector<double> coefficients;
   };
   std::vector<temporary_dcdriftparam> temp_dcdriftparams;

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
      temp_dcdriftparams.push_back({detectoridentifier, approxOrder, coefficients});
   } // while(std::getline(ifs, line))

   std::cout << "[FilterTimeFrameSliceByTrack::ParseDetectorConfig_DCDriftParam] Parsed " << temp_dcdriftparams.size() << " drift parameter entries from file: " << filename << std::endl;

   // to be implemented
   return true;
}



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

