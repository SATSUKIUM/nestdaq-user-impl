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
      std::cout << "\t" << "-> found." << std::endl;
      feaddritem = fChMap->getFEAddrItem(dopeKey_DETtoFE);
      feaddritem.decode();
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
      std::cout << "\t" << "-> found." << std::endl;
      feaddritem = fChMap->getFEAddrItem(dopeKey_DETtoFE);
      feaddritem.decode();
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
   // kldc 2 U 95
   std::cout << "\n\t" << "checking kldc 2 U 95 DETIdItem search..." << std::endl;
   const std::string_view test_detectorname_kldc2u95 = "kldc";
   const std::string_view test_planename_kldc2u95 = "U";
   const uint8_t test_segment_kldc2u95 = 2;
   const uint16_t test_channelnumber_kldc2u95 = 95;
   const std::string_view test_channelname_kldc2u95 = "0";
   _FOUND_DETtoFE = fChMap->getDopeKey_DETtoFE(test_detectorname_kldc2u95, test_planename_kldc2u95, test_segment_kldc2u95, test_channelname_kldc2u95, test_channelnumber_kldc2u95, dopeKey_DETtoFE);
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
            std::cout << "\t" << "-> found detector configuration for kldc 2 U 95." << std::endl;
         }
         else{
            std::cout << "\t" << "-> not found detector configuration for kldc 2 U 95." << std::endl;
         }
      }
   }
   else{
      std::cout << "\t" << "-> not found." << std::endl;
   }
   #endif
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

      // Create GeomItemDC and set its properties
      std::unique_ptr<chmap::GeomItemDC> geomitemdc = std::make_unique<chmap::GeomItemDC>();
/* CLASS DEFINITION
    class GeomItemDC : public GeomItem {
        public:
            virtual ~GeomItemDC() = default;
            void SetWireGeometry(double centerWireNumber_, double wirePitch_, double offset_) {
                centerWireNumber = centerWireNumber_;
                wirePitch = wirePitch_;
                offset = offset_;
            }
        private:
            double centerWireNumber; // もし1.0なら、中心のワイヤーは1番ワイヤー。0.5なら、中心のワイヤーは1番と2番の間にある。
            double wirePitch; // [mm] 測定軸方向のワイヤ間隔
            double offset; // [mm] 測定軸方向のワイヤのオフセット(微調整のため)
    }; // class GeomItemDC
*/
      geomitemdc->SetGlobalPosition(x, y, z);
      geomitemdc->SetResolution(resolution, resolution, resolution);
      geomitemdc->SetRotationAngles(tiltAngle, rotationAngle1, rotationAngle2);
      geomitemdc->SetWireGeometry(wireCenterNumber, wirePitch, offset);

      int missing_count_FEAddrItem = 0;
      int missing_count_DETIdItem = 0;
      // Register KLDC
      if(DetectorName == "kldc"){
         for(int i=0+32; i<128-32; ++i){
            int ChannelNumber = i;
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

