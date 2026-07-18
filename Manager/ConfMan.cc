/*
  ConfMan.cc

  2024/11 K.Shirotori
*/

#include "ConfMan.hh"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

#include "GeomMan.hh"
#include "DCTdcCalibMan.hh"
#include "DCDriftParamMan.hh"
#include "TrTdcCalibMan.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"

const std::string defGeomFile="Geom.param";

ConfMan * ConfMan::confManager_ = 0;

ConfMan::ConfMan( const std::string & filename )
   : ConfFileName_(filename),
     GeomFileName_(defGeomFile),
     DCTdcCalibManager_(0),
     DCDriftParamManager_(0),
     BDCTRangeLow_(0), BDCTRangeHigh_(0),
     BDCTRangeTOT_(0), BDCTRMode_(0),
     KLDCTRangeLow_(0), KLDCTRangeHigh_(0),
     KLDCTRangeTOT_(0), KLDCTRMode_(0),
     TrTdcCalibManager_(0),
     BFTTRangeLow_(0), BFTTRangeHigh_(0),
     BFTTRangeTOT_(0),
     HodoParamManager_(0),
     HodoPHCManager_(0)
{
  static const std::string funcname="[ConfMan::ConfMan]";
  if( confManager_ ){
    std::cerr << funcname << ": constring twice" << std::endl;
    std::exit(-1);
  }
  confManager_=this;
}

ConfMan::~ConfMan( )
{
  confManager_=0;
}

const int BufSize=300;

bool ConfMan::Initialize()
{
  static const std::string funcname="[ConfMan::Initialize]";

  char buf[BufSize], buf1[BufSize];
  double val1, val2, val3;
  int id;

  FILE *fp;
  if((fp=fopen(ConfFileName_.c_str(),"r"))==0){
    std::cerr << funcname << ": file open fail" << std::endl;
    exit(-1);
  }

  while( fgets( buf, BufSize, fp ) != 0 ){
    if( buf[0]!='#' ){ 
      // Geometry
      if( sscanf(buf,"GEO: %s",buf1)==1 ){
        GeomFileName_=buf1;
      }
      else if( sscanf(buf,"DCTDC: %s",buf1)==1 ){
        DCTdcCalibFileName_=buf1;
      }
      else if( sscanf(buf,"DCDRFT: %s",buf1)==1 ){
        DCDriftParamFileName_=buf1;
      }
      else if( sscanf(buf,"KLDCTRMODE: %d", &id )==1 ){
         KLDCTRMode_=id;
      }
      else if( sscanf(buf,"KLDCTRROUTINE: %d", &id )==1 ){
         KLDCTRRoutine_=id;
      }
      else if( sscanf(buf,"BDCTRANGELOW: %lf", &val1 )==1 ){
         BDCTRangeLow_=val1;
      }
      else if( sscanf(buf,"BDCTRANGEHIGH: %lf", &val1 )==1 ){
         BDCTRangeHigh_=val1;
      }
      else if( sscanf(buf,"BDCTRANGETOT: %lf", &val1 )==1 ){
         BDCTRangeTOT_=val1;
      }
      else if( sscanf(buf,"BDCTRMODE: %d", &id )==1 ){
         BDCTRMode_=id;
      }
      else if( sscanf(buf,"KLDCTRANGELOW: %lf", &val1 )==1 ){
         KLDCTRangeLow_=val1;
      }
      else if( sscanf(buf,"KLDCTRANGEHIGH: %lf", &val1 )==1 ){
         KLDCTRangeHigh_=val1;
      }
      else if( sscanf(buf,"KLDCTRANGETOT: %lf", &val1 )==1 ){
         KLDCTRangeTOT_=val1;
      }
      else if( sscanf(buf,"TRTDC: %s",buf1)==1 ){
        TrTdcCalibFileName_=buf1;
      }
      else if( sscanf(buf,"BFTTRANGELOW: %lf", &val1 )==1 ){
         BFTTRangeLow_=val1;
      }
      else if( sscanf(buf,"BFTTRANGEHIGH: %lf", &val1 )==1 ){
         BFTTRangeHigh_=val1;
      }
      else if( sscanf(buf,"BFTTRANGETOT: %lf", &val1 )==1 ){
         BFTTRangeTOT_=val1;
      }
      else if( sscanf(buf,"HDPRM: %s",buf1)==1 ){
        HodoParamFileName_=buf1;
      }
      else if( sscanf(buf,"HDPRC: %s",buf1)==1 ){
        HodoPHCFileName_=buf1;
      }
      else {
        std::cout << funcname << ": un-recognized record\n"
        	  << buf << std::endl;
      }
    }
  }

  fclose(fp);

  std::cout << "------------------" << ConfFileName_ << "------" << std::endl;
  std::cout << "Geom. Param.:  " << GeomFileName_ << std::endl;
  std::cout << "--------------" << std::endl;
  std::cout << "DC TDC Calib. Param.:  " << DCTdcCalibFileName_ << std::endl;
  std::cout << "DC Drift Time Param.:  " << DCDriftParamFileName_ << std::endl;
  std::cout << "BDC TDC cut range:  " << BDCTRangeLow_ << " - " <<  BDCTRangeHigh_ << std::endl;
  std::cout << "BDC TOT cut range:  " << "TOT > " <<  BDCTRangeTOT_ << std::endl;
  std::cout << "BDC Traking mode:  " <<  BDCTRMode_ << std::endl;
  std::cout << "KLDC TDC cut range:  " << KLDCTRangeLow_ << " - " <<  KLDCTRangeHigh_ << std::endl;
  std::cout << "KLDC TOT cut range:  " << "TOT > " <<  KLDCTRangeTOT_ << std::endl;
  std::cout << "KLDC Traking mode:  " <<  KLDCTRMode_ << std::endl;
  std::cout << "KLDC Tracking routine:  " << KLDCTRRoutine_ << std::endl;
  std::cout << "--------------" << std::endl;
  std::cout << "BFT TDC cut range:  " << BFTTRangeLow_ << " - " <<  BFTTRangeHigh_ << std::endl;
  std::cout << "BFT TOT cut range:  " << "TOT > " <<  BFTTRangeTOT_ << std::endl;
  std::cout << "--------------" << std::endl;
  std::cout << "Hodoscope TDC Calib. Param.:  " <<  HodoParamFileName_ << std::endl;
  std::cout << "Hodo PHC Param.:  " << HodoPHCFileName_ << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;

  InitializeParameterFiles();
  InitializeHistograms();

  return true;
}

bool ConfMan::InitializeParameterFiles( void )
{
   GeomManager_ = & GeomMan::GetInstance();
   if( GeomFileName_!="" )
      GeomManager_->Initialize(GeomFileName_);
   else
      GeomManager_->Initialize();
   
   if( DCTdcCalibFileName_!="" )
      DCTdcCalibManager_ = new DCTdcCalibMan(DCTdcCalibFileName_);
   if(DCTdcCalibManager_) DCTdcCalibManager_->Initialize();
   
   if( DCDriftParamFileName_!="" )
      DCDriftParamManager_ = new DCDriftParamMan(DCDriftParamFileName_);
   if(DCDriftParamManager_) DCDriftParamManager_->Initialize();

   if( TrTdcCalibFileName_!="" )
      TrTdcCalibManager_ = new TrTdcCalibMan(TrTdcCalibFileName_);
   if(TrTdcCalibManager_) TrTdcCalibManager_->Initialize();

   if( HodoParamFileName_!="" )
      HodoParamManager_ = new HodoParamMan(HodoParamFileName_);
   if(HodoParamManager_) HodoParamManager_->Initialize();
   
   if( HodoPHCFileName_!="" )
      HodoPHCManager_ = new HodoPHCMan(HodoPHCFileName_);
   if(HodoPHCManager_) HodoPHCManager_->Initialize();
   
   return true;
}
