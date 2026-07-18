/*
  ConfMan.hh

  2024/11 K.Shirotori
*/

#ifndef Confman_h
#define Confman_h 1

#include <string>

class VEvent;

class GeomMan;
class DCTdcCalibMan;
class DCDriftParamMan;
class TrTdcCalibMan;
class HodoParamMan;
class HodoPHCMan;

class ConfMan
{

private:
   //Conf
   std::string ConfFileName_;
   static ConfMan *confManager_;
   
   //Geometry files  
   std::string GeomFileName_;
   GeomMan *GeomManager_;

   //DC
   std::string      DCTdcCalibFileName_;
   DCTdcCalibMan   *DCTdcCalibManager_;
   std::string      DCDriftParamFileName_;
   DCDriftParamMan *DCDriftParamManager_;

   int BDCTRangeLow_, BDCTRangeHigh_, BDCTRangeTOT_;
   int BDCTRMode_;

   int KLDCTRangeLow_, KLDCTRangeHigh_, KLDCTRangeTOT_;
   int KLDCTRMode_;
   int KLDCTRRoutine_;

   //Tr
   std::string      TrTdcCalibFileName_;
   TrTdcCalibMan   *TrTdcCalibManager_;

   int BFTTRangeLow_, BFTTRangeHigh_, BFTTRangeTOT_;

   //Hodo
   std::string   HodoParamFileName_;
   HodoParamMan *HodoParamManager_;
   std::string   HodoPHCFileName_;
   HodoPHCMan   *HodoPHCManager_;
   
public:
  ConfMan( const std::string & filename );
   ~ConfMan();
   
   static ConfMan *GetConfManager( void ) { return confManager_; }
   bool Initialize( void );
   void SetFileName( const std::string & filename ) { ConfFileName_=filename; }
   VEvent* EventAllocator();
   
   //Geometry files  
   GeomMan *GetGeomManager( void ) { return GeomManager_; }

   //DC
   DCTdcCalibMan   *GetDCTdcCalibManager( void ) { return DCTdcCalibManager_; }
   DCDriftParamMan *GetDCDriftParamManager( void ) { return DCDriftParamManager_; } 

   double BDCTRangeLow( void ) const { return BDCTRangeLow_; }
   double BDCTRangeHigh( void ) const { return BDCTRangeHigh_; }
   double BDCTRangeTOT( void ) const { return BDCTRangeTOT_; }
   int BDCTRMode( void ) const { return BDCTRMode_; }

   double KLDCTRangeLow( void ) const { return KLDCTRangeLow_; }
   double KLDCTRangeHigh( void ) const { return KLDCTRangeHigh_; }
   double KLDCTRangeTOT( void ) const { return KLDCTRangeTOT_; }
   int KLDCTRMode( void ) const { return KLDCTRMode_; }
   int KLDCTRRoutine( void ) const { return KLDCTRRoutine_; }

   //Tr
   TrTdcCalibMan   *GetTrTdcCalibManager( void ) { return TrTdcCalibManager_; }

   double BFTTRangeLow( void ) const { return BFTTRangeLow_; }
   double BFTTRangeHigh( void ) const { return BFTTRangeHigh_; }
   double BFTTRangeTOT( void ) const { return BFTTRangeTOT_; }
   
   //Hodo
   HodoParamMan *GetHodoParamManager( void ) { return HodoParamManager_;}
   HodoPHCMan   *GetHodoPHCManager( void ) { return HodoPHCManager_;}
   
private:
   ConfMan( const ConfMan & );
   ConfMan & operator = ( const ConfMan & );
   
public:
   bool InitializeParameterFiles();
   bool InitializeHistograms();
};

#endif
