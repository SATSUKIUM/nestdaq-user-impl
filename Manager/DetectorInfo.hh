/*
  DetectorID.hh
*/

#ifndef DetectorInfo_h
#define DetectorInfo_h 1

//Detector ID
//DCs FTs
const int DetIdBDC1  = 101;
const int DetIdBDC2  = 102;
const int DetIdKLDC1 = 103;
const int DetIdKLDC2 = 104;
const int DetIdBFT1  = 105;
const int DetIdBFT2  = 106;
const int DetIdSFT1  = 107;
const int DetIdSFT2  = 108;

//Timing counters
const int DetIdUTOF = 3;
const int DetIdDTOF = 4;
const int DetIdLTOF = 5;
const int DetIdT0   = 6;
const int DetIdT0r  = 7;
const int DetIdBref = 8;
const int DetIdT1   = 9;
const int DetIdBHT  = 10;

//DC/FT Number of Layers
const int NumOfLayersBDC  = 8;
const int NumOfLayersKLDC = 8;
const int NumOfLayersBFT  = 6;
const int NumOfLayersSFT  = 6;

const int PlOffsBDC  = 100;
const int PlOffsKLDC = 200;
const int PlOffsBFT  = 300;
const int PlOffsSFT  = 400;

//Timing counterSegments
const int NumOfSegUTOF = 1;
const int NumOfSegDTOF = 3;
const int NumOfSegLTOF = 6;
const int NumOfSegT0   = 8;
const int NumOfSegT0r  = 1;
const int NumOfSegBref = 2;
const int NumOfSegT1   = 1;
const int NumOfSegBHT  = 1;

//Number of Wires
//Drift chamber
const int NumOfWireBDCX  = 116;
const int NumOfWireBDCUV = 109;
const int NumOfWireKLDC  = 128;

//Fiber tracker
const int NumOfFiberBFT   = 256;
const int NumOfFiberSFTX  = 1152;
const int NumOfFiberSFTUV = 1024;

#endif
