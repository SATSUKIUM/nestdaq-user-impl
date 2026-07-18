/*
  DCParameters.hh

  2024/05  K.Shirotori
*/

#ifndef DCParameters_h 
#define DCParameters_h 1 

struct DCPairPlaneInfo 
{
  bool flag;
  int id1, id2;
  double CellSize;
};

struct DCComPlaneInfo
{
  int id1, id2;
  DCComPlaneInfo(int i1, int i2) : id1(i1), id2(i2) {}
  const int nBrotherPlanes = 2;
  int GetNLayers() const { return nBrotherPlanes; };
};

extern const DCPairPlaneInfo PPInfoBDC[], PPInfoKLDC[];
extern const int NPPInfoBDC, NPPInfoKLDC;
extern const int NDimKLDC;
extern const DCComPlaneInfo CoPlaneInfoKLDC[];

#ifdef DefStatic
const DCPairPlaneInfo PPInfoBDC[] = {
  { true,  1,  2, 20.82 }, { false,  3,  4, 20.82 },
  { true,  5,  6, 20.82 }, { false,  7,  8, 20.82 }, 
};

const DCPairPlaneInfo PPInfoKLDC[] = {
  { true,  1,  2, 9.007 }, { true,  3,  4, 9.007 },
  { true,  5,  6, 9.007 }, { true,  7,  8, 9.007 }, 
};

const DCComPlaneInfo CoPlaneInfoKLDC[] = { // PPInfoKLDC[0] と PPInfoKLDC[3]が同じ測定軸を持つ組み合わせであるという解釈をして使う。
  {0, 3}, {1, 2}, // 0: KLDC1 VV', 1: KLDC2 V'V, 2: KLDC1 UU', 3: KLDC2 U'U. Reference below. (from ConfMan call for initialize GeomMan with filename "param/param_DetGeom001".)
};
// ================================
// Reference
/*
# KLDC
 201 KLDC1-v-1  00.00    00.00   499.00  -45.00   00.00    00.00   499.00   0.150   64.00  9.007  -0.02
 202 KLDC1-v-2  00.00    00.00   506.80  -45.00   00.00    00.00   506.80   0.150   64.50  9.007   0.04
 203 KLDC1-u-1  00.00    00.00   523.20   45.00   00.00    00.00   523.20   0.150   64.50  9.007  -0.05
 204 KLDC1-u-2  00.00    00.00   531.00   45.00   00.00    00.00   531.00   0.150   64.00  9.007   0.01
 205 KLDC2-u-2  00.00    00.00   626.80   45.00   00.00    00.00   626.80   0.150   64.50  9.007   0.02
 206 KLDC2-u-1  00.00    00.00   619.00   45.00   00.00    00.00   619.00   0.150   64.00  9.007  -0.05
 207 KLDC2-v-2  00.00    00.00   651.00  -45.00   00.00    00.00   651.00   0.150   64.00  9.007   0.02
 208 KLDC2-v-1  00.00    00.00   643.20  -45.00   00.00    00.00   643.20   0.150   64.50  9.007  -0.01
#
*/


const int NPPInfoBDC  = sizeof(PPInfoBDC)/sizeof(DCPairPlaneInfo);
const int NPPInfoKLDC = sizeof(PPInfoKLDC)/sizeof(DCPairPlaneInfo);
const int NDimKLDC = sizeof(CoPlaneInfoKLDC)/sizeof(DCComPlaneInfo);

#endif

const int MinNumOfHitsBDC  = 6;
const int MinNumOfHitsKLDC = 3; // もともと6だったけど、フィッティング方法を変えたので3に

// DL Ranges
const double MinDLBDC[9] = {
   0.0,
  -0.3, -0.3, -0.3, -0.3,
  -0.3, -0.3, -0.3, -0.3
}; 

const double MaxDLBDC[9] = {
    0.0,
   10.7, 10.7, 10.7, 10.7,
   10.7, 10.7, 10.7, 10.7
};

const double MinDLKLDC[9] = {
   0.0,
  -0.3, -0.3, -0.3, -0.3,
  -0.3, -0.3, -0.3, -0.3
}; 

const double MaxDLKLDC[9] = {
   0.0,
   4.8, 4.8, 4.8, 4.8,
   4.8, 4.8, 4.8, 4.8
};

#endif
