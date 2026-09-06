#ifndef DCLTrackHit_h
#define DCLTrackHit_h 1

#include "FilterTimeFrameSliceByTrack.h"

using nestdaq::DCHit;

class DCLTrackHit{
public:
    DCLTrackHit(DCHit* parent_, int nth_, double w_, int leftright) : parent(parent_), nth(nth_), w(w_), leftright(leftright) {
        parent->RegisterHits(this);
    };
    ~DCLTrackHit() = default;

private:
    DCHit* parent;
    int nth;
    double w; // 測定軸方向の座標, ヒット位置
    int leftright; // -1, 1

public:
    double GetWirePosition() const { return parent->GetWirePos(); };
    double GetWireAngle() const { return parent->GetWireAngle(); };
    double GetDriftLength() const { return parent->GetDriftLength(nth); };

    double GetGlobalZ() const{ return parent->GetGlobalZ(); };



} // class nestdaq::DCLTrackHit

#endif // DCLTrackHit_h