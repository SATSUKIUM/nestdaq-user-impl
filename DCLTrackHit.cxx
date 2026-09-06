#include "DCLTrackHit.h"
#include "DCHit.h"

using nestdaq::DCLTrackHit;

DCLTrackHit::DCLTrackHit(const DCHit* parent_, int nth_, double w_, int leftright) : parent(parent_), nth(nth_), w(w_), leftright(leftright) {
    parent->RegisterHits(this);
};

double DCLTrackHit::GetWirePosition() const { return parent->GetWirePos(); };
double DCLTrackHit::GetWireAngle() const { return parent->GetWireAngle(); };
double DCLTrackHit::GetDriftLength() const { return parent->GetDriftLength(nth); };
double DCLTrackHit::GetGlobalZ() const{ return parent->GetGlobalZ(); };