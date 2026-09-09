#include "DCLTrackHit.h"
#include "DCHit.h"
#include "FilterTimeFrameSliceByTrackDebugger.h"

using nestdaq::DCLTrackHit;

DCLTrackHit::DCLTrackHit(DCHit* actualhit, int nth_, double w_, int leftright) : actualhit(actualhit), nth(nth_), w(w_), leftright(leftright) {
    actualhit->RegisterHits(this);
};

// ヒットごとに変わらない情報
double DCLTrackHit::GetWirePosition() const { return actualhit->GetWirePos(); };
double DCLTrackHit::GetWireAngle() const { return actualhit->GetWireAngle(); };
double DCLTrackHit::GetResolution() const { return actualhit->GetResolution(); };
double DCLTrackHit::GetGlobalZ() const{
    return actualhit->GetGlobalZ();
};

// ヒットごとに変わる情報
double DCLTrackHit::GetDriftLength() const { return actualhit->GetDriftLength(nth); };

void DCLTrackHit::clearFlag() { actualhit->clearFlag(nth); };
void DCLTrackHit::setFlag() { actualhit->setFlag(nth); };
bool DCLTrackHit::showFlag() const { return actualhit->showFlag(nth); };