#include "DCLTrackHit.h"

using nestdaq::DCLTrackHit;

DCLTrackHit::DCLTrackHit(const DCHit* parent_, int nth_, double w_, int leftright) : parent(parent_), nth(nth_), w(w_), leftright(leftright) {
    parent->RegisterHits(this);
};