/*
  DCPairHitCluster.cc

  2018/12  K.Shirotori
*/

#include "DCPairHitCluster.h"

DCPairHitCluster::DCPairHitCluster( DCLTrackHit *hitA, DCLTrackHit *hitB )
  : hitA_(hitA), hitB_(hitB), nhits_(0)
{
  if(hitA_) ++nhits_;
  if(hitB_) ++nhits_;
}

DCPairHitCluster::~DCPairHitCluster()
{
}
