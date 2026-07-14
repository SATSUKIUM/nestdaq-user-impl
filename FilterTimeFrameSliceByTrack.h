/**
 * @file FilterTimeFrameSliceByTrack.h
 * @brief Abstract base class for filtering time frame slice
 * @date Created : 2026-07-14 13:41:20 JST
 *       Last Modified : 2026-07-14 13:41:20 JST
 *
 * @author Shinsuke OTA <ota@rcnp.osaka-u.ac.jp>
 * @author Shunichi KASHIMA <kashima@rcnp.osaka-u.ac.jp>
 *
 */


#ifndef NESTDAQ_TIMEFRAMESLICERBYTRACK_H
#define NESTDAQ_TIMEFRAMESLICERBYTRACK_H

#include "fairmq/Device.h"
#include "KTimer.cxx"
#include "SubTimeFrameHeader.h"
#include "TimeFrameHeader.h"
#include "FilterHeader.h"
#include "HeartbeatFrameHeader.h"
#include "FrameContainer.h"
#include "FilterTimeFrameSliceABC.h"

namespace nestdaq {
   class FilterTimeFrameSliceByTrack;
}


class nestdaq::FilterTimeFrameSliceByTrack : public nestdaq::FilterTimeFrameSliceABC {
public:
   FilterTimeFrameSliceByTrack();
   virtual ~FilterTimeFrameSliceByTrack() override = default;

   virtual bool ProcessSlice(TTF& ) override;

};


#endif  // NESTDAQ_TIMEFRAMESLICERBYTRACK_H
