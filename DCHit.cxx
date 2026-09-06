#include "DCHit.h"
#include "DCConstants.h"
#include "FilterTimeFrameSliceByTrackDebugger.h"

using nestdaq::DCHit;

bool DCHit::CalcDriftTimes(double standardTime, const DCTimeRange& DCTimeRange){
    const std::string_view funcname = "[FilterTimeFrameSliceByTrack::DCHit::CalcDriftTimes] ";
    #if CHECK_COUT_MAKEPAIRPLANEHITCLUSTER
    std::cout << funcname << "TDCs.size(): " << TDCs.size() << ", TOTs.size(): " << TOTs.size() << std::endl;
    #endif
    int dctdc_min = DCTimeRange.lower_bound;
    int dctdc_max = DCTimeRange.upper_bound;
    int dctot_min = DCTimeRange.tot_min;

    if(detid == nullptr){
        return false;
    }
    else{
        if(detid->detconf == nullptr){
            return false;
        }
        else{
            const chmap::CalibrationItem_DCTdcCalib* calibitem_dctdccalib = dynamic_cast<const chmap::CalibrationItem_DCTdcCalib*>(detid->detconf->membername_calib_dctdccalib.get());
            if(calibitem_dctdccalib == nullptr){
            return false;
            }
            double offset = calibitem_dctdccalib->GetOffset();
            double scale = calibitem_dctdccalib->GetScale();
            #if CHECK_COUT_DRIFTTIME
            std::cout << funcname << "offset: " << offset << ", scale: " << scale << std::endl;
            #endif

            if(TDCs.size() == 0 || TOTs.size() == 0){
            return false;
            }
            for(size_t i=0; i<TDCs.size(); ++i){ // already ensured that TDCs.size() == TOTs.size()
            uint32_t tdc = TDCs[i];
            uint32_t tot = TOTs[i];
            #if !FILEOUT_DRIFTTIME
            if((tot > dctot_min) && (tdc - standardTime >= dctdc_min) && (tdc - standardTime <= dctdc_max)){
            #else
            if(tot > dctot_min){
            #endif
                double driftTime = scale * (tdc - standardTime) + offset;
                #if CHECK_COUT_DRIFTTIME
                std::cout << "\tdriftTime = scale * (tdc - standardTime) + offset = " << scale << " * (" << tdc << " - " << standardTime << ") + " << offset << " = " << driftTime << std::endl;
                #endif
                DriftTimes.push_back(driftTime);         
                #if FILEOUT_DRIFTTIME
                int plane = static_cast<int>(detid->plane);
                int segment = static_cast<int>(detid->segment);
                int channel_number = static_cast<int>(detid->channel_number);
                gDebugFile << segment << " " << plane << " " << channel_number << "  " << driftTime << std::endl;
                #endif
            }
            } // for(const auto& tdc : tdcs)
        } // if(detid->detconf == nullptr)
    } // if(detid == nullptr)

    #if CHECK_COUT_MAKEPAIRPLANEHITCLUSTER
    std::cout << "\tDriftTimes.size(): " << DriftTimes.size() << std::endl;
    #endif
    return this->CalcDriftLengths();
} // bool nestdaq::FilterTimeFrameSliceByTrack::DCHit::CalcDriftTimes(double standardTime)

bool DCHit::IsValidDriftLength(double min, double max, double dl){
    if(dl < min || dl > max){
        return false;
    }
    else{
        return true;
    }
} // bool nestdaq::FilterTimeFrameSliceByTrack::DCHit::IsValidDriftLength(double min, double max)

bool DCHit::CalcDriftLengths(){
    if(detid == nullptr){
        return false;
    }
    else{
        // catching max, min drift length values from DCConstants.h
        const chmap::DETIdItem* detid = this->GetDETIdItem();
        int iplane = 4 * static_cast<int>(detid->segment) + static_cast<int>(detid->plane);
        double min_drift_length = DCConstants::fMinDLKLDC[iplane];
        double max_drift_length = DCConstants::fMaxDLKLDC[iplane];

        if(detid->detconf == nullptr){
            return false;
        }
        else{
            const chmap::CalibrationItem_DCDriftLength* calibitem_dcdriftlen = dynamic_cast<const chmap::CalibrationItem_DCDriftLength*>(detid->detconf->membername_calib_dcdriftlen.get());
            if(calibitem_dcdriftlen == nullptr){
            return false;
            }
            for(const auto& driftTime : DriftTimes){
                double driftLength = calibitem_dcdriftlen->GetDriftLength(driftTime);
                if(this->IsValidDriftLength(min_drift_length, max_drift_length, driftLength)){
                    DriftLengths.push_back(driftLength);
                }
            } // for(const auto& driftTime : driftTimes)
        } // if(detid->detconf == nullptr)
    } // if(detid == nullptr)
    #if CHECK_COUT_MAKEPAIRPLANEHITCLUSTER
    std::cout << "\tDriftLengths.size(): " << DriftLengths.size() << std::endl;
    #endif
   return true;
} // bool nestdaq::FilterTimeFrameSliceByTrack::DCHit::CalcDriftLengths()

double DCHit::GetGlobalZ() const{
   if(detid == nullptr){
      return 0.0;
   }
   else{
      if(detid->detconf == nullptr){
         return 0.0;
      }
      else{
         const chmap::GeomItemDC* geomitemdc = dynamic_cast<const chmap::GeomItemDC*>(detid->detconf->membername_geom.get());
         if(geomitemdc == nullptr){
            return 0.0;
         }
         else{
            return geomitemdc->GetGlobalZ();
         } // if(geomitemdc == nullptr)
      } // if(detid->detconf == nullptr)
   }
} // double nestdaq::FilterTimeFrameSliceByTrack::DCHit::GetGlobalZ() const

// void DCHit::SetStatusDLRange(double min, double max){ // input unit: mm
//     for(size_t i=0; i<DriftLengths.size(); ++i){
//         double dl = DriftLengths[i];
//         if(dl < DCConstants::DRIFTLENGTH_MIN || dl > DCConstants::DRIFTLENGTH_MAX){
//             DriftLengths[i] = -1.0; // mark as invalid
//         }
//     } // void nestdaq::FilterTimeFrameSliceByTrack::DCHit::SetStatusDLRange()