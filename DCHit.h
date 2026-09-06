#ifndef DCHit_h
#define DCHit_h 1

#include "FilterTimeFrameSliceByTrack.h"

namespace nestdaq {
    class DCLTrackHit;
    class FilterTimeFrameSliceByTrack;
    class DCHit {
    public:
    DCHit(){};
    DCHit(double wirePos, double wireAngle, const chmap::DETIdItem* detid) {
        this->wirePos = wirePos;
        this->wireAngle = wireAngle;
        this->detid = detid;
    };
    DCHit(double wirePos, double wireAngle, const chmap::DETIdItem* detid, uint32_t tdc, uint32_t tot) {
        this->wirePos = wirePos;
        this->wireAngle = wireAngle;
        this->detid = detid;
        this->TDCs.push_back(tdc);
        this->TOTs.push_back(tot);
    };
    ~DCHit() = default;

    private:
    double wirePos;
    double wireAngle;
    
    std::vector<uint32_t> TDCs;
    std::vector<uint32_t> TOTs;
    std::vector<double> DriftTimes;
    std::vector<double> DriftLengths;
    const chmap::DETIdItem* detid;

    mutable std::vector< DCLTrackHit* > Cont_; // このヒットがどのDCLTrackHitに属するかを登録する

    public:
    void AddHit(uint32_t tdc, uint32_t tot){
        TDCs.push_back(tdc);
        TOTs.push_back(tot);
        return;
    }
    const chmap::DETIdItem* GetDETIdItem() const { return detid; };
    int Clear(){
        int n = TDCs.size();
        TDCs.clear();
        DriftTimes.clear();
        DriftLengths.clear();
        return n;
    };
    bool CalcDriftTimes(double standardTime, const DCTimeRange& DCTimeRange);
    bool CalcDriftLengths();

    double GetWirePos() const { return wirePos; };
    double GetWireAngle() const { return wireAngle; };
    double GetDriftLength(int nth) const { return DriftLengths[nth]; };
    double GetGlobalZ() const;

    void RegisterHits( DCLTrackHit* hit) const { // このヒットがどのDCLTrackHitに属するかを登録する
        Cont_.push_back(hit);
    }
    }; // class nestdaq::FilterTimeFrameSliceByTrack::DCHit
} // namespace nestdaq
#endif // DCHit_h