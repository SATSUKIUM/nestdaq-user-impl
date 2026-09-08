#ifndef DCHIT_H_
#define DCHIT_H_ 1

#include "FilterTimeFrameSliceByTrack.h"
#include "DCTimeRange.h"

namespace nestdaq {
    class DCLTrackHit;
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
        
        std::vector<uint32_t> TDCs; // unit: ns, TDCs.size() == TOTs.size()
        std::vector<uint32_t> TOTs; // unit: ns
        std::vector<double> DriftTimes; // unit: ns

        std::vector<double> DriftLengths; // unit: mm, DriftLengths.size() == DriftTimes.size()
        std::vector<bool> IsBelongToGoodTrack; // トラックをヒット数やChiSquareでソートしたのちに、上位のトラックに属するヒットはtrueにされる
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
        double GetResolution() const;
        int GetNumMultiHit() const { return DriftLengths.size(); };
        bool IsValidDriftLength(double min, double max, double dl); // input unit: mm
        bool RangeCheck(int nth) const { return (nth >= 0 && nth < DriftLengths.size()); };

        void RegisterHits( DCLTrackHit* hit) const {
            Cont_.push_back(hit);
        }
        void clearFlag(int nth) { IsBelongToGoodTrack[nth] = false; }
        void setFlag(int nth) { IsBelongToGoodTrack[nth] = true; }
        bool showFlag(int nth) const { return IsBelongToGoodTrack[nth]; }
    }; // class nestdaq::FilterTimeFrameSliceByTrack::DCHit
} // namespace nestdaq
#endif // DCHIT_H_