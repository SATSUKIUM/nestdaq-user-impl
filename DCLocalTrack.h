#ifndef DCLOCALTRACK_H_
#define DCLOCALTRACK_H_ 1

#include <vector>
#include "DCLTrackHit.h"

namespace nestdaq{
    class DCLTrackHit;
    class DCLocalTrack{
        public:
            DCLocalTrack() = default;
            ~DCLocalTrack() = default;
        private:
            std::vector<DCLTrackHit*> dclthits;

        public:
            void AddHit(DCLTrackHit* hit){
                dclthits.push_back(hit);
            }
            bool DoFit();
            bool AngleCorrection();
            std::size_t GetNHits() const {
                return dclthits.size();
            }
            DCLTrackHit* GetHit(std::size_t i) const;

            double GetX0() const { return x0; }
            double GetY0() const { return y0; }
            double GetU0() const { return u0; }
            double GetV0() const { return v0; }
            double GetChiSqr() const { return chisqr; }
            bool GetStatus() const { return status; }

            double CalcX(double z) const { return x0 + u0 * z; }
            double CalcY(double z) const { return y0 + v0 * z; }

            // トラックに属するヒットのフラグを管理する
            void SetFlags(){
                for(auto hit : dclthits){
                    hit->setFlag();
                } // set the flag for each hit
            }
            void ClearFlags(){
                for(auto hit : dclthits){
                    hit->clearFlag();
                } // clear the flag for each hit
            }
            bool showFlag(int idclthit) const { return dclthits[idclthit]->showFlag(); }
            int GetNumOfTrueFlags() const {
                int count = 0;
                for(auto hit : dclthits){
                    if(hit->showFlag() == true){
                        ++count;
                    }
                }
                return count;
            }
            void CalcHitPositions();

        private: // fit info
            bool status; // fit?
            double x0{0.0}, y0{0.0}, u0{0.0}, v0{0.0}; // position, position, slope, slope
            double chisqr{0.0};

        private: // constants
            static constexpr int ReservedNumOfHits = 16;

    }; // class nestdaq::DCLocalTrack

    struct DCLTrackComp 
    : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
    {
    bool operator()( const DCLocalTrack * const p1, 
            const DCLocalTrack * const p2 ) const
    {
        int n1=p1->GetNHits(), n2=p2->GetNHits();
        double chi1=p1->GetChiSqr(),chi2=p2->GetChiSqr();
        if( (n1>n2+1) ){
        return true;
        }
        else if( (n2>n1+1)  ){
        return false;
        }
        else{
        return (chi1<=chi2);
        }
    }
    };

    struct DCLTrackComp1 
    : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
    {
    bool operator()( const DCLocalTrack * const left, 
            const DCLocalTrack * const right ) const
    {
        int n1=left->GetNHits(), n2=right->GetNHits();
        if(n1>n2) return true;
        else if(n2>n1) return false;
        else
        return (left->GetChiSqr()) < (right->GetChiSqr());
    }
    };

    struct DCLTrackComp2 
    : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
    {
    bool operator()( const DCLocalTrack * const p1, 
            const DCLocalTrack * const p2 ) const
    {
        int n1=p1->GetNHits(), n2=p2->GetNHits();
        if(n1<n2) return true;
        else if(n2<n1) return false;
        else
        return (p1->GetChiSqr())<=(p2->GetChiSqr());
    }
    };

    struct DCLTrackComp3 
    : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
    {
    bool operator()( const DCLocalTrack * const p1, 
            const DCLocalTrack * const p2 ) const
    {
        int n1=p1->GetNHits(), n2=p2->GetNHits();
        double chi1=p1->GetChiSqr(),chi2=p2->GetChiSqr();
        double a1=fabs(1.-chi1),a2=fabs(1.-chi2);
        if(a1<a2) return true;
        else if(a2<a1) return false;
        else
        return (n1<=n2);
    }
    };

    struct DCLTrackComp4 
    : public std::binary_function <DCLocalTrack *, DCLocalTrack *, bool>
    {
    bool operator()( const DCLocalTrack * const p1, 
            const DCLocalTrack * const p2 ) const
    {
        int n1=p1->GetNHits(), n2=p2->GetNHits();
        double chi1=p1->GetChiSqr(),chi2=p2->GetChiSqr();
        //if( (n1>n2+1) ){
        //    if( (n1>n2+1) && (fabs(chi1-chi2)<5.) ){
        if( (n1>n2+1) && (fabs(chi1-chi2)<2.) ){
        return true;
        }
        //else if( (n2>n1+1)  ){
        //    else if( (n2>n1+1) && (fabs(chi1-chi2)<5.) ){
        else if( (n2>n1+1) && (fabs(chi1-chi2)<2.) ){
        return false;
        }
        else{
        return (chi1<=chi2);
        }
    }
    };

} // namespace nestdaq

#endif // DCLOCALTRACK_H_