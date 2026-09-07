#ifndef DCLTrackHit_h
#define DCLTrackHit_h 1

namespace nestdaq{
    class DCHit;
    class DCLTrackHit{
    public:
        DCLTrackHit(const DCHit* parent_, int nth_, double w_, int leftright);
        ~DCLTrackHit() = default;

    private:
        const DCHit* parent;
        int nth;
        double w; // 測定軸方向の座標, ヒット位置
        int leftright; // -1, 1

    private:
        // position and slope calculated with the track fit result
        double x0{0.0}, y0{0.0}; // position, position
    public:
        void SetCalPosition(double x0_, double y0_){
            x0 = x0_;
            y0 = y0_;
        }

    public:
        double GetWirePosition() const;
        double GetWireAngle() const;
        double GetDriftLength() const;
        double GetResolution() const;
        int GetLeftRight() const { return leftright; };

        double GetGlobalZ() const;



    }; // class nestdaq::DCLTrackHit
} // namespace nestdaq

#endif // DCLTrackHit_h