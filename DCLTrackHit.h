#ifndef DCLTrackHit_h
#define DCLTrackHit_h 1

namespace nestdaq{
    class DCHit;
    class DCLTrackHit{
    public:
        DCLTrackHit(DCHit* acutualhit, int nth_, double w_, int leftright);
        ~DCLTrackHit() = default;

    private:
        DCHit* actualhit{nullptr}; // pointer to the actual DCHit object
        int nth{0}; // トラックに属する「ヒット」がacutualhitのTDC/TOTの何番目かを示すインデックス
        double w{0.0}; // 測定軸方向の座標, ヒット位置
        int leftright{0}; // -1, 1

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

        // トラックに属するかどうかのフラグを管理
        void clearFlag() { actualhit->clearFlag(nth); };
        void setFlag() { actualhit->setFlag(nth); };
        bool showFlag() const { return actualhit->showFlag(nth); };



    }; // class nestdaq::DCLTrackHit
} // namespace nestdaq

#endif // DCLTrackHit_h