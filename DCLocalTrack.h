#ifndef DCLOCALTRACK_H_
#define DCLOCALTRACK_H_ 1

#include <vector>

namespace nestdaq{
    class DCLTrackHit;
    class DCLocalTrack{
        public:
            DCLocalTrack() = default;
            ~DCLocalTrack() = default;
        private:
            std::vector<DCLTrackHit*> hits;
    };

} // namespace nestdaq

#endif // DCLOCALTRACK_H_