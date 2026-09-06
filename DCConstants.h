#ifndef DCCONSTANTS_H
#define DCCONSTANTS_H 1

namespace DCConstants{
    // max, min value of drift length for KLDC, unit: mm
    const double fMinDLKLDC[8] = {
    -0.3, -0.3, -0.3, -0.3,
    -0.3, -0.3, -0.3, -0.3
    }; 
    const double fMaxDLKLDC[8] = {
        4.8, 4.8, 4.8, 4.8,
        4.8, 4.8, 4.8, 4.8
    };
} // namespace DCConstants

#endif // DCCONSTANTS_H