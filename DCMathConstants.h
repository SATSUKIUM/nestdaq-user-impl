#ifndef DCMATHCONSTANTS_H
#define DCMATHCONSTANTS_H 1

namespace nestdaq::DCMath{
    static constexpr double Deg2Rad = acos(-1.)/180.; // angle [rad] = angle [deg] * Deg2Rad
    static constexpr double Rad2Deg = 180./acos(-1.); // angle [deg] = angle [rad] * Rad2Deg

    const double Infinity = std::numeric_limits<double>::infinity();
    const double TINY = std::numeric_limits<double>::epsilon();
} // namespace nestdaq::DCMath

#endif // DCMATHCONSTANTS_H