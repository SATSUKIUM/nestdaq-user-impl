#ifndef DCMATH_H
#define DCMATH_H 1

namespace DCMath{
    static constexpr double Deg2Rad = acos(-1.)/180.; // angle [rad] = angle [deg] * Deg2Rad
    static constexpr double Rad2Deg = 180./acos(-1.); // angle [deg] = angle [rad] * Rad2Deg
}

#endif // DCMATH_H