#include "../include/Position.h"
#include "../include/SAT_Const.h"

#include <cmath>

Matrix Position(double lon, double lat, double alt) {
    double R_equ = Constants::R_Earth;
    double f = Constants::f_Earth;

    double e2 = f * (2.0 - f);     // Square of eccentricity
    double CosLat = cos(lat);   // (Co)sine of geodetic latitude
    double SinLat = sin(lat);

    // Position vector
    double N = R_equ / std::sqrt(1.0 - e2 * SinLat * SinLat);

    Matrix r(1,3);

    r(1,1) = (N + alt) * CosLat * cos(lon);
    r(1,2) = (N + alt) * CosLat * sin(lon);
    r(1,3) = ((1.0 - e2) * N + alt) * SinLat;

    return r;
}
