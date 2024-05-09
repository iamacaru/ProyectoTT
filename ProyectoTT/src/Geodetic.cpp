#include "../include/Geodetic.h"
#include "../include/SAT_Const.h"

#include <cmath>
#include <stdexcept>
#include <limits>

using namespace std;

Matrix Geodetic(const Matrix& r) {

    const double R_equ = Constants::R_Earth;
    const double f = Constants::f_Earth;

    const double epsRequ = numeric_limits<double>::epsilon() * R_equ;   // Convergence criterion
    const double e2 = f * (2.0 - f);                                    // Cuadrado de la excentricidad

    double X = r(1, 1);                                            // Cartesian coordinates
    double Y = r(1, 2);
    double Z = r(1, 3);
    double rho2 = X * X + Y * Y;                                        // Square of distance from z-axis

    // Check validity of input data
    if (r.norm() == 0.0) {
        printf("ERROR: Invalid input in Geodetic constructor");
        exit(EXIT_FAILURE);
    }

    // Iteration
    double dZ = e2 * Z;

    double ZdZ;
    double Nh;
    double SinPhi;
    double N;
    double dZ_new;

    while (true) {
        ZdZ = Z + dZ;
        Nh = sqrt(rho2 + ZdZ * ZdZ);
        SinPhi = ZdZ / Nh;                                              // Sine of geodetic latitude
        N = R_equ / sqrt(1.0 - e2 * SinPhi * SinPhi);
        dZ_new = N * e2 * SinPhi;

        if (fabs(dZ - dZ_new) < epsRequ) {
            break;
        }

        dZ = dZ_new;
    }

    // Longitude, latitude, altitude
    double lon = atan2(Y, X);
    double lat = atan2(ZdZ, sqrt(rho2));
    double h = Nh - N;

    double resultValues[] = {lon, lat, h};

    return {1, 3, resultValues, 3};
}
