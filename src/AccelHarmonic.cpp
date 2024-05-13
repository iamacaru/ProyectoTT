#include "../include/AccelHarmonic.h"
#include "../include/Legendre.h"
#include "../include/Global.h"

#include <cmath>

Matrix AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max) {
    double r_ref = 6378.1363e3;   // Earth's radius [m]; GGM03S
    double gm = 398600.4415e9;    // [m^3/s^2]; GGM03S

    // Body-fixed position
    Matrix r_bf = E * r;

    // Auxiliary quantities
    double d = r_bf.norm();       // distance
    double latgc = asin(r_bf(3, 1) / d);
    double lon = atan2(r_bf(2, 1), r_bf(1, 1));

    Matrix pnm(n_max+1,m_max+1);
    Matrix dpnm(n_max+1,m_max+1);
    Legendre(n_max,m_max,latgc,pnm,dpnm);

    double dUdr = 0.0;
    double dUdlatgc = 0.0;
    double dUdlon = 0.0;
    double q1 = 0.0, q2 = 0.0, q3 = 0.0;

    for (int n = 0; n <= n_max; ++n) {
        double b1 = (-gm / (d * d)) * pow(r_ref / d, n) * (n + 1);
        double b2 = (gm / d) * pow(r_ref / d, n);
        double b3 = (gm / d) * pow(r_ref / d, n);

        for (int m = 0; m <= m_max; ++m) {
            q1 += pnm(n + 1, m + 1) * ((*Global::Cnm)(n + 1, m + 1) * cos(m * lon) + (*Global::Snm)(n + 1, m + 1) * sin(m * lon));
            q2 += dpnm(n + 1, m + 1) * ((*Global::Cnm)(n + 1, m + 1) * cos(m * lon) + (*Global::Snm)(n + 1, m + 1) * sin(m * lon));
            q3 += m * pnm(n + 1, m + 1) * ((*Global::Snm)(n + 1, m + 1) * cos(m * lon) - (*Global::Cnm)(n + 1, m + 1) * sin(m * lon));
        }

        dUdr += q1 * b1;
        dUdlatgc += q2 * b2;
        dUdlon += q3 * b3;

        q1 = 0.0; q2 = 0.0; q3 = 0.0;
    }

    // Body-fixed acceleration
    double r2xy = r_bf(1, 1) * r_bf(1, 1) + r_bf(2, 1) * r_bf(2, 1);

    double ax = (1.0 / d * dUdr - r_bf(3, 1) / (d * d * sqrt(r2xy)) * dUdlatgc) * r_bf(1, 1) - (1.0 / r2xy * dUdlon) * r_bf(2, 1);
    double ay = (1.0 / d * dUdr - r_bf(3, 1) / (d * d * sqrt(r2xy)) * dUdlatgc) * r_bf(2, 1) + (1.0 / r2xy * dUdlon) * r_bf(1, 1);
    double az = 1.0 / d * dUdr * r_bf(3, 1) + sqrt(r2xy) / (d * d) * dUdlatgc;

    Matrix a_bf(3, 1);
    a_bf(1, 1) = ax;
    a_bf(2, 1) = ay;
    a_bf(3, 1) = az;

    // Inertial acceleration
    Matrix a = E.traspuesta() * a_bf;

    Matrix result(3,1);
    result(1,1) = a(1,1);
    result(2,1) = a(2,1);
    result(3,1) = a(3,1);

    return result;
}