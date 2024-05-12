#include "../include/NutMatrix.h"
#include "../include/MeanObliquity.h"
#include "../include/NutAngles.h"
#include "../include/R_x.h"
#include "../include/R_z.h"

Matrix NutMatrix(double Mjd_TT) {
    // Mean obliquity of the ecliptic
    double eps = MeanObliquity(Mjd_TT);

    // Nutation in longitude and obliquity
    Matrix angles = NutAngles(Mjd_TT);
    double dpsi = angles(1, 1);
    double deps = angles(1, 2);

    // Transformation from mean to true equator and equinox
    return R_x(-eps-deps) * R_z(-dpsi) * R_x(+eps);
}