#include "../include/Gast.h"
#include "../include/EqnEquinox.h"
#include "../include/Gmst.h"

#include <cmath>

double gast(double Mjd_UT1) {
    double gstime = fmod(gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2 * M_PI);
    return gstime;
}