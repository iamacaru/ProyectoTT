#include "../include/EqnEquinox.h"
#include "../include/MeanObliquity.h"
#include "../include/NutAngles.h"

#include <cmath>

double EqnEquinox(double Mjd_TT) {
    // Nutation in longitude and obliquity
    Matrix nutAngles = NutAngles(Mjd_TT);
    double dpsi = nutAngles(1, 1);

    // Equation of the equinoxes
    return dpsi * cos(MeanObliquity(Mjd_TT));
}
