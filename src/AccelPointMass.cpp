#include "../include/AccelPointMass.h"

#include <cmath>

Matrix AccelPointMass(Matrix& r, Matrix& s, double GM) {

    // Relative position vector of satellite w.r.t. point mass
    Matrix d = r - s;

    double norm_d = d.norm();
    double norm_s = s.norm();

    // Acceleration
    Matrix a = -GM * ( d/pow(norm_d, 3) + s/pow(norm_s, 3));

    return a;
}
