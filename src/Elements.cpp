#include "../include/Elements.h"
#include "../include/SAT_Const.h"

#include <cmath>

void elements(Matrix& y, double& p, double& a, double& e, double& i, double& Omega,double& omega, double& M) {
    Matrix r = y.subMatrix(1,3, 1);                                                  // Position
    Matrix v = y.subMatrix(4,6, 1);                                                  // Velocity

    Matrix h = r.cross(v);                                                                   // Areal velocity
    double magh = h.norm();
    p = magh * magh / (Constants::GM_Earth);
    double H = h.norm();

    Omega = atan2 (h(1,1), -h(1,2));                                                // Long. ascend. node
    Omega = fmod(Omega,Constants::pi2);
    i = atan2( sqrt(h(1,1) * h(1,1)+h(1,2) * h(1,2)), h(1,3));    // Inclination
    double u = atan2(r(1,3)*H, -r(1,1) * h(1,2) + r(1,2) * h(1,1));  // Arg. of latitude

    double R  = r.norm();                                                                           // Distance

    a = 1 / (2 / R - v.dot(v) / (Constants::GM_Earth));                                      // Semi-major axis

    double eCosE = 1 - R / a;                                                                       // e*cos(E)
    double eSinE = r.dot(v) / sqrt((Constants::GM_Earth) * a);                            // e*sin(E)

    double e2 = eCosE * eCosE + eSinE * eSinE;
    e  = sqrt(e2);                                                                               // Eccentricity
    double E  = atan2(eSinE,eCosE);                                                           // Eccentric anomaly

    M  = fmod(E-eSinE,(Constants::pi2));                                                      // Mean anomaly

    double nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);                                       // True anomaly

    omega = fmod(u-nu, (Constants::pi2));                                                     // Arg. of perihelion
}