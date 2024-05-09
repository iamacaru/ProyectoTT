#include "../include/AzElPa.h"
#include "../include/SAT_Const.h"

#include <cmath>

void  AzElPa(const Matrix& s, double& Az, double& El, Matrix& dAds, Matrix& dEds) {
    double rho = sqrt(s(1,1) * s(1,1) + s(1,2) * s(1,2));

    // Angles
    Az = atan2(s(1,1), s(1,2));

    if (Az < 0.0) {
        Az += Constants::pi2;
    }

    El = atan(s(1,3) / rho);

    // Partials
    dAds(1, 1) = s(1, 2)/(rho*rho);
    dAds(1, 2) = -s(1, 1)/(rho*rho);
    dAds(1, 3) = 0.0;

    dEds(1, 1) = -s(1, 1)*s(1, 3)/rho/s.dot(s);
    dEds(1, 2) = -s(1, 2)*s(1, 3)/rho/s.dot(s);
    dEds(1, 3) = rho/s.dot(s);
}