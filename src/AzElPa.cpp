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

    double aux = 0.0;

    for (int j = 1; j <= s.columns(); j++) {
        aux += s(1, j) * s(1, j);
    }

    dEds(1, 1) = -s(1, 1)*s(1, 3)/rho/aux;
    dEds(1, 2) = -s(1, 2)*s(1, 3)/rho/aux;
    dEds(1, 3) = rho/aux;
}