#include "../include/Cheb3D.h"

#include <stdexcept>

Matrix Cheb3D(double t, int N, double Ta, double Tb, Matrix&  Cx, Matrix& Cy, Matrix& Cz) {
    // Check validity
    if ((t < Ta) || (t > Tb)) {
        printf("ERROR: Time out of range in Cheb3D");
        exit(EXIT_FAILURE);
    }

    // Clenshaw algorithm
    double tau = (2 * t - Ta - Tb) / (Tb - Ta);

    Matrix f1(1, 3);
    Matrix f2(1, 3);

    Matrix old_f1(1, 3);
    double values[] = {0,0,0};
    for (int i = N - 1; i >= 1; i--) {
        old_f1 = f1;
        values[0] = Cx(1,i);
        values[1] = Cy(1,i);
        values[2] = Cz(1,i);
        f1 = 2 * tau * f1 - f2 + Matrix (1,3, values, 3);
        f2 = old_f1;
    }

    values[0] = Cx(1,1);
    values[1] = Cy(1,1);
    values[2] = Cz(1,1);
    Matrix ChebApp = tau * f1 - f2 + Matrix (1,3, values, 3);

    return ChebApp;
}
