#include "../include/G_AccelHarmonic.h"
#include "../include/AccelHarmonic.h"

Matrix G_AccelHarmonic(Matrix& r, Matrix& U, int n_max, int m_max) {
    double d = 1.0;   // Position increment [m]

    Matrix G(3, 3);
    Matrix dr(3, 1);

    Matrix aux1(r.rows(), r.columns());
    Matrix aux2(r.rows(), r.columns());

    // Gradient
    for (int i = 1; i <= 3; ++i) {
        // Set offset in i-th component of the position vector
        dr(1,1) = 0; dr(2,1) = 0; dr(3,1) = 0;
        dr(i, 1) = d;
        // Acceleration difference
        aux1 = r+dr/2;
        aux2 = r-dr/2;
        Matrix da = AccelHarmonic(aux1, U, n_max, m_max) - AccelHarmonic(aux2, U, n_max, m_max);
        // Derivative with respect to i-th axis
        for (int j = 1; j <= 3; ++j) {
            G(j, i) = da(j, i) / d;
        }
    }

    return G;
}
