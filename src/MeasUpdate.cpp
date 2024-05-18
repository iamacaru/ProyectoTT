#include "MeasUpdate.h"

void MeasUpdate(Matrix& x, Matrix& z, Matrix& g, Matrix& s, Matrix& G, Matrix& P, int n, Matrix& K) {
    int m = z.rows() * z.columns();
    Matrix Inv_W(m, m);

    for (int i = 1; i <= m; ++i) {
        Inv_W(i, i) = s(1, i) * s(1, i);     // Inverse weight (measurement covariance)
    }

    // Kalman gain
    K = P * G.traspuesta() * (Inv_W + G * P * G.traspuesta()).inversa()(1, 1);


    // State update
    x = x + K * (z - g)(1, 1);

    // Covariance update
    P = (Matrix::identidad(n) - K * G) * P;
}