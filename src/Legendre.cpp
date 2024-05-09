#include "../include/Legendre.h"

#include <cmath>

void Legendre(int n, int m, double fi, Matrix& pnm, Matrix& dpnm) {
    pnm = Matrix(n+1,m+1);
    dpnm = Matrix(n+1,m+1);

    pnm(1,1) = 1.0;
    dpnm(1,1) = 0.0;
    pnm(2,2) = sqrt(3.0) * cos(fi);
    dpnm(2,2) = -sqrt(3.0) * sin(fi);

    for (int i = 2; i <= n; i++) {
        pnm(i+1, i+1) = sqrt((2 * i + 1) / (2.0 * i)) * cos(fi) * pnm(i,i);
        dpnm(i+1, i+1) = sqrt((2 * i + 1) / (2.0 * i)) * (cos(fi) * dpnm(i,i) - sin(fi) * pnm(i,i));
    }

    for (int i = 1; i <= n; i++) {
        pnm(i+1,i) = sqrt(2 * i + 1) * sin(fi) * pnm(i,i);
        dpnm(i+1,i) = sqrt(2 * i + 1) * (cos(fi) * pnm(i,i) + sin(fi) * dpnm(i,i));
    }

    int j = 0;
    int k = 2;
    while (true) {
        for (int i = k; i <= n; i++) {
            pnm(i+1,j+1) = sqrt((2.0 * i + 1) / ((i - j) * (i + j))) * ((sqrt(2.0 * i - 1.0) * sin(fi) * pnm(i,j + 1)) -
            (sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0)) * pnm(i - 1,j + 1)));
        }
        j++;
        k++;
        if (j > m) {
            break;
        }
    }

    j = 0;
    k = 2;
    while (true) {
        for (int i = k; i <= n; i++) {
            dpnm(i + 1,j + 1) = sqrt((2.0 * i + 1.0) / ((i - j) * (i + j))) * (sqrt(2 * i - 1) * sin(fi) * dpnm(i,j + 1) +
            sqrt(2 * i - 1) * cos(fi) * pnm(i,j + 1) - (sqrt(((i + j - 1.0) * (i - j - 1)) / (2.0 * i - 3.0)) * dpnm(i - 1,j + 1)));
        }
        j++;
        k++;
        if (j > m) {
            break;
        }
    }
}