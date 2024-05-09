#include "../include/EccAnom.h"

#include <cmath>
#include <stdexcept>
#include <limits>

using namespace std;

double EccAnom(double M, double e) {
    const int maxit = 15;
    int i = 0;

    // Starting value
    M = fmod(M, 2.0 * M_PI);

    double E;
    if (e < 0.8) {
        E = M;
    } else {
        E = M_PI;
    }

    double f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));


    // Iteration
    while (fabs(f) > 10e2 * numeric_limits<double>::epsilon()) {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i++;

        if (i == maxit) {
            printf("ERROR: Convergence problems in EccAnom");
            exit(EXIT_FAILURE);
        }
    }

    return E;
}
