#include "../include/Unit.h"

Matrix unit(const Matrix& vec) {

    double magv = vec.norm();

    Matrix outvec(vec.rows(), vec.columns());

    const double SMALL = 1e-6;

    if (magv > SMALL) {
        for (int j = 1; j < 4; j++) {
            outvec(1, j) = vec(1, j) / magv;
        }
    } else {
        for (int j = 1; j < 4; j++) {
            outvec(1, j) = 0.0;
        }
    }

    return outvec;
}
