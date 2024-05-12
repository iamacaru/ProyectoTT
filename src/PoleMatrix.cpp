#include "../include/PoleMatrix.h"
#include "../include/R_x.h"
#include "../include/R_y.h"

Matrix PoleMatrix(double xp, double yp) {
    return R_y(-xp) * R_x(-yp);
}