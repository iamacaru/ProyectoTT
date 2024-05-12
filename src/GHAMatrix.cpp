#include "../include/GHAMatrix.h"
#include "../include/Gast.h"
#include "../include/R_z.h"

Matrix GHAMatrix(double Mjd_UT1) {
    return R_z(gast(Mjd_UT1));
}
