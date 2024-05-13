#include "../include/VarEqn.h"
#include "../include/Global.h"
#include "../include/IERS.h"
#include "../include/Timediff.h"
#include "../include/SAT_Const.h"
#include "../include/PrecMatrix.h"
#include "../include/NutMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/GHAMatrix.h"

/*--------------------------------------------------------------------------

  Initial Orbit Determination using Gauss and Extended Kalman Filter methods

 References:
   O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
   Applications", Springer Verlag, Heidelberg, 2000.

   D. Vallado, "Fundamentals of Astrodynamics and Applications",
   4th Edition, 2013.

   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003.

 Last modified:   2020/03/16   Meysam Mahooti
--------------------------------------------------------------------------*/

Matrix VarEqn(double x, Matrix& yPhi) {

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_UT1;
    char interp = 'l';

    IERS(*Global::eopdata, Global::AuxParam::Mjd_UTC, interp, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
    Matrix timediff1 = timediff(UT1_UTC, TAI_UTC);
    UT1_TAI = timediff1(1, 1);
    UTC_GPS = timediff1(1, 2);
    UT1_GPS = timediff1(1, 3);
    TT_UTC = timediff1(1, 4);
    GPS_UTC = timediff1(1, 5);
    Mjd_UT1 = Global::AuxParam::Mjd_TT + (UT1_UTC - TT_UTC) / 86400;

    // Transformation matrix
    Matrix P = PrecMatrix(Constants::MJD_J2000, Global::AuxParam::Mjd_TT + x / 86400);
    Matrix N = NutMatrix(Global::AuxParam::Mjd_TT + x/86400);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole, y_pole) * GHAMatrix(Mjd_UT1) * T;

    // State vector components
    double valuesR[] = {yPhi(1,1), yPhi(1,2), yPhi(1,3)};
    Matrix r(1,3, valuesR, 3);
    double valuesV[] = {yPhi(1,4), yPhi(1,5), yPhi(1,6)};
    Matrix v(1,3, valuesV, 3);
    Matrix Phi(6, 6);

    double vectorDouble[1000];

    // Asignar valores al array (por ejemplo, del 1 al 1000)
    for (int i = 0; i < 1000; ++i) {
        vectorDouble[i] = i + 1;
    }
    Matrix prueba{1000,1,vectorDouble,1000};

    // State transition matrix
    for (int j = 1; j <= 6; ++j) {
        for (int i = 1; i <= 6; ++i) {
            Phi(i,j) = prueba(6*j+i,1);
        }
    }







    return {3,3};
}