#include "../include/Accel.h"
#include "../include/Global.h"
#include "../include/IERS.h"
#include "../include/Timediff.h"
#include "../include/SAT_Const.h"
#include "../include/PrecMatrix.h"
#include "../include/NutMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/GHAMatrix.h"
#include "../include/Mjday_TDB.h"
#include "../include/JPL_Eph_DE430.h"
#include "../include/AccelHarmonic.h"
#include "../include/AccelPointMass.h"

Matrix Accel(double x, Matrix& Y) {

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_UT1, Mjd_TT, MJD_TDB;
    char interp = 'l';

    IERS(*Global::eopdata, Global::AuxParam::Mjd_UTC + x / 86400, interp, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
    Matrix timediff1 = timediff(UT1_UTC, TAI_UTC);
    UT1_TAI = timediff1(1, 1);
    UTC_GPS = timediff1(1, 2);
    UT1_GPS = timediff1(1, 3);
    TT_UTC = timediff1(1, 4);
    GPS_UTC = timediff1(1, 5);
    Mjd_UT1 = Global::AuxParam::Mjd_UTC + x/86400 + UT1_UTC / 86400;
    Mjd_TT = Global::AuxParam::Mjd_UTC + x/86400 + TT_UTC/86400;

    Matrix P = PrecMatrix(Constants::MJD_J2000, Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole, y_pole) * GHAMatrix(Mjd_UT1) * T;

    MJD_TDB = Mjday_TDB(Mjd_TT);

    Matrix r_Mercury(3,3), r_Venus(3,3), r_Earth(3,3), r_Mars(3,3), r_Jupiter(3,3),
            r_Saturn(3,3), r_Uranus(3,3), r_Neptune(3,3), r_Pluto(3,3), r_Moon(3,3),
            r_Sun(3,3);

    JPL_Eph_DE430(r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune,
                  r_Pluto, r_Moon, r_Sun, MJD_TDB);

    // Acceleration due to harmonic gravity field
    double values1[] = {Y(1, 1), Y(2, 1), Y(3, 1)};
    Matrix aux(3, 1, values1, 3);
    Matrix a = AccelHarmonic(aux, E, Global::AuxParam::n, Global::AuxParam::m);

    // Luni-solar perturbations
    if (Global::AuxParam::sun) {
        a = a + AccelPointMass(aux, r_Sun, Constants::GM_Sun);
    }

    if (Global::AuxParam::moon) {
        a = a + AccelPointMass(aux, r_Moon, Constants::GM_Moon);
    }

    // Planetary perturbations
    if (Global::AuxParam::planets) {
        a = a + AccelPointMass(aux, r_Mercury, Constants::GM_Mercury);
        a = a + AccelPointMass(aux, r_Venus, Constants::GM_Venus);
        a = a + AccelPointMass(aux, r_Mars, Constants::GM_Mars);
        a = a + AccelPointMass(aux, r_Jupiter, Constants::GM_Jupiter);
        a = a + AccelPointMass(aux, r_Saturn, Constants::GM_Saturn);
        a = a + AccelPointMass(aux, r_Uranus, Constants::GM_Uranus);
        a = a + AccelPointMass(aux, r_Neptune, Constants::GM_Neptune);
        a = a + AccelPointMass(aux, r_Pluto, Constants::GM_Pluto);
    }

    double values2[] = {Y(4, 1), Y(5, 1), Y(6, 1), a(1, 1), a(2, 1), a(3, 1)};
    Matrix dY(6, 1, values2, 6);

    return dY;
}