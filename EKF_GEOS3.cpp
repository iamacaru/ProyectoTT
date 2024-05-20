#include "./include/Matrix.h"
#include "./include/R_z.h"
#include "./include/Global.h"
#include "./include/Timediff.h"
#include "./include/AzElPa.h"
#include "./include/Mjday.h"
#include "./include/Position.h"
#include "./include/IERS.h"
#include "./include/SAT_Const.h"
#include "./include/Gmst.h"
#include "./include/TimeUpdate.h"
#include "./include/Accel.h"
#include "./include/LTC.h"
#include "./include/MeasUpdate.h"
#include "./include/DEInteg.h"
#include "./include/VarEqn.h"

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

/*--------------------------------------------------------------------------

  Initial Orbit Determination using Gauss and Extended Kalman Filter methods

 References:
   O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
   Applications", Springer Verlag, Heidelberg, 2000.

   D. Vallado, "Fundamentals of Astrodynamics and Applications",
   4th Edition, 2013.

   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003.

--------------------------------------------------------------------------*/

int main() {
    Global::Pc();
    Global::GGM03S();
    Global::eop19620101(21413);
    extern double Mjd0;

    int nobs = 46;
    Matrix obs(nobs, 4);

    #ifdef _WIN32
        // CLion en Windows
        ifstream file("../data/GEOS3.txt");
    #else
        // Linux en bash
        ifstream file("./data/GEOS3.txt");
    #endif

    if (!file.is_open()) {
        std::cerr << "Error al abrir el fichero." << std::endl;
        return 1;
    }

    string line;
    int Y, M, D, h, m, s;
    double az, el, Dist;

    for (int i = 1; i <= nobs; ++i) {
        if (!getline(file, line)) {
            break;
        }

        Y = stoi(line.substr(0, 4));
        M = stoi(line.substr(5, 2));
        D = stoi(line.substr(8, 2));
        h = stoi(line.substr(12, 2));
        m = stoi(line.substr(15, 2));
        s = stoi(line.substr(18, 6));
        az = stod(line.substr(25, 8));
        el = stod(line.substr(35, 7));
        Dist = stod(line.substr(44, 10));

        obs(i,1) = Mjday(Y, M, D, h, m, s);
        obs(i,2) = Constants::Rad * az;
        obs(i,3) = Constants::Rad * el;
        obs(i,4) = 1e3 * Dist;
    }

    file.close();

    double sigma_range = 92.5;                   // [m]
    double sigma_az = 0.0224 * Constants::Rad;   // [rad]
    double sigma_el = 0.0139 * Constants::Rad;   // [rad]

    // Kaena Point station
    double lat = Constants::Rad * 21.5748;       // [rad]
    double lon = Constants::Rad * (-158.2706);   // [rad]
    double alt = 300.20;                         // [m]

    Matrix Rs = Position(lon, lat, alt).traspuesta();

    double Mjd1 = obs(1, 1);
    double Mjd2 = obs(9, 1);
    double Mjd3 = obs(18, 1);

    double valuesR2[] = {6.221397628578685e+06, 2.867713779657379e+06, 3.006155985099489e+06};
    Matrix r2(1, 3, valuesR2, 3);
    double valuesV2[] = {4.645047251618071e+03, -2.752215915882051e+03, -7.507999409870331e+03};
    Matrix v2(1, 3, valuesV2, 3);
    // [r2,v2] = anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3), Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);
    // [r2,v2] = anglesdr(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3), Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);

    Matrix Y0_apr = r2.concatenar(v2).traspuesta();

    Mjd0 = Mjday(1995, 1, 29, 02, 38, 0);

    double Mjd_UTC = obs(9,1);

    Global::AuxParam::Mjd_UTC = Mjd_UTC;
    Global::AuxParam::n = 20;
    Global::AuxParam::m = 20;
    Global::AuxParam::sun = 1;
    Global::AuxParam::moon = 1;
    Global::AuxParam::planets = 1;
    double auxMjd0 = Mjd0;

    int n_eqn  = 6;

    Matrix y = DEInteg(Accel, 0, -(obs(9, 1) - Mjd0) * 86400.0, 1e-13, 1e-6, 6, Y0_apr);

    Matrix P(6, 6);

    for (int i = 1; i <= 3; i++) {
        P(i, i) = 1e8;
    }
    for (int i = 4; i <= 6; i++) {
        P(i, i) = 1e3;
    }

    Matrix LT = LTC(lon,lat);

    Matrix yPhi(42, 1);
    Matrix Phi(6, 6);

    // Measurement loop
    double t = 0.0, t_old;
    Matrix Y_old(1, 1), U(1, 1), r(1, 1), S(1, 1), dAds(1,3), dEds(1,3),
    dAdY(1, 1), K(1, 1), aux1(1, 1), aux2(1, 1), aux3(1, 1), dEdY(1, 1),
    dDds(1, 1), dDdY(1, 1), timediff1(1, 1);

    char interp = 'l';
    double x_pole = 0.0, y_pole = 0.0, UT1_UTC = 0.0, LOD = 0.0, dpsi = 0.0, deps = 0.0, dx_pole = 0.0, dy_pole = 0.0,
    TAI_UTC = 0.0, UT1_TAI = 0.0, UTC_GPS = 0.0, UT1_GPS = 0.0, TT_UTC = 0.0, GPS_UTC = 0.0, Mjd_TT = 0.0, Mjd_UT1 = 0.0,
    theta = 0.0, Azim = 0.0, Elev = 0.0;

    for (int i = 1; i <= nobs; i++) {
        // Previous step
        t_old = t;
        Y_old = y;

        // Time increment and propagation
        Mjd_UTC = obs(i, 1);                      // Modified Julian Date
        t = (Mjd_UTC - auxMjd0) * 86400.0;              // Time since epoch [s]

        IERS(*Global::eopdata, Mjd_UTC, interp, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
        timediff1 = timediff(UT1_UTC, TAI_UTC);
        UT1_TAI = timediff1(1, 1); UTC_GPS = timediff1(1, 2); UT1_GPS = timediff1(1, 3); TT_UTC = timediff1(1, 4); GPS_UTC = timediff1(1, 5);
        Mjd_TT = Mjd_UTC + TT_UTC / 86400;
        Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC) / 86400.0;
        Global::AuxParam::Mjd_UTC = Mjd_UTC;
        Global::AuxParam::Mjd_TT = Mjd_TT;

        for (int ii = 1; ii <= 6; ii++) {
            yPhi(ii, 1) = Y_old(ii, 1);
            for (int j = 1; j <= 6; j++) {
                if (ii == j) {
                    yPhi(6 * j + ii, 1) = 1;
                } else {
                    yPhi(6 * j + ii, 1) = 0;
                }
            }
        }

        yPhi = DEInteg (VarEqn, 0, t-t_old, 1e-13, 1e-6, 42, yPhi);

        // Extract state transition matrices
        for (int j = 1; j <= 6; j++) {
            for (int k = 1; k <= 6; k++) {
                Phi(k, j) = yPhi(6 * j + k,1);
            }
        }

        y = DEInteg (Accel, 0, t-t_old, 1e-13, 1e-6,6, Y_old);

        // Topocentric coordinates
        theta = gmst(Mjd_UT1);                    // Earth rotation
        U = R_z(theta);
        r = y.traspuesta().subMatrix(1, 3, 1).traspuesta();
        S = LT * (U * r - Rs);                          // Topocentric position [m]

        // Time update
        P = TimeUpdate(P, Phi);


        // Azimuth and partials
        AzElPa(S.traspuesta(), Azim, Elev, dAds, dEds);    // Azimuth, Elevation
        dAdY = (dAds * LT * U).concatenar(Matrix (1, 3));

        // Measurement update
        aux1(1, 1) = obs(i, 2);
        aux2(1, 1) = Azim;
        aux3(1, 1) = sigma_az;
        MeasUpdate(y, aux1, aux2, aux3, dAdY, P, 6, K);

        // Elevation and partials
        r = y.traspuesta().subMatrix(1, 3, 1).traspuesta();
        S = LT * (U * r - Rs);                          // Topocentric position [m]
        AzElPa(S.traspuesta(), Azim, Elev, dAds, dEds);     // Azimuth, Elevation
        dEdY = (dEds * LT * U).concatenar(Matrix (1 , 3));

        // Measurement update
        aux1(1, 1) = obs(i, 3);
        aux2(1, 1) = Elev;
        aux3(1, 1) = sigma_el;
        MeasUpdate(y, aux1, aux2, aux3, dEdY, P, 6, K);

        // Range and partials
        r = y.traspuesta().subMatrix(1, 3, 1).traspuesta();
        S = LT * (U * r - Rs);                          // Topocentric position [m]
        Dist = S.norm();
        dDds = (S / Dist).traspuesta();         // Range
        dDdY = (dDds * LT * U).concatenar(Matrix (1, 3));

        // Measurement update
        aux1(1, 1) = obs(i, 4);
        aux2(1, 1) = Dist;
        aux3(1, 1) = sigma_range;
        MeasUpdate(y, aux1, aux2, aux3, dDdY, P, 6, K);
    }

    IERS(*Global::eopdata, obs(46,1), interp, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
    timediff1 = timediff(UT1_UTC, TAI_UTC);
    UT1_TAI = timediff1(1, 1); UTC_GPS = timediff1(1, 2); UT1_GPS = timediff1(1, 3); TT_UTC = timediff1(1, 4); GPS_UTC = timediff1(1, 5);
    Mjd_TT = Mjd_UTC + TT_UTC / 86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC) / 86400.0;
    Global::AuxParam::Mjd_UTC = Mjd_UTC;
    Global::AuxParam::Mjd_TT = Mjd_TT;

    Matrix Y0 = DEInteg(Accel, 0, -(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,y);

    double valuesY_true[] = {5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3};
    Matrix Y_true(6, 1, valuesY_true, 6);

    cout << fixed << setprecision(1);

    cout << "Error of Position Estimation" << endl;
    cout << "dX" << setw(10) << Y0(1, 1) - Y_true(1, 1) << " [m]" << endl;
    cout << "dY" << setw(10) << Y0(2, 1) - Y_true(2, 1) << " [m]" << endl;
    cout << "dZ" << setw(10) << Y0(3, 1) - Y_true(3, 1) << " [m]" << endl;

    cout << endl;

    cout << "Error of Velocity Estimation" << endl;
    cout << "dVx" << setw(8) << Y0(4, 1) - Y_true(4, 1) << " [m/s]" << endl;
    cout << "dVy" << setw(8) << Y0(5, 1) - Y_true(5, 1) << " [m/s]" << endl;
    cout << "dVz" << setw(8) << Y0(6, 1) - Y_true(6, 1) << " [m/s]" << endl;

    return 0;
};
