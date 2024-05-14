#include "../include/JPL_Eph_DE430.h"
#include "../include/Global.h"
#include "../include/Cheb3D.h"

void JPL_Eph_DE430(Matrix& r_Mercury, Matrix& r_Venus, Matrix& r_Earth, Matrix& r_Mars, Matrix& r_Jupiter,
                     Matrix& r_Saturn, Matrix& r_Uranus, Matrix& r_Neptune, Matrix& r_Pluto, Matrix& r_Moon,
                     Matrix& r_Sun, double Mjd_TDB) {
    extern double Mjd0;

    double JD = Mjd_TDB + 2400000.5;
    int i = 1;
    while (i <= Global::PC->rows()) {
        if ((*Global::PC)(i,1) <= JD && JD <= (*Global::PC)(i,2)) {
            break;
        }
        i++;
    }
    Matrix PCtemp = (*Global::PC).subMatrix(1, (*Global::PC).columns(), i);

    double t1 = PCtemp(1,1)-2400000.5;                                                                   // MJD at start of interval

    double dt = Mjd_TDB - t1;

    int numElementos = (270 - 231) / 13 + 1;
    Matrix temp(1, numElementos);
    i = 1;
    for (int k = 231; k <= 270; k += 13) {
        temp(1, i) = k;
        i++;
    }

    Matrix Cx_Earth = PCtemp.subMatrix((int)temp(1,1), (int)temp(1,2) - 1, 1);
    Matrix Cy_Earth = PCtemp.subMatrix((int)temp(1,2), (int)temp(1,3) - 1, 1);
    Matrix Cz_Earth = PCtemp.subMatrix((int)temp(1,3), (int)temp(1,4) - 1, 1);
    temp = temp + 39;
    Matrix Cx = PCtemp.subMatrix((int)temp(1,1), (int)temp(1,2) - 1, 1);
    Matrix Cy = PCtemp.subMatrix((int)temp(1,2), (int)temp(1,3) - 1, 1);
    Matrix Cz = PCtemp.subMatrix((int)temp(1,3), (int)temp(1,4) - 1, 1);
    Cx_Earth = Cx_Earth.concatenar(Cx);
    Cy_Earth = Cy_Earth.concatenar(Cy);
    Cz_Earth = Cz_Earth.concatenar(Cz);
    int j = 0;
    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    }
    else if (16 < dt && dt <= 32) {
        j=1;
        Mjd0 = t1+16*j;
    }
    r_Earth = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, Cx_Earth.subMatrix(13*j+1, 13*j+13, 1), Cy_Earth.subMatrix(13*j+1, 13*j+13, 1), Cz_Earth.subMatrix(13*j+1, 13*j+13, 1));







}