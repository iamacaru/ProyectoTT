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
    Matrix Cx_Earth = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
    Matrix Cy_Earth = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
    Matrix Cz_Earth = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1, 4) - 1, 1);
    temp = temp + 39;
    Matrix Cx = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
    Matrix Cy = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
    Matrix Cz = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1, 4) - 1, 1);
    Cx_Earth = Cx_Earth.concatenar(Cx);
    Cy_Earth = Cy_Earth.concatenar(Cy);
    Cz_Earth = Cz_Earth.concatenar(Cz);
    int j = 0;
    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    }
    else if (16 < dt && dt <= 32) {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }
    Matrix aux1 = Cx_Earth.subMatrix(13 * j + 1, 13 * j + 13, 1);
    Matrix aux2 = Cy_Earth.subMatrix(13 * j + 1, 13 * j + 13, 1);
    Matrix aux3 = Cz_Earth.subMatrix(13 * j + 1, 13 * j + 13, 1);
    r_Earth = 1e3 * Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, aux1, aux2, aux3).traspuesta();

    i = 1;
    for (int k = 441; k <= 480; k += 13) {
        temp(1, i) = k;
        i++;
    }
    Matrix Cx_Moon = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
    Matrix Cy_Moon = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
    Matrix Cz_Moon = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1, 4) - 1, 1);
    for (i = 1; i <= 7; i++) {
        temp = temp + 39;
        Cx = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
        Cy = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
        Cz = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1, 4) - 1, 1);
        Cx_Moon = Cx_Moon.concatenar(Cx);
        Cy_Moon = Cy_Moon.concatenar(Cy);
        Cz_Moon = Cz_Moon.concatenar(Cz);
    }
    if (0 <= dt && dt<=4) {
        j = 0;
        Mjd0 = t1;
    } else if (4 < dt && dt<=8) {
        j = 1;
        Mjd0 = t1 + 4 * j;
    } else if (8 < dt && dt <= 12) {
        j = 2;
        Mjd0 = t1 + 4 * j;
    } else if (12 < dt && dt <= 16) {
        j = 3;
        Mjd0 = t1 + 4 * j;
    } else if (16 < dt && dt <= 20) {
        j = 4;
        Mjd0 = t1 + 4 * j;
    } else if (20 < dt && dt <= 24) {
        j = 5;
        Mjd0 = t1 + 4 * j;
    } else if (24 < dt && dt <= 28) {
        j = 6;
        Mjd0 = t1 + 4 * j;
    } else if (28 < dt && dt <= 32) {
        j = 7;
        Mjd0 = t1 + 4 * j;
    }
    aux1 = Cx_Moon.subMatrix(13 * j + 1, 13 * j + 13, 1);
    aux2 = Cy_Moon.subMatrix(13 * j + 1, 13 * j + 13, 1);
    aux3 = Cz_Moon.subMatrix(13 * j + 1, 13 * j + 13, 1);
    r_Moon = 1e3 * Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4, aux1, aux2, aux3).traspuesta();

    i = 1;
    for (int k = 753; k <= 786; k += 11) {
        temp(1, i) = k;
        i++;
    }
    Matrix Cx_Sun = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
    Matrix Cy_Sun = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
    Matrix Cz_Sun = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1, 4) - 1, 1);
    temp = temp+33;
    Cx = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
    Cy = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
    Cz = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1, 4) - 1, 1);
    Cx_Sun = Cx_Sun.concatenar(Cx);
    Cy_Sun = Cy_Sun.concatenar(Cy);
    Cz_Sun = Cz_Sun.concatenar(Cz);
    if (0 <= dt && dt <= 16){
        j = 0;
        Mjd0 = t1;
    } else if (16 < dt && dt <= 32) {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }
    aux1 = Cx_Sun.subMatrix(11 * j + 1, 11 * j + 11, 1);
    aux2 = Cy_Sun.subMatrix(11 * j + 1, 11 * j + 11, 1);
    aux3 = Cz_Sun.subMatrix(11 * j + 1, 11 * j + 11, 1);
    r_Sun = 1e3 * Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+16, aux1, aux2, aux3).traspuesta();

    i = 1;
    for (int k = 3; k <= 45; k += 14) {
        temp(1, i) = k;
        i++;
    }
    Matrix Cx_Mercury = PCtemp.subMatrix((int) temp(1, 1), (int)temp(1, 2) - 1, 1);
    Matrix Cy_Mercury = PCtemp.subMatrix((int) temp(1, 2), (int)temp(1, 3) - 1, 1);
    Matrix Cz_Mercury = PCtemp.subMatrix((int) temp(1, 3), (int)temp(1, 4) - 1, 1);
    for (i = 1; i <= 3; i++) {
        temp = temp+42;
        Cx = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
        Cy = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
        Cz = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1, 4) - 1, 1);
        Cx_Mercury = Cx_Mercury.concatenar(Cx);
        Cy_Mercury = Cy_Mercury.concatenar(Cy);
        Cz_Mercury = Cz_Mercury.concatenar(Cz);
    }
    if (0 <= dt && dt <= 8) {
        j = 0;
        Mjd0 = t1;
    } else if (8 < dt && dt <= 16) {
        j = 1;
        Mjd0 = t1 + 8 * j;
    } else if (16 < dt && dt <= 24) {
        j = 2;
        Mjd0 = t1 + 8 * j;
    } else if (24 < dt && dt <= 32) {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }
    aux1 = Cx_Mercury.subMatrix(14 * j + 1, 14 * j + 14, 1);
    aux2 = Cy_Mercury.subMatrix(14 * j + 1, 14 * j + 14, 1);
    aux3 = Cz_Mercury.subMatrix(14 * j + 1, 14 * j + 14, 1);
    r_Mercury = 1e3 * Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8, aux1, aux2, aux3).traspuesta();

    i = 1;
    for (int k = 171; k <= 201; k += 10) {
        temp(1, i) = k;
        i++;
    }
    Matrix Cx_Venus = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1,2) - 1, 1);
    Matrix Cy_Venus = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1,3) - 1, 1);
    Matrix Cz_Venus = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1,4) - 1, 1);
    temp = temp+30;
    Cx = PCtemp.subMatrix((int)temp(1,1), (int)temp(1,2) - 1, 1);
    Cy = PCtemp.subMatrix((int)temp(1,2), (int)temp(1,3) - 1, 1);
    Cz = PCtemp.subMatrix((int)temp(1,3), (int)temp(1,4) - 1, 1);
    Cx_Venus = Cx_Venus.concatenar(Cx);
    Cy_Venus = Cy_Venus.concatenar(Cy);
    Cz_Venus = Cz_Venus.concatenar(Cz);
    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    } else if (16 < dt && dt <= 32) {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }
    aux1 = Cx_Venus.subMatrix(10 * j + 1, 10 * j + 10, 1);
    aux2 = Cy_Venus.subMatrix(10 * j + 1, 10 * j + 10, 1);
    aux3 = Cz_Venus.subMatrix(10 * j + 1, 10 * j + 10, 1);
    r_Venus = 1e3 * Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16, aux1, aux2, aux3).traspuesta();

    i = 1;
    for (int k = 309; k <= 342; k += 11) {
        temp(1, i) = k;
        i++;
    }
    Matrix Cx_Mars = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
    Matrix Cy_Mars = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
    Matrix Cz_Mars = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1, 4) - 1, 1);
    j = 0;
    Mjd0 = t1;
    aux1 = Cx_Mars.subMatrix(11 * j + 1, 11 * j + 11, 1);
    aux2 = Cy_Mars.subMatrix(11 * j + 1, 11 * j + 11, 1);
    aux3 = Cz_Mars.subMatrix(11 * j + 1, 11 * j + 11, 1);
    r_Mars = 1e3 * Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32, aux1, aux2, aux3).traspuesta();

    i = 1;
    for (int k = 342; k <= 366; k += 8) {
        temp(1, i) = k;
        i++;
    }
    Matrix Cx_Jupiter = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
    Matrix Cy_Jupiter = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
    Matrix Cz_Jupiter = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1, 4) - 1, 1);
    j = 0;
    Mjd0 = t1;
    aux1 = Cx_Jupiter.subMatrix(8 * j + 1, 8 * j + 8, 1);
    aux2 = Cy_Jupiter.subMatrix(8 * j + 1, 8 * j + 8, 1);
    aux3 = Cz_Jupiter.subMatrix(8 * j + 1, 8 * j + 8, 1);
    r_Jupiter = 1e3 * Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0+32, aux1, aux2, aux3).traspuesta();

    i = 1;
    for (int k = 366; k <= 387; k += 7) {
        temp(1, i) = k;
        i++;
    }
    Matrix Cx_Saturn = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
    Matrix Cy_Saturn = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
    Matrix Cz_Saturn = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1, 4) - 1, 1);
    j = 0;
    Mjd0 = t1;
    aux1 = Cx_Saturn.subMatrix(7 * j + 1, 7 * j + 7, 1);
    aux2 = Cy_Saturn.subMatrix(7 * j + 1, 7 * j + 7, 1);
    aux3 = Cz_Saturn.subMatrix(7 * j + 1, 7 * j + 7, 1);
    r_Saturn = 1e3 * Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32, aux1, aux2, aux3).traspuesta();

    i = 1;
    for (int k = 387; k <= 405; k += 6) {
        temp(1, i) = k;
        i++;
    }
    Matrix Cx_Uranus = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
    Matrix Cy_Uranus = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
    Matrix Cz_Uranus = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1, 4) - 1, 1);
    j = 0;
    Mjd0 = t1;
    aux1 = Cx_Uranus.subMatrix(6 * j + 1, 6 * j + 6, 1);
    aux2 = Cy_Uranus.subMatrix(6 * j + 1, 6 * j + 6, 1);
    aux3 = Cz_Uranus.subMatrix(6 * j + 1, 6 * j + 6, 1);
    r_Uranus = 1e3 * Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, aux1, aux2, aux3).traspuesta();

    i = 1;
    for (int k = 405; k <= 423; k += 6) {
        temp(1, i) = k;
        i++;
    }
    Matrix Cx_Neptune = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
    Matrix Cy_Neptune = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
    Matrix Cz_Neptune = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1, 4) - 1, 1);
    j = 0;
    Mjd0 = t1;
    aux1 = Cx_Neptune.subMatrix(6 * j + 1, 6 * j + 6, 1);
    aux2 = Cy_Neptune.subMatrix(6 * j + 1, 6 * j + 6, 1);
    aux3 = Cz_Neptune.subMatrix(6 * j + 1, 6 * j + 6, 1);
    r_Neptune = 1e3 * Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, aux1, aux2, aux3).traspuesta();

    i = 1;
    for (int k = 423; k <= 441; k += 6) {
        temp(1, i) = k;
        i++;
    }
    Matrix Cx_Pluto = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
    Matrix Cy_Pluto = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
    Matrix Cz_Pluto = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1, 4) - 1, 1);
    j = 0;
    Mjd0 = t1;
    aux1 = Cx_Pluto.subMatrix(6 * j + 1, 6 * j + 6, 1);
    aux2 = Cy_Pluto.subMatrix(6 * j + 1, 6 * j + 6, 1);
    aux3 = Cz_Pluto.subMatrix(6 * j + 1, 6 * j + 6, 1);
    r_Pluto = 1e3*Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, aux1, aux2, aux3).traspuesta();

    i = 1;
    for (int k = 819; k <= 839; k += 10) {
        temp(1, i) = k;
        i++;
    }
    Matrix Cx_Nutations = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
    Matrix Cy_Nutations = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
    for (i = 1; i <= 3; i++) {
        temp = temp+20;
        Cx = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
        Cy = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
        Cx_Nutations = Cx_Nutations.concatenar(Cx);
        Cy_Nutations = Cy_Nutations.concatenar(Cy);
    }
    if (0 <= dt && dt <= 8) {
        j = 0;
        Mjd0 = t1;
    } else if (8 < dt && dt <= 16) {
        j = 1;
        Mjd0 = t1 + 8 * j;
    } else if (16 < dt && dt <= 24) {
        j = 2;
        Mjd0 = t1 + 8 * j;
    } else if (24 < dt && dt <= 32) {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }
    aux1 = Cx_Nutations.subMatrix(10 * j + 1, 10 * j + 10, 1);
    aux2 = Cy_Nutations.subMatrix(10 * j + 1, 10 * j + 10, 1);
    aux3 = Matrix(10,1);
    Matrix Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, aux1, aux2, aux3).traspuesta();

    i = 1;
    for (int k = 899; k <= 929; k += 10) {
        temp(1, i) = k;
        i++;
    }
    Matrix Cx_Librations = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
    Matrix Cy_Librations = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
    Matrix Cz_Librations = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1, 4) - 1, 1);
    for (i = 1; i <= 3; i++) {
        temp = temp+30;
        Cx = PCtemp.subMatrix((int)temp(1, 1), (int)temp(1, 2) - 1, 1);
        Cy = PCtemp.subMatrix((int)temp(1, 2), (int)temp(1, 3) - 1, 1);
        Cz = PCtemp.subMatrix((int)temp(1, 3), (int)temp(1, 4) - 1, 1);
        Cx_Librations = Cx_Librations.concatenar(Cx);
        Cy_Librations = Cy_Librations.concatenar(Cy);
        Cz_Librations = Cz_Librations.concatenar(Cz);
    }
    if (0 <= dt && dt <= 8) {
        j = 0;
        Mjd0 = t1;
    } else if (8 < dt && dt <= 16) {
        j = 1;
        Mjd0 = t1 + 8 * j;
    } else if (16 < dt && dt <= 24) {
        j = 2;
        Mjd0 = t1 + 8 * j;
    } else if (24 < dt && dt <= 32) {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }
    aux1 = Cx_Librations.subMatrix(10 * j + 1, 10 * j + 10, 1);
    aux2 = Cy_Librations.subMatrix(10 * j + 1, 10 * j + 10, 1);
    aux3 = Cz_Librations.subMatrix(10 * j + 1, 10 * j + 10, 1);
    Matrix Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, aux1, aux2, aux3).traspuesta();

    double EMRAT = 81.30056907419062;       // DE430
    double EMRAT1 = 1 / (1 + EMRAT);
    Matrix aux4(r_Earth.rows(), r_Earth.columns());
    r_Earth = r_Earth - EMRAT1 * r_Moon;
    r_Mercury = aux4 - r_Earth + r_Mercury;
    r_Venus = aux4 - r_Earth + r_Venus;
    r_Mars = aux4 - r_Earth + r_Mars;
    r_Jupiter = aux4 - r_Earth + r_Jupiter;
    r_Saturn = aux4 - r_Earth + r_Saturn;
    r_Uranus = aux4 - r_Earth + r_Uranus;
    r_Neptune = aux4 - r_Earth + r_Neptune;
    r_Pluto = aux4 - r_Earth + r_Pluto;
    r_Sun = aux4 - r_Earth + r_Sun;
}