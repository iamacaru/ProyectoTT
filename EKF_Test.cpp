#include "./include/Matrix.h"
#include "./include/R_x.h"
#include "./include/R_y.h"
#include "./include/R_z.h"
#include "./include/Sign_.h"
#include "./include/Global.h"
#include "./include/Timediff.h"
#include "./include/Unit.h"
#include "./include/AccelPointMass.h"
#include "./include/AzElPa.h"
#include "./include/Cheb3D.h"
#include "./include/EccAnom.h"
#include "./include/Frac.h"
#include "./include/Geodetic.h"
#include "./include/MeanObliquity.h"
#include "./include/Mjday.h"
#include "./include/Mjday_TDB.h"
#include "./include/NutAngles.h"
#include "./include/Position.h"
#include "./include/IERS.h"
#include "./include/Global.h"
#include "./include/Legendre.h"
#include "./include/SAT_Const.h"


#include <iostream>


#define TOL_ 10e-14

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

using namespace std;

int proMat_01()
{
    double v1[] = {1.0, 2.0, 3.0, 4.0};
    double v2[] = {1.0, 0.0, 0.0, 1.0};
    Matrix m1(2, 2, v1, 4);
    Matrix m2(2, 2, v2, 4);
    Matrix sol(2, 2);
    
    sol = m1 * m2;

    m1.print();
    m2.print();
    sol.print();

    _assert(sol(1,1) == m1(1,1) && sol(1,2) == m1(1,2) && sol(2,1) == m1(2,1) && sol(2,2) == m1(2,2));
    
    return 0;
}

int R_x_01() {
    double alpha = 1.0;
    Matrix sol(3,3);

    sol = R_x(alpha);

    _assert(fabs(sol(1,1) - 1) < TOL_ && fabs(sol(1,2) - 0) < TOL_ && fabs(sol(1,3) - 0 ) < TOL_ &&
            fabs(sol(2,1) - 0) < TOL_ && fabs(sol(2,2) - 0.540302305868140) < TOL_ && fabs(sol(2,3) - 0.841470984807897) < TOL_ &&
            fabs(sol(3,1) - 0) < TOL_ && fabs(sol(3,2) - -0.841470984807897) < TOL_ && fabs(sol(3,3) - 0.540302305868140) < TOL_);

    return 0;
}

int R_y_01() {
    double alpha = 1.0;
    Matrix sol(3,3);

    sol = R_y(alpha);

    _assert(fabs(sol(1,1) - 0.540302305868140) < TOL_ && fabs(sol(1,2) - 0) < TOL_ && fabs(sol(1,3) - -0.841470984807897) < TOL_ &&
            fabs(sol(2,1) - 0) < TOL_ && fabs(sol(2,2) - 1) < TOL_ && fabs(sol(2,3) - 0) < TOL_ &&
            fabs(sol(3,1) - 0.841470984807897) < TOL_ && fabs(sol(3,2) - 0) < TOL_ && fabs(sol(3,3) - 0.540302305868140) < TOL_);

    return 0;
}

int R_z_01() {
    double alpha = 1.0;
    Matrix sol(3,3);

    sol = R_z(alpha);

    _assert(fabs(sol(1,1) - 0.540302305868140) < TOL_ && fabs(sol(1,2) - 0.841470984807897) < TOL_ && fabs(sol(1,3) - 0) < TOL_ &&
            fabs(sol(2,1) - -0.841470984807897) < TOL_ && fabs(sol(2,2) - 0.540302305868140) < TOL_ && fabs(sol(2,3) - 0) < TOL_ &&
            fabs(sol(3,1) - 0) < TOL_ && fabs(sol(3,2) - 0) < TOL_ && fabs(sol(3,3) - 1) < TOL_);

    return 0;
}

int sign_01() {
    _assert(sign_(10, 5) == 10);
    _assert(sign_(10, -5) == -10);
    _assert(sign_(-10, 5) == 10);
    _assert(sign_(-10, -5) == -10);

    return 0;
}

int timediff_01() {
    Matrix sol = timediff(3, 9);

    _assert(fabs(sol(1,1) - -6) < TOL_);
    _assert(fabs(sol(1,2) - 10) < TOL_);
    _assert(fabs(sol(1,3) - 13) < TOL_);
    _assert(fabs(sol(1,4) - 41.1840) < TOL_);
    _assert(fabs(sol(1,5) - -10) < TOL_);

    return 0;
}

int unit_01() {
    double values[] = {3.0, 4.0, 0.0};
    Matrix sol = unit(Matrix(1, 3, values, 3));

    _assert(fabs(sol(1,1) - 0.6) < TOL_);
    _assert(fabs(sol(1,2) - 0.8) < TOL_);
    _assert(fabs(sol(1,3) - 0.0) < TOL_);

    return 0;
}

int AccelPointMass_01() {

    double valuesr[] = {5.0, 10.0, 15.0};
    Matrix r(1, 3, valuesr, 3);

    double valuess[] = {25.0, 30.0, 35.0};
    Matrix s(1, 3, valuess, 3);

    double GM = 500;

    Matrix sol = AccelPointMass(r, s, GM);

    _assert(fabs(sol(1,1) - 0.153884194958199) < TOL_);
    _assert(fabs(sol(1,2) - 0.136548511517370) < TOL_);
    _assert(fabs(sol(1,3) - 0.119212828076542) < TOL_);

    return 0;
}

int AzElPa_01() {

    double values[] = {5.0, 10.0, 15.0};
    Matrix s(1, 3, values, 3);

    double Az;
    double El;
    Matrix dAds(1,3);
    Matrix dEds(1,3);
    AzElPa(s, Az, El, dAds, dEds);
    _assert(fabs(Az - 0.463647609000806) < TOL_);
    _assert(fabs(El - 0.930274014115472) < TOL_);
    _assert(fabs(dAds(1,1) - 0.08) < TOL_);
    _assert(fabs(dAds(1,2) - -0.04) < TOL_);
    _assert(fabs(dAds(1,3) - 0.0) < TOL_);
    _assert(fabs(dEds(1,1) - -0.019166296949998) < TOL_);
    _assert(fabs(dEds(1,2) - -0.038332593899996) < TOL_);
    _assert(fabs(dEds(1,3) - 0.031943828249997) < TOL_);

    return 0;
}

int Cheb3D_01() {
    double values1[] = {1.0,2.0,3.0};
    Matrix Cx(1,3, values1, 3);
    double values2[] = {4.0,5.0,6.0};
    Matrix Cy(1,3, values2, 3);
    double values3[] = {7.0,8.0,9.0};
    Matrix Cz(1,3, values3, 3);

    Matrix sol(1,3);
    sol = Cheb3D(5.0, 2.0, 4.0, 6.0, Cx, Cy, Cz);
    _assert(fabs(sol(1,1) - 1.0) < TOL_);
    _assert(fabs(sol(1,2) - 4.0) < TOL_);
    _assert(fabs(sol(1,3) - 7.0) < TOL_);

    return 0;
}

int EccAnom_01() {
    _assert(fabs(EccAnom(5, 10) - 3.311277978895311) < TOL_);
    _assert(fabs(EccAnom(42, 75) - 3.156847125616838) < TOL_);
    _assert(fabs(EccAnom(17, 64) - 3.161471430843077) < TOL_);

    return 0;
}

int Frac_01() {
    _assert(fabs(Frac(5.75) - 0.75) < TOL_);
    _assert(fabs(Frac(42.541515151) - 0.541515151) < TOL_);
    _assert(fabs(Frac(-5.75) - 0.25) < TOL_);
    _assert(fabs(Frac(-2.2) - 0.8) < TOL_);

    return 0;
}

int Geodetic_01() {
    double values[] = {1, 2, 3};
    Matrix r(1,3,values,3);

    Matrix sol = Geodetic(r);

    _assert(fabs(sol(1,1) - 1.107148717794090) < TOL_);
    _assert(fabs(sol(1,2) - 1.570744136243924) < TOL_);
    _assert(fabs(sol(1,3) - -6.356748616533795e+06) < TOL_);

    return 0;
}

int MeanObliquity_01() {
    _assert(fabs(MeanObliquity(45165) - 0.409132446157106) < TOL_);
    _assert(fabs(MeanObliquity(7) - 0.409413026695083) < TOL_);
    _assert(fabs(MeanObliquity(13857) - 0.409326980794572) < TOL_);

    return 0;
}

int Mjday_01() {
    _assert(fabs(Mjday(1939,5,20) - 29403) < TOL_);
    _assert(fabs(Mjday(2024,4,17,16,59,42) - 6.041770812499989e+04) < TOL_);
    _assert(fabs(Mjday(4,2,29,00,59,59) - -6.774359583449075e+05) < TOL_);
    _assert(fabs(Mjday(0,0,0,0,0,1) - -6.789869999884260e+05) < TOL_);

    return 0;
}

int Mjday_TDB_01() {
    _assert(fabs(Mjday_TDB(55555.0) - 5.555499999999645e+04) < TOL_);
    _assert(fabs(Mjday_TDB(99999.0) - 9.999899999998482e+04) < TOL_);
    _assert(fabs(Mjday_TDB(50.0) - 50.000000001233296) < TOL_);

    return 0;
}

int NutAngles_01() {
    Matrix sol = NutAngles(55555.0);

    _assert(fabs(sol(1,1) - 8.481779033666464e-05) < TOL_);
    _assert(fabs(sol(1,2) - -7.785278092977075e-08) < TOL_);

    return 0;
}

int Position_01() {
    Matrix sol = Position(5.0, 0.75, 200.0);

    _assert(fabs(sol(1,1) - 1.32590301129570906e+06) < TOL_);
    _assert(fabs(sol(1,2) - -4.4822350265126815e+06) < TOL_);
    _assert(fabs(sol(1,3) - 4.32534871154858637e+06) < TOL_);

    return 0;
}

int IERS_01() {
    Global::eop19620101(21413);
    double Mjd_UTC = 37665.0;
    char interp = 'l';

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;

    IERS(*Global::eopdata, Mjd_UTC, interp, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);

    _assert(fabs(x_pole - -6.157133750091107e-08) < TOL_);
    _assert(fabs(y_pole - 1.032653140763312e-06) < TOL_);
    _assert(fabs(UT1_UTC - 0.032633800000000) < TOL_);
    _assert(fabs(LOD - 0.001723000000000) < TOL_);
    _assert(fabs(dpsi - 3.115461196177989e-07) < TOL_);
    _assert(fabs(deps - 2.941364603291555e-08) < TOL_);
    _assert(fabs(dx_pole - 0) < TOL_);
    _assert(fabs(dy_pole - 0) < TOL_);
    _assert(fabs(TAI_UTC - 2) < TOL_);

    return 0;
}

int Legendre_01() {
    int n = 2;
    int m = 3;
    double fi = Constants::pi/4;

    Matrix pnm(n+1,m+1);
    Matrix dpnm(n+1,m+1);
    Legendre(n,m,fi,pnm,dpnm);

    _assert(fabs(pnm(1,1) - 1.0) < TOL_);
    _assert(fabs(pnm(1,2) - 0) < TOL_);
    _assert(fabs(pnm(1,3) - 0) < TOL_);
    _assert(fabs(pnm(1,4) - 0) < TOL_);
    _assert(fabs(pnm(2,1) - 1.224744871391589) < TOL_);
    _assert(fabs(pnm(2,2) - 1.224744871391589) < TOL_);
    _assert(fabs(pnm(2,3) - 0) < TOL_);
    _assert(fabs(pnm(2,4) - 0) < TOL_);
    _assert(fabs(pnm(3,1) - 0.559016994374947) < TOL_);
    _assert(fabs(pnm(3,2) - 1.936491673103709) < TOL_);
    _assert(fabs(pnm(3,3) - 0.968245836551854) < TOL_);
    _assert(fabs(pnm(3,4) - 0) < TOL_);

    _assert(fabs(dpnm(1,1) - 0) < TOL_);
    _assert(fabs(dpnm(1,2) - 0) < TOL_);
    _assert(fabs(dpnm(1,3) - 0) < TOL_);
    _assert(fabs(dpnm(1,4) - 0) < TOL_);
    _assert(fabs(dpnm(2,1) - 1.224744871391589) < TOL_);
    _assert(fabs(dpnm(2,2) - -1.224744871391589) < TOL_);
    _assert(fabs(dpnm(2,3) - 0) < TOL_);
    _assert(fabs(dpnm(2,4) - 0) < TOL_);
    _assert(fabs(dpnm(3,1) - 3.354101966249685) < TOL_);
    _assert(fabs(dpnm(3,2) - 0.000000000000001) < TOL_);
    _assert(fabs(dpnm(3,3) - -1.936491673103709) < TOL_);
    _assert(fabs(dpnm(3,4) - 0) < TOL_);

    return 0;
}



int all_tests()
{
    _verify(proMat_01);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(sign_01);
    _verify(timediff_01);
    _verify(unit_01);
    _verify(AccelPointMass_01);
    _verify(AzElPa_01);
    _verify(Cheb3D_01);
    _verify(EccAnom_01);
    _verify(Frac_01);
    _verify(Geodetic_01);
    _verify(MeanObliquity_01);
    _verify(Mjday_01);
    _verify(Mjday_TDB_01);
    _verify(NutAngles_01);
    _verify(Position_01);
    _verify(IERS_01);
    _verify(Legendre_01);

    return 0;
}


int main()
{
    Global::eop19620101(6);
    Global::eopdata->print();

    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}

