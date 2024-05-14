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
#include "./include/Legendre.h"
#include "./include/SAT_Const.h"
#include "./include/EqnEquinox.h"
#include "./include/Gmst.h"
#include "./include/Gast.h"
#include "./include/GHAMatrix.h"
#include "./include/NutMatrix.h"
#include "./include/PoleMatrix.h"
#include "./include/PrecMatrix.h"
#include "./include/TimeUpdate.h"
#include "./include/VarEqn.h"
#include "./include/AccelHarmonic.h"
#include "./include/G_AccelHarmonic.h"
#include "./include/VarEqn.h"
#include "./include/Elements.h"
#include "./include/JPL_Eph_DE430.h"

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

    double values1[3] = {1.0,2.0,3.0};
    Matrix Cx(1,3, values1, 3);
    double values2[3] = {4.0,5.0,6.0};
    Matrix Cy(1,3, values2, 3);
    double values3[3] = {7.0,8.0,9.0};
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

int EqnEquinox_01() {
    _assert(fabs(EqnEquinox(2451545.0) - -6.70392825380816e-05) < TOL_);
    _assert(fabs(EqnEquinox(2451550.5) - -6.57426418935655e-05) < TOL_);
    _assert(fabs(EqnEquinox(2451560.75) - -6.39866793613151e-05) < TOL_);
    _assert(fabs(EqnEquinox(2451575.3) - -6.30997404226282e-05) < TOL_);

    return 0;
}

int Gmst_01() {
    _assert(fabs(gmst(582003) - 3.93686610296254) < TOL_);
    _assert(fabs(gmst(54417) - 0.902905604417803) < TOL_);
    _assert(fabs(gmst(61597) - 5.03843040108836) < TOL_);
    _assert(fabs(gmst(39457) - 1.15973849105372) < TOL_);

    return 0;
}

int Gast_01() {
    _assert(fabs(gast(52135.4) - 1.86532951201256) < TOL_);
    _assert(fabs(gast(42876.6) - 0.924514810241508) < TOL_);
    _assert(fabs(gast(67234.1) - 2.10962351700529) < TOL_);
    _assert(fabs(gast(39851.5) - 4.80461613750777) < TOL_);

    return 0;
}

int GHAMatrix_01() {
    Matrix sol(3,3);

    sol = GHAMatrix(241575.6);

    _assert(fabs(sol(1,1) - 0.498164074615664) < TOL_ && fabs(sol(1,2) - 0.867082784261295) < TOL_ && fabs(sol(1,3) - 0) < TOL_ &&
            fabs(sol(2,1) - -0.867082784261295 ) < TOL_ && fabs(sol(2,2) - 0.498164074615664) < TOL_ && fabs(sol(2,3) - 0) < TOL_ &&
            fabs(sol(3,1) - 0) < TOL_ && fabs(sol(3,2) - 0) < TOL_ && fabs(sol(3,3) - 1) < TOL_);

    return 0;
}

int NutMatrix_01() {
    Matrix sol(3,3);

    sol = NutMatrix(2524686.6);

    _assert(fabs(sol(1,1) - 0.999999998843278) < TOL_ && fabs(sol(1,2) - -4.43678484217162e-05) < TOL_ && fabs(sol(1,3) - -1.85725364151702e-05) < TOL_ &&
            fabs(sol(2,1) - 4.43685056896323e-05) < TOL_ && fabs(sol(2,2) - 0.99999999838948) < TOL_ && fabs(sol(2,3) - 3.53903222812257e-05) < TOL_ &&
            fabs(sol(3,1) - 1.85709661928042e-05) < TOL_ && fabs(sol(3,2) - -3.5391146275876e-05) < TOL_ && fabs(sol(3,3) - 0.999999999201293) < TOL_);

    return 0;
}

int PoleMatrix_01() {
    Matrix sol(3,3);

    sol = PoleMatrix(15.4, 6.8);

    _assert(fabs(sol(1,1) - -0.95295291688718) < TOL_ && fabs(sol(1,2) - 0.149774827043247) < TOL_ && fabs(sol(1,3) - 0.263530338633677) < TOL_ &&
            fabs(sol(2,1) - 0) < TOL_ && fabs(sol(2,2) - 0.869397490349825) < TOL_ && fabs(sol(2,3) - -0.494113351138608) < TOL_ &&
            fabs(sol(3,1) - -0.303118356745702 ) < TOL_ && fabs(sol(3,2) - -0.470866759240436) < TOL_ && fabs(sol(3,3) - -0.82849487436326) < TOL_);

    return 0;
}

int PrecMatrix_01() {
    Matrix sol(3,3);

    sol = PrecMatrix(2358142.6, 2386475.2);

    _assert(fabs(sol(1,1) - 0.999810970866054) < TOL_ && fabs(sol(1,2) - -0.0180081535228055) < TOL_ && fabs(sol(1,3) - -0.0073300029043494) < TOL_ &&
            fabs(sol(2,1) - 0.0180081534586217) < TOL_ && fabs(sol(2,2) - 0.999837837877244) < TOL_ && fabs(sol(2,3) - -6.60149034239962e-05) < TOL_ &&
            fabs(sol(3,1) - 0.00733000306203428) < TOL_ && fabs(sol(3,2) - -6.59973924696928e-05) < TOL_ && fabs(sol(3,3) - 0.999973132988809) < TOL_);

    return 0;
}

int TimeUpdate_01() {
    double values1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    Matrix P(3,3, values1, 9);
    double values2[] = {9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
    Matrix Phi(3,3, values2, 9);

    Matrix sol1(3,3);
    sol1 = TimeUpdate(P, Phi, 40);
    Matrix sol2(3,3);
    sol2 = TimeUpdate(P, Phi);
    _assert(fabs(sol1(1,1) - 2728) < TOL_ && fabs(sol1(1,2) - 1702) < TOL_ && fabs(sol1(1,3) - 676) < TOL_ &&
            fabs(sol1(2,1) - 1666) < TOL_ && fabs(sol1(2,2) - 1045) < TOL_ && fabs(sol1(2,3) - 424) < TOL_ &&
            fabs(sol1(3,1) - 604) < TOL_ && fabs(sol1(3,2) - 388) < TOL_ && fabs(sol1(3,3) - 172) < TOL_);

    _assert(fabs(sol2(1,1) - 2688) < TOL_ && fabs(sol2(1,2) - 1662) < TOL_ && fabs(sol2(1,3) - 636) < TOL_ &&
            fabs(sol2(2,1) - 1626) < TOL_ && fabs(sol2(2,2) - 1005) < TOL_ && fabs(sol2(2,3) - 384) < TOL_ &&
            fabs(sol2(3,1) - 564) < TOL_ && fabs(sol2(3,2) - 348) < TOL_ && fabs(sol2(3,3) - 132) < TOL_);

    return 0;
}

int AccelHarmonic_01() {
    Global::GGM03S();

    double values1[] = {5542555.93722861, 3213514.86734920, 3990892.97587686};
    Matrix r(3,1, values1, 3);
    double values2[] = {-0.976675972331716, 0.214718082511189, -0.000436019054674645, -0.214718043811152,
                        -0.976676068937815, -0.000134261271504216, -0.000454677699074514, -3.75085994087200e-05, 0.999999895930642};
    Matrix E(3,3, values2, 9);
    int n_max = 20;
    int m_max = 20;

    Matrix sol = AccelHarmonic(r, E, n_max, m_max);

    _assert(fabs(sol(1,1) - -5.13483678540849) < TOL_);
    _assert(fabs(sol(2,1) - -2.97717622353621) < TOL_);
    _assert(fabs(sol(3,1) - -3.70591776714203) < TOL_);

    return 0;
}

int G_AccelHarmonic_01() {
    Global::GGM03S();

    double values1[] = {5542555.93722861, 3213514.86734920, 3990892.97587686};
    Matrix r(3,1, values1, 3);
    double values2[] = {-0.976675972331716, 0.214718082511189, -0.000436019054674645, -0.214718043811152,
                        -0.976676068937815, -0.000134261271504216, -0.000454677699074514, -3.75085994087200e-05, 0.999999895930642};
    Matrix U(3,3, values2, 9);
    int n_max = 20;
    int m_max = 20;

    Matrix sol = G_AccelHarmonic(r, U, n_max, m_max);

    _assert(fabs(sol(1,1) - 5.7003203313144e-07) < TOL_ && fabs(sol(1,2) - 8.67651595015673e-07) < TOL_ && fabs(sol(1,3) - 1.0816935356317e-06) < TOL_ &&
            fabs(sol(2,1) - 8.6765159101887e-07) < TOL_ && fabs(sol(2,2) - -4.23359105550247e-07) < TOL_ && fabs(sol(2,3) - 6.27183703194589e-07) < TOL_ &&
            fabs(sol(3,1) - 1.08169353651988e-06) < TOL_ && fabs(sol(3,2) - 6.27183704526857e-07) < TOL_ && fabs(sol(3,3) - -1.46672928913461e-07) < TOL_);

    // EN MATLAB G_AccelHarmonic([5542555.93722861; 3213514.86734920; 3990892.97587686],[-0.976675972331716, 0.214718082511189, -0.000436019054674645; -0.214718043811152, -0.976676068937815, -0.000134261271504216; -0.000454677699074514, -3.75085994087200e-05, 0.999999895930642],20,20)

    return 0;
}

int VarEqn_01() {
    Global::GGM03S();
    Global::eop19620101(21413);
    Global::AuxParam::Mjd_UTC = 49746.1101504630;
    Global::AuxParam::n = 20;
    Global::AuxParam::m = 20;

    double values[] = {5542555.93722861, 3213514.86734920, 3990892.97587686, 5394.06842166353, -2365.21337882342,
                        -7061.84554200298, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                        0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                        0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
    Matrix yPhi(42,1, values, 42);
    double x = 0.0;

    VarEqn(x, yPhi);

    return 0;
}

int Elements_01() {
    double values[] = {1,2,3,4,5,6};
    Matrix y(1,6,values,6);
    double p, a, e, i, Omega, omega, M;

    elements(y,p, a, e, i, Omega, omega, M);

    _assert(fabs(p - 1.35474011564823e-13) < TOL_);
    _assert(fabs(a - 1.87082869338765) < TOL_);
    _assert(fabs(e - 0.999999999999964) < TOL_);
    _assert(fabs(i - 1.99133066207886) < TOL_);
    _assert(fabs(Omega - 3.6052402625906) < TOL_);
    _assert(fabs(omega - 5.21086941752228) < TOL_);
    _assert(fabs(M - 3.14159030993265) < TOL_);

    return 0;
}

int JPL_Eph_DE430_01() {
    Global::Pc();
    Matrix r_Mercury(3,3), r_Venus(3,3), r_Earth(3,3), r_Mars(3,3), r_Jupiter(3,3),
    r_Saturn(3,3), r_Uranus(3,3), r_Neptune(3,3), r_Pluto(3,3), r_Moon(3,3), r_Sun(3,3);
    double Mjd_TDB = 4.974611199287850e+04;

    JPL_Eph_DE430(r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune,
                  r_Pluto, r_Moon, r_Sun, Mjd_TDB);

    r_Mercury.print();
    r_Venus.print();
    r_Earth.print();
    r_Mars.print();
    r_Jupiter.print();
    r_Saturn.print();
    r_Uranus.print();
    r_Neptune.print();
    r_Pluto.print();
    r_Moon.print();
    r_Sun.print();

    cout << r_Mercury(1,1) - 8.377549589569574*1.0e+10 << endl;
    cout << r_Mercury(2,1) - -6.529112491344618*1.0e+10 << endl;
    cout << r_Mercury(3,1) - -2.339131210124083*1.0e+10 << endl;
    cout << r_Mercury(1,1) - 83775495895.6957 << endl;
    cout << r_Mercury(2,1) - -65291124913.4462 << endl;
    cout << r_Mercury(3,1) - -23391312101.2408 << endl;
    cout << r_Venus(1,1) - -0.152296655739533*1.0e+11 << endl;
    cout << r_Venus(2,1) - -1.101349926375630*1.0e+11 << endl;
    cout << r_Venus(3,1) - -0.410218036256266*1.0e+11 << endl;


    _assert(fabs(r_Mercury(1,1) - 8.377549589569574*1.0e+10) < TOL_ && fabs(r_Mercury(2,1) - -6.529112491344618*1.0e+10) < TOL_ && fabs(r_Mercury(3,1) - -2.339131210124083*1.0e+10) < TOL_);
    _assert(fabs(r_Venus(1,1) - -0.152296655739533*1.0e+11) < TOL_ && fabs(r_Venus(2,1) - -1.101349926375630*1.0e+11) < TOL_ && fabs(r_Venus(3,1) - -0.410218036256266*1.0e+11) < TOL_);
    _assert(fabs(r_Earth(1,1) - -0.924709612291923*1.0e+11) < TOL_ && fabs(r_Earth(2,1) - 1.063949183894933*1.0e+11) < TOL_ && fabs(r_Earth(3,1) - 0.461301399094037*1.0e+11) < TOL_);
    _assert(fabs(r_Mars(1,1) - -8.827841300824323*1.0e+10) < TOL_ && fabs(r_Mars(2,1) - 4.696476977829840*1.0e+10) < TOL_ && fabs(r_Mars(3,1) - 2.907102650267133*1.0e+10) < TOL_);
    _assert(fabs(r_Jupiter(1,1) - -2.983859364660936*1.0e+11) < TOL_ && fabs(r_Jupiter(2,1) - -7.544982589107291*1.0e+11) < TOL_ && fabs(r_Jupiter(3,1) - -3.144105185682280*1.0e+11) < TOL_);
    _assert(fabs(r_Saturn(1,1) - 1.482033999505117*1.0e+12) < TOL_ && fabs(r_Saturn(2,1) - -0.453872894236397*1.0e+12) < TOL_ && fabs(r_Saturn(3,1) - -0.249402247811211*1.0e+12) < TOL_);
    _assert(fabs(r_Uranus(1,1) - 1.412367984017928*1.0e+12) < TOL_ && fabs(r_Uranus(2,1) - -2.511355045786904*1.0e+12) < TOL_ && fabs(r_Uranus(3,1) - -1.118108651902601*1.0e+12) < TOL_);
    _assert(fabs(r_Neptune(1,1) - 1.871250770052329*1.0e+12) < TOL_ && fabs(r_Neptune(2,1) - -3.928976313605399*1.0e+12) < TOL_ && fabs(r_Neptune(3,1) - -1.655020476718540*1.0e+12) < TOL_);
    _assert(fabs(r_Pluto(1,1) - -2.171414794259537*1.0e+12) < TOL_ && fabs(r_Pluto(2,1) - -3.915433128334976*1.0e+12) < TOL_ && fabs(r_Pluto(3,1) - -0.552716250355845*1.0e+12) < TOL_);
    _assert(fabs(r_Moon(1,1) - 0.893833723127192*1.0e+08) < TOL_ && fabs(r_Moon(2,1) - -3.366038321179463*1.0e+08) < TOL_ && fabs(r_Moon(3,1) - -1.146487877519952*1.0e+08) < TOL_);
    _assert(fabs(r_Sun(1,1) - 0.922982517284766*1.0e+11) < TOL_ && fabs(r_Sun(2,1) - -1.053751960790543*1.0e+11) < TOL_ && fabs(r_Sun(3,1) - -0.456863672263533*1.0e+11) < TOL_);

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
    _verify(EqnEquinox_01);
    _verify(Gmst_01);
    _verify(Gast_01);
    _verify(GHAMatrix_01);
    _verify(NutMatrix_01);
    _verify(PoleMatrix_01);
    _verify(PrecMatrix_01);
    _verify(TimeUpdate_01);
    _verify(AccelHarmonic_01);
    _verify(G_AccelHarmonic_01);
    _verify(VarEqn_01);
    _verify(Elements_01);
    //_verify(JPL_Eph_DE430_01);

    return 0;
}


int main()
{

    extern double Mjd0;
    Mjd0 = Mjday(1995,1,29,02,38,0);

    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}

