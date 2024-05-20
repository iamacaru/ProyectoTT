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
#include "./include/JPL_Eph_DE430.h"
#include "./include/Accel.h"
#include "./include/LTC.h"
#include "./include/MeasUpdate.h"
#include "./include/DEInteg.h"

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
    Matrix sol = R_x(alpha);

    double values[] = {1, 0, 0,
                       0, 0.54030230586814, 0.841470984807897,
                       0, -0.841470984807897, 0.54030230586814};
    Matrix comp(3, 3, values, 9);

    _assert(sol.isEqual(comp, TOL_));

    return 0;
}

int R_y_01() {
    double alpha = 1.0;
    Matrix sol = R_y(alpha);

    double values[] = {0.54030230586814, 0, -0.841470984807897,
                       0, 1, 0,
                       0.841470984807897, 0, 0.54030230586814};
    Matrix comp(3, 3, values, 9);

    _assert(sol.isEqual(comp, TOL_));

    return 0;
}

int R_z_01() {
    double alpha = 1.0;
    Matrix sol = R_z(alpha);

    double values[] = {0.54030230586814, 0.841470984807897, 0,
                       -0.841470984807897, 0.54030230586814,
                       0, 0, 0, 1};
    Matrix comp(3, 3, values, 9);

    _assert(sol.isEqual(comp, TOL_));

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
    double UT1_UTC = 3;
    double TAI_UTC = 9;
    Matrix sol = timediff(UT1_UTC, TAI_UTC);

    double values[] = {-6, 10, 13,  41.184, -10};
    Matrix comp(1, 5, values, 5);

    _assert(sol.isEqual(comp, TOL_));

    return 0;
}

int unit_01() {
    double values1[] = {3.0, 4.0, 0.0};
    Matrix vec(1, 3, values1, 3);
    Matrix sol = unit(vec);

    double values2[] = {0.6, 0.8, 0.0};
    Matrix comp(1, 3, values2, 3);

    _assert(sol.isEqual(comp, TOL_));

    return 0;
}

int AccelPointMass_01() {
    double values1[] = {5.0, 10.0, 15.0};
    Matrix r(1, 3, values1, 3);
    double values2[] = {25.0, 30.0, 35.0};
    Matrix s(1, 3, values2, 3);
    double GM = 500;
    Matrix sol = AccelPointMass(r, s, GM);

    double values3[] = {0.153884194958199, 0.136548511517370, 0.119212828076542};
    Matrix comp(1, 3, values3, 3);

    _assert(sol.isEqual(comp, TOL_));

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

    double values1[] = {0.08, -0.04, 0.0};
    Matrix comp1(1, 3, values1, 3);
    double values2[] = {-0.019166296949998, -0.038332593899996, 0.031943828249997};
    Matrix comp2(1, 3, values2, 3);

    _assert(fabs(Az - 0.463647609000806) < TOL_);
    _assert(fabs(El - 0.930274014115472) < TOL_);
    _assert(dAds.isEqual(comp1, TOL_));
    _assert(dEds.isEqual(comp2, TOL_));

    return 0;
}

int Cheb3D_01() {
    double values1[] = {1.0, 2.0, 3.0};
    Matrix Cx(1,3, values1, 3);
    double values2[] = {4.0, 5.0, 6.0};
    Matrix Cy(1,3, values2, 3);
    double values3[] = {7.0, 8.0, 9.0};
    Matrix Cz(1,3, values3, 3);
    Matrix sol = Cheb3D(5.0, 2.0, 4.0, 6.0, Cx, Cy, Cz);

    double values4[] = {1.0, 4.0, 7.0};
    Matrix comp(1, 3, values4, 3);

    _assert(sol.isEqual(comp, TOL_));

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
    double values1[] = {1, 2, 3};
    Matrix r(1,3,values1,3);
    Matrix sol = Geodetic(r);

    double values2[] = {1.107148717794090, 1.570744136243924, -6.356748616533795e+06};
    Matrix comp(1, 3, values2, 3);

    _assert(sol.isEqual(comp, TOL_));

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
    double Mjd_TT = 55555.0;
    Matrix sol = NutAngles(Mjd_TT);

    double values[] = {8.481779033666464e-05, -7.785278092977075e-08};
    Matrix comp(1, 2, values, 2);

    _assert(sol.isEqual(comp, TOL_));

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
    Matrix pnm(n+1, m+1);
    Matrix dpnm(n+1, m+1);
    Legendre(n, m, fi, pnm, dpnm);

    double values1[] = {1.0, 0.0, 0.0, 0.0,
                        1.224744871391589, 1.224744871391589, 0.0, 0.0,
                        0.559016994374947, 1.936491673103709, 0.968245836551854, 0.0};
    Matrix comp1(3, 4, values1, 12);
    double values2[] = {0, 0, 0, 0,
                        1.224744871391589, -1.224744871391589, 0, 0,
                        3.354101966249685, 0.000000000000001, -1.936491673103709, 0};
    Matrix comp2(3, 4, values2, 12);

    _assert(pnm.isEqual(comp1, TOL_));
    _assert(dpnm.isEqual(comp2, TOL_));

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
    double Mjd_UT1 = 241575.6;
    Matrix sol = GHAMatrix(Mjd_UT1);

    double values[] {0.498164074615664, 0.867082784261295, 0,
                     -0.867082784261295, 0.498164074615664, 0,
                     0, 0, 1};
    Matrix comp(3, 3, values, 9);

    _assert(sol.isEqual(comp, TOL_));

    return 0;
}

int NutMatrix_01() {
    double Mjd_TT = 2524686.6;
    Matrix sol = NutMatrix(Mjd_TT);

    double values[] {0.999999998843278, -4.43678484217162e-05, -1.85725364151702e-05,
                     4.43685056896323e-05, 0.99999999838948, 3.53903222812257e-05,
                     1.85709661928042e-05, -3.5391146275876e-05, 0.999999999201293};
    Matrix comp(3, 3, values, 9);

    _assert(sol.isEqual(comp, TOL_));

    return 0;
}

int PoleMatrix_01() {
    double xp = 15.4;
    double yp = 6.8;
    Matrix sol = PoleMatrix(xp, yp);

    double values[] {-0.95295291688718, 0.149774827043247, .263530338633677,
                     0, 0.869397490349825, -0.494113351138608,
                     -0.303118356745702, -0.470866759240436, -0.82849487436326};
    Matrix comp(3, 3, values, 9);

    _assert(sol.isEqual(comp, TOL_));

    return 0;
}

int PrecMatrix_01() {
    double Mjd_1 = 2358142.6;
    double Mjd_2 = 2386475.2;
    Matrix sol = PrecMatrix(Mjd_1, Mjd_2);

    double values[] {0.999810970866054, -0.0180081535228055, -0.0073300029043494,
                     0.0180081534586217, 0.999837837877244, -6.60149034239962e-05,
                     0.00733000306203428, -6.59973924696928e-05, 0.999973132988809};
    Matrix comp(3, 3, values, 9);

    _assert(sol.isEqual(comp, TOL_));

    return 0;
}

int TimeUpdate_01() {
    double values1[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    Matrix P(3,3, values1, 9);
    double values2[] = {9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
    Matrix Phi(3,3, values2, 9);
    double Qdt = 40;
    Matrix sol1= TimeUpdate(P, Phi, Qdt);
    Matrix sol2 = TimeUpdate(P, Phi);

    double values3[] = {2728, 1702, 676,
                      1666, 1045, 424,
                      604, 388, 172};
    Matrix comp1(3, 3, values3, 9);
    double values4[] = {2688, 1662, 636,
                        1626, 1005, 384,
                        564, 348, 132};
    Matrix comp2(3, 3, values4, 9);

    _assert(sol1.isEqual(comp1, TOL_));
    _assert(sol2.isEqual(comp2, TOL_));

    return 0;
}

int AccelHarmonic_01() {
    double values1[] = {5542555.93722861, 3213514.86734920, 3990892.97587686};
    Matrix r(3,1, values1, 3);
    double values2[] = {-0.976675972331716, 0.214718082511189, -0.000436019054674645, -0.214718043811152,
                        -0.976676068937815, -0.000134261271504216, -0.000454677699074514, -3.75085994087200e-05, 0.999999895930642};
    Matrix E(3,3, values2, 9);
    int n_max = 20;
    int m_max = 20;
    Matrix sol = AccelHarmonic(r, E, n_max, m_max);

    double values3[] = {-5.13483678540849, -2.97717622353621, -3.70591776714203};
    Matrix comp(3, 1, values3, 9);

    _assert(sol.isEqual(comp, TOL_));

    return 0;
}

int G_AccelHarmonic_01() {
    double values1[] = {5542555.93722861, 3213514.86734920, 3990892.97587686};
    Matrix r(3,1, values1, 3);
    double values2[] = {-0.976675972331716, 0.214718082511189, -0.000436019054674645, -0.214718043811152,
                        -0.976676068937815, -0.000134261271504216, -0.000454677699074514, -3.75085994087200e-05, 0.999999895930642};
    Matrix U(3,3, values2, 9);
    int n_max = 20;
    int m_max = 20;
    Matrix sol = G_AccelHarmonic(r, U, n_max, m_max);

    double values3[] = {0.057003203757233 * 1.0e-05, 0.086765159323932 * 1.0e-05, 0.108169353563170 * 1.0e-05,
                        0.086765159412749*1.0e-05,  -0.042335910688251*1.0e-05,   0.062718370319459*1.0e-05,
                        0.108169354007259 * 1.0e-05, 0.062718370319459 * 1.0e-05, -0.014667292891346 * 1.0e-05};
    Matrix comp(3,3,values3,9);

    _assert(sol.isEqual(comp, TOL_));

    return 0;
}

int VarEqn_01() {

    Global::AuxParam::Mjd_UTC = 49746.1101504630;
    Global::AuxParam::Mjd_TT = 49746.1108586111;
    Global::AuxParam::n = 20;
    Global::AuxParam::m = 20;

    double x = 0.0;
    double values1[] = {5542555.93722861, 3213514.8673492, 3990892.97587685, 5394.06842166351, -2365.21337882342,
                         -7061.84554200295, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                         0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                         1, 0, 0, 0, 0, 0, 0, 1};
    Matrix yPhi(42,1, values1, 42);

    Matrix sol = VarEqn(x, yPhi);

    double values2[] = {5394.06842166351, -2365.21337882342, -7061.84554200295, -5.1348367854085, -2.97717622353621, -3.70591776714204,
                        0, 0, 0, 5.70032035795975e-07, 8.67651593239316e-07, 1.08169354007259e-06, 0, 0, 0, 8.67651590574781e-07,
                        -4.23359106882515e-07, 6.27183702306411e-07, 0, 0, 0, 1.08169353651988e-06, 6.27183704082768e-07,
                        -1.46672928913461e-07, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
    Matrix comp(42, 1, values2, 42);

    _assert(sol.isEqual(comp, TOL_));

    return 0;
}

int JPL_Eph_DE430_01() {

    Matrix r_Mercury(3,3), r_Venus(3,3), r_Earth(3,3), r_Mars(3,3), r_Jupiter(3,3),
    r_Saturn(3,3), r_Uranus(3,3), r_Neptune(3,3), r_Pluto(3,3), r_Moon(3,3),
    r_Sun(3,3);
    double Mjd_TDB = 4.974611199287850e+04;
    JPL_Eph_DE430(r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune,
                  r_Pluto, r_Moon, r_Sun, Mjd_TDB);

    double values1[] = {83775495895.6957, -65291124913.4462, -23391312101.2408};
    Matrix comp1(3, 1, values1, 3);
    double values2[] = {-15229665573.9533, -110134992637.563, -41021803625.6266};
    Matrix comp2(3, 1, values2, 3);
    double values3[] = {-92470961229.1923, 106394918389.493, 46130139909.4037};
    Matrix comp3(3, 1, values3, 3);
    double values4[] = {-88278413008.2432, 46964769778.2984, 29071026502.6713};
    Matrix comp4(3, 1, values4, 3);
    double values5[] = {-298385936466.094, -754498258910.729, -314410518568.228};
    Matrix comp5(3, 1, values5, 3);
    double values6[] = {1482033999505.12, -453872894236.397, -249402247811.211};
    Matrix comp6(3, 1, values6, 3);
    double values7[] = {1412367984017.93, -2511355045786.9, -1118108651902.6};
    Matrix comp7(3, 1, values7, 3);
    double values8[] = {1871250770052.33, -3928976313605.4, -1655020476718.54};
    Matrix comp8(3, 1, values8, 3);
    double values9[] = {-2171414794259.54, -3915433128334.98, -552716250355.845};
    Matrix comp9(3, 1, values9, 3);
    double values10[] = {89383372.3127192, -336603832.117946, -114648787.751995};
    Matrix comp10(3, 1, values10, 3);
    double values11[] = {92298251728.4766, -105375196079.054, -45686367226.3533};
    Matrix comp11(3, 1, values11, 3);

    _assert(r_Mercury.isEqual(comp1, 10e-4));
    _assert(r_Venus.isEqual(comp2, 10e-4));
    _assert(r_Earth.isEqual(comp3, 10e-4));
    _assert(r_Mars.isEqual(comp4, 10e-4));
    _assert(r_Jupiter.isEqual(comp5, 10e-4));
    _assert(r_Saturn.isEqual(comp6, 10e-3));
    _assert(r_Uranus.isEqual(comp7, 10e-3));
    _assert(r_Neptune.isEqual(comp8, 10e-3));
    _assert(r_Pluto.isEqual(comp9, 10e-3));
    _assert(r_Moon.isEqual(comp10, 10e-5));
    _assert(r_Sun.isEqual(comp11, 10e-4));

    return 0;
}

int Accel_01() {
    Global::AuxParam::Mjd_UTC = 4.974611128472211e+04;
    Global::AuxParam::n = 20;
    Global::AuxParam::m = 20;
    Global::AuxParam::sun = 1;
    Global::AuxParam::moon = 1;
    Global::AuxParam::planets = 1;

    double x = 0.0;
    double values1[] = {6221397.62857869, 2867713.77965738, 3006155.98509949,
                        4645.04725161807,-2752.21591588205, -7507.99940987033};
    Matrix Y(6,1, values1, 6);
    Matrix sol = Accel(x, Y);

    double values2[] = {4645.04725161807, -2752.21591588205, -7507.99940987033,
                        -5.92414951314006, -2.73076669788113, -2.86933570556259};
    Matrix comp(6, 1, values2, 6);

    _assert(sol.isEqual(comp, TOL_));

    return 0;
}

int LTC_01() {
    double lon = -2.76234307910694;
    double lat = 0.376551295459273;
    Matrix sol = LTC(lon, lat);

    double values[] = {0.370223471399199, -0.928942722252092, 0,
                       0.341586711932422, 0.136136938528208, 0.929938305587722,
                       -0.863859421119156, -0.344284987681776, 0.367715580035218};
    Matrix comp(3, 3, values, 9);

    _assert(sol.isEqual(comp, TOL_));

    return 0;
}

int MeasUpdate_01() {
    double values1[] = {5738566.5776918, 3123975.34092958, 3727114.48156063,
                        5199.63329181126, -2474.43881044665, -7195.16752553894};
    Matrix x(6, 1, values1, 6);
    Matrix z(1, 1);
    z(1, 1) = 1.0559084894933;
    Matrix g(1, 1);
    g(1, 1) = 1.05892995381513;
    Matrix s(1, 1);
    s(1, 1) = 0.00039095375244673;
    double values2[] = {9.59123748602943e-08, 2.16050345227544e-07, -3.27382770920699e-07, 0, 0, 0};
    Matrix G(1, 6, values2, 6);
    double values3[] = {101453348.207917, 120429.109355752, 148186.145010685, 39372.9209771494, 3284.21674871861, 4014.15727499737,
                        120429.109355752, 101309543.076737, 84141.6477319108, 3284.34773346912, 35369.9224485894, 2255.66799443683,
                        148186.145010685, 84141.6477319108, 101344434.103716, 4014.41933457261, 2255.72532464628, 36274.7873567659,
                        39372.9209771494, 3284.34773346912, 4014.41933457261, 1001.21615369033, 1.32096249010467, 1.60455480925268,
                        3284.21674871861, 35369.9224485894, 2255.72532464628, 1.32096249010467, 999.576829597177, 0.892927374761907,
                        4014.15727499737, 2255.66799443683, 36274.7873567659, 1.60455480925268, 0.892927374761907, 999.924178045209};
    Matrix P(6, 6, values3, 36);
    int n = 6;
    Matrix K(6, 1);
    MeasUpdate(x, z, g, s, G, P, n, K);

    double values4[] = {582691.206462721, 1312775.30841773, -1989454.89979559,
                        190.367307357728, 433.242659422225, -660.433799251448};
    Matrix comp1(6, 1, values4, 6);
    double values5[] = {5736805.99700083, 3120008.83717262, 3733125.54856025,
                        5199.05810378403, -2475.74783768479, -7193.17204837758};
    Matrix comp2(6, 1, values5, 6);
    double values6[] = {95796502.307684, -12624173.0727019, 19462086.3185293, 37524.809129534, -921.76222363725, 10425.7388953049,
                        -12624173.0727019, 72596566.3325406, 43597431.3212774, -879.359516331218, 25894.0537871216, 16700.6535128224,
                        19462086.3185293, 43597431.3212774, 35401902.436512, 10324.3398033036, 16615.9994838993, 14384.0288752816,
                        37524.809129534, -879.359516331219, 10324.3398033036, 1000.61236892024, -0.053145927681848, 3.69924152397372,
                        -921.76222363725, 25894.0537871216, 16615.9994838993, -0.0531459276818478, 996.449599436078, 5.6600675709344,
                        10425.7388953049, 16700.6535128224, 14384.0288752816, 3.69924152397372, 5.6600675709344, 992.657163953104};
    Matrix comp3(6, 6, values6, 36);

    _assert(K.isEqual(comp1, 10e-9));
    _assert(x.isEqual(comp2, 10e-9));
    _assert(P.isEqual(comp3, 10e-7));

    return 0;
}

int DEInteg_01() {
    Global::AuxParam::Mjd_UTC = 4.974611128472211e+04;
    Global::AuxParam::n = 20;
    Global::AuxParam::m = 20;
    Global::AuxParam::sun = 1;
    Global::AuxParam::moon = 1;
    Global::AuxParam::planets = 1;

    double t = 0;
    double tout = -1.349999919533730e+02;
    double relerr = 1.00000000000000e-13;
    double abserr = 1.00000000000000e-06;
    int n_eqn = 6;
    double values1[] = {6221397.62857869, 2867713.77965738, 3006155.98509949,
                       4645.04725161807, -2752.21591588205, -7507.99940987033};
    Matrix y(6, 1, values1, 6);
    Matrix sol = DEInteg(Accel, t, tout, relerr, abserr, n_eqn, y);

    double values2[] = {5542555.93722861, 3213514.86734920, 3990892.97587686,
                        5394.06842166353, -2365.21337882342, -7061.84554200298};
    Matrix comp(6, 1, values2, 6);


    _assert(sol.isEqual(comp, 10e-9));

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
    _verify(JPL_Eph_DE430_01);
    _verify(Accel_01);
    _verify(LTC_01);
    _verify(MeasUpdate_01);
    _verify(DEInteg_01);

    return 0;
}


int main()
{
    Global::Pc();
    Global::GGM03S();
    Global::eop19620101(21413);

    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
