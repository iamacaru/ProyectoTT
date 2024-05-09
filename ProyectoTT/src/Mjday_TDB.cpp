#include "../include/Mjday_TDB.h"

#include <cmath>

double Mjday_TDB(double Mjd_TT) {
    // Compute Julian Centureis of TT
    double T_TT = (Mjd_TT - 51544.5) / 36525.0;

    // Compute Modified Julian Date of TDB
    double Mjd_TDB = Mjd_TT + (0.001658 * std::sin(628.3076 * T_TT + 6.2401) +
                              0.000022 * std::sin(575.3385 * T_TT + 4.2970) +
                              0.000014 * std::sin(1256.6152 * T_TT + 6.1969) +
                              0.000005 * std::sin(606.9777 * T_TT + 4.0212) +
                              0.000005 * std::sin(52.9691 * T_TT + 0.4444) +
                              0.000002 * std::sin(21.3299 * T_TT + 5.5431) +
                              0.000010 * std::sin(628.3076 * T_TT + 4.2490)) / 86400.0;

    return Mjd_TDB;
}
