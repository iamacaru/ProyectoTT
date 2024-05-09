#include "../include/MeanObliquity.h"
#include "../include/SAT_Const.h"

double MeanObliquity(double Mjd_TT) {
    double T = (Mjd_TT - Constants::MJD_J2000) / 36525.0;

    return Constants::Rad * (84381.448/3600-(46.8150+(0.00059-0.001813*T)*T)*T/3600);
}
