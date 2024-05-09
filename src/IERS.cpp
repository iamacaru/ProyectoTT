#include "../include/IERS.h"
#include "../include/SAT_Const.h"

#include <cmath>

void IERS(const Matrix& eop, double Mjd_UTC, char interp, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole,
          double& dy_pole, double& TAI_UTC) {

    if (interp == 'l') {
        double mjd = floor(Mjd_UTC);

        int i = 1;
        while (i <= eop.rows()) {
            if (mjd == eop(i, 4)) {
                break;
            }
            i++;
        }

        double values1[] = {eop(i, 1), eop(i, 2), eop(i, 3), eop(i, 4),
                           eop(i, 5), eop(i, 6), eop(i, 7), eop(i, 8),
                           eop(i, 9), eop(i, 10), eop(i, 11), eop(i, 12),
                           eop(i, 13), };
        Matrix preeop(1, 13, values1, 13);
        double values2[] = {eop(i+1, 1), eop(i+1, 2), eop(i+1, 3), eop(i+1, 4),
                           eop(i+1, 5), eop(i+1, 6), eop(i+1, 7), eop(i+1, 8),
                           eop(i+1, 9), eop(i+1, 10), eop(i+1, 11), eop(i+1, 12),
                           eop(i+1, 13), };
        Matrix nexteop(1, 13, values2, 13);

        double mfme = 1440.0 * (Mjd_UTC - floor(Mjd_UTC));
        double fixf = mfme / 1440.0;

        // Establecimiento de parámetros de rotación de la Tierra (IERS)
        x_pole  = (preeop(1,5) + (nexteop(1,5) - preeop(1,5)) * fixf) / Constants::Arcs;
        y_pole  = (preeop(1,6) + (nexteop(1,6) - preeop(1,6)) * fixf) / Constants::Arcs;
        UT1_UTC = preeop(1,7) + (nexteop(1,7) - preeop(1,7)) * fixf;
        LOD = preeop(1,8) + (nexteop(1,8) - preeop(1,8)) * fixf;
        dpsi = (preeop(1,9) + (nexteop(1,9) - preeop(1,9)) * fixf) / Constants::Arcs;
        deps = (preeop(1,10) + (nexteop(1,10) - preeop(1,10)) * fixf) / Constants::Arcs;
        dx_pole = (preeop(1,11) + (nexteop(1,11) - preeop(1,11)) * fixf) / Constants::Arcs;
        dy_pole = (preeop(1,12) + (nexteop(1,12) - preeop(1,12)) * fixf) / Constants::Arcs;
        TAI_UTC = preeop(1,13);

    } else if (interp == 'n'){
        double mjd = floor(Mjd_UTC);

        int i = 0;
        while (i <= eop.rows()) {
            if (mjd == eop(i, 4)) {
                break;
            }
            i++;
        }

        double values1[] = {eop(i, 1), eop(i, 2), eop(i, 3), eop(i, 4),
                            eop(i, 5), eop(i, 6), eop(i, 7), eop(i, 8),
                            eop(i, 9), eop(i, 10), eop(i, 11), eop(i, 12),
                            eop(i, 13), };
        Matrix eopAux(1, 13, values1, 13);

        x_pole = eopAux(i, 5) / Constants::Arcs;
        y_pole = eopAux(i, 6) / Constants::Arcs;
        UT1_UTC = eopAux(i, 7);
        LOD = eopAux(i, 8);
        dpsi = eopAux(i, 9) / Constants::Arcs;
        deps = eopAux(i, 10) / Constants::Arcs;
        dx_pole = eopAux(i, 11) / Constants::Arcs;
        dy_pole = eopAux(i, 12) / Constants::Arcs;
        TAI_UTC = eopAux(i, 13);
    }
}
