#ifndef PROYECTO_IERS_H
#define PROYECTO_IERS_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 IERS: Management of IERS time and polar motion data

--------------------------------------------------------------------------*/

void IERS(const Matrix& eop, double Mjd_UTC, char interp, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD,
          double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC);

#endif