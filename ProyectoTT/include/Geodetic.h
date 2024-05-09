#ifndef PROYECTO_GEODETIC_H
#define PROYECTO_GEODETIC_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 Geodetic.m

 Purpose:
   geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
   from given position vector (r [m])

--------------------------------------------------------------------------*/

Matrix Geodetic(const Matrix& r);

#endif