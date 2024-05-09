#ifndef PROYECTO_POSITION_H
#define PROYECTO_POSITION_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 Position.m

 Purpose:
   Position vector (r [m]) from geodetic coordinates (Longitude [rad],
   latitude [rad], altitude [m])

--------------------------------------------------------------------------*/

Matrix Position(double lon, double lat, double alt);

#endif