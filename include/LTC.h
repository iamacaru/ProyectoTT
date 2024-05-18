#ifndef PROYECTO_LTC_H
#define PROYECTO_LTC_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 LTC.m

 Purpose:
   Transformation from Greenwich meridian system to
   local tangent coordinates

 Inputs:
   lon      -Geodetic East longitude [rad]
   lat      -Geodetic latitude [rad]

 Output:
   M        -Rotation matrix from the Earth equator and Greenwich meridian
             to the local tangent (East-North-Zenith) coordinate system

--------------------------------------------------------------------------*/

Matrix  LTC(double lon, double lat);

#endif