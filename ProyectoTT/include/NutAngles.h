#ifndef PROYECTO_NUTANGLES_H
#define PROYECTO_NUTANGLES_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 NutAngles.m

 Purpose:
   Nutation in longitude and obliquity

 Input:
   Mjd_TT     Modified Julian Date (Terrestrial Time)

 Output:
   dpsi,deps  Nutation Angles

--------------------------------------------------------------------------*/

Matrix NutAngles(double Mjd_TT);

#endif