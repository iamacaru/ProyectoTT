#ifndef PROYECTO_GHAMATRIX_H
#define PROYECTO_GHAMATRIX_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 GHAMatrix.m

 Purpose:
   Transformation from true equator and equinox to Earth equator and
   Greenwich meridian system

 Input:
   Mjd_UT1   Modified Julian Date UT1

 Output:
   GHAmat    Greenwich Hour Angle matrix

--------------------------------------------------------------------------*/

Matrix GHAMatrix(double Mjd_UT1);

#endif