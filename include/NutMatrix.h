#ifndef PROYECTO_NUTMATRIX_H
#define PROYECTO_NUTMATRIX_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 NutMatrix.m

 Purpose:
   Transformation from mean to true equator and equinox

 Input:
   Mjd_TT    Modified Julian Date (Terrestrial Time)

 Output:
   NutMat    Nutation matrix

--------------------------------------------------------------------------*/

Matrix NutMatrix(double Mjd_TT);

#endif