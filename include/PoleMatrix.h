#ifndef PROYECTO_POLEMATRIX_H
#define PROYECTO_POLEMATRIX_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 PoleMatrix.m

 Purpose:
   Transformation from pseudo Earth-fixed to Earth-fixed coordinates
   for a given date

 Input:
   Pole coordinte(xp,yp)

 Output:
   PoleMat   Pole matrix

--------------------------------------------------------------------------*/

Matrix PoleMatrix(double xp, double yp);

#endif