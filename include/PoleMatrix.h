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

/*!
 * @file PoleMatrix.h
 * @brief Transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date
 *
 * @param xp Pole coordinte xp
 * @param yp Pole coordinte yp
 * @return Pole matrix
 */

Matrix PoleMatrix(double xp, double yp);

#endif