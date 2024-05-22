#ifndef PROYECTO_PRECMATRIX_H
#define PROYECTO_PRECMATRIX_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 PrecMatrix.m

 Purpose:

   Precession transformation of equatorial coordinates

 Input:
   Mjd_1     Epoch given (Modified Julian Date TT)
   MjD_2     Epoch to precess to (Modified Julian Date TT)

 Output:
   PrecMat   Precession transformation matrix

--------------------------------------------------------------------------*/

/*!
 * @file PrecMatrix.h
 * @brief Precession transformation of equatorial coordinates
 *
 * @param Mjd_1 Epoch given (Modified Julian Date TT)
 * @param MjD_2 Epoch to precess to (Modified Julian Date TT)
 * @return Precession transformation matrix
 */

Matrix PrecMatrix(double Mjd_1, double Mjd_2);

#endif