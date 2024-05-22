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

/*!
 * @file NutMatrix.h
 * @brief Transformation from mean to true equator and equinox
 *
 * @param Mjd_TT    Modified Julian Date (Terrestrial Time)
 * @return Nutation matrix
 */

Matrix NutMatrix(double Mjd_TT);

#endif