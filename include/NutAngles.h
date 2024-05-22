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

/*!
 * @file NutAngles.h
 * @brief Computes the nutation in longitude and obliquity
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return A Matrix object containing the nutation in longitude and obliquity in radians
 */

Matrix NutAngles(double Mjd_TT);

#endif