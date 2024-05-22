#ifndef PROYECTO_GEODETIC_H
#define PROYECTO_GEODETIC_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 Geodetic.m

 Purpose:
   geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
   from given position vector (r [m])

--------------------------------------------------------------------------*/

/*!
 * @file Geodetic.h
 * @brief geodetic coordinates (Longitude [rad], latitude [rad], altitude [m]) from given position vector (r [m])
 *
 * @param r Position vector in meters.
 * @return A Matrix containing the geodetic coordinates: [Longitude, Latitude, Altitude].
 */

Matrix Geodetic(const Matrix& r);

#endif