#ifndef PROYECTO_POSITION_H
#define PROYECTO_POSITION_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 Position.m

 Purpose:
   Position vector (r [m]) from geodetic coordinates (Longitude [rad],
   latitude [rad], altitude [m])

--------------------------------------------------------------------------*/

/*!
 * @file Position.h
 * @brief Position vector (r [m]) from geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
 *
 * @param lon Longitude in radians.
 * @param lat Latitude in radians.
 * @param alt Altitude in meters.
 * @return A Matrix object containing the position vector (r [m]) in ECEF coordinates.
 */

Matrix Position(double lon, double lat, double alt);

#endif