#ifndef PROYECTO_AZELPA_H
#define PROYECTO_AZELPA_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 Purpose:
  Computes azimuth, elevation and partials from local tangent coordinates

 Input:
   s      Topocentric local tangent coordinates (East-North-Zenith frame)

 Outputs:
   A      Azimuth [rad]
   E      Elevation [rad]
   dAds   Partials of azimuth w.r.t. s
   dEds   Partials of elevation w.r.t. s

--------------------------------------------------------------------------*/

/*!
 * @file AzElPa.h
 * @brief Computes azimuth, elevation and partials from local tangent coordinates
 *
 * @param s Topocentric local tangent coordinates (East-North-Zenith frame)
 * @param A Azimuth [rad]
 * @param E Elevation [rad]
 * @param dAds Partials of azimuth w.r.t. s
 * @param dEds Partials of elevation w.r.t. s
 */

void  AzElPa(const Matrix& s, double& Az, double& El, Matrix& dAds, Matrix& dEds);

#endif