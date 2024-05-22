#ifndef PROYECTO_IERS_H
#define PROYECTO_IERS_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 IERS: Management of IERS time and polar motion data

--------------------------------------------------------------------------*/

/*!
 * @file IERS.h
 * @brief Computes various Earth orientation parameters using IERS data
 *
 * @param eop Earth Orientation Parameters matrix containing data from IERS
 * @param Mjd_UTC Modified Julian Date in Coordinated Universal Time (UTC)
 * @param interp Interpolation flag ('N' for no interpolation, 'L' for linear interpolation)
 * @param x_pole Output parameter for polar motion in the x-direction (arcseconds)
 * @param y_pole Output parameter for polar motion in the y-direction (arcseconds)
 * @param UT1_UTC Output parameter for the difference between Universal Time (UT1) and Coordinated Universal Time (UTC) (seconds)
 * @param LOD Output parameter for the Length of Day (seconds)
 * @param dpsi Output parameter for the nutation correction in longitude (arcseconds)
 * @param deps Output parameter for the nutation correction in obliquity (arcseconds)
 * @param dx_pole Output parameter for the rate of change of polar motion in the x-direction (arcseconds/day)
 * @param dy_pole Output parameter for the rate of change of polar motion in the y-direction (arcseconds/day)
 * @param TAI_UTC Output parameter for the difference between International Atomic Time (TAI) and UTC (seconds)
 */

void IERS(const Matrix& eop, double Mjd_UTC, char interp, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD,
          double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC);

#endif