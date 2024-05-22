#ifndef PROYECTO_JPL_EPH_DE430_H
#define PROYECTO_JPL_EPH_DE430_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 JPL_Eph_DE430: Computes the sun, moon, and nine major planets' equatorial
                position using JPL Ephemerides

 Inputs:
   Mjd_TDB         Modified julian date of TDB

 Output:
   r_Earth(solar system barycenter (SSB)),r_Mars,r_Mercury,r_Venus,
   r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,
   r_Sun(geocentric equatorial position ([m]) referred to the
   International Celestial Reference Frame (ICRF))

 Notes: Light-time is already taken into account

--------------------------------------------------------------------------*/

/*!
 * @file JPL_Eph_DE430.h
 * @brief Computes the sun, moon, and nine major planets' equatorial position using JPL Ephemerides
 *
 * @param r_Mercury Position vector of Mercury (output)
 * @param r_Venus Position vector of Venus (output)
 * @param r_Earth Position vector of Earth (solar system barycenter) (output)
 * @param r_Mars Position vector of Mars (output)
 * @param r_Jupiter Position vector of Jupiter (output)
 * @param r_Saturn Position vector of Saturn (output)
 * @param r_Uranus Position vector of Uranus (output)
 * @param r_Neptune Position vector of Neptune (output)
 * @param r_Pluto Position vector of Pluto (output)
 * @param r_Moon Position vector of Moon (output)
 * @param r_Sun Position vector of Sun (output)
 * @param Mjd_TDB Modified julian date of TDB
 * @return
 */

void JPL_Eph_DE430(Matrix& r_Mercury, Matrix& r_Venus, Matrix& r_Earth, Matrix& r_Mars, Matrix& r_Jupiter,
                     Matrix& r_Saturn, Matrix& r_Uranus, Matrix& r_Neptune, Matrix& r_Pluto, Matrix& r_Moon,
                     Matrix& r_Sun, double Mjd_TDB);

#endif