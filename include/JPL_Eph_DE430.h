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

void JPL_Eph_DE430(Matrix& r_Mercury, Matrix& r_Venus, Matrix& r_Earth, Matrix& r_Mars, Matrix& r_Jupiter,
                     Matrix& r_Saturn, Matrix& r_Uranus, Matrix& r_Neptune, Matrix& r_Pluto, Matrix& r_Moon,
                     Matrix& r_Sun, double Mjd_TDB);

#endif