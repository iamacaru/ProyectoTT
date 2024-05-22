#ifndef PROYECTO_EQNEQUINOX_H
#define PROYECTO_EQNEQUINOX_H

/*--------------------------------------------------------------------------

 EqnEquinox.m

 Purpose:
   Computation of the equation of the equinoxes

 Input:
   Mjd_TT    Modified Julian Date (Terrestrial Time)

 Output:
    EqE      Equation of the equinoxes

 Notes:
   The equation of the equinoxes dpsi*cos(eps) is the right ascension of
   the mean equinox referred to the true equator and equinox and is equal
   to the difference between apparent and mean sidereal time.

--------------------------------------------------------------------------*/

/*!
 * @file EqnEquinox.h
 * @brief Computation of the equation of the equinoxes
 * @details The equation of the equinoxes dpsi*cos(eps) is the right ascension of the mean equinox referred to the true equator and equinox and is equal to the difference between apparent and mean sidereal time.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return Equation of the equinoxes
 */

double EqnEquinox(double Mjd_TT);

#endif