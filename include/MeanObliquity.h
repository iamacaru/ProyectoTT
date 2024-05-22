#ifndef PROYECTO_MEANOBLIQUITY_H
#define PROYECTO_MEANOBLIQUITY_H

/*--------------------------------------------------------------------------

 MeanObliquity.m

 Purpose:
   Computes the mean obliquity of the ecliptic

 Input:
   Mjd_TT    Modified Julian Date (Terrestrial Time)

 Output:
   MOblq     Mean obliquity of the ecliptic [rad]

--------------------------------------------------------------------------*/

/*!
 * @file MeanObliquity.h
 * @brief Computes the mean obliquity of the ecliptic
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return Mean obliquity of the ecliptic [rad]
 */


double MeanObliquity(double Mjd_TT);

#endif