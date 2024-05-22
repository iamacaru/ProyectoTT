#ifndef PROYECTO_GAST_H
#define PROYECTO_GAST_H

/*--------------------------------------------------------------------------

 GAST.m

 Purpose:
   Greenwich Apparent Sidereal Time

 Input:
   Mjd_UT1   Modified Julian Date UT1

 Output:
   gstime    GAST in [rad]

--------------------------------------------------------------------------*/

/*!
 * @file Gast.h
 * @brief Greenwich Apparent Sidereal Time
 *
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return GAST in [rad]
 */

double gast(double Mjd_UT1);

#endif