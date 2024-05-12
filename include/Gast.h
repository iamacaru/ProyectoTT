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

double gast(double Mjd_UT1);

#endif