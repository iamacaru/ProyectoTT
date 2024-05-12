#ifndef PROYECTO_GMST_H
#define PROYECTO_GMST_H

/*--------------------------------------------------------------------------

 Purpose:
   Greenwich Mean Sidereal Time

 Input:
   Mjd_UT1    Modified Julian Date UT1

 Output:
   gmstime	   GMST in [rad]

--------------------------------------------------------------------------*/

double gmst(double Mjd_UT1);

#endif