#ifndef PROYECTO_ACCEL_H
#define PROYECTO_ACCEL_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 Accel.m

 Purpose:
   Computes the acceleration of an Earth orbiting satellite due to
    - the Earth's harmonic gravity field,
    - the gravitational perturbations of the Sun and Moon
    - the solar radiation pressure and
    - the atmospheric drag

 Inputs:
   Mjd_TT      Terrestrial Time (Modified Julian Date)
   Y           Satellite state vector in the ICRF/EME2000 system

 Output:
   dY		    Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system

--------------------------------------------------------------------------*/

Matrix Accel(double x, Matrix& Y);

#endif