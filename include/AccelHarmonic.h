#ifndef PROYECTO_ACCELHARMONIC_H
#define PROYECTO_ACCELHARMONIC_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 AccelHarmonic.m

 Purpose:
   Computes the acceleration due to the harmonic gravity field of the
   central body

 Inputs:
   r           Satellite position vector in the inertial system
   E           Transformation matrix to body-fixed system
   n_max       Maximum degree
   m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)

 Output:
   a           Acceleration (a=d^2r/dt^2)

 Last modified:   2015/08/12   M. Mahooti

--------------------------------------------------------------------------*/

Matrix AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max);

#endif