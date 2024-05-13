#ifndef PROYECTO_G_ACCELHARMONIC_H
#define PROYECTO_G_ACCELHARMONIC_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 G_AccelHarmonic.m

 Purpose:
   Computes the gradient of the Earth's harmonic gravity field

 Inputs:
   r           Satellite position vector in the true-of-date system
   U           Transformation matrix to body-fixed system
   n           Gravity model degree
   m 			Gravity model order

 Output:
   G    		Gradient (G=da/dr) in the true-of-date system

--------------------------------------------------------------------------*/

Matrix G_AccelHarmonic(Matrix& r, Matrix& U, int n_max, int m_max);

#endif