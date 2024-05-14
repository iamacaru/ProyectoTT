
#ifndef PROYECTO_CHEB3D_H
#define PROYECTO_CHEB3D_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 Chebyshev approximation of 3-dimensional vectors

 Inputs:
     N       Number of coefficients
     Ta      Begin interval
     Tb      End interval
     Cx      Coefficients of Chebyshev polyomial (x-coordinate)
     Cy      Coefficients of Chebyshev polyomial (y-coordinate)
     Cz      Coefficients of Chebyshev polyomial (z-coordinate)

%--------------------------------------------------------------------------*/

Matrix Cheb3D(double t, int N, double Ta, double Tb, Matrix& Cx, Matrix& Cy, Matrix& Cz);

#endif