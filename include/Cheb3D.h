
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

/*!
 * @file Cheb3D.h
 * @brief Chebyshev approximation of 3-dimensional vectors
 *
 * @param t Value in the interval [Ta, Tb] where the approximation is evaluated
 * @param N Number of coefficients
 * @param Ta Begin interval
 * @param Tb End interval
 * @param Cx Coefficients of Chebyshev polyomial (x-coordinate)
 * @param Cy Coefficients of Chebyshev polyomial (y-coordinate)
 * @param Cz Coefficients of Chebyshev polyomial (z-coordinate)
 * @return Approximated 3-dimensional vector at position `t`
 */

Matrix Cheb3D(double t, int N, double Ta, double Tb, Matrix& Cx, Matrix& Cy, Matrix& Cz);

#endif