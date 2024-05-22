#ifndef PROYECTO_DEINTEG_H
#define PROYECTO_DEINTEG_H

#include "Matrix.h"

/*----------------------------------------------------------------------------

 Purpose:
   Numerical integration methods for ordinaray differential equations

   This module provides implemenation of the variable order variable
   stepsize multistep method of Shampine & Gordon.

 Reference:

   Shampine, Gordon: "Computer solution of Ordinary Differential Equations",
   Freeman and Comp., San Francisco (1975).

----------------------------------------------------------------------------*/

/*!
 * @file DEInteg.h
 * @brief Numerical integration methods for ordinaray differential equations
 * This module provides implemenation of the variable order variable stepsize multistep method of Shampine & Gordon.
 *
 * @param func Function that computes the derivatives of the system state vector.
 * @param t Initial time.
 * @param tout Time at which the solution is desired.
 * @param relerr Relative error tolerance.
 * @param abserr Absolute error tolerance.
 * @param n_eqn Number of equations in the system.
 * @param y State vector at the initial time `t`, updated to the state vector at `tout`.
 * @return State vector at the time `tout`.
 */

Matrix DEInteg(Matrix (*func)(double, Matrix&), double t, double tout, double relerr, double abserr, int n_eqn, Matrix& y);

#endif