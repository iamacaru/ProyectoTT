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

Matrix DEInteg(Matrix (*func)(double, Matrix&), double t, double tout, double relerr, double abserr, int n_eqn, Matrix& y);

#endif