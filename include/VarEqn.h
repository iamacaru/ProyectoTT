#ifndef PROYECTO_VAREQN_H
#define PROYECTO_VAREQN_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

  Initial Orbit Determination using Gauss and Extended Kalman Filter methods

 References:
   O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
   Applications", Springer Verlag, Heidelberg, 2000.

   D. Vallado, "Fundamentals of Astrodynamics and Applications",
   4th Edition, 2013.

   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003.

 Last modified:   2020/03/16   Meysam Mahooti
--------------------------------------------------------------------------*/

Matrix VarEqn(double x, Matrix& yPhi);

#endif