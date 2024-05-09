#ifndef PROYECTO_UNIT_H
#define PROYECTO_UNIT_H

#include "Matrix.h"


/*--------------------------------------------------------------------------

  this function calculates a unit vector given the original vector. if a
  zero vector is input, the vector is set to zero.

  input:
    vec         - vector

  output:
    outvec      - unit vector

--------------------------------------------------------------------------*/

Matrix unit(const Matrix& vec);

#endif