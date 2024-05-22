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

/*!
 * @file Unit.h
 * @brief This function calculates a unit vector given the original vector. If a zero vector is input, the vector is set to zero
 *
 * @param vec vector
 * @return unit vector
 */

Matrix unit(const Matrix& vec);

#endif