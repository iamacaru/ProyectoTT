#ifndef PROYECTO_R_X_H
#define PROYECTO_R_X_H

#include "Matrix.h"

/*--------------------------------------------------------------------------
  input:
    angle       - angle of rotation [rad]

  output:
    rotmat      - vector result
--------------------------------------------------------------------------*/

/*!
 * @file R_x.h
 * @brief Computes the rotation matrix for a rotation about the x-axis
 *
 * @param alpha angle of rotation [rad]
 * @return vector result
 */

Matrix R_x(double alpha);

#endif
