#ifndef PROYECTO_R_Y_H
#define PROYECTO_R_Y_H

#include "Matrix.h"

/*--------------------------------------------------------------------------
  input:
    angle       - angle of rotation [rad]

  output:
    rotmat      - vector result
--------------------------------------------------------------------------*/

/*!
 * @file R_y.h
 * @brief Computes the rotation matrix for a rotation about the y-axis
 *
 * @param alpha angle of rotation [rad]
 * @return vector result
 */

Matrix R_y(double alpha);

#endif
