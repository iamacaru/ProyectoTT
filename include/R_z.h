#ifndef PROYECTO_R_Z_H
#define PROYECTO_R_Z_H

#include "Matrix.h"

/*--------------------------------------------------------------------------
  input:
    angle       - angle of rotation [rad]

  output:
    rotmat      - vector result
--------------------------------------------------------------------------*/

/*!
 * @file R_z.h
 * @brief Computes the rotation matrix for a rotation about the z-axis
 *
 * @param alpha angle of rotation [rad]
 * @return vector result
 */

Matrix R_z(double alpha);

#endif
