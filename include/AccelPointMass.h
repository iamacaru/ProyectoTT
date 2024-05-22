#ifndef PROYECTO_ACCELPOINTMASS_H
#define PROYECTO_ACCELPOINTMASS_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 AccelPointMass: Computes the perturbational acceleration due to a point
				  mass

 Inputs:
   r           Satellite position vector
   s           Point mass position vector
   GM          Gravitational coefficient of point mass

 Output:
   a    		Acceleration (a=d^2r/dt^2)

--------------------------------------------------------------------------*/

/*!
 * @file AccelPointMass.h
 * @brief Computes the perturbational acceleration due to a point mass
 *
 * @param r Satellite position vector
 * @param s Point mass position vector
 * @param GM Gravitational coefficient of point mass
 * @return Acceleration (a=d^2r/dt^2)
 */

Matrix AccelPointMass(Matrix& r, Matrix& s, double GM);

#endif