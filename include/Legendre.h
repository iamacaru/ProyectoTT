#ifndef PROYECTO_LEGENDRE_H
#define PROYECTO_LEGENDRE_H

#include "Matrix.h"

/*!
 * @file Legendre.h
 * @brief Computes the associated Legendre polynomials and their derivatives
 *
 * @param n Maximum degree of the Legendre polynomials
 * @param m Maximum order of the Legendre polynomials
 * @param fi Geocentric latitude [rad]
 * @param pnm Matrix of associated Legendre polynomials (output)
 * @param dpnm Matrix of derivatives of associated Legendre polynomials (output)
 */

void Legendre(int n, int m, double fi, Matrix& pnm, Matrix& dpnm);

#endif