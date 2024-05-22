#ifndef PROYECTO_MEASUPDATE_H
#define PROYECTO_MEASUPDATE_H

#include "Matrix.h"

/*!
 * @file MeasUpdate.h
 * @brief Performs the measurement update step for a Kalman Filter
 *
 * @param x State vector
 * @param z Measurement vector
 * @param g Measurement function
 * @param s Covariance of measurement noise
 * @param G Sensitivity matrix
 * @param P Error covariance matrix
 * @param n Number of measurements
 * @param K Kalman gain (output)
 */

void MeasUpdate(Matrix& x, Matrix& z, Matrix& g, Matrix& s, Matrix& G, Matrix& P, int n, Matrix& K);

#endif