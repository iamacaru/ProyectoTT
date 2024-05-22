#ifndef PROYECTO_TIMEUPDATE_H
#define PROYECTO_TIMEUPDATE_H

#include "Matrix.h"

/*!
 * @file TimeUpdate.h
 * @brief Performs the time update step in a Kalman filter.
 *
 * @param P The process covariance matrix
 * @param Phi The state transition matrix
 * @param Qdt The process noise covariance
 * @return The updated process covariance matrix
 */

Matrix TimeUpdate(Matrix&  P, Matrix&  Phi, double Qdt = 0.0);

#endif