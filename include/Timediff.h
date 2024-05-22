#ifndef PROYECTO_TIMEDIFF_H
#define PROYECTO_TIMEDIFF_H

#include "Matrix.h"

/*--------------------------------------------------------------------------

 Time differences [s]

--------------------------------------------------------------------------*/

/*!
 * @file Timediff.h
 * @brief Computes various time differences
 *
 * @param UT1_UTC Difference between UT1 and UTC time [s]
 * @param TAI_UTC Difference between TAI and UTC time [s]
 * @return A matrix containing the computed time differences
 */

Matrix timediff(double UT1_UTC, double TAI_UTC);

#endif
