#ifndef PROYECTO_ECCANOM_H
#define PROYECTO_ECCANOM_H

/*--------------------------------------------------------------------------

 Purpose:
   Computes the eccentric anomaly for elliptic orbits

 Inputs:
   M         Mean anomaly in [rad]
   e         Eccentricity of the orbit [0,1]

 Output:
             Eccentric anomaly in [rad]

--------------------------------------------------------------------------*/

/*!
 * @file EccAnom.h
 * @brief Computes the eccentric anomaly for elliptic orbits
 *
 * @param M Mean anomaly in [rad]
 * @param e Eccentricity of the orbit [0,1]
 * @return Eccentric anomaly in [rad]
 */

double EccAnom(double M, double e);

#endif