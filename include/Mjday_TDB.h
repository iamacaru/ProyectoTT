#ifndef PROYECTO_MJDAY_TDB_H
#define PROYECTO_MJDAY_TDB_H

/*--------------------------------------------------------------------------

 Mjday_TDB: Computes the Modified Julian Date for barycentric dynamical
            time

  Inputs:
    Mjd_TT      - Modified julian date (TT)

  Output:
    Mjd_TDB     - Modified julian date (TDB)

 Reference:
 Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill;
 New York; 3rd edition(2007).

 Last modified:   2015/08/12   M. Mahooti

--------------------------------------------------------------------------*/

/*!
 * @file Mjday_TDB.h
 * @brief Computes the Modified Julian Date for barycentric dynamical time
 *
 * @param Mjd_TT Modified julian date (TT)
 * @return Mjd_TDB Modified julian date (TDB)
 */

double Mjday_TDB(double Mjd_TT);

#endif