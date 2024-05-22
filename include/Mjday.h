#ifndef PROYECTO_MJDAY_H
#define PROYECTO_MJDAY_H

/*--------------------------------------------------------------------------
  inputs:
    year        - year
    mon         - month
    day         - day
    hr          - universal time hour
    min         - universal time min
    sec         - universal time sec

  output:
    Mjd         - Modified julian date
--------------------------------------------------------------------------*/

/*!
 * @file Mjday.h
 * @brief Computes the Modified Julian Date (MJD) from calendar date and time
 *
 * @param year year
 * @param mon month
 * @param day day
 * @param hr universal time hour
 * @param min universal time min
 * @param sec universal time sec
 * @return Modified julian date
 */

double Mjday(int yr, int mon, int day, int hr = 0, int min = 0, double sec = 0.0);

#endif