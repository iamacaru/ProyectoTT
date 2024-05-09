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

double Mjday(int yr, int mon, int day, int hr = 0, int min = 0, double sec = 0.0);

#endif