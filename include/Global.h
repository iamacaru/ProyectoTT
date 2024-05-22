#ifndef PROJECT_GLOBAL_H
#define PROJECT_GLOBAL_H

#include "Matrix.h"
#include <cstdio>
#include <cstdlib>

/*!
 * @file Global.h
 * @brief Class containing global variables and functions.
 */
class Global {
public:
    //! Pointer to Earth Orientation Parameters data.
    static Matrix *eopdata;

    /*!
     * @brief Function to initialize Earth Orientation Parameters for a specific date.
     *
     * @param f Flag indicating the date for which Earth Orientation Parameters are needed.
     */
    static void eop19620101(int f);

    //! Coefficients of the Geopotential Model.
    static Matrix *Cnm;
    static Matrix *Snm;

    //! Function to load the GGM03S gravity model coefficients.
    static void GGM03S();

    //! Pointer to the planetary constants data.
    static Matrix *PC;

    //! Function to load Planetary Constants.
    static void Pc();

    //! Modified Julian Date at epoch.
    double Mjd0;

    //! Structure containing auxiliary parameters.
    struct AuxParam {
        //! Modified Julian Date in UTC.
        static double Mjd_UTC;
        //! Modified Julian Date in Terrestrial Time.
        static double Mjd_TT;
        //! Maximum degree of the gravitational model.
        static int n;
        //! Maximum order of the gravitational model.
        static int m;
        //! Flag indicating whether to include Sun in the gravity model.
        static int sun;
        //! Flag indicating whether to include Moon in the gravity model.
        static int moon;
        //! Flag indicating whether to include other planets in the gravity model.
        static int planets;
    };
};

#endif
