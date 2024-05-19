#ifndef PROYECTO_GLOBAL_H
#define PROYECTO_GLOBAL_H

#include "Matrix.h"
#include <cstdio>
#include <cstdlib>

class Global {
public:
    static Matrix *eopdata;
    static void eop19620101 (int f);

    static Matrix *Cnm;
    static Matrix *Snm;
    static void GGM03S ();

    static Matrix *PC;
    static void Pc ();

    double Mjd0;

    struct AuxParam {
        static double Mjd_UTC;
        static double Mjd_TT;
        static int n;
        static int m;
        static int sun;
        static int moon;
        static int planets;
    };
};

#endif