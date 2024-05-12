#ifndef PROYECTO_TIMEUPDATE_H
#define PROYECTO_TIMEUPDATE_H

#include "Matrix.h"

Matrix TimeUpdate(Matrix&  P, Matrix&  Phi, double Qdt = 0.0);

#endif