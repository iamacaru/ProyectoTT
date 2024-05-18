cmake_minimum_required(VERSION 3.23)
project(Proyecto)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

add_executable(Proyecto
        include/Matrix.h
        include/R_x.h
        src/Matrix.cpp
        src/R_x.cpp
        EKF_Test.cpp include/R_y.h src/R_y.cpp include/R_z.h src/R_z.cpp include/SAT_Const.h include/Global.h src/Global.cpp include/Sign_.h src/Sign_.cpp src/Timediff.cpp include/Timediff.h
        include/Unit.h
        src/Unit.cpp
        include/AccelPointMass.h
        src/AccelPointMass.cpp
        include/AzElPa.h
        src/AzElPa.cpp
        include/Cheb3D.h
        src/Cheb3D.cpp
        include/EccAnom.h
        src/EccAnom.cpp
        include/Frac.h
        src/Frac.cpp
        include/Geodetic.h
        src/Geodetic.cpp
        include/Legendre.h
        src/Legendre.cpp
        include/MeanObliquity.h
        src/MeanObliquity.cpp
        include/Mjday.h
        src/Mjday.cpp
        include/Mjday_TDB.h
        src/Mjday_TDB.cpp
        include/NutAngles.h
        src/NutAngles.cpp
        include/Position.h
        src/Position.cpp
        include/IERS.h
        src/IERS.cpp
        EKF_GEOS3.cpp
        include/EqnEquinox.h
        src/EqnEquinox.cpp
        include/Gmst.h
        src/Gmst.cpp
        include/Gast.h
        src/Gast.cpp
        include/GHAMatrix.h
        src/GHAMatrix.cpp
        include/NutMatrix.h
        src/NutMatrix.cpp
        include/PoleMatrix.h
        src/PoleMatrix.cpp
        include/PrecMatrix.h
        src/PrecMatrix.cpp
        include/TimeUpdate.h
        src/TimeUpdate.cpp
        include/VarEqn.h
        src/VarEqn.cpp
        include/AccelHarmonic.h
        src/AccelHarmonic.cpp
        include/G_AccelHarmonic.h
        src/G_AccelHarmonic.cpp
        include/JPL_Eph_DE430.h
        src/JPL_Eph_DE430.cpp
        include/Accel.h src/Accel.cpp
        include/LTC.h
        src/LTC.cpp
        include/MeasUpdate.h
        src/MeasUpdate.cpp
        include/DEInteg.h
        src/DEInteg.cpp)
