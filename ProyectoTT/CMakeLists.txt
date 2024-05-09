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
        src/EKF_GEOS3.cpp
)
