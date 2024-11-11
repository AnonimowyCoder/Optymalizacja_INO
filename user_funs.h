#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix df1(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1R(matrix, matrix = NAN, matrix = NAN);

matrix ff2T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff2R(matrix x, matrix ud1=NAN, matrix ud2=NAN);
matrix df2();