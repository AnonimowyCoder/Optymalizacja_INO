#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix ff1R(matrix, matrix = NAN, matrix = NAN);
matrix df1(double, matrix, matrix = NAN, matrix = NAN);
matrix df2(double, matrix, matrix = NAN, matrix = NAN);
matrix ff2T(matrix, matrix = NAN, matrix = NAN);
matrix ff2R(matrix, matrix = NAN, matrix = NAN);
matrix fun3(matrix, matrix =NAN, matrix =NAN);
matrix fR3(matrix, matrix ud1=NAN, matrix ud2 = NAN);
matrix df3(double, matrix, matrix ud1=NAN, matrix ud2=NAN);
matrix fun4(matrix x, matrix ud1, matrix ud2);
matrix grad4(matrix x, matrix ud1, matrix ud2);
matrix hesj4(matrix x, matrix ud1, matrix ud2);

matrix fT4(matrix x, matrix ud1, matrix ud2);

matrix fR4(matrix x, matrix ud1, matrix ud2);
matrix gf(matrix x, matrix ud1, matrix ud2);
matrix ff5T(matrix x, matrix ud1=NAN, matrix ud2=NAN);