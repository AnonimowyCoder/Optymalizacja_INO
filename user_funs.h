#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix ff2T(matrix, matrix, matrix);

matrix ff3T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff3R(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix df3(double t, matrix Y, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff4T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix gf4T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix hf4T(matrix x, matrix ud1, matrix ud2);

matrix ff5T_1(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff5T_2(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff5T_comb(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix trivial(matrix x, matrix ud1=NAN, matrix ud2=NAN);