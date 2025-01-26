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
matrix df2(double t, matrix Y, matrix ud1 = NAN, matrix ud2= NAN);

matrix ff3T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff3R(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix df3(double t, matrix Y, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff4T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix gf4T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix hf4T(matrix x, matrix ud1, matrix ud2);
matrix ff4R(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix gf4(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff5T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff5T_1(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff5T_2(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff5T_comb(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix trivial(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff5R(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

//LAB 6
matrix ff6T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix df6(double t, matrix Y, matrix ud1 =NAN, matrix ud2=NAN);
matrix ff6R(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);