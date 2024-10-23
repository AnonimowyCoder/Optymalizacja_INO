#include"user_funs.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = pow(x, 2);
	return y;
}
matrix ff1L(matrix x, matrix ud1, matrix ud2) {
	double PI = 3.14;
	matrix y = -cos(0.1 * x(0)) * exp(-pow(0.1 * x(0) - 2 * PI, 2)) + 0.002 * pow(0.1 * x(0), 2);
	return y;

}
matrix ff1R(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	matrix Y0 = matrix(3, new double[3] { 5,1,10 });
	matrix* Y = solve_ode(df1R, 0, 1, 2000, Y0, ud1, x);
	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (int i = 0; i < n; i++) {
		if (max < Y[1](i, 2))
			max = Y[1](i, 2);
	}

	y = abs(max - 50);
	return y;
}

matrix df1R(double t, matrix Y, matrix ud1, matrix ud2) {
	double a = 0.98;
	double b = 0.63;
	double g = 9.81;

	double P_A = 0.5;
	double V_A = 5.0;
	double T_A = 90;

	double P_B = 1;
	double D_B = 0.00365665;

	double F_in = 0.010;
	double T_in = 20;

	double outflow_A,outflow_B;

	if (Y(0) > 0) {
		outflow_A = a * b * m2d(ud2) * sqrt(2.0 * g * Y(1) / P_B);
	}
	else {
		outflow_A = 0.0;
	}
	if (Y(1) > 0) {
		outflow_B = a * b * D_B * sqrt(2.0 * g * Y(1) / P_B);
	}
	else {
		outflow_B = 0.0;
	}

	matrix dY=matrix(3, 1);
	dY(0) = -outflow_A;
	dY(1) = (outflow_A+F_in-outflow_B);
	dY(2)= (F_in / Y(1)) * (T_in - Y(2)) + (outflow_A / Y(1)) * (T_A - Y(2));
	return dY;
}