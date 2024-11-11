#include"user_funs.h"
#include <cmath>

#define M_PI 3.14159265358979323846

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

matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
	double arg = m2d(x);
	matrix result = -cos(0.1 * arg) * exp(-pow(0.1 * arg - 2 * 3.14, 2)) + 0.002 * pow(0.1 * arg, 2);
	return result;
}

// t = czas, Y = wektor stanu, zawiera aktualne wartosci: Y(0) = Va, Y(1)= Vb, Y(2) = Tb
matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
	//Dane z konspektu
	double a = 0.98; //lepkosc
	double b = 0.63; //zwezenie strumienia cieczy
	double g = 9.81; //przyspieszenie ziemskie

	double tA = 90.0; //temperatura w A

	double Pa = 0.5;
	double Pb = 1.0;
	double Db = 0.00365665;

	double fOutA = 0, fOutB = 0;	//Szybkosc wyplywu wody
	double fInB = 0.01, tInB = 20.0;	//Woda wplywajaca do zbiornika 0.01m3/s, 20C
	matrix dY(3, 1);	// Wektor jest 3x1 : dVA, dVB, dTB

	if (Y(0) > 0) {
		fOutA = a * b * m2d(ud2) * sqrt(2.0 * g * Y(0) / Pa);	//szybkosc wyplywu z A 
	}
	else {
		fOutA = 0.0;
	}
	if (Y(1) > 0) {
		fOutB = a * b * Db * sqrt(2.0 * g * Y(1) / Pb);	//szybkosc wyplywu z B
	}
	else {
		fOutB = 0.0;
	}

	dY(0) = -fOutA;	 //woda wyplywa z VA wiec na -
	dY(1) = (fOutA + fInB - fOutB); //z A wplywa i z B wyplywa
	dY(2) = (fInB / Y(1)) * (tInB - Y(2)) + (fOutA / Y(1)) * (tA - Y(2));	//zmiana temperatury w B
	return dY;
}


matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(3, new double[3] {5, 1, 20});	//Poczatkowy stan ukladu Va=5, Vb=1,Tb=20
	matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, ud1, x); // rozniczka

	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (int i = 1; i < n; ++i)
		if (max < Y[1](i, 2)) max = Y[1](i, 2);

	y = abs(max - 50); //wartoc bezwzgledna miedzy max w zbiorniku a 50
	return y;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2)
{
	double x1 = x(0);
	double x2 = x(1);
	return pow(x1, 2) + pow(x2, 2) - cos(2.5 * M_PI * x1) - cos(2.5 * M_PI * x2) + 2;
}

//matrix ff2R(matrix x, matrix ud1, matrix ud2) {}

//matrix df2(){}