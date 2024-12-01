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

//LAB 1

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

//LAB 2

matrix ff2T(matrix x, matrix ud1, matrix ud2)
{
	double x1 = x(0);
	double x2 = x(1);
	return pow(x1, 2) + pow(x2, 2) - cos(2.5 * M_PI * x1) - cos(2.5 * M_PI * x2) + 2;
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2) {
	double b = 0.5; //wspolczynnik tarcia [N*m*s]
	double Mr = 1; //masa ramienia robota [kg]
	double Mc = 5; //Masa ciezarka [kg]
	double l = 1;  // dlugosc ramienia robota [m]
	double I = ((Mr * pow(l, 2)) / 3.0) + (Mc * pow(l, 2)); //Obliczenie momentu bezwladnosci [kg*m^2]

	double alfa_zad = M_PI; // wartosc docelowa kata
	double omega_zad = 0.0;  // wartosc docelowa predkosci katowej

	double k1 = ud1(0, 0); //wektor ud1 przechowuje wspolczynniki wzmocnienia
	double k2 = ud1(1, 0);

	double alfa = Y(0, 0);
	double omega = Y(1, 0);

	double M = k1 * (alfa_zad - alfa) + k2 * (omega_zad - omega);

	matrix dY(2, 1);	//wektor rozniczki
	dY(0, 0) = omega;
	dY(1, 0) = (M - b * omega) / I;
	return dY;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2) {
	double t0 = 0.0; //czas poczatkowy
	double dt = 0.1; //krok czasowy
	double t_end = 100; //czas koncowy
	matrix Y0(2, 1); //warunki poczatkowe


	matrix* diff_sol = solve_ode(df2, t0, dt, t_end, Y0, x); //rozwiazwywanie ukladu RR

	double alfa_zad = M_PI;
	double omega_zad = 0.0;

	double alfa_t = 0;
	double omega_t = 0;
	double M = 0;
	double Q = 0.0;
	for (int i = 0; i < get_len(diff_sol[0]); i++) {
		alfa_t = diff_sol[1](i, 0);
		omega_t = diff_sol[1](i, 1);
		M = x(0, 0) * (alfa_zad - alfa_t) + x(1, 0) * (omega_zad - omega_t);
		Q += dt * (10 * pow(alfa_zad - alfa_t, 2) + pow(omega_zad - omega_t, 2) + M * M);
	}

	return Q;

}

//LAB 3
matrix ff3T(matrix x, matrix ud1, matrix ud2)
{
	//ud2(0) = c - współczynnik c>0
	double x1 = x(0);
	double x2 = x(1);
	double x1_PI = x1 / M_PI;
	double x2_PI = x2 / M_PI;
	double up = sin(M_PI * sqrt(pow(x1_PI, 2) + pow(x2_PI, 2)));
	double down = M_PI * sqrt(pow(x1_PI, 2) + pow(x2_PI, 2));
	matrix y = up / down;
	
	//zewn S(x)
	if (ud2(1) > 1) {
		if (1 > x(0)) 
			y = y + ud2(0) * pow(-x(0) + 1, 2); //g1
		if (1 > x(1)) 
			y = y + ud2(0) * pow(-x(1) + 1, 2); //g2
		if (norm(x) > ud1(0)) 
			y = y + ud2(0) * pow(norm(x) - ud1(0), 2); //g3 //ud1(0) = parametr a
		return y;
	}
	//wewn S(x)
	if (1 > x(0)) y = 1e10;
	else y = y - ud2(0) / (1 - x(0)); //c/1-g1

	if (1 > x(1)) y = 1e10;
	else y = y - ud2(0) / (1 - x(1)); //c / 1-g2

	if (norm(x) > ud1(0)) y = 1e10;
	else y = y - ud2(0) / (norm(x) - ud1(0)); //ud1(0) = parametr a // c/ norm(x) - a

	return y;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2)
{
	//Dane startowe
	double t0 = 0.0;
	double tend = 7.0;
	double dt = 0.01;
	matrix Y0(4, new double[4] {
		x(0),	//V0x =V0x
		0,		//V0y = 0
		0,		//x0 = 0
		100,	//y0 = 100 m
		});

	matrix omega(1, new double[1]{ x(1)});
	matrix* diff_solution = solve_ode(df3, t0,dt,tend, Y0,ud1, omega);

	int n = get_len(diff_solution[0]);
	int maxval1, maxval2;
	cout << diff_solution[1] << endl;

	//Implementacja Funkcji Kary zewnetrznej	S(x) = E pow(max(0, gi()),2);
	
	for(int i = 0; i < n; i++) {
		//Sprawdzenie warunku czy przy y 50, x e <4.5 ;5.5>
		/*if(abs(diff_solution[1](i,1) -50) )
		if(abs(diff_solution[1](i,1)) )*/

	}

	
	return NULL;
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2)
{
	//Y(0) = V0x
	//Y(1) = V0y
	//Y(2) = x0
	//Y(3) = y0
	//ud2(0) = omega

	//Definicje zmiennych
	double g = 9.81; //g [m/s]
	double m = 0.6; //masa 0.6 [kg]
	double r = 0.12; // promien 0.12 [m]
	double ro = 1.2; //gestosc powietrza [kg/m3]
	double C = 0.47; // Wspolczynnik oporu kuli

	double S = M_PI * pow(r,2); //Przekruj kuli [m2]

	//Sily oporu powietrza
	double Dx = 0.5 * C * ro * S * Y(0) * abs(Y(0));
	double Dy = 0.5 * C * ro * S * Y(1) * abs(Y(1));

	//sila Magnusa
	double FMx = ro * Y(1) * m2d(ud2(0)) * M_PI * pow(r, 3);
	double FMy = ro * Y(0) * m2d(ud2(0)) * M_PI * pow(r, 3);

	//wektor rozniczki
	matrix dY(4,1);

	dY(0) = Y(2);					//x0
	dY(1) = Y(3);					//y0

	dY(2) = (-Dx -FMx)/m;			//ax
	dY(3) = (-Dy - FMy - m*g)/m;	//ay
	return dY;
}
