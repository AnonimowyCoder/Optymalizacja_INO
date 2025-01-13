#include"user_funs.h"
#include <cmath>
#include "opt_alg.h"

#define M_PI 3.14159265358979323846 
//define M_PI 3.14

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


	// Zewnętrzny przypadek S(x)
	if (ud2(1) > 1) {
		if (x1 < 1)
			y = y + ud2(0) * pow(1 - x1, 2); // g1
		if (x2 < 1)
			y = y + ud2(0) * pow(1 - x2, 2); // g2
		if (norm(x) > ud1(0))
			y = y + ud2(0) * pow(norm(x) - ud1(0), 2); // g3 (ud1(0) - parametr a)
		return y;
	}

	//Wewnętrzny przypadek S(x)
	if (x1 < 1)
		return 1e10;
	y = y - ud2(0) / (1 - x1); // c / (1 - g1)

	if (x2 < 1)
		return 1e10;
	y = y - ud2(0) / (1 - x2); // c / (1 - g2)

	if (norm(x) > ud1(0))
		return 1e10;
	y = y - ud2(0) / (norm(x) - ud1(0)); // c / (norm(x) - a)

	return y;
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2) {
	double g = 9.81; //g [m/s]
	double m = 0.6; //masa 0.6 [kg]
	double r = 0.12; // promien 0.12 [m]
	double ro = 1.2; //gestosc powietrza [kg/m3]
	double C = 0.47; // Wspolczynnik oporu kuli

	double S = M_PI * pow(r, 2);

	double Dx = 0.5 * C * ro * S * Y(1) * abs(Y(1));
	double Dy = 0.5 * C * ro * S * Y(3) * abs(Y(3));

	double Fmx = M_PI * ro * Y(3) * m2d(ud2) * pow(r, 3);
	double Fmy = M_PI * ro * Y(1) * m2d(ud2) * pow(r, 3);

	matrix dY(4, 1);
	dY(0) = Y(1);						//x0 = v0x
	dY(1) = (-Dx - Fmx) / m;			//V0x = ax
	dY(2) = Y(3);						//y0 = V0y
	dY(3) = (-m * g - Dy - Fmy) / m;	//V0y = ay

	return dY;

}

matrix ff3R(matrix x, matrix ud1, matrix ud2) {
	matrix y;

	double t0 = 0.0;
	double tend = 7.0;
	double dt = 0.01;

	//Dane wejsciowe
	matrix Y0(4, new double[4] {
		0,      // x0 = 0
		x(0),   // V0x
		100,    // y0 = 100 m
		0       // V0y
		});


	matrix* Y = solve_ode(df3, t0, dt, tend, Y0, ud1, x(1));
	int n = get_len(Y[0]);

	int punkt50 = -1;
	int punkt0 = -1;
	double wartosc50 = 1e10;
	double wartosc0 = 1e10;

	for (int i = 0; i < n - 1; ++i) {
		if (abs(Y[1](i, 2) - 50) < wartosc50) {
			wartosc50 = abs(Y[1](i, 2) - 50);
			punkt50 = i;
		}
		if (abs(Y[1](i, 2)) < wartosc0) {
			wartosc0 = abs(Y[1](i, 2));
			punkt0 = i;
		}
	}


	y = -Y[1](punkt0, 0);

	if (abs(x(0)) - 10 > 0)
		y = y + ud2() * pow(abs(x(0)) - 10, 2); //WARUNEK KARA DLA X

	if (abs(x(1)) - 15 > 0)
		y = y + ud2() * pow(abs(x(1)) - 15, 2); //KARA DLA OMEGA

	if (abs(Y[1](punkt50, 0) - 5) - 0.5 > 0)
		y = y + ud2() * pow(abs(Y[1](punkt50, 0) - 5) - 0.5, 2); //KARA DLA y=50 

	ofstream out("symulacja_real.csv");
	out << "Czas, X , Y " << endl;

	for (int i = 0; i < n; ++i)
		out << Y[0](i) << "," << Y[1](i, 0) << "," << Y[1](i, 2) << endl;

	return y;
}


// LAB 4
matrix ff4T(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	if (isnan(ud2(0, 0))) {	//jeśli wyznaczamy po prostu wartość funkcji dwuwymiarowej

		y = pow((x(0) + 2 * x(1) - 7), 2) + pow((2 * x(0) + x(1) - 5), 2);

	}
	else {		//jeśli wyznaczamy wartość funkcji w jakimś kierunku
		y = ff4T(ud2[0] + x * ud2[1], 0, ud1);
	}
	return y;
}

matrix gf4T(matrix x, matrix ud1, matrix ud2) {
	matrix g(2, 1);
	g(0) = -34 + 10 * x(0) + 8 * x(1);
	g(1) = -38 + 8 * x(0) + 10 * x(1);
	return g;
}

matrix hf4T(matrix x, matrix ud1, matrix ud2) {
	matrix h(2, 2);
	h(0, 0) = 10; h(0, 1) = 8;
	h(1, 0) = 8;  h(1, 1) = 10;
	return h;
}

//matrix gf4(matrix theta, matrix ud1, matrix ud2) {
//	int m = 100;
//	int n = get_dim(theta);  // Rozmiar wektora theta
//	matrix g(n, 1);          // Gradient
//	static matrix X(n, m), Y(1, m);
//
//	if (solution::g_calls == 1) {
//		ifstream in("XData.txt");
//		in >> X;  // Wczytaj dane uczące do macierzy X
//		in.close();
//		in.open("YData.txt");
//		in >> Y;  // Wczytaj etykiety do macierzy Y
//		in.close();
//	}
//
//	double h;  // Hipoteza
//	for (int j = 0; j < n; ++j) {
//		for (int i = 0; i < m; ++i) {
//			h = m2d(trans(theta) * X[i]);  // Hipoteza: h_theta(x_i)
//			h = 1 / (1 + exp(-h));         // Funkcja sigmoid
//			g(j) += X(j, i) * (h - Y(0, i));  // Oblicz gradient
//		}
//		g(j) /= m;  // Normalizuj przez liczbę przykładów
//	}
//	return g;
//}

//matrix ff4R(matrix theta, matrix ud1, matrix ud2) {
//	matrix y;           // Wartość funkcji kosztu
//	int m = 100;        // Liczba przykładów
//	int n = get_len(theta);  // Liczba parametrów
//	static matrix X(n, m), Y(1, m);
//
//	if (solution::f_calls == 1) {
//		ifstream in("XData.txt");
//		in >> X;  // Wczytaj dane uczące do macierzy X
//		in.close();
//		in.open("YData.txt");
//		in >> Y;  // Wczytaj etykiety do macierzy Y
//		in.close();
//	}
//
//	double h;  // Hipoteza
//	y = 0;     // Inicjalizacja funkcji kosztu
//	for (int i = 0; i < m; i++) {
//		h = m2d(trans(theta) * X[i]);  // Hipoteza: h_theta(x_i)
//		h = 1.0 / (1.0 + exp(-h));     // Funkcja sigmoid
//		y = Y(0, i) * log(h) + (1 - Y(0, i)) * log(1 - h);  // Funkcja kosztu
//	}
//	y = y/ m;  // Normalizuj przez liczbę przykładów
//	return y;
//}

matrix ff4R(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	int m = 100;
	int n = get_len(x);
	static matrix X(n, m), Y(1, m);
	if (solution::f_calls == 1) {
		ifstream in("XData.txt");
		in >> X;
		in.close();
		in.open("YData.txt");
		in >> Y;
		in.close();
	}
	int P = 0;
	double h;
	y = 0;
	for (int i = 0; i < m; i++) {
		h = m2d(trans(x) * X[i]);
		h = 1.0 / (1.0 + exp(-h));
		y = y - Y(0, i) * log(h) - (1 - Y(0, i)) * log(1 - h);
	}
	y = y / m;
	return y;
}

matrix gf4(matrix x, matrix ud1, matrix ud2) {
	int m = 100;
	int n = get_len(x);
	matrix g(n, 1);
	static matrix X(n, m), Y(1, m);
	if (solution::g_calls == 1) {
		ifstream in("XData.txt");
		in >> X;
		in.close();
		in.open("YData.txt");
		in >> Y;
		in.close();
	}

	double h;
	for (int j = 0; j < n; ++j) {
		for (int i = 0; i < m; ++i) {
			h = m2d(trans(x) * X[i]);
			h = 1 / (1 + exp(-h));
			g(j) = g(j) + X(j, i) * (h - Y(0, i));
		}
		g(j) = g(j) / m;
	}
	return g;
}


//matrix trivial(matrix x, matrix ud1, matrix ud2) {
//	matrix y;
//
//	if (isnan(ud2(0, 0)))
//		y = pow(x(0), 2) + pow(x(1), 2);
//	else
//		y = trivial(ud2[0] + x * ud2[1], 0, NAN);
//	return y;
//}

matrix ff5T(matrix x, matrix ud1, matrix ud2)
{
	matrix y(2, 1); // ud1(0) = a	
	y(0) = ud1(0) * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2));
	y(1) = (1/ ud1(0)) * (pow(x(0) + 2, 2) + pow(x(1) + 2, 2));
	return y;
}

matrix ff5T_1(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	double a = ud1(0, 0); //przypisanie wartosci a
	y = a * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2));
	//cout << "y1 = " << y << endl;
	return y;
}
//ud1 = a
matrix ff5T_2(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	double a = ud1(0, 0); //przypisanie wartosci a
	y = (1 / a) * (pow(x(0) + 2, 2) + pow(x(1) + 2, 2));
	//cout << "y2 = " << y << endl;

	return y;
}

//zamiana problemu wielokryterialnego na jednokryterialny
// 
//ud1(1,0) = w
//ud1(0,0) = a;		ud2 = wyznaczanie minimum po kierunku
matrix ff5T_comb(matrix x, matrix ud1, matrix ud2) {
	double w = ud1(1, 0);
	matrix y;
	if (isnan(ud2(0, 0))) {	//jeśli wyznaczamy po prostu wartość funkcji dwuwymiarowej
		//cout << "Wyznaczam wartosc funkcji w ogolnie" << endl;

		y = w * ff5T_1(x, ud1, ud2) + (1 - w) * ff5T_2(x, ud1, ud2);
	}
	else {		//jeśli wyznaczamy wartość funkcji w jakimś kierunku
		//cout << "Wyznaczam wartosc funkcji w kierunku" << endl;
		y = ff5T_comb(ud2[0] + x * ud2[1], ud1, NAN);
	}
	//cout << "zwracana wartosc y_comb to: " << y << endl;
	return y;
}




matrix ff5R(matrix x, matrix ud1, matrix ud2) {
	// Stałe materiałowe i geometryczne
	const double P = 1000.0;    // Siła [N]
	const double E = 207e9;     // Moduł Younga [Pa]
	const double rho = 7800.0;  // Gęstość [kg/m^3]
	matrix y;

	if (isnan(ud2(0, 0))) {
		double l = x(0);  // Długość [mm]
		double d = x(1);  // Średnica [mm]

		// Konwersja na metry dla obliczeń
		double l_m = l * 0.001;   // [m]
		double d_m = d * 0.001;   // [m]

		y = matrix(3, 1);

		// Masa [kg]
		double mass = rho * (M_PI * pow(d_m, 2) / 4.0) * l_m;
		y(0) = mass;

		// Ugięcie [mm]
		double deflection = (64.0 * P * pow(l_m, 3)) / (3.0 * E * M_PI * pow(d_m, 4));
		y(1) = deflection * 1000.0;  // Konwersja na [mm]

		// Naprężenie [Pa]
		double stress = (32.0 * P * l_m) / (M_PI * pow(d_m, 3));
		y(2) = stress;

	}
	else {
		// Transformacja zmiennych dla metody optymalizacji
		matrix xt = ud2[0] + x * ud2[1];
		matrix yt = ff5R(xt, ud1);

		// Funkcja celu z normalizacją
		y = ud1 * (yt(0) - 0.12) / (15.3 - 0.12) +
			(1 - ud1) * (yt(1) - 0.042) / (3.2 - 0.042);

		// Funkcja kary
		const double c = 1e10;

		// Ograniczenia na długość [mm]
		if (xt(0) < 200) {
			y = y + c * pow(200 - xt(0), 2);
		}
		if (xt(0) > 1000) {
			y = y + c * pow(xt(0) - 1000, 2);
		}

		// Ograniczenia na średnicę [mm]
		if (xt(1) < 10) {
			y = y + c * pow(10 - xt(1), 2);
		}
		if (xt(1) > 50) {
			y = y + c * pow(xt(1) - 50, 2);
		}

		// Ograniczenie na ugięcie [mm]
		if (yt(1) > 5) {
			y = y + c * pow(yt(1) - 5, 2);
		}

		// Ograniczenie na naprężenie [Pa]
		if (yt(2) > 300e6) {
			y = y + c * pow(yt(2) - 300e6, 2);
		}
	}

	return y;
}
