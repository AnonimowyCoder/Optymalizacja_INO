#include"user_funs.h"
#define M_PI 3.1415926

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2] { m2d(x), 0.5 });
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
	double I = m * pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;
	return dY;
}
matrix ff2T(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	//z = pow(x, 2) + pow(y, 2) - cos(2.5 * PI * x) - cos(2.5 * PI * y) + 2;
	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
	return y;
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
	const double g = 9.81;   // Przyspieszenie grawitacyjne [m/s²]
	const double m = 0.6;    // Masa [kg]
	const double r = 0.12;   // Promień [m]
	const double ro = 1.2;   // Gęstość powietrza [kg/m³]
	const double C = 0.47;   // Współczynnik oporu

	double S = M_PI * pow(r, 2); // Pole przekroju poprzecznego

	// Siły aerodynamiczne
	double Dx = 0.5 * C * ro * S * Y(1) * abs(Y(1));
	double Dy = 0.5 * C * ro * S * Y(3) * abs(Y(3));

	// Siły Magnusa
	double Fmx = M_PI * ro * Y(3) * m2d(ud2) * pow(r, 3);
	double Fmy = M_PI * ro * Y(1) * m2d(ud2) * pow(r, 3);

	matrix dY(4, 1);
	dY(0) = Y(1);                      // x(t)
	dY(1) = (-Dx - Fmx) / m;           // v_x(t)
	dY(2) = Y(3);                      // y(t)
	dY(3) = (-m * g - Dy - Fmy) / m;   // v_y(t)

	return dY;
}


matrix ff3R(matrix x, matrix ud1, matrix ud2) {
	matrix y;

	double t0 = 0.0;
	double tend = 7.0;
	double dt = 0.01;

	// Dane początkowe
	matrix Y0(4, new double[4] {
		0,      // x0 = 0
			x(0),   // V0x
			100,    // y0 = 100 m
			0       // V0y
		});

	matrix* Y = solve_ode(df3, t0, dt, tend, Y0, ud1, x(1));
	int n = get_len(Y[0]);

	// Ustalenie punktów dotyczących ograniczeń
	int punkt50 = -1;
	int punkt0 = -1;
	double wartosc50 = 1e10;
	double wartosc0 = 1e10;

	for (int i = 0; i < n - 1; ++i) {
		double odleglosc_50 = abs(Y[1](i, 2) - 50);
		double odleglosc_0 = abs(Y[1](i, 2));

		if (odleglosc_50 < wartosc50) {
			wartosc50 = odleglosc_50;
			punkt50 = i;
		}
		if (odleglosc_0 < wartosc0) {
			wartosc0 = odleglosc_0;
			punkt0 = i;
		}
	}

	// Położenie końcowe
	y = -1 * Y[1](punkt0, 0);

	// Warunki kar
	if (abs(x(0)) - 10 > 0) {
		double kara_x = ud2() * pow(abs(x(0)) - 10, 2);
		y = y + kara_x;  // Kara dla x
	}

	if (abs(x(1)) - 15 > 0) {
		double kara_omega = ud2() * pow(abs(x(1)) - 15, 2);
		y = y + kara_omega;  // Kara dla omega
	}

	if (abs(Y[1](punkt50, 0) - 5) - 0.5 > 0) {
		double kara_y50 = ud2() * pow(abs(Y[1](punkt50, 0) - 5) - 0.5, 2);
		y = y + kara_y50;  // Kara dla y=50
	}

	// Zapis wyników do pliku
	ofstream out("symulacja_real.csv");
	out << "Czas,X,Y\n";

	for (int i = 0; i < n; ++i) {
		out << Y[0](i) << "," << Y[1](i, 0) << "," << Y[1](i, 2) << "\n";
	}

	return y;
}

matrix ff4T(matrix x, matrix ud1, matrix ud2) {
	matrix y;

	if (isnan(ud2(0, 0)))
		y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
	else
		y = ff4T(ud2[0] + x * ud2[1], 0, NAN);
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

//ud1 ma strukture:
//ud1(0,0) = a
//ud1(1,0) = w


//ud1 = a
matrix ff5T_1(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	double a = ud1(0,0); //przypisanie wartosci a
	y = a * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2));
	//cout << "y1 = " << y << endl;
	return y;
}
//ud1 = a
matrix ff5T_2(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	double a = ud1(0,0); //przypisanie wartosci a
	y = (1/a) * (pow(x(0) + 2, 2) + pow(x(1) + 2, 2));
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
matrix trivial(matrix x, matrix ud1, matrix ud2) {
	matrix y;

	if (isnan(ud2(0, 0)))
		y = pow(x(0), 2) + pow(x(1), 2);
	else
		y = trivial(ud2[0] + x * ud2[1], 0, NAN);
	return y;
}
