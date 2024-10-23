/*********************************************
Kod stanowi uzupełnienie materiałów do ćwiczeń
w ramach przedmiotu metody optymalizacji.
Kod udostępniony na licencji CC BY-SA 3.0
Autor: dr inż. Łukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#define _USE_MATH_DEFINES
#include"opt_alg.h"
#include"math.h"
#include<vector>
#include<cstdio>

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{

	matrix ud1(1, 1, -100), ud2(1, 1, 100);
	double* eskp = new double[2];
	double x0 = rand() % 201 - 100, d = 1.0, alfa =5.29;
	int Nmax = 10000;
	int l = 0;
	double eps = 1e-2;
	double gamma = 1e-200;
	solution fibo;
	solution lagi;
	//ofstream expToFile("./eksp.csv");
	//ofstream fibToFile("./fib.csv");
	//ofstream lafToFile("./lag.csv");
	ofstream slagToFile("./slag.csv");
	ofstream sfimToFile("./sfib.csv");

	/*for (int i = 0; i < 100; i++) {
		eskp = expansion(ff1L, x0, d, alfa, Nmax);
		expToFile << x0 << ";" << eskp[0] << ";" << eskp[1] << ";" << solution::f_calls << endl;
		solution::clear_calls();
		x0 = rand() % 201 - 100;

		fibo = fib(ff1L, eskp[0], eskp[1], eps, NULL, NULL);
		fibToFile << m2d(fibo.x) << ";" << m2d(fibo.y) << ";" << solution::f_calls << endl;
		solution::clear_calls();

		lagi = lag(ff1L, eskp[0], eskp[1], eps, gamma, Nmax, NULL, NULL);
		lafToFile << m2d(lagi.x) << ";" << m2d(lagi.y) << ";" << solution::f_calls << endl;
		solution::clear_calls();
	}
	expToFile.close();
	fibToFile.close();
	lafToFile.close();*/

	//zad1_bezpara
	/*fibo = fib(ff1L, -100, 100, eps);
	lagi = lag(ff1L, -100, 100, eps, gamma, Nmax);
	cout << m2d(fibo.x) << ";" << m2d(fibo.y) << ";" << solution::f_calls << endl;
	cout << m2d(lagi.x) << ";" << m2d(lagi.y) << ";" << solution::f_calls << endl;*/
	//symulacje
	double epss = 0.00000000001;
	double gammas = 0.000001;
	fibo = fib(ff1R, 0.0001, 0.01,epss,ud1,ud2);
	//cout << m2d(fibo.x) << ";" << m2d(fibo.y) << ";" << solution::f_calls << endl;
	lagi= lag(ff1R, 0.0001, 0.01, epss,gammas, Nmax, ud1, ud2);
	//cout << m2d(lagi.x) << ";" << m2d(lagi.y) << ";" << solution::f_calls << endl;
	matrix Y0 = matrix(3, new double[3] {5, 1, 20});
	matrix* Yf = solve_ode(df1R, 0, 1, 2000, Y0, ud1, fibo.x(0));
	sfimToFile << Yf[1] << endl;
	matrix* Yl= solve_ode(df1R, 0, 1, 2000, Y0, ud1, lagi.x(0));
	slagToFile << Yl[1] << endl;
	sfimToFile.close();
	slagToFile.close();
}