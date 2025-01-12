/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"

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
		lab5();
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
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{

}

void lab2()
{
	double epsilon = 1e-6;
	int Nmax = 10000;
	double s = 0.05;
	double alpha = 0.8;
	double start_point_array[2] = { 0.4,0.4 };
	matrix start_point(2, start_point_array);
	cout << "start point: " << start_point << endl;
	solution opt;
	opt = HJ(ff2T, start_point, s, alpha, epsilon, Nmax);
	cout << "Policzomne" << endl << endl;
	cout << "Watosc x: " << endl << opt.x << endl;
	cout << "Wartosc y: " << opt.y << endl;
	//cout << opt << endl;
}

void lab3()
{
	matrix x0 = matrix(2, 1);
	srand(time(NULL));

	//Plik csv
	ofstream testfunkc("testfunc.csv");
	if (!testfunkc.is_open())throw string("PLIK CSV NIE OTWARTY");

	// Nag��wki w pliku CSV
	testfunkc << "Parametr a,Lp.,x1(0),x2(0),"
		<< "ZE_x1,ZE_x2,r,y,Liczba wywolan funkcji celu,"
		<< "WE_x1,WE_x2,r,y,Liczba wywolan funkcji celu\n";

	//Inicjalizacja metody Neldera-Meada
	//double s = 0.5; //dlugosc synopsu
	//double alpha = 1.0; //wspolczynnik odbicia
	//double beta = 0.5; //zawezenie <1
	//double gamma = 2.0; // ekspansjia >1
	//double delta = 0.5; //redukcja <1
	const double epsilon = 1e-3; //dokladnosc
	const int Nmax = 100000; // liczba wywolan

	//Inicjal funkcji kary
	matrix parametr_a[3] = { 4.0, 4.4934, 5.0 }; //wysylany do ud1
	solution wynik_zewn, wynik_wewn;
	int wybrany = 0;
	double c = 2.0;

	//100 optymalizacji:
	for (int i = 0; i < 100; i++) {
		//break;
		do // g(x1, x2) <= a
		{
			x0(0) = (((double)rand() / RAND_MAX) * 5) + 1;
			x0(1) = (((double)rand() / RAND_MAX) * 5) + 1;

		} while (norm(x0) > parametr_a[wybrany]);

		cout << "Wylosowane pkt startowe:\n";
		cout << x0 << endl;

		wynik_zewn = pen(ff3T, x0, c, 2, epsilon, Nmax, parametr_a[wybrany]);
		cout << wynik_zewn;

		testfunkc << m2d(parametr_a[wybrany]) << "," << i + 1 << "," << x0(0) << "," << x0(1) << ","
			<< wynik_zewn.x(0) << "," << wynik_zewn.x(1) << "," << norm(wynik_zewn.x) << "," << wynik_zewn.y(0) << "," << solution::f_calls << ",";// << "\n";
		solution::clear_calls();

		wynik_wewn = pen(ff3T, x0, c, 0.5, epsilon, Nmax, parametr_a[wybrany]);
		testfunkc << wynik_wewn.x(0) << "," << wynik_wewn.x(1) << "," << norm(wynik_wewn.x) << "," << wynik_wewn.y(0) << "," << solution::f_calls << "," << "\n";
		solution::clear_calls();
	}
	testfunkc.close();

	//matrix rozniczka = df3();
	//cout << "wynik rozniczki\n";

	//FUNKCJA RZECZYWISTA:
	solution::clear_calls();
	x0(0) = 0;		//V0x						
	x0(1) = 0;		//Omega [rad/s]				
	solution wynik_real = pen(ff3R, x0, c, 2, epsilon, Nmax);
	cout << wynik_real;

}
void lab4()
{
	double tab_test[] = { 10.0,7.0 };
	//punkt startowy 
	matrix x0 = matrix(2, tab_test);

	const double step = 0.1;
	const double epsilon = 1e-4; //dokladnosc
	const int Nmax = 10000; // liczba wywolan

	solution test;
	test = SD(ff4T, gf4T, x0,-1, epsilon, Nmax);
	cout << "##################################" << endl;
	cout << "Wyniki: " << endl;
	cout << test.x << endl;
	cout << test.y;
}

void lab5()
{
	double tab_test[] = { 2.0,1.0 };
	matrix x0 = matrix(2, tab_test);
	//ustawienie macierzy parametrów
	matrix ud1(2,1,0.0);
	ud1.set_row(1.0, 0);
	ud1.set_row(0.5, 1);
	
	const double epsilon = 1e-4; //dokladnosc
	const int Nmax = 10000; // liczba wywolan

	solution test;
	test = Powell(ff5T_comb, x0, epsilon, Nmax,ud1);
	//test = Powell(ff4T, x0, epsilon, Nmax, ud1);
	//test = Powell(trivial, x0, epsilon, Nmax, ud1);

	cout << "##################################" << endl;
	cout << "Wyniki: " << endl;
	cout << test.x << endl;
	cout << test.y;



}

void lab6()
{

}