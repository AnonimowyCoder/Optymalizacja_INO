/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"
#include <vector>

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

/*
void lab1() {
	
	ofstream file("lab1_ftest_7_5.csv"); // Otwieranie pliku CSV

	// Generator liczb losowych
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> distribution(-100.0, 100.0);

	// Sprawdzenie, czy plik został poprawnie otwarty
	if (!file.is_open()) {
		cout << "file error" << endl;
		return;
	}

	// Nagłówki w pliku CSV, dostosowane do struktury tabeli w Excelu
	file << "Współczynnik ekspansji,Lp.,x(0),a,b,Liczba wywołań funkcji celu,"
		<< "Fib_x*,Fib_y*,Fib_Liczba wywołań funkcji celu,Fib_Minimum,"
		<< "Lag_x*,Lag_y*,Lag_Liczba wywołań funkcji celu,Lag_Minimum\n";



	double wspolczynnik_ekspansji = 6.5; //współczynnika ekspansji

	for (int i = 0; i < 100; i++) {
		
		double startPoint = distribution(gen); // Losowanie punktu startowego

		double* part = expansion(ff1T, startPoint, 0.5, wspolczynnik_ekspansji, 2000); // Ekspansja wstępna, zawężenie przedziału

		solution minimumFib = fib(ff1T, part[0], part[1], 0.0000001); //Metoda Fibonacciego

		solution minimumLag = lag(ff1T, part[0], part[1], 0.00001, 0.0000001, 1000); //Metody Lagrange'a

		string fib_minimum_type = (minimumFib.x(0) < 50) ? "lokalne" : "globalne"; 		// Określenie, czy minimum jest lokalne czy globalne

		string lag_minimum_type = (minimumLag.x(0) < 50) ? "lokalne" : "globalne";

		// Zapis wyników do pliku CSV zgodnie z szablonem
		file << wspolczynnik_ekspansji << ","   // Współczynnik ekspansji
			<< i + 1 << ","                    // Lp.
			<< startPoint << ","               // Punkt startowy x(0)
			<< part[0] << ","                  // Dolny zakres przedziału (a)
			<< part[1] << ","                  // Górny zakres przedziału (b)
			<< solution::f_calls << ","        // Liczba wywołań funkcji celu dla ekspansji
			<< minimumFib.x(0) << ","          // Wynik optymalizacji Fibonacciego (x*)
			<< minimumFib.y(0) << ","          // Wartość funkcji celu dla Fibonacciego (y*)
			<< minimumFib.f_calls << ","       // Liczba wywołań funkcji celu dla Fibonacciego
			<< fib_minimum_type << ","         // Typ minimum (lokalne/globalne) dla Fibonacciego
			<< minimumLag.x(0) << ","          // Wynik optymalizacji Lagrange'a (x*)
			<< minimumLag.y(0) << ","          // Wartość funkcji celu dla Lagrange'a (y*)
			<< minimumLag.f_calls << ","       // Liczba wywołań funkcji celu dla Lagrange'a
			<< lag_minimum_type << "\n";       // Typ minimum (lokalne/globalne) dla Lagrange'a
	}

	file.close(); // Zamknięcie pliku po zakończeniu
}
*/

void lab1()
{

	matrix ud1(1, 1, -100), ud2(1, 1, 100);
	double* eskp = new double[2];
	double x0 = rand() % 201 - 100, d = 1.0, alfa = 5.29;
	int Nmax = 10000;
	int l = 0;
	solution fibonachi;
	solution lagrange;

	ofstream lagrangeToFile("./lagrange_rozwiazanie.csv");
	ofstream fibonachiToFile("./fibonachi_rozwiazanie.csv");

	
	//symulacje

	fibonachi = fib(ff1R, 0.0001, 0.01, 0.0000001, ud1, ud2);
	//cout << m2d(fibonachi.x) << ";" << m2d(fibonachi.y) << ";" << solution::f_calls << endl;
	cout << fibonachi.x(0) << endl << fibonachi.y(0) << endl << fibonachi.f_calls << endl;

	printf("\n");

	lagrange = lag(ff1R, 0.0001, 0.01, 0.00001, 0.0000001, Nmax, ud1, ud2);
	//cout << m2d(lagrange.x) << ";" << m2d(lagrange.y) << ";" << solution::f_calls << endl;
	cout << lagrange.x(0) << endl << lagrange.y(0) << endl << lagrange.f_calls << endl;

	matrix Y0 = matrix(3, new double[3] {5, 1, 20});

	matrix* Yf = solve_ode(df1, 0, 1, 2000, Y0, ud1, fibonachi.x(0));
	fibonachiToFile << Yf[1] << endl;

	matrix* Yl = solve_ode(df1, 0, 1, 2000, Y0, ud1, lagrange.x(0));
	lagrangeToFile << Yl[1] << endl;

	fibonachiToFile.close();
	lagrangeToFile.close();
}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
