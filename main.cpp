/*********************************************
Kod stanowi uzupełnienie materiałów do ćwiczeń
w ramach przedmiotu metody optymalizacji.
Kod udostępniony na licencji CC BY-SA 3.0
Autor: dr inż. Łukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"
#include <fstream>

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
	return 0;
}

void lab0()
{
	////Funkcja testowa
	//double epsilon = 1e-2;
	//int Nmax = 10000;
	//matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	//solution opt;
	//a(0) = -1;
	//a(1) = 2;
	//opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	//cout << opt << endl << endl;
	//solution::clear_calls();

	////Wahadlo
	//Nmax = 1000;
	//epsilon = 1e-2;
	//lb = 0;
	//ub = 5;
	//double teta_opt = 1;
	//opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	//cout << opt << endl << endl;
	//solution::clear_calls();

	////Zapis symulacji do pliku csv
	//matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });
	//matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	//ofstream Sout("symulacja_lab0.csv");
	//Sout << hcat(Y[0], Y[1]);
	//Sout.close();
	//Y[0].~matrix();
	//Y[1].~matrix();
}

//void lab1() {
//	ofstream file("lab1_tab3.csv");
//	random_device rd;
//	mt19937 gen(rd());
//	uniform_real_distribution<double> distribution(-100.0, 100.0);
//	if (!file.is_open()) {
//		cout << "file error" << endl;
//		return;
//	}
//	for (int i = 0; i < 100; i++) {
//		double startPoint = distribution(gen);
//		double* part = expansion(ff1T, startPoint, 0.5, 1.5, 1000);
//		file << startPoint << ":" << part[0] << ":" << part[1] << ":" << solution::f_calls << ":";
//
//		solution minimumFib = fib(ff1T, part[0], part[1], 0.0000001);
//		file << minimumFib.x(0) << ":" << minimumFib.y(0) << ":" << minimumFib.f_calls << ":";
//		if (minimumFib.x(0) < 50)file << "lokalne";
//		else file << "globalne";
//		file << ":";
//
//		solution minimumLag = lag(ff1T, part[0], part[1], 0.00001, 0.0000001, 1000);
//		file << minimumLag.x(0) << ":" << minimumLag.y(0) << ":" << minimumLag.f_calls << ":";
//		if (minimumLag.x(0) < 50)file << "lokalne";
//		else file << "globalne";
//		file << endl;
//
//		matrix ud1(1, 1, -100), ud2(1, 1, 100);
//		solution fib_result = lag(ff1R, 0.0001, 0.01, 0.0000001, 0.0000001, 1000);
//		cout << fib_result.x(0) << endl;
//		cout << fib_result.y(0) << endl;
//		matrix Y0 = matrix(3, new double[3]{ 5, 1, 10 });
//		matrix* Y = solve_ode(df1, 0, 1, 1000, Y0, ud1, fib_result.x(0));
//		file << Y[1] << endl;
//		cout << solution::f_calls << endl;
//	}
//	file.close();
//}


void lab1() {
	// Otwieranie pliku CSV do zapisu wyników
	ofstream file("pojednyczaIteracja.csv");

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

	double wspolczynnik_ekspansji = 1.5; // Przykład wartości współczynnika ekspansji, można dostosować

	// Pętla dla 100 optymalizacji
	//for (int i = 0; i < 100; i++) {
		// Losowanie punktu startowego
		double startPoint = distribution(gen);

		// Ekspansja wstępna, zawężenie przedziału
		//double* part = expansion(ff1T, startPoint, 0.5, wspolczynnik_ekspansji, 2000);

		// Wynik metody Fibonacciego
		solution minimumFib = fib(ff1T, -100, 100, 0.0000001);

		// Wynik metody Lagrange'a
		solution minimumLag = lag(ff1T, -100,100, 0.00001, 0.0000001, 1000);

		// Określenie, czy minimum jest lokalne czy globalne (Fibonacciego)
		string fib_minimum_type = (minimumFib.x(0) < 50) ? "lokalne" : "globalne";

		// Określenie, czy minimum jest lokalne czy globalne (Lagrange'a)
		string lag_minimum_type = (minimumLag.x(0) < 50) ? "lokalne" : "globalne";

		// Zapis wyników do pliku CSV zgodnie z szablonem
		file << wspolczynnik_ekspansji << ","   // Współczynnik ekspansji
			//<< i + 1 << ","                    // Lp.
			<< startPoint << ","               // Punkt startowy x(0)
			//<< part[0] << ","                  // Dolny zakres przedziału (a)
			//<< part[1] << ","                  // Górny zakres przedziału (b)
			<< solution::f_calls << ","        // Liczba wywołań funkcji celu dla ekspansji
			<< minimumFib.x(0) << ","          // Wynik optymalizacji Fibonacciego (x*)
			<< minimumFib.y(0) << ","          // Wartość funkcji celu dla Fibonacciego (y*)
			<< minimumFib.f_calls << ","       // Liczba wywołań funkcji celu dla Fibonacciego
			<< fib_minimum_type << ","         // Typ minimum (lokalne/globalne) dla Fibonacciego
			<< minimumLag.x(0) << ","          // Wynik optymalizacji Lagrange'a (x*)
			<< minimumLag.y(0) << ","          // Wartość funkcji celu dla Lagrange'a (y*)
			<< minimumLag.f_calls << ","       // Liczba wywołań funkcji celu dla Lagrange'a
			<< lag_minimum_type << "\n";       // Typ minimum (lokalne/globalne) dla Lagrange'a
	//}

	//wspolczynnik_ekspansji = 3.0;
	//for (int i = 0; i < 100; i++) {
	//	// Losowanie punktu startowego
	//	double startPoint = distribution(gen);

	//	// Ekspansja wstępna, zawężenie przedziału
	//	double* part = expansion(ff1T, startPoint, 0.5, wspolczynnik_ekspansji, 2000);

	//	// Wynik metody Fibonacciego
	//	solution minimumFib = fib(ff1T, part[0], part[1], 0.0000001);

	//	// Wynik metody Lagrange'a
	//	solution minimumLag = lag(ff1T, part[0], part[1], 0.00001, 0.0000001, 1000);

	//	// Określenie, czy minimum jest lokalne czy globalne (Fibonacciego)
	//	string fib_minimum_type = (minimumFib.x(0) < 50) ? "lokalne" : "globalne";

	//	// Określenie, czy minimum jest lokalne czy globalne (Lagrange'a)
	//	string lag_minimum_type = (minimumLag.x(0) < 50) ? "lokalne" : "globalne";

	//	// Zapis wyników do pliku CSV zgodnie z szablonem
	//	file << wspolczynnik_ekspansji << ","   // Współczynnik ekspansji
	//		<< i + 1 << ","                    // Lp.
	//		<< startPoint << ","               // Punkt startowy x(0)
	//		<< part[0] << ","                  // Dolny zakres przedziału (a)
	//		<< part[1] << ","                  // Górny zakres przedziału (b)
	//		<< solution::f_calls << ","        // Liczba wywołań funkcji celu dla ekspansji
	//		<< minimumFib.x(0) << ","          // Wynik optymalizacji Fibonacciego (x*)
	//		<< minimumFib.y(0) << ","          // Wartość funkcji celu dla Fibonacciego (y*)
	//		<< minimumFib.f_calls << ","       // Liczba wywołań funkcji celu dla Fibonacciego
	//		<< fib_minimum_type << ","         // Typ minimum (lokalne/globalne) dla Fibonacciego
	//		<< minimumLag.x(0) << ","          // Wynik optymalizacji Lagrange'a (x*)
	//		<< minimumLag.y(0) << ","          // Wartość funkcji celu dla Lagrange'a (y*)
	//		<< minimumLag.f_calls << ","       // Liczba wywołań funkcji celu dla Lagrange'a
	//		<< lag_minimum_type << "\n";       // Typ minimum (lokalne/globalne) dla Lagrange'a
	//}

	//wspolczynnik_ekspansji = 7.5;
	//for (int i = 0; i < 100; i++) {
	//	// Losowanie punktu startowego
	//	double startPoint = distribution(gen);

	//	// Ekspansja wstępna, zawężenie przedziału
	//	double* part = expansion(ff1T, startPoint, 0.5, wspolczynnik_ekspansji, 2000);

	//	// Wynik metody Fibonacciego
	//	solution minimumFib = fib(ff1T, part[0], part[1], 0.0000001);

	//	// Wynik metody Lagrange'a
	//	solution minimumLag = lag(ff1T, part[0], part[1], 0.00001, 0.0000001, 1000);

	//	// Określenie, czy minimum jest lokalne czy globalne (Fibonacciego)
	//	string fib_minimum_type = (minimumFib.x(0) < 50) ? "lokalne" : "globalne";

	//	// Określenie, czy minimum jest lokalne czy globalne (Lagrange'a)
	//	string lag_minimum_type = (minimumLag.x(0) < 50) ? "lokalne" : "globalne";

	//	// Zapis wyników do pliku CSV zgodnie z szablonem
	//	file << wspolczynnik_ekspansji << ","   // Współczynnik ekspansji
	//		<< i + 1 << ","                    // Lp.
	//		<< startPoint << ","               // Punkt startowy x(0)
	//		<< part[0] << ","                  // Dolny zakres przedziału (a)
	//		<< part[1] << ","                  // Górny zakres przedziału (b)
	//		<< solution::f_calls << ","        // Liczba wywołań funkcji celu dla ekspansji
	//		<< minimumFib.x(0) << ","          // Wynik optymalizacji Fibonacciego (x*)
	//		<< minimumFib.y(0) << ","          // Wartość funkcji celu dla Fibonacciego (y*)
	//		<< minimumFib.f_calls << ","       // Liczba wywołań funkcji celu dla Fibonacciego
	//		<< fib_minimum_type << ","         // Typ minimum (lokalne/globalne) dla Fibonacciego
	//		<< minimumLag.x(0) << ","          // Wynik optymalizacji Lagrange'a (x*)
	//		<< minimumLag.y(0) << ","          // Wartość funkcji celu dla Lagrange'a (y*)
	//		<< minimumLag.f_calls << ","       // Liczba wywołań funkcji celu dla Lagrange'a
	//		<< lag_minimum_type << "\n";       // Typ minimum (lokalne/globalne) dla Lagrange'a
	//}

	// Zamknięcie pliku po zakończeniu
	file.close();
}

void lab2()
{
	//srand(time(nullptr));
	//ofstream Sout("1.csv");
	//double s = 0.6, alphaHJ = 0.5, alphaR = 1.3, beta = 0.5, epsilon = 1e-3,Nmax=1000;
	//matrix x0(2, 1, 0.5);
	//matrix s0(2, 1, s);
	//matrix res(100, 2);
	//matrix resHJ(100, 4);
	//matrix resR(100, 4);
	//for (int i = 0; i < 1; i++){
	//	x0(0) = -1.0 + (2.0 * rand() / RAND_MAX);
	//	x0(1) = -1.0 + (2.0 * rand() / RAND_MAX);
	//	res(i, 0) = x0(0);
	//	res(i, 1) = x0(1);
	//	cout << endl << "hj" << endl;
	//	solution hj = HJ(ff2T, x0, s, alphaHJ, epsilon, Nmax);
	//	resHJ(i, 0) = hj.x(0);
	//	resHJ(i, 1) = hj.x(1);
	//	resHJ(i, 2) = hj.y(0);
	//	resHJ(i, 3) = solution::f_calls;
	//	cout << endl << "rosen" << endl;
	//	solution rosen = Rosen(ff2T, x0, s0, alphaR, beta, epsilon, Nmax);
	//	resR(i, 0) = rosen.x(0);
	//	resR(i, 1) = rosen.x(1);
	//	resR(i, 2) = rosen.y(0);
	//	resR(i, 3) = solution::f_calls;
	//}
	//return;
	//matrix excel = hcat(res, resHJ);
	//excel = hcat(excel, resR);
	////Sout << excel;
	//cout << HJ(ff2T, x0, s, alphaHJ, epsilon, Nmax)<<endl;
	//cout<< Rosen(ff2T, x0, s0, alphaR, beta, epsilon, Nmax)<<endl;

	//cout << "Problem rzeczywisty:" << endl;
	//x0 = 10 * rand_mat(2, 1);
	//cout << x0 << endl << endl;
	//cout << HJ(ff2R, x0, s, alphaHJ, epsilon, Nmax) << endl;
	//cout << Rosen(ff2R, x0, s0, alphaR, beta, epsilon, Nmax) << endl;
	//solution::clear_calls();
}
void lab3()
{
	//matrix x0 = matrix(2, 1, 1.0);

	//double c0 = 2;
	//matrix a[3] = { 4, 4.4934, 5 };

	//const double epsilon = 1e-3;
	//const int Nmax = 10000;
	//
	//double pkt1[100] = { 0 };
	//double pkt2[100] = { 0 };


	////for (int i = 0; i < 100; i++) {
	////	do
	////		x0 = 5 * rand_mat(2, 1) + 1; // generowanie punktow znajdujacych sie w obszarze funkcji
	////	while (norm(x0) > a[ktory]);

	////	pkt1[i] = x0(0);
	////	pkt2[i] = x0(1);
	////}
	////// Wygenerowane wartosci sa dluzsze niz wyswietla
	////// Prowadzi to do niedokladnosci - wyswietlany wynik moze sie roznic od tego dla dokladnie takiej wartosci

	/////**/
	////for (int i = 0; i < 100; i++) {
	////	x0(0) = pkt1[i];
	////	x0(1) = pkt2[i];

	////	cout << x0(0) << ";" << x0(1) << ";";

	////	solution zew_opt = pen(fun3, x0, c0, 2, epsilon, Nmax, a[ktory]);
	////	cout << zew_opt.x(0) << ";" << zew_opt.x(1) << ";" << norm(zew_opt.x) << ";" << zew_opt.y[0] << solution::f_calls << ";";
	////	solution::clear_calls();

	////	solution wew_opt = pen(fun3, x0, c0, 0.5, epsilon, Nmax, a[ktory]);
	////	cout << wew_opt.x(0) << ";" << wew_opt.x(1) << ";" << norm(wew_opt.x) << ";" << wew_opt.y[0] << solution::f_calls << endl;
	////	solution::clear_calls();
	////}


	//// Problem rzeczywisty
	//
	//matrix x1 = matrix(2, 1);
	//x1(0) = 1.; // wybieram punkty dowolnie
	//x1(1) = -3.;

	//cout << pen(fR3, x1, c0, 2, epsilon, Nmax) << endl;
}

void lab4()
{
	//matrix x0 = matrix(2, 1, 0.0);
	//const double epsilon = 1e-3;
	//const int Nmax = 10000;

	//double h;

	//double pkt1[100] = { 0 };
	//double pkt2[100] = { 0 };

	//// cout << SD(fT4, grad4, x0, h, epsilon, Nmax) << endl;
	//// cout << CG(fT4, grad4, x0, h, epsilon, Nmax) << endl;
	//// cout << Newton(fT4, grad4, hesj4, x0, h, epsilon, Nmax) << endl;


	//
	////for (int i = 0; i < 100; i++) {

	////	x0 = 20 * rand_mat(2, 1) - 10; // generowanie punktow znajdujacych sie w obszarze funkcji
	////	pkt1[i] = x0(0);
	////	pkt2[i] = x0(1);
	////}

	////int counter;
	////counter = 1;
	//////counter = 100;
	////for(int i = 0; i < counter; i++)
	////{
	////	//x0(0) = pkt1[i];
	////	//x0(1) = pkt2[i];

	////	x0(0) = 5.5275;
	////	x0(1) = -1.29;

	////	cout << x0(0) << " " << x0(1) << " ";


	////	//h = 0.05;
	////	//h = 0.12;
	////	h = -1;

	////	solution SD_sol = SD(fT4, grad4, x0, h, epsilon, Nmax);
	////	cout << SD_sol.x(0) << " " << SD_sol.x(1) << " " << SD_sol.y[0] << SD_sol.f_calls << " " << SD_sol.g_calls << " ";
	////	solution::clear_calls();

	////	solution CG_sol = CG(fT4, grad4, x0, h, epsilon, Nmax);
	////	cout << CG_sol.x(0) << " " << CG_sol.x(1) << " " << CG_sol.y[0] << CG_sol.f_calls << " " << CG_sol.g_calls << " ";
	////	solution::clear_calls();

	////	solution Nw_sol = Newton(fT4, grad4, hesj4, x0, h, epsilon, Nmax);
	////	cout << Nw_sol.x(0) << " " << Nw_sol.x(1) << " " << Nw_sol.y[0] << Nw_sol.f_calls << " " << Nw_sol.g_calls << " " << Nw_sol.H_calls << " ";
	////	solution::clear_calls();



	////	cout << endl;
	////}
	//


	///*
	//x0(0) = 6.2231;
	//x0(1) = -0.238002;

	//cout << Newton(fT4, grad4, hesj4, x0, -1, epsilon, Nmax);
	//*/



	//matrix x1(3, 1, 0.0); // wyliczone z wcześniejszych

	////h = 0.01;
	////h = 0.001;
	//h = 0.0001;

	//solution Real_sol = CG(fR4, gf, x1, h, 0.000001, Nmax);

	//int m = 100;
	//static matrix X(3, m), Y(1, m);
	//ifstream in("XData.txt");
	//in >> X;
	//in.close();
	//in.open("YData.txt");
	//in >> Y;
	//in.close();

	//double P = 0.0;

	//for (int i = 0; i < 100; i++) {
	//	double h = 1.0 / (1 + exp(-(trans(Real_sol.x) * X[i])()));
	//	if (lroundf(h) == Y(0, i))
	//		h = 1;
	//	else
	//		h = 0;

	//	P += h;
	//}
	//P /= m;

	//cout << Real_sol.x(0, 0) << " " << Real_sol.x(1, 0) << " " << Real_sol.x(2, 0) << " " << Real_sol.y(0, 0) << " " << P << " " << Real_sol.g_calls << endl;
}

void lab5()
{
	////int N = 2, Nmax = 10000, mi = 20, lambda = 40;
	////double epsilon = 1e-3;

	////matrix limits(2, 2), sigma0(2, 1);
	////limits(0, 0) = limits(1, 0) = -5;
	////limits(0, 1) = limits(1, 1) = 5;
	////sigma0(0) = sigma0(1) = 1;

	////solution solEvo;	

	////string s = "NIE";
	////sigma0(0) = sigma0(1) = 100; // sigma - 0.01 / 0.1 / 1 / 10 / 100
	////for (int i = 0; i < 100; ++i) {
	////	solution::clear_calls();
	////	solEvo = EA(ff5T, N, limits, mi, lambda, sigma0, epsilon, Nmax);

	////	if (solEvo.y < 0.001)
	////		s = "TAK";

	////	cout << solEvo.x(0) << ";" << solEvo.x(1) << ";" << solEvo.y << solEvo.f_calls << ";" << s << endl;
	////	s = "NIE";
	////}
	//int N = 2, Nmax = 10000, mi = 20, lambda = 40;
	//double epsilon = 1e-3;

	//matrix limits(2, 2), sigma0(2, 1);
	//limits(0, 0) = limits(1, 0) = -5;
	//limits(0, 1) = limits(1, 1) = 5;
	//sigma0(0) = sigma0(1) = 1;

	//solution solEvo;
	////cout << solEvo.x(0) << " " << solEvo.x(1) << " " << solEvo.y << " " << solEvo.f_calls << endl;

	//
	//string s;
	////sigma0(0) = sigma0(1) = 0.01;
	////sigma0(0) = sigma0(1) = 0.1;
	////sigma0(0) = sigma0(1) = 1;
	////sigma0(0) = sigma0(1) = 10;
	//sigma0(0) = sigma0(1) = 100;

	/*for (int i = 0; i < 100; ++i) {
		solution::clear_calls();
		s = "NIE";

		solEvo = EA(ff5T, N, limits, mi, lambda, sigma0, epsilon, Nmax);

		if (solEvo.y < 0.001)
			s = "TAK";

		cout << solEvo.x(0) << " " << solEvo.x(1) << " " << solEvo.y << " " << solEvo.f_calls << " " << s << endl;
	}*/
	
}

void lab6()
{

}
