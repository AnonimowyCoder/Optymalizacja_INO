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
		lab6();
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

	double s = 0.1; // dl kroku HJ
	matrix s0(2, 1, s);
	double alphaHJ = 0.5; //alpha HJ
	double epsilon = 0.000001; //dokladnosc wyniku -> HJ.y < eps;
	double Nmax = 10000; //Max wykonan
	double alphaR = 1.13;
	double beta = 0.5;
	matrix x(2, 1);
	matrix wynik_test_HJ(100,4);
	//Generowanie ziarno
	srand(time(NULL));

	//Plik csv
	ofstream testfunkc("testfunc.csv");
	if (!testfunkc.is_open())throw string("PLIK CSV NIE OTWARTY");

	// Nagłówki w pliku CSV
	testfunkc << "Dl korku,Lp.,x1(0),x2(0),"
		<< "HJ_x1,HJ_x2,y,Liczba wywolan funkcji celu,Minimum globalne,"
		<< "R_x1,R_x2,y,Liczba wywolan funkcji celu,Minimum globalne\n";

	//100 optymalizacji:
	//for (int i = 0; i < 100; i++) {
	//	//generowanie punktow startowych <-1,1>
	//	x(0) = (((double)rand() / RAND_MAX) * 2) - 1;
	//	x(1) = (((double)rand() / RAND_MAX) * 2) - 1;

	//	//cout << "HJ-test start" << endl;
	//	//cout << "x1: " << x(0) << "\nx2: " << x(1) << endl;
	//	solution HJ_test = HJ(ff2T, x, s, alphaHJ, epsilon, Nmax);
	//	//cout << "HJ test done" << endl;

	//	//sprawdzenie czy globalne czy lokalne
	//	string HJ_czy = (HJ_test.x(0) < 0.1 && HJ_test.x(0) > -0.1 && HJ_test.x(1) < 0.1 && HJ_test.x(1) > -0.1) ? "TAK" : "NIE" ; 		// Określenie, czy minimum jest lokalne czy globalne

	//	//cout << "Osiagnieta f(x1,x2): " << HJ_test.y << "\n";
	//	//cout << " x1 = " << HJ_test.x(0)
	//	//	<< " x2 = " << HJ_test.x(1) << "\n";
	//	//cout << "Minimum " << HJ_czy << endl;

	//	testfunkc << s << "," << i + 1 << "," << x(0) << "," << x(1) << ","
	//		<< HJ_test.x(0) << "," << HJ_test.x(1) << "," << HJ_test.y(0) << ","<< solution::f_calls << "," << HJ_czy << "\n";
	//	solution::clear_calls();
	//}
	
	

	for (int i = 0; i < 100; i++) {
		//generowanie punktow startowych <-1,1>
		x(0) = (((double)rand() / RAND_MAX) * 2) - 1;
		x(1) = (((double)rand() / RAND_MAX) * 2) - 1;

		//cout << "R-test start" << endl;
		//cout << "x1: " << x(0) << "\nx2: " << x(1) << endl;
		solution R_test = Rosen(ff2T, x,s0, alphaR, beta, epsilon, Nmax);
		//cout << "R test done" << endl;

		//sprawdzenie czy globalne czy lokalne
		string R_czy = (R_test.x(0) < 0.1 && R_test.x(0) > -0.1 && R_test.x(1) < 0.1 && R_test.x(1) > -0.1) ? "TAK" : "NIE"; 		// Określenie, czy minimum jest lokalne czy globalne

		/*	cout << "Osiagnieta f(x1,x2): " << HJ_test.y << "\n";
			cout << " x1 = " << HJ_test.x(0)
				<< " x2 = " << HJ_test.x(1) << "\n";
			cout << "Minimum " << HJ_czy << endl;*/
		cout << R_test << endl;

		testfunkc << s << "," << i + 1 << "," << x(0) << "," << x(1) << ","
			<< R_test.x(0) << "," << R_test.x(1) << "," << R_test.y(0) << "," << solution::f_calls << "," << R_czy << "\n";

	}

	testfunkc.close();

}

void lab3()
{
	matrix x0 = matrix(2, 1);
	srand(time(NULL));

	//Plik csv
	ofstream testfunkc("testfunc.csv");
	if (!testfunkc.is_open())throw string("PLIK CSV NIE OTWARTY");

	// Nagłówki w pliku CSV
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
	matrix parametr_a[3] = {4.0, 4.4934, 5.0}; //wysylany do ud1
	solution wynik_zewn, wynik_wewn;
	int wybrany = 0;
	double c = 2.0;

	//100 optymalizacji:
	//for (int i = 0; i < 100; i++) {
	//	//break;
	//	do // g(x1, x2) <= a
	//	{
	//		x0(0) = (((double)rand() / RAND_MAX) * 5) + 1;
	//		x0(1) = (((double)rand() / RAND_MAX) * 5) + 1;

	//	} while (norm(x0) > parametr_a[wybrany]);

	//	cout << "Wylosowane pkt startowe:\n";
	//	cout << x0 << endl;
	//	
	//	wynik_zewn = pen(ff3T, x0, c, 2, epsilon, Nmax, parametr_a[wybrany]);
	//	cout << wynik_zewn;

	//	testfunkc << m2d(parametr_a[wybrany]) << "," << i + 1 << "," << x0(0) << "," << x0(1) << ","
	//		<< wynik_zewn.x(0) << "," << wynik_zewn.x(1) << "," << norm(wynik_zewn.x) << "," << wynik_zewn.y(0) << "," << solution::f_calls << ",";// << "\n";
	//	solution::clear_calls();

	//	wynik_wewn = pen(ff3T, x0, c, 0.5, epsilon, Nmax, parametr_a[wybrany]);
	//	testfunkc	<< wynik_wewn.x(0) << "," << wynik_wewn.x(1) << "," << norm(wynik_wewn.x) << "," << wynik_wewn.y(0) << "," << solution::f_calls << "," << "\n";
	//	solution::clear_calls();
	//}
	//testfunkc.close();

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
	
	matrix theta(3, 1, 0.0); 
	double h;
	int Nmax = 1e4;          

	//h = 0.01;
	//h = 0.001;
	h = 0.0001;

	solution Real_sol = CG(ff4R, gf4, theta, h, 0.000001, Nmax);

	const int m = 100;      
	static matrix X(3, m), Y(1, m);

	ifstream in("XData.txt");
	in >> X;
	in.close();

	in.open("YData.txt");
	in >> Y;
	in.close();

	double P = 0.0; 

	for (int i = 0; i < m; ++i) {
		double h = 1.0 / (1 + exp(-(trans(Real_sol.x) * X[i])()));
		P += (lroundf(h) == Y(0, i)) ? 1.0 : 0.0;
	}

	P /= m;

	cout << "Theta: " << Real_sol.x(0, 0) << " " << Real_sol.x(1, 0) << " " << Real_sol.x(2, 0)
		<< "\nKoszt: " << Real_sol.y(0, 0)
		<< "\nDokładność: " << P
		<< "\nLiczba wywołań gradientu: " << Real_sol.g_calls
		<< endl;

}

void lab5()
{
	//Testy

	//double tab_test[] = { 2.01822,1.95171};
	//matrix x0 = matrix(2, tab_test);
	////ustawienie macierzy parametrów
	//matrix ud1(2, 1, 0.0);
	//ud1.set_row(1.0, 0);
	//ud1.set_row(0.5, 1);

	////cout << ud1(0, 0); //a =1
	////cout << ud1(1, 0); //w =0.5

	//const double epsilon = 1e-4; //dokladnosc
	//const int Nmax = 10000; // liczba wywolan

	//solution test;
	//test = Powell(ff5T_comb, x0, epsilon, Nmax, ud1);
	////test = Powell(ff4T, x0, epsilon, Nmax, ud1);
	////test = Powell(trivial, x0, epsilon, Nmax, ud1);

	//cout << "##################################" << endl;
	//cout << "Wyniki: " << endl;
	//cout << test.x << endl;
	//cout << "WYNIKI: " << endl;
	//cout << test.y; 

	//cout << ff5T_1(test.x, ud1) << endl;
	//cout << ff5T_2(test.x, ud1) << endl;

	//cout << "X.1 = " << test.x(0) << endl;
	//cout << "X.2 = " << test.x(1) << endl;


	//------------------------------------
	//Problem testowy

	//ofstream output("lab5.csv");
	//if (!output.is_open()) throw string("Nie udało się otworzyć pliku lab5.csv");
	//output << "step, w, x1(0), x2(0), x1*, x2*, f1*,f2*, f_calls\n";
	//double w = 0.0;
	//std::random_device rd;
	//std::mt19937 gen(rd());
	//std::uniform_real_distribution<> distrib(-10.0, 10.0);
	//const double epsilon = 1e-6; //dokladnosc
	//const int Nmax = 10000; // liczba wywolan
	//int step = 1;
	//while (w < 1.01) {
	//	solution::clear_calls();
	//	double x1start = distrib(gen);
	//	double x2start = distrib(gen);
	//	double tab_test[] = { x1start, x2start };
	//	matrix x0 = matrix(2, tab_test);
	//	double a = 1;
	//	matrix ud1(2, 1, 0.0);
	//	ud1.set_row(a, 0);
	//	ud1.set_row(w, 1);
	//	solution result = Powell(ff5T_comb, x0, epsilon, Nmax, ud1);
	//	matrix f1 = ff5T_1(result.x, ud1);
	//	matrix f2 = ff5T_2(result.x, ud1);
	//	output << step << "," << w << "," << x1start << "," << x2start << ',' << result.x(0) << "," << result.x(1) << "," << m2d(f1) << "," << m2d(f2) << "," << solution::f_calls << ",";

	//	solution::clear_calls();
	//	a = 10;
	//	ud1.set_row(a, 0);
	//	ud1.set_row(w, 1);
	//	result = Powell(ff5T_comb, x0, epsilon, Nmax, ud1);
	//	f1 = ff5T_1(result.x, ud1);
	//	f2 = ff5T_2(result.x, ud1);
	//	output << result.x(0) << "," << result.x(1) << "," << m2d(f1) << "," << m2d(f2) << "," << solution::f_calls << ",";

	//	solution::clear_calls();
	//	a = 100;
	//	ud1.set_row(a, 0);
	//	ud1.set_row(w, 1);
	//	result = Powell(ff5T_comb, x0, epsilon, Nmax, ud1);
	//	f1 = ff5T_1(result.x, ud1);
	//	f2 = ff5T_2(result.x, ud1);
	//	output << result.x(0) << "," << result.x(1) << "," << m2d(f1) << "," << m2d(f2) << "," << solution::f_calls << "\n";

	//	step++;
	//	w += 0.01;
	//}

	//----------------------------------
	//Rzeczywisty Problem

	matrix x0(2, 1); // Macierz startowa
	matrix ud1(1);   // Dodatkowe parametry
	double epsilon = 1e-3; // Dokładność
	int Nmax = 30000; // Maksymalna liczba iteracji

	std::ofstream Sout("symulacja.csv"); // Plik wyjściowy
	if (!Sout.is_open()) throw string("Nie udało się otworzyć pliku symulacja.csv");

	// Nagłówki kolumn
	Sout << "w,l(0) [mm],d(0) [mm],l* [mm],d* [mm],masa* [kg],ugiecie* [mm],Liczba wywolan funkcji celu\n";

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> distribI(200.0, 1000.0);
	std::uniform_real_distribution<> distribD(10.0, 50.0);


	// Pętla po różnych wartościach w
	for (double w = 0.0; w <= 1.01; w += 0.01) {
		// Losowe wartości początkowe dla x0
		x0(0) = distribI(gen); // l <200, 1000> [mm]
		x0(1) = distribD(gen); // d <10, 50> [mm]


		// Wyświetlenie wartości początkowych
		std::cout << "x0: " << x0 << std::endl;

		// Optymalizacja metodą Powella
		solution result1 = Powell(ff5R, x0, epsilon, Nmax, w);

		// Zapis wyników do pliku
		Sout << w << ","                  // w
			<< x0(0) << ","              // l(0) [mm]
			<< x0(1) << ","              // d(0) [mm]
			<< result1.x(0) << ","       // l* [mm]
			<< result1.x(1) << ","       // d* [mm]
			<< result1.y(0) << ","       // masa* [kg]
			<< result1.y(1) << ","       // ugięcie* [mm]
			<< solution::f_calls << "\n"; // Liczba wywołań funkcji celu

		// Wyczyszczenie liczby wywołań funkcji celu
		solution::clear_calls();
	}

	// Zamknięcie pliku
	Sout.close();

	//std::ofstream sim("symulacja.csv"); // Plik wyjściowy
	//if (!sim.is_open()) throw string("Nie udało się otworzyć pliku symulacja.csv");

	//// Nagłówki kolumn
	//sim << "w,l(0) [mm],d(0) [mm],l* [mm],d* [mm],masa* [kg],ugiecie* [mm],Liczba wywolan funkcji celu\n";
	//double epsilon = 1e-3; // Dokładność
	//int Nmax = 30000; // Maksymalna liczba iteracji


	//std::random_device rd;
	//std::mt19937 gen(rd());
	//std::uniform_real_distribution<> distribI(200.0, 1000.0);
	//std::uniform_real_distribution<> distribD(10.0, 50.0);
	//matrix x0 = matrix(2, 1);
	//matrix ud1(1);
	//double w = 0.0;
	//while (w < 1.01) {
	//	solution::clear_calls();
	//	x0(0) = distribI(gen); // l <200, 1000> [mm]
	//	x0(1) = distribD(gen); // d <10, 50> [mm]
	//	solution result = Powell(ff5R, x0, epsilon, Nmax, ud1);
	//	sim << w << "," << x0(0) << "," << x0(1) << "," << result.x(0) << "," << result.x(1) << "," << result.y(0) << "," << result.y(1) << "," << solution::f_calls << "\n";
	//	w += 0.01;
	//}
}

matrix load_data(const std::string& filename) {
	std::ifstream file(filename);
	if (!file.is_open()) throw std::runtime_error("Nie można otworzyć pliku");

	matrix data;
	file >> data; // Wczytanie całej macierzy z pliku
	return data;
}

matrix fileToMatrix(int rows, int cols, std::string filename) {
	matrix result = matrix(rows, cols);

	ifstream file(filename);
	if (!file.is_open()) {
		throw string("Nie mozna otworzyæ pliku: " + filename);
	}

	int i = 0;
	std::string line;
	while (std::getline(file, line))
	{
		int j = 0;
		std::istringstream lineStream(line);
		std::string value;
		while (lineStream >> value) {
			result(i, j) = std::stod(value);
			++j;
		}
		++i;
	}

	return result;
}

matrix loadExperimentalData(const std::string& filename)
{
	std::ifstream file(filename);
	if (!file.is_open())
	{
		std::cerr << "Nie można otworzyć pliku: " << filename << std::endl;
		return matrix(); // Zwraca pustą macierz
	}

	std::vector<double> x1_vals;
	std::vector<double> x2_vals;
	std::string line;

	while (std::getline(file, line))
	{
		std::replace(line.begin(), line.end(), ',', '.'); // Zamiana ',' na '.'
		std::replace(line.begin(), line.end(), ';', ' '); // Zamiana ';' na ' '

		std::istringstream iss(line);
		double x1, x2;
		if (iss >> x1 >> x2)
		{
			x1_vals.push_back(x1);
			x2_vals.push_back(x2);
		}
	}

	file.close();

	// Tworzenie macierzy na podstawie wczytanych danych
	int rows = x1_vals.size();
	matrix data_matrix(rows, 2); // Tworzymy macierz o wymiarach rows x 2

	for (int i = 0; i < rows; ++i)
	{
		data_matrix(i, 0) = x1_vals[i]; // Pierwsza kolumna
		data_matrix(i, 1) = x2_vals[i]; // Druga kolumna
	}

	return data_matrix;
}

void lab6()
{
	//double sigma_tab[] = { 0.01, 0.1, 1, 10, 100 };
	//int N = 2;

	//matrix lb(N, 1);
	//lb(0) = -5;
	//lb(1) = -5;

	//matrix ub(N, 1);
	//ub(0) = 5;
	//ub(1) = 5;

	//int mi = 20;
	//int lambda = 40;
	//double epsilon = 1e-3;
	//int Nmax = 10000;

	//ofstream results("optimization_results.csv");
	//if (!results.is_open()) throw string("Nie można otworzyć pliku CSV");

	//results << "x1*,x2*,y*,f_calls,minimum globalne\n";

	//for (int t = 0; t < 5; t++)
	//{
	//    for (int i = 0; i < 100; i++)
	//    {
	//        solution result = EA(ff6T, N, lb, ub, mi, lambda, sigma_tab[t], epsilon, Nmax);
	//        results << result.x(0) << "," << result.x(1) << "," << result.y(0) << "," << solution::f_calls << "," << (result.y(0) < epsilon ? "TAK" : "NIE") << "\n";
	//        solution::clear_calls();
	//    }
	//}

	int N = 2;
	matrix lb(N, 1), ub(N, 1);
	for (int i = 0; i < N; i++)
	{
		lb(i, 0) = -10.0;
		ub(i, 0) = 10.0;
	}

	// Parametry algorytmu
	int mi = 5;
	int lambda = 10;
	double eps = 1e-3;
	int Nmax = 10000;

	N = 2; // Liczba zmiennych decyzyjnych: b1 i b2
	matrix lb_real(N, 1), ub_real(N, 1);
	lb_real(0, 0) = 0.1; // Dolna granica dla b1
	lb_real(1, 0) = 0.1; // Dolna granica dla b2
	ub_real(0, 0) = 3.0; // Górna granica dla b1
	ub_real(1, 0) = 3.0; // Górna granica dla b2

	double sigmaValues[] = { 0.01, 0.1, 1, 10, 100 };
	// Parametry algorytmu
	mi = 5;
	lambda = 10;
	eps = 1e-2; // 1e-1 fast calculations and ok results, 1e-2 better result but 3 mins calculation
	Nmax = 10000;

	cout << "\nOptymalizacja problemu rzeczywistego:\n";
	// OPTYMALIZACJA
	matrix sigma0_real(N, 1, sigmaValues[0]);

	solution::clear_calls();

	matrix ud1_real;
	matrix ud2_real = loadExperimentalData("polozenia.txt"); // Experimantal data

	solution bestSolReal = EA(ff6R, N, lb_real, ub_real, mi, lambda, sigma0_real, eps, Nmax, ud1_real, ud2_real);

	double b1 = bestSolReal.x(0, 0);
	double b2 = bestSolReal.x(1, 0);
	double J = bestSolReal.y(0, 0);
	int calls = solution::f_calls;
	cout << "b1*, b2*, J(b*), Liczba wywolan funkcji celu\n";
	cout << b1 << ", " << b2 << ", " << J << ", " << calls << "\n";

}
