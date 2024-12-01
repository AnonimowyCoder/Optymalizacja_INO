#include"opt_alg.h"
#include <iomanip>


// lab 0
solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}


// lab 1
double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0,0 };
		solution X0(x0), X1(x0 + d);
		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);
		if (X0.y == X1.y) {
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);
			return p;
		}
		if (X1.y > X0.y) {
			d = -d;
			X1.x = X0.x + d;
			X1.fit_fun(ff, ud1, ud2);
			if (X1.y >= X0.y) {
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x - d);
				return p;
			}
		}
		
		return p;

	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2) {
	try {
		solution::clear_calls();
		solution Xopt;
		vector<double> fi = { 1,1 };
		int k = 1;
		while (fi[k] < (b - a) / epsilon) {
			fi.push_back(fi[k] + fi[k - 1]);
			k++;
		}
		solution Ai(a),
			Bi(b),
			Ci(b - (fi[k - 1] / fi[k] * (b - a))),
			Di(a + b - Ci.x(0));
		for (int i = 0; i < k - 2; i++) {
			Ci.fit_fun(ff);
			Di.fit_fun(ff);
			if (Ci.y(0) < Di.y(0))Bi = Di;
			else Ai = Ci;
			Ci.x = Bi.x(0) - fi[k - i - 2] / fi[k - i - 1] * (Bi.x(0) - Ai.x(0));
			Di.x = Ai.x(0) + Bi.x(0) - Ci.x(0);
		}
		Xopt = Ci;
		Xopt.fit_fun(ff);
		return Xopt;
	}
	catch (string ex_info) {
		throw ("solution fib(...):\n" + ex_info);
	}
}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2) {
	try {
		solution::clear_calls();
		solution Ai = a,
			Bi = b,
			Ci = (a + b) / 2,
			Xopt = Ci,
			Di = 0;
		int i = 0;
		double l, m;
		do {
			Ai.fit_fun(ff, ud1, ud2); Bi.fit_fun(ff, ud1, ud2); Ci.fit_fun(ff, ud1, ud2);
			l = Ai.y(0) * (pow(Bi.x(0), 2) - pow(Ci.x(0), 2))
				+ Bi.y(0) * (pow(Ci.x(0), 2) - pow(Ai.x(0), 2))
				+ Ci.y(0) * (pow(Ai.x(0), 2) - pow(Bi.x(0), 2));
			m = Ai.y(0) * (Bi.x(0) - Ci.x(0))
				+ Bi.y(0) * (Ci.x(0) - Ai.x(0))
				+ Ci.y(0) * (Ai.x(0) - Bi.x(0));
			if (m <= 0) {
				Xopt.flag = 0;
				break;
			}
			Di = l / 2 / m;
			Di.fit_fun(ff, ud1, ud2);
			if (Ai.x(0) < Di.x(0) && Di.x(0) < Ci.x(0)) {
				if (Di.y(0) < Ci.y(0)) {
					Bi = Ci;
					Ci = Di;
				}
				else Ai = Di;
			}
			else {
				if (Ci.x(0) < Di.x(0) && Di.x(0) < Bi.x(0)) {
					if (Di.y(0) < Ci.y(0)) {
						Ai = Ci;
						Ci = Di;
					}
					else Bi = Di;
				}
				else {
					Xopt.flag = 0;
					break;
				}
			}
			if (Xopt.f_calls > Nmax) {
				Xopt.flag = 0;
				break;
			}
			Xopt = Di;
		} while (Bi.x(0) - Ai.x(0) >= epsilon || abs(Di.x(0) - Xopt.x(0)) >= gamma);
		return Xopt;
	}
	catch (string ex_info) {
		throw ("solution lag(...):\n" + ex_info);
	}
}


// lab 2 x0 - punkt startowy, s - dl kroku, 0 < alpha < 1 -> dl zmniejszania kroku, epslon - dokladnosc, Nmax - l. wyowlan
solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{

		solution::clear_calls();
		solution Xopt = x0;
		solution xB = x0;
		solution x_B;

		do
		{
			//cout << "START\n";
			//cout << Xopt << endl;
			cout << Xopt.x(0) << ";" << Xopt.x(1) << endl;
			xB = Xopt;
			xB.fit_fun(ff, ud1, ud2);
			Xopt = HJ_trial(ff, xB, s, ud1, ud2);
			if (Xopt.y < xB.y) {
				do {
					x_B = xB;
					xB = Xopt;
					Xopt = 2 * xB.x - x_B.x;
					Xopt.fit_fun(ff, ud1, ud2);
					Xopt = HJ_trial(ff, Xopt, s, ud1, ud2);
					if (solution::f_calls > Nmax) {
						Xopt.flag = 0;
						cout << Xopt << endl;
						cout << xB << endl;
						throw string("Liczba wywolan przekracza Nmax");
					}
				} while (Xopt.y < xB.y);
				Xopt = xB;
			}
			else
				s = alpha * s;
			if (solution::f_calls > Nmax) {
				Xopt.flag = 0;
				cout << Xopt << endl;
				cout << xB << endl;
				throw string("Liczba wywolan przekracza Nmax");
			}
		} while (s >= epsilon);

		Xopt = xB; 
		Xopt.flag = 1;
		cout << "Xopt\n" << Xopt << endl;
		//cout << "xB\n" << xB << endl;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}


solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		solution tmp(XB);
		int n = get_dim(XB);  //n = wymiar problemu
		matrix jednostkowa = ident_mat(n); //tworzenie macierzy jendostkowej e(j)
		for (int i = 0; i < n; i++) {	//szukanie minimum w otoczeniu o promieniu s

			tmp = XB.x + s * jednostkowa[i]; // uzywamy i tego wektora z macierzy jednostkowej.
			tmp.fit_fun(ff, ud1, ud2);

			if (tmp.y < XB.y) //Jezeli y jest mniejszy po kroku w przod
				XB = tmp;
			else {
				tmp = XB.x - s * jednostkowa[i]; //Jezeli y jest mniejszy po kroku w tyl
				tmp.fit_fun(ff, ud1, ud2);
				if (tmp.y < XB.y)
					XB = tmp;
			}
		}
		return XB; //zwraca znaleziony pkt i wartosc
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		//inicjalizacja podstawowych zmiennych;
		solution::clear_calls();
		solution Xopt(x0);
		solution xB;

		int n = get_dim(x0); //liczba wymiarow
		matrix jendostkowa = ident_mat(n); //macierz jednostkowa; wersory
		matrix lambda(n, 1); //wartosci przesuniecia
		matrix p(n, 1);	// liczba niepowodzen
		matrix s(s0); // przesuniecia
		//matrix Q(n, n);

		double dlugosc_kroku; //uzywany do zakonczenia petli do...while
		Xopt.fit_fun(ff, ud1, ud2);

		//Rozpoczecie poszukiwan
		do //until:		max (dlugosc_kroku) < epsilon
		{
			//cout << Xopt.x(0) << ";" << Xopt.x(1) << endl;
			for (int j = 0; j < n; j++) { // j - wymiar
				xB = Xopt.x + s(j) * jendostkowa[j]; //jednostkowa[j] = wektor macierzy jednostkowej;
				xB.fit_fun(ff, ud1, ud2);
				if (xB.y < Xopt.y) {
					Xopt = xB; //Xopt = xB.x + s(j) + jendostkowa[j]; tylko bez niepotrzebnego powtorzenia fit_fun();
					lambda(j) = lambda(j) + s(j);
					s(j) = alpha * s(j);
				}
				else{
					s(j) = -beta * s(j); 
					p(j)++; 
				}
			}// wyniki labda(j), s(j) i p(j) sa juz dla kolejnej iteracji i+1;

			//Czy zmieniac kierunek poszukiwan
			bool nowy_kierunek = true;
			for (int j = 0; j < n; j++) {
				if (lambda(j) != 0 && p(j) != 0)
					continue;
				else {
					nowy_kierunek = false;
					break; 
				}
			}
			//Zmiana kierunku
			if (nowy_kierunek) {
				//Zmiana kierunkow bazy d(j)
				//Definicja Q^(i)
				matrix Q(n, n);
				for (int j = 0; j < n; j++) { // Wypelnienie Macierzy Q do skosu wartosciami z wektora lambda(j)
					for (int k = 0; k <= j; k++)
						Q(j, k) = lambda(j);
				}
				Q = Q * jendostkowa;
				//Definicja wektora v
				for (int j = 0; j < n; j++) {
					// Ustalenie pierwszego wektora
					if (j == 0) {
						matrix v = Q[0];
						jendostkowa.set_col(v / norm(v), j);
					}
					else {
						matrix v = Q[j];
						for (int k = 0; k < j; k++) {
							double suma = (trans(v) * jendostkowa[k])(0, 0);//liczenie sumy
							v = v - suma * jendostkowa[k];
						}
						//normalizacja i zapisanie do j-tej kolumny
						jendostkowa.set_col(v / norm(v), j);
					}
				}
				//nowe tablice do poszukiwan
				lambda = matrix(n, 1);
				p = matrix(n, 1);
				s = s0;
			}

			if (solution::f_calls > Nmax) {
				Xopt.flag = 0;
				cout << Xopt<< endl;
				throw string("Liczba wywolan przekracza Nmax");
			}

			//Znalezienie max kroku w tablicy wektorow dlugosci kroku
			dlugosc_kroku = abs(s(0)); //j(|s(j)|)
			for (int j = 1; j < n; ++j) {
				if (dlugosc_kroku < abs(s(j)))
					dlugosc_kroku = abs(s(j));
			}
		} while (dlugosc_kroku >= epsilon);

		Xopt.flag = 1;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

//lab 3
solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		solution result(x0);
		matrix c_as_ud2(2, new double[2]{c, dc});
		double alpha = 1,
			beta = 0.5,
			gamma = 2,
			delta = 0.5,
			s = 0.5;
		
		while(true) {
			Xopt = sym_NM(ff, result.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c_as_ud2); //ud1 to parametr a

			c_as_ud2(0) = dc * c_as_ud2(0);

			if (solution::f_calls > Nmax) {
				Xopt.flag = 0; //0 nie powiodlo sie
				throw string("Przekracza Nmax");
			}
			if (norm(Xopt.x - result.x) < epsilon) {
				Xopt.flag = 1; // 1 -> powiodlo sie
				return Xopt;
			}
			result = Xopt;
		}
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int n = get_dim(x0);

		matrix jednostkowa = ident_mat(n);
		solution* p = new solution[n + 1];
		p[0].x = x0;
		p[0].fit_fun(ff, ud1, ud2); // Obliczenie wartości funkcji celu dla punktu początkowego 

		//Obliczenie funkcji celu dla pozostałych punktów podczas inicjalizacji
		for (int i = 1; i < n + 1; ++i) {
			p[i].x = p[0].x + s * jednostkowa[i - 1];  // Macierz jednostkowa [i-1]
			p[i].fit_fun(ff, ud1, ud2);  // Obliczenie funkcji celu dla każdego wierzchołka
		}

		solution p_odb, p_z, p_e;
		int max, min;

		while (true) {
			max = min = 0;

			//Porównanie funkcji celu
			for (int i = 1; i < n + 1; ++i) {
				if (p[i].y < p[min].y)
					min = i;
				if (p[i].y > p[max].y)
					max = i;
			}

			matrix p_suma(n, 1);
			for (int i = 0; i < n + 1; ++i) { 
				if (i != max)
					p_suma = p_suma + p[i].x; 
			}
			p_suma = p_suma / n;
			p_odb.x = p_suma + alpha * (p_suma - p[max].x); // Punkt odbicia
			p_odb.fit_fun(ff, ud1, ud2); 

			if (p_odb.y < p[min].y) {
				p_e.x = p_suma + gamma * (p_odb.x - p_suma); // Punkt ekspansji
				p_e.fit_fun(ff, ud1, ud2); 
				if (p_e.y < p_odb.y) {
					p[max] = p_e;
				}
				else {
					p[max] = p_odb;
				}
			}
			else {
				if (p[min].y <= p_odb.y && p_odb.y < p[max].y) {
					p[max] = p_odb;
				}
				else {
					p_z.x = p_suma + beta * (p[max].x - p_suma); // Punkt kontrakcji
					p_z.fit_fun(ff, ud1, ud2);
					if (p_z.y >= p[max].y) {
						// Redukcja sympleksu
						for (int i = 0; i < n + 1; ++i) {
							if (i != min) {
								p[i].x = p[min].x + delta * (p[i].x - p[min].x);
								p[i].fit_fun(ff, ud1, ud2); 
							}
						}
					}
					else {
						p[max] = p_z;
					}
				}
			}

			if (solution::f_calls > Nmax) {
				p[min].flag = 0;
				cout << p[min];
				throw string("Liczba wywolan przekracza Nmax");

			}
			// Kryterium zakończenia
			double max_dist = norm(p[min].x - p[0].x); 
			for (int i = 1; i < n + 1; ++i) { 
				double dist = norm(p[min].x - p[i].x);
				if (dist > max_dist)
					max_dist = dist;
			}
			if (max_dist < epsilon) {
				solution result = p[min];
				delete[] p; 
				result.flag = 1;
				return result;
			}
		}
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}



//lab 4
solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
