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
double* expansion(matrix(*ff)(matrix, matrix, matrix), double _x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };
		solution x0(_x0);
		solution x1(_x0 + d);

		std::vector<solution> x_val;
		x_val.push_back(x0.x);
		x_val.push_back(x1.x);

		x_val[0].fit_fun(ff, ud1, ud2);
		x_val[1].fit_fun(ff, ud1, ud2);

		if (x_val[1].y(0) == x_val[0].y(0))
		{
			p[0] = x_val[0].x(0); p[1] = x_val[1].x(0);
			return p;
		}
		if (x_val[1].y(0) > x_val[0].y(0))
		{
			d = -d;
			x_val[1].x = _x0 + d;

			x_val[1].fit_fun(ff, ud1, ud2);
			if (x_val[1].y(0) >= x_val[0].y(0))
			{
				p[0] = x_val[1].x(0); p[1] = x_val[0].x(0) - d;
				return p;
			}
		}

		int i = 0;

		while (true)
		{
			if (solution::f_calls > Nmax)
			{
				delete[] p;
				return NULL;
			}

			i++;

			x_val.push_back(_x0 + pow(alpha, i) * d);
			x_val[i + 1].fit_fun(ff, ud1, ud2);
			if (x_val[i].y(0) <= x_val[i + 1].y(0))
			{
				break;
			}
		}


		if (d > 0)
		{
			p[0] = x_val[i - 1].x(0); p[1] = x_val[i + 1].x(0);
			return p;
		}
		else
		{
			p[0] = x_val[i + 1].x(0); p[1] = x_val[i - 1].x(0);
			return p;
		}
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
// dla h0 = -1 mamy zmiennokrokow¹ metodê
// wykorzystuje metodê ekspancji i z³otego podzia³u (golden)
solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		int i = 0;
		int wym = get_len(x0);	//wymiar problemu
		solution X_0(x0), X_i, gold;
		double* zakres;
		matrix ud_min(wym, 2); //macierz do poszukiwañ w jednym wymiarze
		do {
			X_0.grad(gf, ud1, ud2); //obliczenie gradientu
			cout << "wartosc X_0.g = " << X_0.g << endl << endl << endl;

			if (h0 == -1) {		//przypadek zmiennokrokowy
				ud_min.set_col(X_0.x, 0);
				ud_min.set_col(-1 * X_0.g, 1);

				zakres = expansion(ff, 0, 1, 1.2, Nmax, ud1, ud_min);	//wyznaczenie zakresu dla metody golden
				cout << "Wyznaczony zakres to: " << zakres[0] << ",  " << zakres[1] << endl;
				gold = golden(ff, zakres[0], zakres[1], epsilon, Nmax, ud1, ud_min);
				cout << "wartosc kroku gold.x = " << gold.x << endl << endl;

				X_i.x = X_0.x - gold.x * X_0.g;			//wyznaczenie kolejnego po³o¿enia x
			}
			else 	//przypadek sta³okrokowy
				X_i.x = X_0.x - h0 * X_0.g;	//wyliczanie kolejnego po³o¿enia x

			i++;

			if (norm(X_i.x - X_0.x) < epsilon) {
				X_i.fit_fun(ff, ud1, ud2);
				return X_i;
			}

			X_0.x = X_i.x;

		} while (i < Nmax);

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
		ofstream FILE("OUTPUT_CG.txt");
		//solution Xopt;
		int n = get_len(x0);
		solution X, X1;
		X.x = x0;
		matrix d(n, 1), P(n, 2);
		solution h;
		double* ab, beta;
		//X.grad(gf);
		//d = -X.g;
		d = -X.grad(gf, ud1, ud2);
		while (true) {
			if (h0 < 0) {
				//cout << "here" << endl;
				//P[0] = X.x;
				//P[1] = d;
				P.set_col(X.x, 0);
				P.set_col(d, 1);
				ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
				h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, P);
				X1.x = X.x + h.x * d;
			}
			else {
				X1.x = X.x + h0 * d;
			}

			//log(X1);

			if (norm(X1.x - X.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax) {
				X1.fit_fun(ff, ud1);
				X1.flag = 0;
				FILE.close();
				return X1;
			}
			X1.grad(gf);
			beta = pow(norm(X1.g), 2) / pow(norm(X.g), 2);
			d = -X1.g + beta * d;
			X = X1;
		}
		//return Xopt;
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
		int i = 0;
		int wym = get_len(x0);	//wymiar problemu
		solution X_0(x0), X_i, gold;
		double* zakres;
		matrix ud_min(wym, 2); //macierz do poszukiwañ w jednym wymiarze
		matrix step(2, 1);
		do {
			step = -inv(X_0.hess(Hf, ud1, ud2)) * X_0.grad(gf, ud1, ud2);

			if (h0 == -1) {		//przypadek zmiennokrokowy
				ud_min.set_col(X_0.x, 0);
				ud_min.set_col(-1 * X_0.g, 1);

				zakres = expansion(ff, 0, 1, 1.2, Nmax, ud1, ud_min);	//wyznaczenie zakresu dla metody golden
				gold = golden(ff, zakres[0], zakres[1], epsilon, Nmax, ud1, ud_min);
				cout << "wartosc kroku gold.x = " << gold.x << endl << endl;
				X_i.x = X_0.x - gold.x * X_0.g;			//wyznaczenie kolejnego po³o¿enia x
			}
			else 	//przypadek sta³okrokowy
				X_i.x = X_0.x + h0 * step;	//wyliczanie kolejnego po³o¿enia x

			i++;

			if (norm(X_i.x - X_0.x) < epsilon) {
				X_i.fit_fun(ff, ud1, ud2);
				return X_i;
			}

			X_0.x = X_i.x;

		} while (i < Nmax);

	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}


// a, b - poczatkowy przedzial poszukiwan
// epsilon - dok³adnoœæ obliczeñ

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int i = 0;
		double alfa = (sqrt(5) - 1) / 2;
		solution a_sol(a);
		solution b_sol(b);
		solution c_sol(b - alfa * (b - a));
		solution d_sol(a + alfa * (b - a));
		c_sol.fit_fun(ff, ud1, ud2);
		d_sol.fit_fun(ff, ud1, ud2);

		while ((b_sol.x - a_sol.x) > epsilon && i < Nmax) {
			if (c_sol.y < d_sol.y) {
				b_sol.x = d_sol.x;
				d_sol.x = c_sol.x;
				c_sol.x = b_sol.x - alfa * (b_sol.x - a_sol.x);
				c_sol.fit_fun(ff, ud1, ud2);
				d_sol.fit_fun(ff, ud1, ud2);
			}
			else {
				a_sol.x = c_sol.x;
				c_sol.x = d_sol.x;
				d_sol.x = a_sol.x + alfa * (b_sol.x - a_sol.x);
				c_sol.fit_fun(ff, ud1, ud2);
				d_sol.fit_fun(ff, ud1, ud2);
			}
			i++;
		}

		if (i >= Nmax) {
			throw string("Przekroczono maksymaln¹ liczbê iteracji Nmax");
		}
		Xopt.x = (a_sol.x + b_sol.x) / 2.0;
		Xopt.fit_fun(ff, ud1, ud2);

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}


//LAB 5
solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		const int n_wym = 2;
		solution Xopt;
		int i = 0;

		//narzędzia do poszukiwań w jednym wymiarze
		double* zakres;		//zakres znaleziony metodą ekspansji
		matrix ud_min(n_wym, 2); //macierz do poszukiwan w jednym wymiarze
		solution gold;		//rozwiązanie optymalne na kierunku (długość kroku)

		//utworzenie wersorów
		matrix tab_dj[n_wym];
		// j=0; j=1
		for (int j = 0; j < n_wym; j++) {
			matrix e_j(n_wym, 1, 0.0);
			e_j.set_row(1.0, j);
			tab_dj[j] = e_j;
			//cout << "utworzone wersor nr " << j << ":  " << tab_dj[j]<<endl;
		}

		// punkty pomocnicze w tablicy
		matrix tab_p[4];	//n_wym+2
		for (int j = 0; j < 4; j++) {
			matrix p(n_wym, 1, 0.0);
			tab_p[j] = p;
		}

		//glowna petla 
		solution X(x0);
		do {
			tab_p[0] = X.x;

			for (int j = 1; j <= n_wym; j++) {
				//optymalizacja po kiernku
				//j=1, j=2
				ud_min.set_col(X.x, 0);		//punkt początkowy
				ud_min.set_col(tab_dj[j - 1], 1);		//kierunek (wykorzystujemy wersory)

				//cout << "(pierwsze)wyswietlam ud_min do obliczen po kierunku : " << endl << ud_min << endl;

				zakres = expansion(ff, 0, 1, 1.2, Nmax, ud1, ud_min);	//wyznaczenie zakresu dla metody golden
				//cout << "(pierwsze)Wyznaczony zakres to: " << zakres[0] << ",  " << zakres[1] << endl;
				gold = golden(ff, zakres[0], zakres[1], epsilon, Nmax, ud1, ud_min);
				//cout << "(pierwsze)wartosc kroku gold.x = " << gold.x << endl << endl;

				//obliczenie punktu pj
				tab_p[j] = tab_p[j - 1] + gold.x() * tab_dj[j - 1];
				//cout << "obliczony punkt p[j] = " << endl << tab_p[j] << endl;
			}
			// sprawdzenie warunku zakończenia
			if (norm(tab_p[n_wym] - X.x) < epsilon) {
				X.fit_fun(ff, ud1, ud2);
				return X;
			}

			//wyznaczenie nowych kierunków
			tab_dj[0] = tab_dj[1];

			//wyznaczony nowy kierunk
			tab_dj[1] = (tab_p[n_wym] - tab_p[0]) * (1 / norm(tab_p[n_wym] - tab_p[0]));

			//wyznaczenie nowej długość kroku

			ud_min.set_col(tab_p[n_wym], 0);		//punkt początkowy
			ud_min.set_col(tab_dj[1], 1);		//nowy kierunek

			//cout << "(drugie)wyswietlam ud_min do obliczen po kierunku: " << endl << ud_min << endl;

			zakres = expansion(ff, 0, 1, 1.2, Nmax, ud1, ud_min);	//wyznaczenie zakresu dla metody golden
			//cout << "(drugie)Wyznaczony zakres to: " << zakres[0] << ",  " << zakres[1] << endl;
			gold = golden(ff, zakres[0], zakres[1], epsilon, Nmax, ud1, ud_min);
			//cout << "(drugie)wartosc kroku gold.x = " << gold.x << endl << endl;
			tab_p[n_wym + 1] = tab_p[n_wym] + gold.x() * tab_dj[1];
			X.x = tab_p[n_wym + 1];
			//cout << endl << "NOWY PUNKT X.x TO: " << endl << X.x << endl;
			i++;
		} while (X.f_calls < Nmax);
		cout << "ZWRACAM Xopt\n";
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}


//LAB 6

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
