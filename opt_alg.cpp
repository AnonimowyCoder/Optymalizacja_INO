#include"opt_alg.h"

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

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

///////////////////////////////
// ff -> fit function (funkcja celu)
// x0 -> punkt startowy
// s -> pocz�tkowa d�ugo�� kroku
// alpha -> wsp�czynnik zmniejszenia kroku
// epsilon -> dok�adno��
// Nmax -> maksymalna liczba wywo�a� funkcji celu
solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt(x0);
		solution XB;
		solution XB_d;
		do {
			XB = Xopt;
			// etap roboczy
			Xopt = HJ_trial(ff, XB, s, ud1, ud2);
			cout << "przed fit fun, po HJ_trial" << endl << "Xopt y = " << Xopt.y << endl;
			XB.fit_fun(ff);
			Xopt.fit_fun(ff);
			cout << "przed if, po fit fun " << endl << "Xopt y = " << Xopt.y << endl;
			if (Xopt.y < XB.y) {
				do {
					XB_d = XB;
					XB = Xopt;
					Xopt.x = (2 * XB.x) - XB_d.x;
					Xopt = HJ_trial(ff, Xopt, s, ud1, ud2);
					Xopt.fit_fun(ff);
					XB.fit_fun(ff);
				} while (Xopt.y < XB.y);
				Xopt.x = XB.x;
			}
			else
				s = alpha * s;
		} while (s > epsilon && Nmax > Xopt.f_calls);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}


// ff -> fit function (funkcja celu)
// XB -> aktualna warto�� x
// s -> d�ugo�� kroku
solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		int dim = get_dim(XB);
		cout << "dim = " << dim << endl;

		for (int j = 0; j < dim; j++) {
			// utworzenie wersora
			matrix e_j(dim, 1, 0.0);
			cout << "tu dziala, iteracja nr: " << j << endl;
			e_j.set_row(1.0, j);
			cout << "XB.x: " << endl << XB.x << endl;
			cout << "e_j: " << endl << e_j << endl;

			//utworzenie dodatkowych solution
			solution XB_plus(XB);
			solution XB_minus(XB);

			//wyliczenie bazowego po�o�enia x
			XB.fit_fun(ff);

			//wyliczenie nowego po�o�enia x
			XB_plus.x = XB.x + (e_j * s);
			cout << "Po operacji dodania wersora, XB_plus.x: " << endl << XB_plus.x << endl;
			XB_plus.fit_fun(ff);

			XB_minus.x = XB.x - (e_j * s);
			XB_minus.fit_fun(ff);

			//poprawa w kierunku e_j
			if (XB_plus.y < XB.y) {
				cout << "blok if (XB_plus.y < XB.y)" << endl;
				XB.x = XB.x + (e_j * s);
			}
			//poprawa w kierunku -e_j
			else {

				if (XB_minus.y < XB.y) {
					cout << "blok if (XB_minus.y < XB.y)" << endl;
					XB.x = XB.x - (e_j * s);
				}
			}
		}
		return XB;
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
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		solution result(x0);
		matrix c_as_ud2(2, new double[2] {c, dc});
		double alpha = 1,
			beta = 0.5,
			gamma = 2,
			delta = 0.5,
			s = 0.5;

		while (true) {
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
		p[0].fit_fun(ff, ud1, ud2); // Obliczenie warto�ci funkcji celu dla punktu pocz�tkowego 

		//Obliczenie funkcji celu dla pozosta�ych punkt�w podczas inicjalizacji
		for (int i = 1; i < n + 1; ++i) {
			p[i].x = p[0].x + s * jednostkowa[i - 1];  // Macierz jednostkowa [i-1]
			p[i].fit_fun(ff, ud1, ud2);  // Obliczenie funkcji celu dla ka�dego wierzcho�ka
		}

		solution p_odb, p_z, p_e;
		int max, min;

		while (true) {
			max = min = 0;

			//Por�wnanie funkcji celu
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
			// Kryterium zako�czenia
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

//Lab_4 
///////////////////////////
// dla h0 = -1 mamy zmiennokrokowa metode
// wykorzystuje metode ekspancji i zlotego podzialu (golden)
solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		int i = 0;
		int wym = get_len(x0);	//wymiar problemu
		solution X_0(x0), X_i, gold;
		double* zakres;
		matrix ud_min(wym, 2); //macierz do poszukiwan w jednym wymiarze
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

				X_i.x = X_0.x - gold.x * X_0.g;			//wyznaczenie kolejnego polozenia x
			}
			else 	//przypadek stalokrokowy
				X_i.x = X_0.x - h0 * X_0.g;	//wyliczanie kolejnego polozenia x

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



solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	try
	{
		int i = 0;
		int wym = get_len(x0);	//wymiar problemu
		solution X_0(x0), X_i, gold;
		double* zakres;
		matrix ud_min(wym, 2); //macierz do poszukiwa� w jednym wymiarze
		matrix step_0(2, 1); //zmienna przechowuj�ca kierunek kroku
		step_0 = -X_0.grad(gf, ud1, ud2); //obliczenie gradientu

		do {

			if (h0 == -1) {		//przypadek zmiennokrokowy
				ud_min.set_col(X_0.x, 0);
				ud_min.set_col(step_0, 1);

				zakres = expansion(ff, 0, 1, 1.2, Nmax, ud1, ud_min);	//wyznaczenie zakresu dla metody golden
				cout << "Wyznaczony zakres to: " << zakres[0] << ",  " << zakres[1] << endl;
				gold = golden(ff, zakres[0], zakres[1], epsilon, Nmax, ud1, ud_min);
				cout << "wartosc kroku gold.x = " << gold.x << endl << endl;

				X_i.x = X_0.x + gold.x * step_0;			//wyznaczenie kolejnego po�o�enia x
			}
			else 	//przypadek sta�okrokowy
				X_i.x = X_0.x + h0 * step_0;	//wyliczanie kolejnego po�o�enia x

			i++;

			if (norm(X_i.x - X_0.x) < epsilon) {
				X_i.fit_fun(ff, ud1, ud2);
				return X_i;
			}

			X_0.x = X_i.x;
			// wyznaczenie nowego kierunku
			matrix grad_i(2, 1);
			grad_i = -X_0.grad(gf, ud1, ud2);

			double beta = pow(norm(grad_i), 2) / pow(norm(step_0), 2);
			step_0 = grad_i + beta * step_0;


		} while (i < Nmax);

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
		matrix ud_min(wym, 2); //macierz do poszukiwa� w jednym wymiarze
		matrix step(2, 1);
		do {
			step = -inv(X_0.hess(Hf, ud1, ud2)) * X_0.grad(gf, ud1, ud2);

			if (h0 == -1) {		//przypadek zmiennokrokowy
				ud_min.set_col(X_0.x, 0);
				ud_min.set_col(-1 * X_0.g, 1);

				zakres = expansion(ff, 0, 1, 1.2, Nmax, ud1, ud_min);	//wyznaczenie zakresu dla metody golden
				gold = golden(ff, zakres[0], zakres[1], epsilon, Nmax, ud1, ud_min);
				cout << "wartosc kroku gold.x = " << gold.x << endl << endl;
				X_i.x = X_0.x - gold.x * X_0.g;			//wyznaczenie kolejnego po�o�enia x
			}
			else 	//przypadek sta�okrokowy
				X_i.x = X_0.x + h0 * step;	//wyliczanie kolejnego po�o�enia x

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
// epsilon - dok�adno�� oblicze�

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		//solution Xopt;
		double alfa = (sqrt(5) - 1) / 2;
		solution A, B, C, D;
		A.x = a;
		B.x = b;
		C.x = B.x - alfa * (B.x - A.x);
		C.fit_fun(ff, ud1, ud2);
		D.x = A.x + alfa * (B.x - A.x);
		D.fit_fun(ff, ud1, ud2);

		while (true) {
			if (C.y < D.y) {
				B = D;
				D = C;
				C.x = B.x - alfa * (B.x - A.x);
				C.fit_fun(ff, ud1, ud2);
			}
			else {
				A = C;
				C = D;
				D.x = A.x + alfa * (B.x - A.x);
				D.fit_fun(ff, ud1, ud2);
			}
			if (B.x - A.x < epsilon || solution::f_calls > Nmax) {
				A.x = (A.x + B.x) / 2;
				A.fit_fun(ff, ud1, ud2);
				A.flag = 0;
				return A;
			}
		}

		//return Xopt;
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
				ud_min.set_col(tab_dj[j-1], 1);		//kierunek (wykorzystujemy wersory)
				
				cout << "(pierwsze)wyswietlam ud_min do obliczen po kierunku : " <<endl<< ud_min << endl;

				zakres = expansion(ff, 0, 1, 1.2, Nmax, ud1, ud_min);	//wyznaczenie zakresu dla metody golden
				cout << "(pierwsze)Wyznaczony zakres to: " << zakres[0] << ",  " << zakres[1] << endl;
				gold = golden(ff, zakres[0], zakres[1], epsilon, Nmax, ud1, ud_min);
				cout << "(pierwsze)wartosc kroku gold.x = " << gold.x << endl << endl;

				//obliczenie punktu pj
				tab_p[j] = tab_p[j - 1] + gold.x() * tab_dj[j - 1];
				cout << "obliczony punkt p[j] = " << endl<<tab_p[j] << endl;
			}
			// sprawdzenie warunku zakończenia
			if (norm(tab_p[n_wym] - X.x) < epsilon) {
				X.fit_fun(ff, ud1, ud2);
				return X;
			}

			//wyznaczenie nowych kierunków
			tab_dj[0] = tab_dj[1];
			
			//wyznaczony nowy kierunk
            tab_dj[1] = (tab_p[n_wym] - tab_p[0])* (1/ norm(tab_p[n_wym] - tab_p[0]));

			//wyznaczenie nowej długość kroku

			ud_min.set_col(tab_p[n_wym], 0);		//punkt początkowy
			ud_min.set_col(tab_dj[1], 1);		//nowy kierunek

			cout << "(drugie)wyswietlam ud_min do obliczen po kierunku: " << endl << ud_min << endl;

			zakres = expansion(ff, 0, 1, 1.2, Nmax, ud1, ud_min);	//wyznaczenie zakresu dla metody golden
			cout << "(drugie)Wyznaczony zakres to: " << zakres[0] << ",  " << zakres[1] << endl;
			gold = golden(ff, zakres[0], zakres[1], epsilon, Nmax, ud1, ud_min);
			cout << "(drugie)wartosc kroku gold.x = " << gold.x << endl << endl;
			tab_p[n_wym + 1] = tab_p[n_wym] + gold.x() * tab_dj[1];
			X.x = tab_p[n_wym + 1];
			cout <<endl<< "NOWY PUNKT X.x TO: " << endl << X.x << endl;
			i++;
		} while (X.f_calls<Nmax);

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