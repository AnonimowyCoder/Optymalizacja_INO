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
double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };
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
// s -> pocz¹tkowa d³ugoœæ kroku
// alpha -> wspó³czynnik zmniejszenia kroku
// epsilon -> dok³adnoœæ
// Nmax -> maksymalna liczba wywo³añ funkcji celu
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
			cout <<"przed fit fun, po HJ_trial"<<endl << "Xopt y = " << Xopt.y << endl;
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
		} while (s > epsilon&&Nmax>Xopt.f_calls);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}


// ff -> fit function (funkcja celu)
// XB -> aktualna wartoœæ x
// s -> d³ugoœæ kroku
solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		int dim = get_dim(XB);
		cout << "dim = " << dim << endl;

		for (int j = 0; j < dim; j++) {
			// utworzenie wersora
			matrix e_j(dim, 1,0.0);
			cout << "tu dziala, iteracja nr: " << j << endl;
			e_j.set_row(1.0, j);
			cout << "XB.x: "<<endl << XB.x << endl;
			cout << "e_j: " <<endl<< e_j << endl;

			//utworzenie dodatkowych solution
			solution XB_plus(XB);
			solution XB_minus(XB);

			//wyliczenie bazowego po³o¿enia x
			XB.fit_fun(ff);

			//wyliczenie nowego po³o¿enia x
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
		p[0].fit_fun(ff, ud1, ud2); // Obliczenie wartoœci funkcji celu dla punktu pocz¹tkowego 

		//Obliczenie funkcji celu dla pozosta³ych punktów podczas inicjalizacji
		for (int i = 1; i < n + 1; ++i) {
			p[i].x = p[0].x + s * jednostkowa[i - 1];  // Macierz jednostkowa [i-1]
			p[i].fit_fun(ff, ud1, ud2);  // Obliczenie funkcji celu dla ka¿dego wierzcho³ka
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
			// Kryterium zakoñczenia
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
			X_0.grad(gf,ud1,ud2); //obliczenie gradientu
			cout << "wartosc X_0.g = " << X_0.g << endl << endl << endl;

			if (h0 == -1) {		//przypadek zmiennokrokowy
				ud_min.set_col(X_0.x, 0);
				ud_min.set_col(-1 * X_0.g, 1);

				zakres = expansion(ff, 0, 1, 1.2, Nmax, ud1, ud_min);	//wyznaczenie zakresu dla metody golden
				cout << "Wyznaczony zakres to: " << zakres[0] << ",  " << zakres[1] << endl;
				gold = golden(ff, zakres[0], zakres[1], epsilon, Nmax, ud1, ud_min);
				cout << "wartosc kroku gold.x = " << gold.x << endl<<endl;
				
				X_i.x = X_0.x - gold.x * X_0.g;			//wyznaczenie kolejnego po³o¿enia x
			}
			else 	//przypadek sta³okrokowy
				X_i.x = X_0.x - h0 * X_0.g;	//wyliczanie kolejnego po³o¿enia x
			
			i++;

			if (norm(X_i.x - X_0.x) < epsilon) {
				X_i.fit_fun(ff,ud1,ud2);
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
		matrix ud_min(wym, 2); //macierz do poszukiwañ w jednym wymiarze
		matrix step_0(2, 1); //zmienna przechowuj¹ca kierunek kroku
		step_0 = -X_0.grad(gf, ud1, ud2); //obliczenie gradientu

		do {

			if (h0 == -1) {		//przypadek zmiennokrokowy
				ud_min.set_col(X_0.x, 0);
				ud_min.set_col(step_0, 1);

				zakres = expansion(ff, 0, 1, 1.2, Nmax, ud1, ud_min);	//wyznaczenie zakresu dla metody golden
				cout << "Wyznaczony zakres to: " << zakres[0] << ",  " << zakres[1] << endl;
				gold = golden(ff, zakres[0], zakres[1], epsilon, Nmax, ud1, ud_min);
				cout << "wartosc kroku gold.x = " << gold.x << endl << endl;

				X_i.x = X_0.x + gold.x * step_0;			//wyznaczenie kolejnego po³o¿enia x
			}
			else 	//przypadek sta³okrokowy
				X_i.x = X_0.x + h0 * step_0;	//wyliczanie kolejnego po³o¿enia x

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
