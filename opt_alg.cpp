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
			cout << "START\n";
			cout << Xopt << endl;
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
				} while (Xopt.y >= xB.y);
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
		} while (s < epsilon);

		Xopt = xB; 
		Xopt.flag = 1;
		cout << "Xopt\n" << Xopt << endl;
		cout << "xB\n" << xB << endl;
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
		solution tmp(XB); //tworzenie tmp soultion 
		int n = get_dim(XB);  //n = wymiar problemu
		matrix jednostkowa = ident_mat(n); //tworzenie macierzy jendostkowej
		for (int i = 0; i < n; i++) {	//szukanie minimum w otoczeniu o promieniu s // moze i =0
			tmp = XB.x + s * jednostkowa(i, i);
			tmp.fit_fun(ff, ud1, ud2);
			if (tmp.y < XB.y) //Jezeli y jest mniejszy po kroku w przod
				XB = tmp;
			else {
				tmp = XB.x - s * jednostkowa(i, i); //Jezeli y jest mniejszy po kroku w tyl
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
		// Inicjalizacja podstawowych zmiennych
		solution::clear_calls();
		int n = get_dim(x0);
		matrix d = ident_mat(n);
		matrix lambda(n, 1), p(n, 1), s(s0);
		solution xB(x0);
		xB.fit_fun(ff, ud1, ud2);

		while (true)
		{
			// Iteracja po wymiarach
			for (int j = 0; j < n; j++)
			{
				solution next_xB(xB.x + s(j) * d[j]);
				next_xB.fit_fun(ff, ud1, ud2);
				if (next_xB.y < xB.y)
				{
					xB = next_xB;
					lambda(j) += s(j);
					s(j) *= alpha;
				}
				else
				{
					s(j) *= -beta;
					p(j)++;
				}
			}

			// Sprawdzenie czy trzeba zmienić kierunek poszukiwań
			bool change = true;
			for (int j = 0; j < n; j++)
			{
				if (lambda(j) == 0 || p(j) == 0)
				{
					change = false;
					break;
				}
			}

			if (change)
			{
				matrix Q(n, n);
				for (int j = 0; j < n; j++)
				{
					for (int k = 0; k <= j; k++)
					{
						Q(j, k) = lambda(j);
					}
				}

				Q = Q * d;
				matrix v = Q[0];
				d.set_col(v / norm(v), 0);
				for (int j = 1; j < n; j++)
				{
					matrix sum(n, 1);
					for (int k = 0; k < j; k++)
					{
						sum = sum + (trans(Q[j]) * d[k]) * d[k];
					}
					v = Q[j] - sum;
					d.set_col(v / norm(v), j);
				}

				s = s0;
				lambda = matrix(n, 1);
				p = matrix(n, 1);
			}

			// Warunek zakończenia poszukiwań
			if (solution::f_calls > Nmax)
			{
				xB.flag = 0;
				return xB;
			}

			double max_step = abs(s(0));
			for (int j = 1; j < n; ++j)
			{
				if (max_step < abs(s(j)))
					max_step = abs(s(j));
			}

			if (max_step < epsilon)
			{
				xB.flag = 1;
				break;
			}
		}

		return xB;
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
		//Tu wpisz kod funkcji

		return Xopt;
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
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

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
