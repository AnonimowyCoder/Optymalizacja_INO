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
			cout << i << "fib: " << Bi.x - Ai.x << endl;
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
			cout << "lag: " << Bi.x - Ai.x << endl;
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



solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2) {
	try {
		solution::clear_calls();
		solution Xopt=x0;
		solution xB;
		solution xZ;
		while (true) {
			cout << Xopt.x(0) << ";" << Xopt.x(1) << endl;
			xB = Xopt;
			xB.fit_fun(ff, ud1, ud2);
			Xopt = HJ_trial(ff, xB, s, ud1, ud2);
			if (Xopt.y < xB.y) {
				while (true) {
					xZ = xB;
					xB = Xopt;
					Xopt = 2 * xB.x - xZ.x;
					Xopt.fit_fun(ff,ud1,ud2);
					Xopt = HJ_trial(ff, Xopt, s, ud1, ud2);
					if (solution::f_calls > Nmax) {
						Xopt.flag = 0;
						cout << "Error" << endl;
						break;
					}
					if (Xopt.y >= xB.y)break;
				}
				Xopt = xB;
			}
			else s = alpha * s;
			if (solution::f_calls > Nmax) {
				cout << "Error" << endl;
				Xopt.flag = 0;
				break;
			}
			if (s < epsilon) {
				Xopt.flag = 1;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info) {
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution X, double s, matrix ud1, matrix ud2){
	try{
		solution xB(X);
		int n = get_dim(X);
		matrix e = ident_mat(n);
		for (int j = 0; j < n; j++) {
			xB = X.x+ s * e[j];
			xB.fit_fun(ff);
			if (xB.y < X.y)X = xB;
			else {
				xB = X.x - s * e[j];
				xB.fit_fun(ff);
				if (xB.y < X.y)X = xB;
			}
		}
		return X;
	}
	catch (string ex_info){
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}


solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2){
	try
	{
		solution::clear_calls();
		int n = get_dim(x0);
		matrix d = ident_mat(n),
			lambda(n, 1),
			p(n, 1),
			s(s0);			
		solution xB(x0);
		xB.fit_fun(ff, ud1, ud2);
		while (true) {
			cout << xB.x(0) << ";" << xB.x(1) << endl;
			for (int j = 0; j < n; j++) {
				solution next_xB(xB.x + s(j) * d[j]);
				next_xB.fit_fun(ff, ud1, ud2);
				if (next_xB.y < xB.y) {
					xB = next_xB;
					lambda(j) += s(j);
					s(j) *= alpha;
				}
				else {
					s(j) *= -beta;
					p(j)++;
				}
			}
			bool change = true;
			for (int j = 0; j < n; j++) {
				if (lambda(j) == 0 || p(j) == 0) {
					change = false;
					break;
				}
			}
			if (change) {
				matrix Q(n, n);
				for (int j = 0; j < n; j++) {
					for (int k = 0; k <= j; k++)Q(j, k) = lambda(j);
				}
				Q = Q * d;
				matrix v(n, 1);
				v = Q[0];
				d.set_col(v / norm(v), 0);
				for (int j = 1; j < n; j++) {
					matrix sum(n, 1);
					for (int k = 0; k < j; k++)sum = sum + (trans(Q[j]) * d[k]) * d[k];
					v = Q[j] - sum;
					d.set_col(v / norm(v), j);
				}
				s = s0;			
				lambda = matrix(n, 1);
				p = matrix(n, 1);		
			}
			if (solution::f_calls > Nmax) {
				xB.flag = 0;
				return xB;
			}
			double max = abs(s(0));
			for (int j = 1; j < n; ++j) {
				if (max < abs(s(j))) max = abs(s(j));
			}
			if (max < epsilon) {
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


solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		//solution Xopt;
		double alpha = 1, beta = 0.5, gamma = 2, delta = 0.5, s = 0.5;
		solution X(x0), X1;
		matrix c0(2, new double[2]{ c,dc });
		while (true) {
			X1 = sym_NM(ff, X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c0);
			if (norm(X1.x - X.x) < epsilon || solution::f_calls > Nmax) {
				X1.flag = 0;
				return X1;
			}
			X = X1;
			c0(0) = c0(0) * dc;
		}
		//return Xopt;
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
		//solution Xopt;
		int n = get_len(x0);
		matrix D = ident_mat(n);
		int N = n + 1;
		solution* S = new solution[N];
		S[0].x = x0;		
		S[0].fit_fun(ff, ud1, ud2);

		for (int i = 1; i < N; ++i) {
			S[i].x = S[0].x + s * D[i - 1]; 
			S[i].fit_fun(ff, ud1, ud2);
		}

		solution PR, PE, PN;
		matrix pc;
		int i_min, i_max;

		while (true) {
			i_min = i_max = 0;
			for (int i = 1; i < N; ++i) {
				if (S[i].y(0) < S[i_min].y(0))
					i_min = i;
				if (S[i].y(0) > S[i_max].y(0))
					i_max = i;
			}

			pc = matrix(n, 1);

			for (int i = 0; i < N; ++i)
				if (i != i_max)
					pc = pc + S[i].x;

			pc = pc / (N - 1);
			PR.x = pc + alpha * (pc - S[i_max].x);
			PR.fit_fun(ff, ud1, ud2);

			if (PR.y(0) < S[i_max].y(0) && S[i_min].y(0) <= PR.y(0))
				S[i_max] = PR;				
			else if (PR.y(0) < S[i_min].y(0))	{
				PE.x = pc + gamma * (PR.x - pc);
				PE.fit_fun(ff, ud1, ud2);
				
				if (PE.y(0) < PR.y(0))
					S[i_max] = PE;
				else
					S[i_max] = PR;
			}
			else {
				PN.x = pc + beta * (S[i_max].x - pc);
				PN.fit_fun(ff, ud1, ud2);
				if (PN.y(0) < S[i_max].y(0))
					S[i_max] = PN;
				else {
					for (int i = 0; i < N; ++i)
						if (i != i_min) {
							S[i].x = delta * (S[i].x + S[i_min].x);
							S[i].fit_fun(ff, ud1, ud2);
						}
				}
			}
			double max_s = norm(S[i_min].x - S[0].x);

			for (int i = 1; i < N; ++i)
				if (max_s < norm(S[i_min].x - S[i].x))
					max_s = norm(S[i_min].x - S[i].x);

			if (max_s < epsilon)
				return S[i_min];
		}
		//return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}
// show_every_iteration
//#define log(X1) cout << X1.x(0) << " " << X1.x(1) << endl;
#define log(X1) FILE << X1.x(0) << " " << X1.x(1) << endl;
//#define log(X1) ;

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		ofstream FILE("OUTPUT_SD.txt");
		//solution Xopt;
		solution X, X1;
		X.x = x0;
		int n = get_len(x0);
		matrix d(n, 1), P(n, 2);
		solution h;
		double* ab;
		while (true) {
			//X.grad(gf);
			//d = -X.g;
			d = -X.grad(gf, ud1, ud2);
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

			log(X1);

			if (norm(X1.x - X.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax) {
				X1.fit_fun(ff, ud1, ud2);
				X1.flag = 0;
				FILE.close();
				return X1;
			}
			X = X1;
		}

		//return Xopt;
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

			log(X1);

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
		ofstream FILE("OUTPUT_Newton.txt");
		//solution Xopt;
		int n = get_len(x0);
		solution X, X1;
		X.x = x0;
		matrix d(n, 1), P(n, 2);
		solution h;
		double* ab;
		while (true) {
			X.grad(gf);
			X.hess(Hf);
			d = -inv(X.H) * X.g;
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

			log(X1);

			if (norm(X1.x - X.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax) {
				X1.fit_fun(ff, ud1);
				X1.flag = 0;
				FILE.close();
				return X1;
			}
			X = X1;
		}

		//return Xopt;
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
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

//solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix limits, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
//{
//	try
//	{
//		//solution Xopt;
//		//Tu wpisz kod funkcji
//		//mi liczebnosc populacji
//		solution* P = new solution[mi + lambda];
//		solution* Pm = new solution[mi];
//		default_random_engine gen;
//		gen.seed(static_cast<unsigned int>(chrono::system_clock::now().time_since_epoch().count()));
//		normal_distribution<double> distr(0.0, 1.0);
//		matrix IFF(mi, 1), temp(N, 2); //IFF macierz z przystosowaniami //temp - kopia osobnikia
//		double r, s, s_IFF;
//		double tau = pow(2 * N, -0.5), tau1 = pow(2 * pow(N, 0.5), -0.5); //tau, tau1 - mutacja
//		int j_min; // najlepsze rozwiazanie
//		for (int i = 0; i < mi; ++i)
//		{
//			P[i].x = matrix(N, 2);
//			for (int j = 0; j < N; ++j)
//			{
//				P[i].x(j, 0) = (limits(j, 1) - limits(j, 0)) * rand_mat(1, 1)() + limits(j, 0);
//				P[i].x(j, 1) = sigma0(j);
//			}
//			P[i].fit_fun(ff, ud1, ud2);
//			if (P[i].y < epsilon)
//				return P[i];
//
//		}
//		while (true)
//		{
//			s_IFF = 0;
//			for (int i = 0; i < mi; ++i)
//			{
//				IFF(i) = 1 / P[i].y();
//				s_IFF += IFF(i);
//			}
//			for (int i = 0; i < lambda; ++i)
//			{
//				r = s_IFF * rand_mat(1, 1)();
//				s = 0;
//				for (int j = 0; j < mi; ++j)
//				{
//					s += IFF(j);
//					if (r <= s)
//					{
//						P[mi + i] = P[j]; // j - wylosowany osobnik
//						break;
//					}
//				}
//			}
//			//mutajca
//			for (int i = 0; i < lambda; ++i)
//			{
//				r = distr(gen);
//				for (int j = 0; j < N; ++j)
//				{
//					P[mi + i].x(j, 1) *= exp(tau1 * r + tau * distr(gen));
//					P[mi + i].x(j, 0) += P[mi + i].x(j, 1) * distr(gen);
//				}
//			}
//			//krzyzowanie
//			for (int i = 0; i < lambda; i += 2)
//			{
//				r = rand_mat(1, 1)();
//				temp = P[mi + i].x;  //jeden z rodzicow
//				P[mi + i].x = r * P[mi + i].x + (1 - r) * P[mi + i + 1].x;  //pierwszy potomek
//				P[mi + i + 1].x = r * P[mi + i + 1].x + (1 - r) * temp;  //drugi potomek
//			}
//			//ocena osobnikow
//			for (int i = 0; i < lambda; ++i)
//			{
//				P[mi + i].fit_fun(ff, ud1, ud2);
//				if (P[mi + i].y < epsilon) 	//ocena rozwiazania
//					return P[mi + i];
//
//			}
//			//wskazanie najelpszych osobnikow
//			for (int i = 0; i < mi; ++i)
//			{
//				j_min = 0;
//				for (int j = 1; j < mi + lambda; ++j)
//					if (P[j_min].y > P[j].y)
//						j_min = j;
//				Pm[i] = P[j_min];
//				P[j_min].y = 1e10;
//			}
//			for (int i = 0; i < mi; ++i)
//				P[i] = Pm[i];  //P[i] najlepsza populacja
//			if (solution::f_calls > Nmax)
//				return P[0];
//
//		}
//		//return Xopt;
//	}
//	catch (string ex_info)
//	{
//		throw ("solution EA(...):\n" + ex_info);
//	}
//}
