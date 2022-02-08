#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

double x_i(const int i, double L, const double delx)
{
	return -0.5 * L + (delx * (i + 1.));
}

void resetA(gsl_matrix* A, const int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			gsl_matrix_set(A, i, j, 0);
		}
	}
}

void setA(gsl_matrix* A, int alfa, const int n, const double N, const double L, const double delx)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			double xi = x_i(i, L, delx);
			double p = 1 + (4. * alfa * xi * xi);
			if (i == j)
				gsl_matrix_set(A, i, j, 2. * N / (p * delx * delx));
			if (i == j - 1)
				gsl_matrix_set(A, i, j, -1. * N / (p * delx * delx));
			if (i == j + 1)
				gsl_matrix_set(A, i, j, -1. * N / (p * delx * delx));
		}
	}
}

void evalf(std::ofstream* name, gsl_vector_complex* eval, int alfa)
{
	for (int i = 0; i < 6; i++)
	{
		gsl_complex cval = gsl_vector_complex_get(eval, i);
		double val = GSL_REAL(cval);
		*name << alfa << "\t" << std::fixed << std::setprecision(6) << sqrt(val);
	}

	*name << "\n";
}

void evecf(std::ofstream* name, gsl_matrix_complex* evec, const int n, const double delx, const double L)
{
	for (int i = 0; i < 6; i++)
	{
		*name << i << " " << i << "\n";
		*name << std::fixed << std::setprecision(6) << -5.0 << " " << 0.0 << "\n";

		for (int j = 0; j < n; j++)
		{
			gsl_complex cv = gsl_matrix_complex_get(evec, j, i);
			*name << std::fixed << std::setprecision(6) << x_i(j, L, delx) << " " << GSL_REAL(cv) << "\n";
		}

		*name << std::fixed << std::setprecision(6) << 5.0 << " " << 0.0 << "\n";
	}
}

int main()
{
	std::ofstream file; 
	std::ofstream file2; 
	std::ofstream file3; 

	file.open("results.txt");
	file2.open("first.txt");
	file3.open("second.txt");

	const double L = 10.0, N = 1.0;
	const int n = 200;
	const double delx = L / (n + 1.0);

	gsl_matrix* A = gsl_matrix_calloc(n, n);
	gsl_vector_complex* eval = gsl_vector_complex_calloc(n); 
	gsl_matrix_complex* evec = gsl_matrix_complex_calloc(n, n); 
	gsl_eigen_nonsymmv_workspace* w = gsl_eigen_nonsymmv_alloc(n);

	for (int alfa=0; alfa <= 100; alfa += 2)
	{
		resetA(A, n);
		setA(A, alfa, n, N, L, delx);
		gsl_eigen_nonsymmv(A, eval, evec, w);
		gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

		evalf(&file, eval, alfa);

		if (alfa == 0)
		{
			evecf(&file2, evec, n, delx, L);
		}
		else if (alfa == 100)
		{
			evecf(&file3, evec, n, delx, L);
		}
	}

	gsl_matrix_free(A);
	gsl_vector_complex_free(eval);
	gsl_matrix_complex_free(evec);
	gsl_eigen_nonsymmv_free(w);

	file.close();
	file2.close();
	file3.close();

	return 0;
}

