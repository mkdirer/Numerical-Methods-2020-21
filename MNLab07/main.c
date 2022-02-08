#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double Wn(double *x_m, double **y, int n, double x) {
	double result = 0.0;

	for (int j  = 0; j <= n; j++) {
		double product = 1.0;

		for (int i = 0; i < j; i++) {
			product *= (x - x_m[i]);
		}

		result += y[j][j] * product;
	}

	return result;
}

int main() {
	double deltaX;
	double xMin = -5.0, xMax = 5.0;

	//FILE * f = fopen("Newton.dat", "w");
	FILE * f = fopen("Czebyszew.dat", "w");
	FILE * f2 = fopen("x_m.dat", "w");
	for (int n = 5; n <= 20; n += 5) {
		double *x = malloc((n + 1) * sizeof(double));
		double **y = malloc((n + 1) * sizeof(double *));
		for (int i = 0; i < n + 1; i++)
			y[i] = malloc((n+1) * sizeof(double));

		deltaX = (xMax - xMin) / n;

		for (int i = 0; i < n + 1; i++) {
			//x[i] = xMin + i * deltaX;
			x[i] = ((xMin - xMax) * cos(M_PI * (2.0 * i + 1.0) / (2.0 * n + 2.0)) + (xMin + xMax)) / 2.0;
			y[i][0] = 1.0 / (1.0 + x[i] * x[i]);
			fprintf(f2, "%10g %10g\n", x[i], y[i][0]); 
		}
		fprintf(f2, "\n"); 

		for(int j = 1; j <= n; j++) 
			for (int i = j; i <= n; i++) 
				y[i][j] = (y[i][j-1] - y[i-1][j-1]) / (x[i] - x[i-j]);


		for (double i = xMin; i <= xMax; i += 0.01) {
			double W = Wn(x, y, n, i);
			fprintf(f, "%10g %10g %10g\n", i, 1/(1 + i * i), W); 
		}

	fprintf(f, "\n"); 

		free(x);
		for (int i = 0; i < n + 1; i++)
			free(y[i]);
		free(y);
	}

	fclose(f);
	fclose(f2);
	return 0;
}

