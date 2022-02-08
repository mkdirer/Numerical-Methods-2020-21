#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double licz_r(double* a, double* b, int n, double x0)
{
	b[n] = 0.0;
	for(int z = n - 1; z>=0; --z)
  {
    b[z] = a[z+1] + x0*b[z+1];
  }
		
	return (double)(a[0] + x0*b[0]);  
}

int main()
{
  //plik wynikowy
  FILE *plik1 = fopen("wyniki.dat", "w+");
  
  //zmienne
  int N=5;
  int I_MAX=30;
	double* a = malloc((N+1) * sizeof(double));
  double* b = malloc((N+1) * sizeof(double));
	double* c = malloc(N * sizeof(double));
  a[0] = 240.0;
  a[1] = -196.0;
  a[2] = -92.0;
  a[3] = 33.0;
  a[4] = 14.0;
  a[5] = 1.0;
	double x0, x1, Rj, Rjp;
	int n;

  //obliczenia
  for(int L=1; L<=N; L++)
  {
    n=N-L+1;
    x0=0.0;

    for(int it=1; it<I_MAX; it++)
    {
      Rj = licz_r(a, b, n, x0);
      Rjp = licz_r(b, c, n-1, x0);
      x1 = x0 - Rj/Rjp;

			printf("m.z. x%d w %d kroku = %10g  Rj = %12g  Rjp = %12g\n", L, it, x1,Rj,Rjp);
      fprintf(plik1, "m.z. x%d w %d kroku = %10g  Rj = %12g  Rjp = %12g\n", L, it, x1,Rj,Rjp);

			if(fabs(x1 - x0) < 1.0E-7)
      {
        printf("\n");
        fprintf(plik1, "\n");
        break;
      }
			x0 = x1; 
		}
		
		for(int i = 0; i<=n-1; i++)
    {
      a[i] = b[i];
    }
	}
	// zwolnienie pamięci
	free(a);
	free(b);
  free(c);
	return 0;
}