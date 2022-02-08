#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Dyrektywy zakladajace, ze te trzy pliki sa skopiowane do aktualnego katalogu. */
#include "nrutil.h"
//#include "nrutil.c" // To mozna usunac, jesli plik jest dodany w poleceniu kompilacji.
//#include "gaussj.c" // To tez mozna usunac, jesli plik jest dodany w poleceniu kompilacji.

/* Dyrektywy dla Taurusa (nie wymagaja kopiowania plikow, ale Taurus musi dzialac...) */
// #include "/opt/NR/numerical_recipes.c/nrutil.h"
// #include "/opt/NR/numerical_recipes.c/nrutil.c"
// #include "/opt/NR/numerical_recipes.c/gaussj.c"

#define N 200 // rozmiar macierzy M: NxN

int main(void)
{
	float **M, **b;
	//	Alokacja macierzy
	M = matrix(1, N, 1, N);
	b = matrix(1, N, 1, 1);

	// 	Wypelnienie macierzy M i wektora b
	for (int i = 1; i <= N; ++i)
	{
		b[i][1] = 0.0;
		for (int j = 1; j <= N; ++j)
			M[i][j] = 0.0;
	}
	b[1][1] = 1.0;
	b[2][1] = 0.0*0.1;
	//b[3][1] = -16.0;

  for(int i=1;i<=N; ++i)
  {
    M[i][i] = 1.0;
    if(i+1<=N) M[i+1][i] = 0.1*0.1 -2.0;
    if(i+2 <=N) M[i+2][i] = 1.0;
  }
	//M[1][1] = 2.0;
	//M[1][2] = 3.0;
	//M[1][3] = 0.0;
	M[2][1] = -1.0;
	//M[2][2] = 1.0;
	//M[2][3] = 4.0;
	//M[3][1] = 3.0;
	//M[3][2] = 11.0;
	//M[3][3] = 1.0;

	//	Rozwiazanie ukladu rownan Mx=b - wywolanie procedury:
	gaussj(M, N, b, 1);

  //Zapisanie wynikow do pliku
  FILE *file = fopen("out.txt", "w");
  if(file == NULL)
  {
    //Nie otrzymalismy pliku 
    printf("Nie mozna otworzyc/utworzyc pliku. \n");
    fclose(file);
    exit(EXIT_FAILURE);
  }else{
    for(int i =1; i<=N; ++i)
    fprintf(file, "%g\t%g\n", i*0.1, b[i][1]);
    fclose(file);
  }
	//	Wypisanie rozwiazania, ktore procedura gaussj(M, N, b, 1); zapisala w wektorze b.
	
  for (int i = 1; i <= N; ++i)
		printf("%g\n", b[i][1]);

	//	Zwolnienie pamieci
	free_matrix(M, 1, N, 1, N);
	free_matrix(b, 1, N, 1, 1);

	return 0;
}
