#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"

#define N 3

void printM(char *text, float ** matrix, int size)
{
  printf("%s", text);
  for(int i=1; i<=size; i++)
  {
    for(int j =1; j<=size; j++)
      printf("%9g ", matrix[i][j]);
    printf("\n");
  }
}

//Zapisywanie do pliku 
void fprintM(FILE *file, char *text, float **matrix, int size)
{
  fprintf(file, "%s" , text);
  for(int i=1; i<= size; i++)
  {
    for(int j=1; j<= size; j++)
    {
      fprintf(file, "%9g" , matrix[i][j]);
    }
    fprintf(file, "\n");
  }
}

int main(void)
{
FILE *file = fopen("odp.dat", "w");
//Alokacja macierzy nxn
float **A = matrix(1,N,1,N);
float **B = matrix(1,N,1,N);
int ite = 1;
//Uzupelnianie macierzy

for(int i=1; i<=N; i++)
{
  for(int j= 1; j<=N; j++)
  {
    A[i][j] = ite;
    B[i][j] = ite;
    ite++;
  }
}
B[1][1] = 1.1;

//Wypisywanie macierzy
printM("macierz A: \n", A, N);
fprintM(file, "macierz A: \n", A, N);
printM("macierz B: \n", B, N);
fprintM(file, "macierz B: \n", B, N);

//Wektory permutacji oraz wypelnienie ich procedura ludcmp i inicjalizacja znaku permutacji 
int *indA = ivector(1, N);
int *indB = ivector(1, N);
float dA;
float dB;
ludcmp(A,N, indA, &dA);
ludcmp(B,N, indB, &dB);

//Macierz L i U dla macierzy A i dB

float ** LA = matrix(1,N,1,N);
float ** LB = matrix(1, N, 1, N);
float ** UA = matrix(1,N,1,N);
float ** UB = matrix(1, N, 1, N);

for(int i = 1; i <= N; i++) {
		for(int j = 1; j <= N; j ++) {
			if (i == j) {
				LA[i][j] = 1;
				LB[i][j] = 1;
				UA[i][j] = A[i][j];
				UB[i][j] = B[i][j];
			}
			else if(i < j) {
				LA[i][j] = 0;
				LB[i][j] = 0;
				UA[i][j] = A[i][j];
				UB[i][j] = B[i][j];
			}
			else {
				LA[i][j] = A[i][j];
				LB[i][j] = B[i][j];
				UA[i][j] = 0;
				UB[i][j] = 0;
			}
			ite++;
		}
	}

printM("\nmacierz L dla A: \n", LA, N);
fprintM(file, "\nmacierz L dla A: \n", LA, N);
printM("\nmacierz L dla B: \n", LB, N);
fprintM(file, "\nmacierz L dla B: \n", LB, N);
printM("\nmacierz U dla A: \n", UA, N);
fprintM(file, "\nmacierz U dla A: \n", UA, N);
printM("\nmacierz U dla B: \n", UB, N);
fprintM(file, "\nmacierz U dla B: \n", UB, N);

//Alokacja i uzupelnienie wektorow wolnych, na których wykonujem procedure
float *a1 = vector(1, N);
float *a2 = vector(1, N);
float *a3 = vector(1, N);
float *b1 = vector(1, N);
float *b2 = vector(1, N);
float *b3 = vector(1, N);

  a1[1] = 1; a1[2] = 0; a1[3] = 0;
	a2[1] = 0; a2[2] = 1; a2[3] = 0;
	a3[1] = 0; a3[2] = 0; a3[3] = 1;
	b1[1] = 1; b1[2] = 0; b1[3] = 0;
	b2[1] = 0; b2[2] = 1; b2[3] = 0;
	b3[1] = 0; b3[2] = 0; b3[3] = 1;
	lubksb(A, N, indA, a1);
	lubksb(A, N, indA, a2);
	lubksb(A, N, indA, a3);
	lubksb(B, N, indB, b1);
	lubksb(B, N, indB, b2);
	lubksb(B, N, indB, b3);
	//Alokacja, uzupelnienie i wypisanie macierzy odwrotnych
	float ** A_inv = matrix(1, N, 1, N);
	float ** B_inv = matrix(1, N, 1, N);

	for(int i = 1; i <= N; i++) {
        	A_inv[i][1] = a1[i];
        	A_inv[i][2] = a2[i];
        	A_inv[i][3] = a3[i];
        	B_inv[i][1] = b1[i];
        	B_inv[i][2] = b2[i];
        	B_inv[i][3] = b3[i];
    	}

  printM("\nMacierz A^-1: \n", A_inv, N);
  fprintM(file, "\nMacierz A^-1: \n", A_inv, N);
  printM("\nMacierz B^-1: \n", B_inv, N);
  fprintM(file, "\nMacierz B^-1: \n", B_inv, N);

  ite =1;

  for(int i=1; i<=N; i++)
{
  for(int j= 1; j<=N; j++)
  {
    A[i][j] = ite;
    B[i][j] = ite;
    ite++;
  }
}
B[1][1] = 1.1;

//Liczenie norm
float maxA =0;
float maxB =0;
float maxA_inv =0;
float maxB_inv =0;

for(int i = 1; i <= N; i++) {
		for(int j = 1; j <= N; j ++) {
			if(fabs(A[i][j]) > maxA)
				maxA = fabs(A[i][j]);
			if(fabs(B[i][j]) > maxB)
				maxB = fabs(B[i][j]);
			if(fabs(A_inv[i][j]) > maxA_inv)
				maxA_inv= fabs(A_inv[i][j]);
			if(fabs(B_inv[i][j]) > maxB_inv)
				maxB_inv= fabs(B_inv[i][j]);
		}
	}

	float kA = maxA * maxA_inv;
	float kB = maxB * maxB_inv;

	printf("\n||A|| = %g\n||A^-1|| = %f = %g\nkA = %f = %g\n", maxA, maxA_inv, maxA_inv, kA, kA);
	fprintf(file, "\n||A|| = %g\n||A^-1|| = %f = %g\nkA = %f = %g\n", maxA, maxA_inv, maxA_inv, kA, kA);
	printf("\n||B|| = %g\n||B^-1|| = %f\nkB = %f\n", maxB, maxB_inv, kB);	
	fprintf(file, "\n||B|| = %g\n||B^-1|| = %f\nkB = %f\n", maxB, maxB_inv, kB);		
	ite = 1;

  for(int i = 1; i <= N; i++) {
		for(int j = 1; j <= N; j ++) {
			A[i][j] = ite;
			B[i][j] = ite;
			ite++;
		}
	}
	B[1][1] = 1.1;

  //Alokowanie macierzy dla iloczynu macierzy 
  float ** A_ilocz = matrix(1,N,1,N);
  float ** B_ilocz = matrix(1,N,1,N);

  for(int i = 1; i <= N; i++)
		for(int j = 1; j <= N; j++)
			A_ilocz[i][j] = B_ilocz[i][j] = 0;
	int k;
	float sumA = 0, sumB = 0;

	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			sumA = sumB = 0;
			for (int k = 1; k <= N; k++) {
				sumA = sumA + A[i][k] * A_inv[k][j];
				sumB = sumB + B[i][k] * B_inv[k][j];
			}
			A_ilocz[i][j] = sumA;
			B_ilocz[i][j] = sumB;
		}
	}

  printM("\n A x A^-1: \n", A_ilocz, N);
  fprintM(file, "\n A x A^-1: \n", A_ilocz, N);
  printM("\n B x B^-1: \n", B_ilocz, N);
  fprintM(file, "\n B x B^-1: \n", B_ilocz, N);

  fclose(file);
	//Zwalnianie pamieci
	free_matrix(A, 1, N, 1, N);
	free_matrix(B, 1, N, 1, N);
	free_matrix(LA, 1, N, 1, N);
	free_matrix(LB, 1, N, 1, N);
	free_matrix(UA, 1, N, 1, N);
	free_matrix(UB, 1, N, 1, N);
	free_matrix(A_inv, 1, N, 1, N);
	free_matrix(B_inv, 1, N, 1, N);
	free_matrix(A_ilocz, 1, N, 1, N);
	free_matrix(B_ilocz, 1, N, 1, N);
	free_vector(a1, 1, N);
	free_vector(a2, 1, N);
	free_vector(a3, 1, N);
	free_vector(b1, 1, N);
	free_vector(b2, 1, N);
	free_vector(b3, 1, N);
	free_ivector(indA, 1, N);
	free_ivector(indB, 1, N);
	
	return 0;
}