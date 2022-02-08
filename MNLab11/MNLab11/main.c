#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include "nrutil.c"
#include "four1.c"
#include <time.h>
#include <math.h>

float frand();
float Max(float a, float b);


int main() {
  FILE *file1 = fopen("k8.dat", "w");
  FILE *file2 = fopen("k10.dat", "w");
  FILE *file3 = fopen("k12.dat", "w");
  FILE *f[] = {file1, file2, file3};
	srand(time(NULL));
	float T = 1.0;
	const int k[] = {8, 10, 12};
    
  for(int i = 0; i < 3; i++){
    FILE* file = f[i];
    const int N = pow(2, k[i]);
    float* g1 = vector(1, 2 * N);
    float* g2 = vector(1, 2 * N);
    float* f = vector(1, 2 * N);
    float* y = vector(1, 2 * N);
    const float dt = (float)(3 * T / N);
    const float sigm = T / 20.0;
        
        
    for(int i = 1; i <= N; i++) {
      float ti = dt * (i - 1);
      f[2 * i - 1] = y[2 * i - 1] = sin((2 * M_PI / T) * ti) +  sin(2 * (2 * M_PI / T) * ti) + sin(3 * (2 * M_PI / T) * ti) + frand();
	    f[2 * i] = y[2 * i] = g1[2 * i] = g2[2 * i] = 0;
	    g1[2 * i - 1] = g2[2 * i - 1] = (1.0 / (sigm * sqrt(2 * M_PI))) * exp(-(ti * ti) / (2 * sigm * sigm));
	}
	    
	four1(f, N, 1);
	four1(g1, N, 1);
	four1(g2, N, -1);

  for (int i = 1; i <= N; i++) {
	  float a1 = f[2 * i - 1];
	  float b1 = f[2 * i];
	  float a2 = g1[2 * i - 1] + g2[2 * i - 1];
	  float b2 = g1[2 * i] + g2[2 * i];
	  f[2 * i - 1] = a1 * a2 - b1 * b2;
	  f[2 * i] = a1 * b2 + a2 * b1;
	}

	four1(f, N, -1);
  float max = fabs(f[1]);

  for(int i = 2; i <= N; i++) {
    max = Max(fabs(max), fabs(f[2 * i - 1]));
  }

  for(int j = 1; j <= N; j++) {
    float ti = dt * (j - 1);
	  fprintf(file, "%12g %12g\n", ti, y[2 * j - 1]);
  }

  fprintf(file, "\n\n");

  for (int j = 1; j <= N; j++) {
	  float ti = dt * (j - 1);
    fprintf(file, "%12g %12gf\n", ti, f[2 * j - 1] * 2.5 / max);
  }
	
	fclose(file);
  free_vector(f, 1, 2 * N);
	free_vector(y, 1, 2 * N);
	free_vector(g1, 1, 2 * N);
	free_vector(g2, 1, 2 * N);
        
  }
	
  return 0;
}

float Max(float a, float b)
{
	return ((a)>(b) ? (a):(b));
}

float frand() {
    return (float)rand()/RAND_MAX - 0.5;
}
