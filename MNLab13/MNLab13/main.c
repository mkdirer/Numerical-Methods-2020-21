#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Numerical_Recipes/nrutil.h"
#include "Numerical_Recipes/nrutil.c"
#include "Numerical_Recipes/gauleg.c"
#include "Numerical_Recipes/gaulag.c"
#include "Numerical_Recipes/gammln.c"
#include "Numerical_Recipes/gauher.c"


float func1(float x)
{
	return 1.0 / ( x * sqrt( x * x - 1.0) );
}
float func2_L(float x)
{
	return log(x) * exp( - (x * x) );
}
float func2_H(float x)
{
	return log(fabs(x)) / 2.0;  
}
float func3(float x)
{
	return sin( 2.0 * x ) * exp( - (3.0 * x) );
}
float func3_L(float x)
{
	return sin( 2.0 * x ) * exp( - (2.0 * x) );
}


int main()
{
	float c1_a = M_PI / 3.0;
	float c2_a = -0.8700577;
	float c3_a = 2.0 / 13.0;
	float integral;
	
	FILE *results = fopen("out.dat", "w");
	if(!results) 
		printf("Error during creating a file: 'out.dat' !");
	
	
//////////////// ad1 ////////////////

	float x1 = 1.0;
	float x2 = 2.0;
	float *x, *w;
	
	
	for(int n = 2; n <= 100; n++)
	{	
			x = vector(1, n);
			w = vector(1, n);

		///////////kwadratura Gaussa-Legandre'a

			gauleg(x1, x2, x, w, n);
			
			integral = 0.0;
			for (int i = 1; i <= n; i++)
				integral += w[i] * func1(x[i]);
			
			fprintf(results,"%d\t%f\n", n, fabs(integral - c1_a) );
			
			free_vector(x, 1, n);
			free_vector(w, 1, n);
	}
	printf("Gauss-Legandre:  c1 =  %f   blad = %f\n", integral, fabs(integral - c1_a) );
	fprintf(results,"\n\n");
	
	
//////////////// ad2 ////////////////
	
	for(int n = 2; n <= 100; n+=2)
	{

	///////////kwadratura Gaussa-Hermite'a
		
		x = vector(1, n);
		w = vector(1, n);

		gauher(x, w, n);
		
		integral = 0.0;
		
		for (int i = 1; i <= n; i++)
			integral += w[i] * func2_H(x[i]);
			
		fprintf(results,"%d\t%f\n", n, fabs(integral - c2_a) );
		
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}
	printf("Gauss-Hermite:   c2 = %f   blad = %f\n", integral, fabs(integral - c2_a) );
	fprintf(results,"\n\n");


	for(int n = 2; n <= 100; n++){
	////////////kwadratura Gaussa-Legandre'a
	
		x1 = 0.0;
		x2 = 5.0;
		x = vector(1, n);
		w = vector(1, n);


		gauleg(x1, x2, x, w, n);
		
		
		integral = 0.0;
		for (int i = 1; i <= n; i++)
			integral += w[i] * func2_L(x[i]);
			
		fprintf(results,"%d\t%f\n", n, fabs(integral - c2_a) );
		

		free_vector(x, 1, n);
		free_vector(w, 1, n);
		
	}
	printf("Gauss-Legandre:  c2 = %f   blad = %f\n", integral, fabs(integral - c2_a) );
	fprintf(results,"\n\n");
	
	
//////////////// ad 3////////////////

	for(int n = 2; n <= 20; n++)
	{
	/////////////kwadratura Gauss-Laguerrre'a
		
		// alokacja	
		float alf = 0.0;
		
		x = vector(1, n);
		w = vector(1, n);


		gaulag(x, w, n, alf);
		
		integral = 0.0;
		for (int i = 1; i <= n; i++)
			integral += w[i] * func3_L(x[i]);
			
		fprintf(results,"%d\t%f\n", n, fabs(integral - c3_a) );
		
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}
	printf("Gauss-Laguerrre: c3 =  %f   blad = %f\n", integral, fabs(integral - c3_a) );

	fclose(results);
	
	return 0;
}