#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


double d_rand(const double min, const double max);
double funkcja(double x, double y);

int main(){
  int N=200;
  double xmin = -10.0;
  double xmax = 10.0;
  double ymin = -10.0;
  double ymax = 10.0;
  double deltax =0.0;
  double deltay=0.0;
  double T=0.0;
  double x0 = 5.0;
  double y0 = 5.0;
  double x[N];
  double y[N]; 
  srand(time(0));

  FILE* w0 = fopen("w0.dat", "w");
  FILE* fileT = fopen("T.dat", "w");

  for(int i=0; i<N;i++)
    x[i] = y[i]=5.0; 

  for(int it=0; it<=20; it++){
    T = 10./pow(2,it);

    for(int k =0; k<100; k++){
      for(int i=0; i<N;i++){
        deltax = d_rand(-1.0, 1.0);
        deltay = d_rand(-1.0, 1.0);

        while( fabs(x[i]+deltax) > xmax || fabs(x[i]+deltax) < xmin )
          deltax = d_rand(-1.0, 1.0);
        while( fabs(y[i]+deltay) > ymax || fabs(y[i]+deltay) < ymin )
          deltay = d_rand(-1.0, 1.0);

        if( funkcja(x[i]+deltax, y[i]+deltay) < funkcja(x[i], y[i]) ){
          x[i] += deltax;
          y[i] += deltay;
         }
         else if( d_rand(0.0, 1.0) < exp( (-1)* (funkcja(x[i]+deltax, y[i]+deltay) - funkcja(x[i],y[i]))/T  ) ){
           x[i] += deltax;
           y[i] += deltay;
         }
      }
      fprintf(w0, "%f\n", funkcja(x[0], y[0]));
    }

    if(it==0 || it==7 || it == 20){
      for(int i=0; i<N; i++)
        fprintf(fileT, "%f\t%f\n", x[i], y[i]);
      fprintf(fileT, "\n\n");
    }

  }

  fclose(w0); 
  fclose(fileT);

  double fun_x=0.0;
  double fun_y = 0.0;
  double fun_min = 0.0;

  fun_x = x[0]; fun_y = y[0];
  fun_min = funkcja(x[0],y[0]);

  for(int i=1; i<N;i++){
    if( funkcja(x[i],y[i]) < fun_min ){
      fun_x = x[i];
      fun_y = y[i];
      fun_min = funkcja(x[i],y[i]);
    }
  }
  printf("wartosc minimalna: %f , w punkcie (%g,%g)\n", fun_min, fun_x,fun_y);

  return 0;
}

/* *************************************** */

double d_rand(const double min, const double max){
  double r = (double)rand()/RAND_MAX; 
  r = r * (max-min)+min;
  return r;
}

/* *************************************** */

double funkcja(double x, double y){
  return sin(x)*sin(y) - exp( (-1)*pow((x+M_PI/2),2) - pow((y-M_PI/2),2) );
}