#include <iostream>
#include <cmath>
#include <fstream>

void printtab(double *tab, const int N)
{
  for(int i=0; i<=N; i++)
  {
    std::cout<<tab[i]<<std::endl;
  }
}

void fprinttab(double *tab, int N, double h, const char* str)
{
  std::fstream file;
  file.open(str, std::ios::out);
  for(int i=0; i<=N; i++)
  {
    //printf("%.2f %g\n", i*h, x_n[i])
    file<<h*i<<"\t"<<tab[i]<<std::endl;
  }
}

int main()
{

	const int N = 2000;
	int it = 0;
	double beta =0.4; //0.0; //0.4; //=0.4;
  double F0 = 0.0;//0.0; //0.0;//0.1;
  double omega = 0.8;//0.8; //=0.8;// =0.8;
	
	

	const char * file = "out2.dat"; //"out1.dat"//"out2.dat" //"out3.dat"
	double h = 0.02, xs[N+1], xn[N+1], b[N+1];
	double d0[N+1], d1[N+1], d2[N+1];
	double sumxn = 1e-5, sumxs = 0;


	double a1 = 1;
	double a2 = (h*h) - 2 - beta * h;
	double a3 = 1 + beta * h;
	d1[0]=d2[0]= 0;
	d0[0]=b[0] = 1;
	d0[1] = 1;
	d1[1] = -1;
	d2[1] = 0;
	b[1] = 0;
	for (int i = 2; i <=N; i++)
	{
		d0[i] = a3;
		d1[i] = a2;
		d2[i] = a1;
		b[i] = h*h*F0*sin(omega*h*i);
	}

  for(int i =0; i<=N; i++)
  {
    xs[i] =0;
    xn[i]=0;
  }

  for(;true;it++)
  {
    if(fabs(sumxn - sumxs)<1e-6)
      break;
    xn[0] = b[0]/d0[0];
    xn[1] = (b[1] - d1[1]*xs[0]) / d0[1];
      for(int i=2; i<=N; i++)
      {
        xn[i] = (b[i]-d1[i]*xs[i-1]-d2[i]*xs[i-2])/d0[i];
      }
      sumxs = sumxn = 0;
      for(int i=0; i<=N; i++)
      {
        sumxs += xs[i] *xs[i];
        sumxn += xn[i] *xn[i];

      }
      for(int i=0; i<=N; i++)
      {
        xs[i]=xn[i];
      }
  }

  printtab(xn,N);
  fprinttab(xn, N, h, file);
  //std::cout<<std::endl<<it;

}