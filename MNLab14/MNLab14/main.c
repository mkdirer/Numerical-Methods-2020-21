#include <stdio.h>
#include "gser.c"
#include "gcf.c"
#include "gammp.c"
#include "nrutil.h"
#include "nrutil.c" 
#include "gammln.c"
#include <math.h>

double gestosc(double o, double sr, double x);
double dystrybuanta(double o, double sr, double x);
double roz_jedno(int a, int c, unsigned long int m);

int main(){
    double srednia1 = 0;
    double srednia2 = 0;
    double val1[10000];
    double val2[10000];
    FILE *file = fopen("U.dat", "w");
    double sig = (1-0)/12.0;
    int tab1[12] = {0};
    int tab2[12] = {0};

    for(int i=0; i<5000; i++){
        double x1 = roz_jedno(123, 1, pow(2, 15));
        double x2 = roz_jedno(123, 1, pow(2, 15));
        fprintf(file, "%f %f\n", x1 , x2 );
        
        for(int j =0; j<12; j++){
            if(x1>=0+j*sig && x1<=0+(j+1)*sig) tab1[j]+=1;
            if(x2>=0+j*sig && x2<=0+(j+1)*sig) tab1[j]+=1;
        }
        val1[i] = x1;
        val1[i+5000] = x2;
        srednia1+=(x1+x2);
    }

    fprintf(file, "\n\n");

    for(int i=0; i<5000; i++){
        double x1 = roz_jedno(69069, 1, pow(2, 32));
        double x2 = roz_jedno(69069, 1, pow(2, 32));
        fprintf(file, "%f %f\n", x1, x2 );

        for(int j =0; j<12; j++){
            if(x1>=0+j*sig && x1<=0+(j+1)*sig) tab2[j]+=1;
            if(x2>=0+j*sig && x2<=0+(j+1)*sig) tab2[j]+=1;
        }
        val2[i] = x1;
        val2[i+5000] = x2;
        srednia2 += (x1+x2);
    }
    fclose(file);


    file =fopen("U_hist.dat", "w");

    for(int i=0; i<12; i++){
        double srodek = (0+i*sig + 0+(i+1)*sig)/2;
        fprintf(file, "%f %f\n", srodek , tab1[i]/10000.0 );
    }
    fprintf(file, "\n\n");
    for(int i=0; i<12; i++){
        double srodek = (0+i*sig + 0+(i+1)*sig)/2;
        fprintf(file, "%f %f\n", srodek , tab2[i]/10000.0 );
    }



    srednia1/=10000.;
    srednia2/=10000.;
    double odchylenie1 = 0;
    double odchylenie2 = 0;
    for(int i=0; i<10000; i++){
        odchylenie1+= (val1[i]-srednia1)*(val1[i]-srednia1);
        odchylenie2+= (val2[i]-srednia2)*(val2[i]-srednia2);
    }
    odchylenie1 = sqrt(odchylenie1/10000);
    odchylenie2 = sqrt(odchylenie2/10000);
    printf("Srednia = %g, odchylenie =%g dla pierwszej wersji\n", srednia1, odchylenie1);
    printf("Srednia = %g, odchylenie =%g dla drugiej wersji\n", srednia2, odchylenie2);
    fclose(file);



    double srednia = 0.2;
    double o = 0.5;
    double srednia_r = 0;
    double odchylenie_r = 0;
    sig = (-(srednia-3*o)+(srednia+3*o))/12.0;
    double N[10000];
     int t[12]={0};
    for(int i=0; i<10000; i++){
        int odrzucono = 1;
        while(odrzucono){
            double u1 = roz_jedno(69069, 1, pow(2, 32))*6*(o)+(srednia-3*o);
            double u2 = roz_jedno(69069, 1, pow(2, 32));
            if(u2<=gestosc(o, srednia, u1)){
                N[i] = u1;
                srednia_r += u1;
                odrzucono = 0;
                 for(int j =0; j<12; j++){
                    if(u1>=srednia-3*o+j*sig && u1<=srednia-3*o+(j+1)*sig) t[j]+=1;                    
                }
            }
        }
    }
    srednia_r/=10000.;
   
    for(int i=0; i<10000; i++){
        odchylenie_r+= (N[i]-srednia_r)*(N[i]-srednia_r);
    }
    odchylenie_r = sqrt(odchylenie_r/10000.);
    double war = odchylenie_r*odchylenie_r;
    printf("Wariancja = %g, srednia = %g, odchylenie =%g  dla rozkladu normalnego\n", war, srednia_r, odchylenie_r);
    double X = 0;
    file =fopen("N_hist.dat", "w");
    for(int i=0; i<12; i++){
        double jmin = srednia-3*o+i*sig ;
        double jmax = srednia-3*o+(i+1)*sig;
        double srodek = (jmin + jmax)/2;
        fprintf(file, "%f %f\n", srodek , t[i]/10000.0 );
        double pj = dystrybuanta(o, srednia, jmax) -dystrybuanta(o, srednia, jmin);
        printf("j: %d => pj: %g\n", i, pj);
        X += ((t[i]-10000*pj)*(t[i]-10000*pj))/(10000*pj);
    }
    printf("\nStatystyka testowa X^2 = %g\n", X);
    fclose(file);
    


    double v = 12 - 2- 1.0;
    double P = gammp(v/2, X/2);
    printf("Poziom ufnosci: %g\n", P);
    double alfa = 1- P;
    printf("Poziom istotnosci: %g\n ", alfa);

    double E =  16.91; 
    if(X<E){
      printf("Hipoteza H_0 nie została odrzucona na poziomie istotności alfa\n");
    }
    else 
    {
      printf("Hipoteza H_0 została odrzucona na poziomie istotności alfa\n");
    }

    return 0;
}

double roz_jedno(int a, int c, unsigned long int m){
    static long int x = 10;
    x = (a*x+c)%m;
    return x/(m+1.0);
}


double gestosc(double o, double srednia, double x){
    return (1/(o*sqrt(2*M_PI))*exp(-(x-srednia)*(x-srednia)/(2*o*o)));
}


double dystrybuanta(double o, double srednia, double x){
    return (1+erf((x-srednia)/(sqrt(2)*o)))/2.0;
}
