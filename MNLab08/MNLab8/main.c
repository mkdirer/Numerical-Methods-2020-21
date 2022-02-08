#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nrutil.h"
//#include "nrutil.c"
//#include "gaussj.c"

void wyzM(float *xw, float *yw, float *m, int n, float alfa, float beta)
{
    int xmax = 5;
    int xmin = -5;
    double h = (xmax - xmin) / (n - 1.);
    int i, j;
    double lambda = h / (h + h);
    double mi = 1 - lambda;

    float **A, **d;
    A = matrix(1, n, 1, n);
    d = matrix(1, n, 1, 1);

    // Macierz A
    for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++)
        {
            if (i == j)
                A[i][j] = 2.;
            else if (i == (j + 1))
                A[i][j] = mi;
            else if (i == (j - 1))
                A[i][j] = lambda;
            else
                A[i][j] = 0.;
        }
    A[1][1] = A[n][n] = 1.;
    A[1][2] = A[n][n - 1] = 0.;

    // Macierz d
    
    for (i = 1; i <=n; i++)
    {
      

        d[i][1] = (6. / (h + h)) * (((yw[i + 1] - yw[i]) / h) - ((yw[i] - yw[i - 1]) / h));
    }

    d[1][1] = alfa;
    d[n][1] = beta;

    // gaussj Am = d
    gaussj(A, n, d, 1);

    //Przypisanie wartosci do macierzy m
    for (i = 1; i <= n; i++)
    {
        m[i] = d[i][1];
    }

    free_matrix(A, 1, n, 1, n);
    free_matrix(d, 1, n, 1, 1);
}

float wyzSx(float *xw, float *yw, float *m, int n, float x)
{
    int i;
    double sx;
    int xmax = 5;
    int xmin = -5;
    double h = (xmax - xmin) / (n - 1.);

    float *A, *B;
    A = vector(1, n);
    B = vector(1, n);

    for (i = 2; i <= n; i++)
    {
        A[i] = (yw[i] - yw[i - 1]) / h - (h / 6) * (m[i] - m[i - 1]);
        B[i] = yw[i - 1] - m[i - 1] * (pow(h, 2) / 6);
    }

    for (i = 2; i <= n; i++)
    {
        if (xw[i - 1] <= x && x <= xw[i])
        {
            sx = m[i - 1] * pow(xw[i] - x, 3) / (6 * h) + m[i] * pow(x - xw[i - 1], 3) / (6 * h) + A[i] * (x - xw[i - 1]) + B[i];
            break;
        }
    }

    free_vector(A, 1, n);
    free_vector(B, 1, n);
    return sx;
}

float fun(float x)
{
    return (1 / (1 + x * x));
}

float fun2(float x)
{
    return cos(2 * x);
}

float pochodna(float x)
{
    float x_step = 0.01;
    float wyn;
    wyn = (fun(x - x_step) - 2 * fun(x) + fun(x + x_step)) / pow(x_step, 2);
    return wyn;
}

int main()
{
    int xmax = 5;
    int xmin = -5;
    double x_step = 0.01;
    int n;
    float alfa, beta;
    double h;

    float *wx;
    float *wy;
    float *m;


    FILE *pochodne, *f1, *f2;
	pochodne  = fopen("pochodne.dat", "w+");
	f1  = fopen("file1.dat", "w+");
	f2  = fopen("file2.dat", "w+");

    /////////////////////zadanie 4///////////////////

    alfa = beta = 0.;

    n = 10;
    wx = vector(1, n);
    wy = vector(1, n);
    m = vector(1, n);
    h = (xmax - xmin) / (n - 1.);

    // funkcja 1.

    for (int l = 1; l <= n; l++)
    {
        wx[l] = xmin + (l - 1) * h;
        wy[l] = pochodna(wx[l]);
    }

    wyzM(wx, wy, m, n, alfa, beta);

    for (double i = xmin, ko=1; i <= xmax; i += h, ko++)
        fprintf(pochodne, "%20g%20g%20g\n", i, m[abs((int)ko)], wyzSx(wx, wy, m, n, i)); //m[i] i f'' numerycznie

    free_vector(wx, 1, n);
    free_vector(wy, 1, n);
    free_vector(m, 1, n);

    ////////////zadanie 5///////////////

    /////////////// n =5 ////////////////////

    n = 5;
    wx = vector(1, n);
    wy = vector(1, n);
    m = vector(1, n);
    h = (xmax - xmin) / (n - 1.);

    /////////////////////// funkcja 1.
    

    for (int l = 1; l <= n; l++)
    {
        wx[l] = xmin + (l - 1) * h;
        wy[l] = fun(wx[l]);
    }
    wyzM(wx, wy, m, n, alfa, beta);

    //uzupenienie m
    for (double x = xmin; x <= xmax; x += x_step)
        fprintf(f1, "%g\t%10g\t%10g\n", x, fun(x), wyzSx(wx, wy, m, n, x));
    fprintf(f1,"\n\n");

    ////////////////////// funkcja 2.

    for (int l = 1; l <= n; l++)
    {
        wx[l] = xmin + (l - 1) * h;
        wy[l] = fun2(wx[l]);
    }
    wyzM(wx, wy, m, n, alfa, beta);

    //uzupenienie m
    for (double x = xmin; x <= xmax; x += x_step)
        fprintf(f2, "%g\t%10g\t%10g\n", x, fun2(x), wyzSx(wx, wy, m, n, x));
    fprintf(f2,"\n\n");

    free_vector(wx, 1, n);
    free_vector(wy, 1, n);
    free_vector(m, 1, n);

    /////////////// n = 8 /////////////////////

    n = 8;
    wx = vector(1, n);
    wy = vector(1, n);
    m = vector(1, n);
    h = (xmax - xmin) / (n - 1.);

    ///////////////////////////////// funkcja 1.

    for (int l = 1; l <= n; l++)
    {
        wx[l] = xmin + (l - 1) * h;
        wy[l] = fun(wx[l]);
    }
    wyzM(wx, wy, m, n, alfa, beta);

    //uzupenienie m
    for (double x = xmin; x <= xmax; x += x_step)
        fprintf(f1, "%g\t%10g\t%10g\n", x, fun(x), wyzSx(wx, wy, m, n, x));
    fprintf(f1,"\n\n");

    ///////////////////////////////// funkcja 2.
    

    for (int l = 1; l <= n; l++)
    {
        wx[l] = xmin + (l - 1) * h;
        wy[l] = fun2(wx[l]);
    }
    wyzM(wx, wy, m, n, alfa, beta);

    //uzupenienie m
    for (double x = xmin; x <= xmax; x += x_step)
        fprintf(f2, "%g\t%10g\t%10g\n", x, fun2(x), wyzSx(wx, wy, m, n, x));
    fprintf(f2,"\n\n");

    free_vector(wx, 1, n);
    free_vector(wy, 1, n);
    free_vector(m, 1, n);

    /////////////// n = 21 ////////////////////

    n = 21;
    wx = vector(1, n);
    wy = vector(1, n);
    m = vector(1, n);
    h = (xmax - xmin) / (n - 1.);

    ////////////////////////////////// funkcja 1.

    for (int l = 1; l <= n; l++)
    {
        wx[l] = xmin + (l - 1) * h;
        wy[l] = fun(wx[l]);
    }
    wyzM(wx, wy, m, n, alfa, beta);

    //uzupenienie m
    for (double x = xmin; x <= xmax; x += x_step)
        fprintf(f1, "%g\t%10g\t%10g\n", x, fun(x), wyzSx(wx, wy, m, n, x));
    fprintf(f1,"\n\n");
    //////////////////////////////////// funkcja 2.

    for (int l = 1; l <= n; l++)
    {
        wx[l] = xmin + (l - 1) * h;
        wy[l] = fun2(wx[l]);
    }
    wyzM(wx, wy, m, n, alfa, beta);

    //uzupenienie m
    for (double x = xmin; x <= xmax; x += x_step)
        fprintf(f2, "%g\t%10g\t%10g\n", x, fun2(x), wyzSx(wx, wy, m, n, x));
    fprintf(f2,"\n\n");

    free_vector(wx, 1, n);
    free_vector(wy, 1, n);
    free_vector(m, 1, n);

    fclose(f1);
    fclose(f2);
    fclose(pochodne);

    return 0;
}

