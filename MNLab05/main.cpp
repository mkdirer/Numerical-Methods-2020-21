#include <iostream>
#include <cmath>
#define N 7
#define IT_MAX 12
#define max(X,Y) ((X)>(Y)? (X):(Y))
#define min(X,Y) ((X)<(Y)? (X):(Y))

double skalar(double *a, double *b, int n)
{
    double x=0.;

    for(int i=0;i<n;i++)
        x+= a[i] * b[i];
    return x;
}

void wektor_liczba(double *v1, double *v2, double val, int n)
{
    for(int i=0;i<n;i++) v2[i] = v1[i]*val;
}

void mnozMacierz(double**a, double **b, double **c, int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            c[i][j]=0.;
            for(int k=0;k<n;k++)
                c[i][j] += a[i][k]*b[k][j];
        }
    }
}


void macierz_wektor(double**A, double*x, double*y, double m, double n)
{   int jmin, jmax;
    for(int i=0;i<n;i++){
        jmin=max(0,i-m);
        jmax=min(i+m,n-1);
        y[i]=0;
        for(int j=jmin;j<=jmax;j++)
            y[i]+=A[i][j]*x[j];
    }
}


void metodaHott(double **W, double *xk, double lambdak)
{
    for(int i=0;i<N; i++)
        for(int j=0;j<N;j++)
            W[i][j] = W[i][j] - lambdak * xk[i] *xk[j];
}

void transM(double **a)
{
    double **pom = new double*[N];
    int i, j;
    for(i=0;i<N;i++)
    {
        pom[i] = new double[N];
        for (j = 0; j < N; j++) {
            pom[i][j] = a[j][i];
        }
    }
    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
            a[i][j]=pom[i][j];
        }
}

void printM(double **m)
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
            printf("%g ", m[i][j]);
        std :: cout << std :: endl;
    }

}

void printW(double *v)
{
    for(int i=0;i<N;i++)
        std :: cout << v[i] << " ";
    std :: cout << std :: endl;
}

int main() {

    FILE *plik1 = fopen("lambda.dat", "w+");
    FILE *plik2;
    FILE *plik3;

    double **A = new double*[N];
    double **X = new double*[N];
    double **W = new double*[N];
    double *L = new double[N];

    for(int i=0;i<N;i++)
    {
        A[i] = new double[N];
        W[i] = new double[N];
        X[i] = new double[N];
    }

    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
        {
            A[i][j] = 1./(sqrt(2. + std::fabs(i-j)));
            W[i][j] = 1./(sqrt(2. + std::fabs(i-j)));
        }


    int Kval = N;
    double *xp1 = new double [N];
    double *xp2 = new double [N];
    double lambda=0.;

    for(int k=0;k<Kval;k++)
    {
        for(int i=0;i<N;i++)
            xp1[i]=1;

        for(int i=1;i<=IT_MAX;i++)
        {
            macierz_wektor(W, xp1, xp2, N, N);
            lambda = skalar(xp2, xp1, N)/skalar(xp1, xp1, N);
            fprintf(plik1, "%d %g \n", i, lambda);
            wektor_liczba(xp2, xp1, 1./sqrt(skalar(xp2, xp2, N)), N);
        }
        fprintf(plik1, "\n \n");
        //wartosci wlasne wpisuje wierszami, by pozniej dokonac transpozycji macierzy
        for(int i=0;i<N;i++)
            X[k][i] = xp1[i];

        L[k] = lambda;
        metodaHott(W, xp1, lambda);

    }

    transM(X);
    printM(X);
    plik2 = fopen("MacierzWynik.dat", "w");

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            fprintf(plik2, "%15g", X[i][j]);
        fprintf(plik2, "\n");
    }

    std :: cout << std::endl;
    printW(L);

    double **D = new double*[N];
    double **pom = new double*[N];
    for(int i=0;i<N;i++)
    {
        D[i] = new double[N];
        pom[i] = new double[N];
    }
    transM(X);

    mnozMacierz(X, A, pom, N);
    transM(X);
    mnozMacierz(pom, X, D, N);

    std :: cout << std::endl;
    printM(D);
    plik3 = fopen("MacierzD.dat", "w");
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(plik3, "%15g", D[i][j]);
        }
        fprintf(plik3, "\n");
    }

    for(int i=0;i<N;i++)
    {
        delete [] D[i];
        delete [] A[i];
        delete [] pom[i];
        delete [] X[i];
        delete [] W[i];
    }

    delete [] D;
    delete [] A;
    delete [] pom;
    delete [] X;
    delete [] W;
    delete [] xp1;
    delete [] xp2;
    delete [] L;

    fclose(plik1);

}