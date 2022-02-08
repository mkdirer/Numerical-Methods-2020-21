#include <iostream>
#include <math.h>

const int n = 8;
const int a = 0;
const int b = 1;


void print (double **array)
{
	for (int i=0;i<=n;i++)
	{
		for (int j=0;j<i+1;j++)
			std::cout<<array[i][j]<<" ";

		std::cout<<std::endl;
	}
}

void print_to_file(double **array, const char *filename)
{
    FILE *file = fopen(filename, "wt");
    for (int i=0;i<=n;i++)
	  {
		  for (int j=0;j<i+1;j++)
        fprintf(file, "%12.10g  ", array[i][j]);

    fprintf(file, "\n");
	  }
    fclose(file);
}

void richardson (double **array)
{
	for (int i=1;i<=n;i++)
	{
		for (int j=1;j<i+1;j++)
		{
			array[i][j] = (pow(4,j)*array[i][j-1]-array[i-1][j-1])/(pow(4,j)-1);
		}
	}
}

double func (double arg)
{
	return (log(pow(arg,3)+3*pow(arg,2) + arg + 0.1)) * (sin(18*arg));
}

double index(int iteration,double hw)
{
	return (a + iteration*hw);
}

int main ()
{

	double **array;
	double sum;
	double hw;
	int N;


	array = new double * [n+1];

	for (int i=0;i<=n;i++)
		array[i] = new double [i+1];


	//********** Simpson *************

	std::cout<<std::endl<<std::endl<<"Simpson:"<<std::endl<<std::endl;

	for (int i=0;i<=n;i++)
	{
		sum = 0.0;	
		N = pow(2,i+1);
		hw = (b-a)/static_cast<double>(N);

		for (int j=0;j<=static_cast<int>(N/2-1.0);j++)
		{
			sum += (hw/3.0)*(func(a+2*j*hw)+4*func(a+ (2*j+1)*hw)+func(a+ (2*j+2)*hw));
		}

		array[i][0] = sum;
	}

	richardson(array);
	print(array);
  print_to_file(array, "Simpson.txt");


	//*************Milne***************
	std::cout<<std::endl<<std::endl<<"Milne:"<<std::endl<<std::endl;

	for (int i=0;i<=n;i++)
	{
		sum = 0.0;	
		N = pow(2,i+2);
		hw = (b-a)/static_cast<double>(N);

		for (int j=0;j<=static_cast<int>(N/4-1.0);j++)
		{
			sum += ((4*hw)/90.0)*(7*func(a+4*j*hw)+32*func(a+(4*j+1)*hw)+12*func(a+(4*j+2)*hw)+32*func(a+(4*j+3)*hw)+7*func(a+(4*j+4)*hw));
		}

		array[i][0] = sum;
	}

	richardson (array);
	print(array);
  print_to_file(array, "Milne.txt");


//****************delete*********************

		for (int i=0;i<=n;i++)
			delete [] array[i];
		
		delete [] array;

}


