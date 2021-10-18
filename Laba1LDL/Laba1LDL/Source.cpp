#include <iostream>
using namespace std;
const int n = 3;
double* Gety(double[n], double**);
double* Getz(double[n], double[n][n]);
double* Getx(double[n], double**);
double** GetL(double[n][n], double[n][n]);
void LDL(double[n][n], double [n]);
int main()
{
	double A[n][n] = { 14, 2, 2,
					   2, 11, 5,
					   2, 5, 11 };
	double b[n] = { -16, -1, -19 };
	LDL(A, b);
}
double** GetL(double A[n][n], double D[n][n])
{
	double A2[n][n] = { 0, 0, 0,
					  0, 0, 0,
					  0, 0, 0 };
	for (int i = 0; i < n; i++)
	{
		A2[i][0] = A[i][0];
	}
	double** L = new double*[n];
	for (int i = 0; i < n; i++)
	{
		L[i] = new double[n];
	}
	for (int i = 0; i < n; i++)
	{
		A2[i][i] = 1;
	}
	for (int j = 0; j < n; j++)
	{
		for (int i = j + 1; i < n; i++)
		{
			double sum = 0;
			for (int k = 0; k <= j - 1; k++)
			{
				sum += A2[i][k] * L[j][k];
			}
			A2[i][j] = A[i][j] - sum;
		}
		double sum = 0;
		for (int k = 0; k <= j - 1; k++)
		{
			sum += A2[j][k] * L[j][k];
		}
		D[j][j] = A[j][j] - sum;
		for (int i = j + 1; i < n; i++)
		{
			L[i][j] = A2[i][j] / D[j][j];
		}
	}
	return L;
}

void LDL(double A[n][n], double b[n])
{
	double D[n][n] = { 0, 0, 0,
					  0, 0, 0,
					  0, 0, 0 };
	double** L = GetL(A, D);
	double* y = Gety(b, L);
	double* z = Getz(y, D);
	double* x = Getx(z, L);
	for (int i = 0; i < n; i++)
	{
		cout << x[i] << ' ';
	}
}
double* Gety(double b[n], double** L)
{
	double* y = new double[n];
	y[0] = b[0];
	for (int i = 1; i < n; i++)
	{
		double sum = 0;
		for (int k = 0; k <= i - 1; k++)
		{
			sum += L[i][k] * y[k];
		}
		y[i] = b[i] - sum;
	}
	return y;
}
double* Getz(double y[n], double D[n][n])
{
	double* z =new double[n];
	for (int i = 0; i < n; i++)
	{
		z[i] = y[i] / D[i][i];
	}
	return z;
}
double* Getx(double z[n], double** L)
{
	double* x=new double[n];
	x[n - 1] = z[n - 1];
	for (int i = n - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int k = i + 1; k < n; k++)
		{
			sum += L[k][i] * x[k];
		}
		x[i] = z[i] - sum;
	}
	return x;
}
