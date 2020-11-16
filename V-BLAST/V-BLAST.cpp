#include<iostream>
#include<cmath>
#include <Eigen/Eigen>
#include<Eigen/Core>
#include<Eigen/Dense>
using namespace std;
using namespace Eigen;
#define nR 3
#define nT 3
int main()
{
	int i, j, k;
	MatrixXcd H(nR, nT);
	typedef Matrix<complex<double>, nR, 1> VectornRcd;
	VectornRcd x;
	MatrixXcd G;
	double lenth = 0;
	int min = 0;
	int LastMin = 1;
	VectornRcd w;
	double cReal[10] = { 0 };
	double cImag[10] = { 0 };
	complex<double> y;
	VectornRcd t;
	MatrixXcd t1, t2;
	double temp;
	int numRows, numCols;
	cout << "Please input the real part of H:" << endl;
	for (i = 0; i < nR; i++)
	{
		for (j = 0; j < nT; j++)
		{
			cin >> temp;
			H.real()(i, j) = temp;
		}
	}
	cout << "Please input the imaginary part of H:" << endl;
	for (i = 0; i < nR; i++)
	{
		for (j = 0; j < nT; j++)
		{
			cin >> temp;
			H.imag()(i, j) = temp;
		}
	}
	cout << "H = " << endl << H << endl;
	cout << "Please input the real part of x:" << endl;
	for (i = 0; i < nR; i++)
	{
		cin >> temp;
		x.real()(i) = temp;
	}
	cout << "Please input the imaginary part of x:" << endl;
	for (i = 0; i < nR; i++)
	{
		cin >> temp;
		x.imag()(i) = temp;
	}
	for (i = 1; i <= nT; i++)
	{
		G = H.adjoint();
		G *= H;
		G = G.inverse();
		G *= H.adjoint();//compute pseudo-inverse
		lenth = G.row(0).norm();
		min = 0;
		for (j = 0; j < nT - i + 1; j++)
		{
			if (G.row(j).norm() < lenth)
			{
				lenth = G.row(j).norm();
				min = j;
			}
		}//compute k_i
		w = G.row(min);//w^T=g_i^(k)
		y = w.dot(x);//y_k_i=w_k_i^T*x_i
		if (min >= LastMin)
		{
			if (abs(y.real()) > abs(y.imag()))
			{
				cReal[min + i - 1] = 1;
			}
			else
			{
				cImag[min + i - 1] = 1;
			}
		}
		else
		{
			if (abs(y.real()) > abs(y.imag()))
			{
				cReal[min] = 1;
			}
			else
			{
				cImag[min] = 1;
			}
		}
		LastMin = min;//c=Q[y_k_i]
		for (j = 0; j < nR; j++)
		{
			t.real()(j) = H.real()(j, min) * cReal[min] - H.imag()(j, min) * cImag[min];
			t.imag()(j) = H.real()(j, min) * cImag[min] + H.imag()(j, min) * cReal[min];
		}//t = h_k_i*c_k_i
		numRows = H.rows();
		numCols = H.cols() - 1;
		if (min < numCols && i != nT)
		{
			H.block(0, min, numRows, numCols - min) = H.block(0, min + 1, numRows, numCols - min);
			H.conservativeResize(numRows, numCols);
		}
		x -= t;
		}//(7) and (8)
	cout << "c = " << endl;
	for (j = 0; j < nT; j++)
	{
		cout << cReal[j] << " + " << cImag[j] << "i" << endl;
	}
}