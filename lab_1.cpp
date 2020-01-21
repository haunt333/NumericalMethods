#include "pch.h"
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double x) {
	return exp(-x*x);
}

double rectangle(int n, double a, double b) { //метод прямоугольников
	double result = 0;
	double h = 0.001;
	//double h = 0.01;

	for (int i = 0; i < n; i++) {
		double x1 = a + i * h;
		result += h * f(x1);
	}

	return result;
}

double trapeze(int n, double a, double b) { //метод трапеции
	double result = 0;
	double h = 0.001;
	//double h = 0.01;

	for (int i = 1; i < n; i++) {
		double x1 = a + i * h;
		double x2 = a + (i - 1) * h;
		result += h * (f(x1) + f(x2)) / 2; 
	}
	
	return result;
}

double simpson(int n, double a, double b) { //метод симпсона
	double result = 0;
	double h = 0.001;
	//double h = 0.01;
	
	for (int i = 1; i < n; i++) {
		double x1 = a + i * h;
		double x2 = a + (i - 1) * h;
		result += h * (f(x1) + 4 * f(x1 / 2 + x2 / 2) + f(x2)) / 6;
	}

	return result;
}

double montecarlo(int n, double a, double b) { //метод монте-карыча
	double result = 0;

	for (int i = 0; i < n; i++) {
		double x = (double)(rand()) / RAND_MAX;
		double y = (double)(rand()) / RAND_MAX;
		if (f(x) > y)
			result += 1;
	}

	return result = result / n * (b - a);
}

void main()
{
	double a = 0, b = 1; // пределы интегрирования
	int n = 1000; //число промежутков разбиения

	cout << "integral = " << setprecision(10) << trapeze(n, a, b) << " by trapeze method\n";
	cout << "integral = " << setprecision(10) << simpson(n, a, b) << " by simpson method\n";
	cout << "integral = " << setprecision(10) << rectangle(n, a, b) << "by rectangle method\n";
	cout << "integral = " << setprecision(10) << montecarlo(n, a, b) << " by monte-carlo method";
}

