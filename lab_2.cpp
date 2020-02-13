#include <iostream>
#include <cmath>

using namespace std;

double E10 = 1e-10;
double E05 = 1e-5;
double E12 = 1e-12;
double c;

class Complex   
{
private:
    double re, im;     

public:
    Complex() {};
    Complex(double r) {
        re = r;
        im = 0;
    }
    Complex(double r, double i) {
        re = r;
        im = i;
    }
    Complex(const Complex& c) {
        re = c.re;
        im = c.im;
    }
    ~Complex() {}
    double abs() {
        return sqrt(re * re + im * im);
    }
    Complex& operator = (Complex& c) {
        re = c.re;
        im = c.im;
        return (*this);
    }
    Complex& operator += (Complex& c) {
        re += c.re;
        im += c.im;
        return *this;
    }
    Complex operator + (const Complex& c) {
        return Complex(re + c.re, im + c.im);
    }
    Complex operator - (const Complex& c) {
        return Complex(re - c.re, im - c.im);
    }
    Complex operator * (const Complex& c) {
        return Complex(re * c.re - im * c.im, re * c.im + im * c.re);
    }
    Complex operator / (const Complex& c)
    {
        Complex temp;

        double r = c.re * c.re + c.im * c.im;
        temp.re = (re * c.re + im * c.im) / r;
        temp.im = (im * c.re - re * c.im) / r;

        return temp;
    }
    friend ostream& operator<< (ostream&, const Complex&);
    friend istream& operator>> (istream&, Complex&);

};

ostream& operator<< (ostream& out, const Complex& c)
{
    out << "(" << c.re << ", " << c.im << ")";
    return out;
}

istream& operator>> (istream& in, Complex& c)
{
    in >> c.re >> c.im;
    return in;
}

double f24(double x) { 
    return log10(x) * log10(x) + 0.75 * log10(x) - 0.25;
}

double g24(double x) {
    return log10(x) * log10(x) + 2 * log10(x) + 1;
}

double f04(double x) {
    return log(2 * x - x * x) + 2 - sqrt(x);
}

double df04(double x) {
    return 2 * (x - 1) / (x * x - 2 * x) - 1 / (2 * sqrt(x));
}

Complex f03(Complex z) {
    Complex tzs(2.7, 0);
    Complex one(1, 0);
    return z*z*z*z - tzs * z*z*z + z - one;
}

Complex df03(Complex z){
    Complex four(4, 0);
    Complex ya_debil(2.025, 0);
    return four * (z - ya_debil) * z*z;
}

double support(double x) {
    return x - f04(x) / df04(x);
}

void bisection_method_f(double a, double b) {
    while ((b - a) >= E10) {
        c = (a + b) / 2.0;
        if (f24(c) == 0.0) {
            break;
        }
        else if (f24(c) * f24(a) < 0) { //проверка на смену знака
            b = c; //двигаемся влево
        }
        else {
            a = c; //двигаемся вправо
        }
    }
    cout << c << endl;
}

void bisection_method_g(double a, double b) {
    while ((b - a) >= E10) {
        c = (a + b) / 2.0;
        if (g24(c) == 0.0) {
            break;
        }
        else if (g24(c) * g24(a) > 0) { //проверка на смену знака
            b = c; //двигаемся влево
        }
        else {
            a = c; //двигаемся вправо
        }
    }
    cout << c << endl;
}

void newton_method(double x, double E) {
    int iteration = 0; double past_x = 0;
    while (fabs(f04(x)) >= E) {
        x -= f04(x) / df04(x); //по формуле
        iteration++;
    }
    cout << x << endl;
    cout << iteration << endl;
}

void ez_newton_method(double x, double E) {
    int iteration = 0; double x0 = x;
    while (fabs(f04(x)) >= E) {
        x -= f04(x) / df04(x0); //по формуле
        iteration++;
    }
    cout << x << endl;
    cout << iteration << endl;
}

void secant_method(double x0, double x1, double E) {
    int iteration = 0; double x = x1;
    while (fabs(f04(x)) >= E) {
        x -= f04(x) * (x - x0) / (f04(x) - f04(x0)); //по формуле
        iteration++;
    }
    cout << x << endl;
    cout << iteration << endl;
}

void newton_method_complex(Complex z, double E) {
    int iteration = 0;
   
    while ((f03(z)).abs() > E) {
        auto k = f03(z) / df03(z);
        auto b =  z - k; //по формуле
        z = b;
        iteration++;
    }
    cout << z << endl;
    cout << iteration << endl;
}


int main()
{
    double a = 0.01, b = 3; bool mod = 1;

    Complex z(-0.1, 0);
    Complex z1(2.4, 0);
    Complex z2(0.4, -0.3);
    Complex z3(0.3, 0.2);

    bisection_method_f(a, b);
    bisection_method_f(a + 1, b);
    bisection_method_g(a, b);

    newton_method(1.3, E12);    
    newton_method(1.3, E05);
    newton_method(0.091, E12);
    newton_method(0.091, E05);

    ez_newton_method(1.6, E12);
    ez_newton_method(1.6, E05);
    ez_newton_method(0.095, E12);
    ez_newton_method(0.095, E05);

    secant_method(1.6, 1.8, E12);
    secant_method(1.6, 1.8, E05);
    secant_method(0.01, 0.1, E12);
    secant_method(0.01, 0.1, E05);

    newton_method_complex(z, E05);
    newton_method_complex(z1, E05);
    newton_method_complex(z2, E05);
    newton_method_complex(z3, E05);

}
