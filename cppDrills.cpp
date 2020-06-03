
// quantumdrills.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Quantum Computing For Computer Science Drills Benjamin Schreyer
// each drill is labeled
#include <iostream>
#include <string>
#include <math.h>
//Drill 1.1.1
struct complex
{
	//real
	double a;
	//imaginary
	double b;
};
//Drill 1.3.1
struct polarComplex
{
	double m;//magnitude
	double t;//angle
};
//Drill 1.1.1
struct complex complexMultiply(struct complex x, struct complex y)
{
	//definition of complex multiplication
	struct complex o;
	o.a = x.a * y.a - x.b * y.b;
	o.b = x.a * y.b + y.a * x.b;
	return o;
}
//Drill 1.1.1
struct complex complexAdd(struct complex x, struct complex y)
{
	//definition of complex addition
	struct complex o;
	o.a = x.a + y.a;
	o.b = x.b + y.b;
	return o;
}
//Drill 1.2.1
struct complex complexSubtract(struct complex x, struct complex y)
{
	//use addition but multiply the second term by (-1,0)
	struct complex n;
	n.a = -1.0;
	n.b = 0.0;
	return complexAdd(x, complexMultiply(n, y));
}
//Drill 1.2.1
struct complex complexDivide(struct complex n, struct complex d)
{
	//definition of complex division
	struct complex o;
	o.a = (n.a * d.a + n.b * d.b) / (d.a * d.a + d.b * d.b);
	o.b = (d.a * d.b - n.a * d.b) / (d.a * d.a + d.b * d.b);
	return o;
};

//Dril 1.2.1
double complexModulus(struct complex x)
{
	//length of a 2d vector same as complex modulus
	return pow(x.a * x.a + x.b * x.b, 0.5);
}

//Drill 1.3.1
struct complex polarToCartesian(struct polarComplex x)
{
	struct complex n;
	n.a = x.m;
	n.b = 0;
	struct complex a;
	a.a = cos(x.t);
	a.b = sin(x.t);
	return complexMultiply(n, a);
}
//Drill 1.3.1
struct polarComplex cartesianToPolar(struct complex x)
{
	polarComplex o;
	o.m = complexModulus(x);
	o.t = atan(x.b / x.a);
	return o;
}

//Drill 1.3.1
std::string complexToPolarString(struct complex x)
{
	return std::to_string(atan(x.b / x.a)) + "ang " + std::to_string(complexModulus(x));
}
//Drill 1.1.1
std::string complexToString(struct complex x)
{
	return std::to_string(x.a) + " + " + std::to_string(x.b) + "i";
}

int main()
{
	//Drill 1.1.1
	std::cout << "Give two complex numbers, first the real magnitude then the complex magniute(4 total floats):\n";
	struct complex w;
	struct complex z;
	std::cin >> w.a;
	std::cin >> w.b;
	std::cin >> z.a;
	std::cin >> z.b;
	std::cout << "\nYour given numbers: " + complexToString(w) + " | " + complexToString(z) + "\n\n";
	std::cout << "Drill 1.1.1\n";
	std::cout << "Multiplied: " + complexToString(complexMultiply(w, z)) + " Added: " + complexToString(complexAdd(w, z) ) + "\n\n";
	//Drill 1.2.1
	std::cout << "Drill 1.2.1\n";
	std::cout << "Subracted: " + complexToString(complexSubtract(w, z)) + " Divided: " + complexToString(complexDivide(w, z)) + "\n";
	std::cout << "Length/Modulus of your two complex numbers: " + std::to_string(complexModulus(w)) + " | " + std::to_string(complexModulus(z)) + "\n\n";
	//Drill 1.3.1
	std::cout << "Drill 1.3.1\n";
	std::cout << "Polar of your complex numbers: " + complexToPolarString(w) + " | " + complexToPolarString(z) + " Testing of polar to complex(will just repeat the numbers you gave):" + complexToString(polarToCartesian(cartesianToPolar(w))) + " | " + complexToString(polarToCartesian(cartesianToPolar(z))) + "\n";
	//prevent window close
	std::cin >> z.b;



}

