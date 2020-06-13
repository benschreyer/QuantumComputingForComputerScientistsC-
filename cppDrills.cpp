// quantumdrills.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Quantum Computing For Computer Science Drills Benjamin Schreyer Summer 2020

#include <iostream>
#include <string>
#include <math.h>
#include <vector>

//Drill 1.1.1
struct  complex
{
	//real
	double a;
	//imaginary
	double b;
};


typedef   complex complex;
typedef std::vector<std::vector<  complex>>  complexMatrix;
//Drill 1.3.1
struct  polarComplex
{
	double m;//magnitude
	double t;//angle
};
typedef struct polarComplex polarComplex;
//Drill 1.1.1
  complex complexMultiply(  complex x,   complex y)
{
	//definition of complex multiplication
	  complex o;
	o.a = x.a * y.a - x.b * y.b;
	o.b = x.a * y.b + y.a * x.b;
	return o;
}
//Drill 1.1.1
  complex complexAdd(  complex x,   complex y)
{
	//definition of complex addition
	  complex o;
	o.a = x.a + y.a;
	o.b = x.b + y.b;
	return o;
}
//Drill useful
  
//Drill 1.2.1
  complex complexSubtract(  complex x,   complex y)
{
	//use addition but multiply the second term by (-1,0)
	  complex n;
	n.a = -1.0;
	n.b = 0.0;
	return complexAdd(x, complexMultiply(n, y));
}
//Drill 1.2.1
  complex complexDivide(  complex n,   complex d)
{
	//definition of complex division
	  complex o;
	o.a = (n.a * d.a + n.b * d.b) / (d.a * d.a + d.b * d.b);
	o.b = (d.a * d.b - n.a * d.b) / (d.a * d.a + d.b * d.b);
	return o;
};

//Dril 1.2.1
double complexModulus(  complex x)
{
	//length of a 2d vector same as complex modulus
	return pow(x.a * x.a + x.b * x.b, 0.5);
}

//Drill 1.3.1
  complex polarToCartesian(  polarComplex x)
{
	  //use trig
	  complex n;
	n.a = x.m;
	n.b = 0;
	  complex a;
	a.a = cos(x.t);
	a.b = sin(x.t);
	return complexMultiply(n, a);
}
//Drill 1.3.1
  polarComplex cartesianToPolar(  complex x)
{
	  //use inverse trig
	polarComplex o;
	o.m = complexModulus(x);
	o.t = atan(x.b / x.a);
	return o;
}

//Drill 1.3.1
std::string complexToPolarString(  complex x)
{
	return std::to_string(atan(x.b / x.a)) + "ang " + std::to_string(complexModulus(x));
}
//Drill 1.1.1
std::string complexToString(  complex x)
{
	return std::to_string(x.a) + " + " + std::to_string(x.b) + "i";
}
//Drill 2.1.1/2.2.1
complexMatrix & complexMatrixScalarMultiplication(float r, complexMatrix & m)
{
	//each element is assignedto itself times r
	for (int i = 0; i < m.size(); i++)
	{
		for (int j = 0; j < m[0].size(); j++)
		{
			m[i][j] = complexMultiply({r,0.0}, m[i][j]);
		}
	}
	return m;
}
//Drill 2.1.1/2.2.1
complexMatrix& printComplexMatrix(complexMatrix& m)
{
	//std::cout << "\n\nOOOGGGG" << m.size() << "\n\n";
	for (int i = 0; i < m.size(); i++)
	{
		for (int j = 0; j < m[0].size(); j++)
		{
			std::cout << "(" + complexToString(m[i][j]) + ") ";
		}
		std::cout << "\n";
	}
	return m;
}
complexMatrix & complexMatrixAdditiveInverse(complexMatrix & m)
{
	//multiply each entry in m by -1
	complexMatrixScalarMultiplication(-1.0, m);

	return m;
}
complexMatrix & complexMatrixAddition(complexMatrix & m, complexMatrix & n, complexMatrix & res)
{
	//each entry of the output is the sum of the corresponding entrys in input matrices

	//printComplexMatrix(m);
//	std::cout << "\n\n";
	//printComplexMatrix(n);
	res.resize(m.size());
	for (int i = 0; i < m.size(); i++)
	{
		res[i].resize(m[0].size());
	}
	for (int i = 0; i < m.size(); i++)
	{
		for (int j = 0; j < m[0].size(); j++)
		{
			res[i][j] = complexAdd(m[i][j],n[i][j]);
		}
	}
	return res;
}

//select row then column [][]
//Drill 2.2.2/2.2.3
complexMatrix& complexMatrixMultiplication(complexMatrix & a, complexMatrix & b, complexMatrix & out)
{
	//row by row column by column matrix multiplication, row is accessed at the first []

	//std::cout << "\n\n";
	//printComplexMatrix(a);
	//std::cout << "\n\n";
//	printComplexMatrix(b);
	out.resize(a.size());
	//std::cout <<"\n\n" << a.size() <<"\n\n" << b[0].size()<<"\n\n";
	for (int i = 0; i < out.size(); i++)
	{
		out[i].resize(b[0].size());
	}
	//std::cout << a.size() << b[0].size();
	for (int i = 0; i < a.size(); i++)
	{
		for (int k = 0; k < b[0].size(); k++)
		{
			complex c = { 0.0,0.0 };
			for (int j = 0; j < a[0].size(); j++)
			{
			//	std::cout << complexToString(b[j][k]) + "*" + complexToString(a[i][j]) + "\n\n\n";
				c = complexAdd(c, complexMultiply(b[j][k], a[i][j]));
			}
			out[i][k] = c;
		}
	}

	return out;
}
//Drill 2.4.1
complexMatrix& complexMatrixConjugate(complexMatrix& a)
{
	//the conjugate of a complex number a + bi  is a - bi
	//std::cout << a[0].size() << "   " << a.size();
	for (int i = 0; i < a.size(); i++)
	{
		for (int j = 0; j < a[0].size(); j++)
		{
			a[i][j] = { a[i][j].a,-1 * a[i][j].b };
		}
	}
	
	
	return a;
}
complexMatrix& copyMatrix(complexMatrix& a,complexMatrix& b)
{
	//copy a matrix
	b.resize(a.size());
	//std::cout << a.size();
	for (int i = 0; i < a.size(); i++)
	{
		b[i].resize(a[0].size());
		for (int j = 0; j < a[0].size(); j++)
		{
			b[i][j] = a[i][j];
		}
	}

	return b;
}
complexMatrix& complexMatrixTranspose(complexMatrix& a)
{
	//entry i,j of the matrix is place at entry j,i of the output
	int s = a[0].size();
	int z = a.size();
	complexMatrix b;
	copyMatrix(a, b);
	a.resize(s);
	for (int i = 0; i < s; i++)
	{
		a[i].resize(z);
	}


	for (int i = 0; i < s; i++)
	{
		
		for (int j = 0; j < z; j++)
		{
			a[i][j] = b[j][i];
		}
	}

	return a;
}
complexMatrix& complexMatrixAdjoint(complexMatrix& a)
{
	//take the transpose and adjoint of a, order doesnt matter
	return complexMatrixTranspose(complexMatrixConjugate(a));
}
complex complexTrace(complexMatrix& m)
{
	//add all values in the matrix where i = j or Kroneker Delta i,j = 1
	complex c = { 0.0,0.0 };
	for (int i = 0; i < m.size() && i < m[0].size(); i++)
	{
		c = complexAdd(c, m[i][i]);
	}

	return c;
}

complex complexMatrixInnerProduct(complexMatrix& a, complexMatrix& b)
{
	//inner product defined for complex matrices, vectors are matrices with width or height 1
	complexMatrix g;
	return complexTrace(complexMatrixMultiplication((a), complexMatrixAdjoint(b),g));
}
//Drill 2.4.2
complex complexMatrixNorm(complexMatrix& a)
{
	//sqrt(a^2 + b^2 + c^2...) makes sense in 3R otherwise just a definition
	complexMatrix d;
	copyMatrix(a, d);
	polarComplex c = cartesianToPolar(complexMatrixInnerProduct(a,d));
	c.t /= 2.0;
	c.m = pow(c.m,0.5);
	return polarToCartesian(c);
}
//Drill 2.4.3
complex complexMatrixDistance(complexMatrix& a, complexMatrix& b)
{
	//subtract matricesand the take norm to get distance
	complexMatrix g;
	complexMatrixAddition(a, complexMatrixAdditiveInverse(b), g);
	
	return complexMatrixNorm(g);
}
//Drill 2.6.1
bool doubleEquality(double d1, double d2)
{
	//probably not neccesary at all just use ==
	if (abs(d1 - d2) < 0.0001)
	{
		return true;
	}
	return false;
}
bool complexEquality(complex a, complex b)
{
	//equality for complex values
	return doubleEquality(a.a,b.a) && doubleEquality(a.b, b.b);
}
bool complexMatrixEquality(complexMatrix& a, complexMatrix& b)
{
	//are all the entrys of two matrices the same?
	if (a.size() == b.size() && b[0].size() == a[0].size())
	{
		for (int i = 0; i < a.size(); i++)
		{
			b[i].resize(a[0].size());
			for (int j = 0; j < a[0].size(); j++)
			{
				if (!complexEquality(b[i][j], a[i][j]))
					return false;
			}
		}
		return true;
	}
	else
		return false;
}
bool complexMatrixIsHermetian(complexMatrix& m)
{
	//is the adjoint of a matrix itself?
	complexMatrix g;
	copyMatrix(m, g);
	std::cout << "\n\n";
	printComplexMatrix(g);
	complexMatrixAdjoint(g);
	std::cout << "\n\n";
	printComplexMatrix(g);
	std::cout << "\n\n";
	printComplexMatrix(m);
	std::cout << "\n\n";
	return complexMatrixEquality(m,g);
}
//Drill 2.6.2
bool complexMatrixIsUnitary(complexMatrix& m)
{
	//matrix times adjoint equal to identity?
	complexMatrix ma;
	copyMatrix(m, ma);
	complexMatrixAdjoint(ma);
	complexMatrix I;
	complexMatrix o1;
	complexMatrix o2;
	I.resize(m.size());
	for (int i = 0; i < m.size(); i++)
	{
		I[i].resize(m[0].size());
	}
	for (int i = 0; i < I.size(); i++)
	{
		
		for (int j = 0; j < I[0].size(); j++)
		{
			if (i == j)
			{
				I[i][j] = { 1.0,0.0 };

			}
			else 
			{
				I[i][j] = { 0.0,0.0 };
			}
		}
	}
	if (complexMatrixEquality(complexMatrixMultiplication(m, ma, o1), I) && complexMatrixEquality(complexMatrixMultiplication(ma, m, o2), I))
	{
		return true;
	}
	return false;
}
complexMatrix& complexMatrixTensorProduct(complexMatrix& a,complexMatrix& b,complexMatrix & out)
{
	//tensor of two matrices, scalar each matrix of b by corresponding value of a
	out.resize(a.size() * b.size());
	for(int i = 0;i < out.size();i++)
	{
		out[i].resize(a[0].size() * b[0].size());
	}
	for (int i = 0; i < out.size(); i++)
	{
		for (int j = 0; j < out[0].size(); j++)
		{
			out[i][j] = complexMultiply(a[i / b.size()][j / b[0].size()], b[i % b.size()][j%b[0].size()]);
		}
	}
	return out;
}
//4.1.1
double ketStateProbability(complexMatrix& ket, int x)
{
	//find the probability of finding the particle at discrete location x
	return pow(complexModulus(ket[x][0]), 2) / pow(complexMatrixNorm(ket).a, 2);
}
complex ketStateTransitionAmplitude(complexMatrix& ketFrom, complexMatrix& ketToo)
{
	//amplitude of transition between two states
	//complexToString();
	complexMatrixAdjoint(ketToo);
	complexMatrix t;
	
	complexMatrixMultiplication(ketToo, ketFrom, t);
	complexMatrixAdjoint(ketToo);
	return complexTrace(t);
}
int main()
{

	//Drill 1.1.1
	std::cout << "Give two complex numbers, first the real magnitude then the complex magniute(4 total floats):\n";
	complex w;
	complex z;
	std::cin >> w.a;
	std::cin >> w.b;
	std::cin >> z.a;
	std::cin >> z.b;
	std::cout << "\nYour given numbers: " + complexToString(w) + " | " + complexToString(z) + "\n\n";
	std::cout << "Drill 1.1.1\n";
	std::cout << "Multiplied: " + complexToString(complexMultiply(w, z)) + " Added: " + complexToString(complexAdd(w, z)) + "\n\n";
	//Drill 1.2.1
	std::cout << "Drill 1.2.1\n";
	std::cout << "Subracted: " + complexToString(complexSubtract(w, z)) + " Divided: " + complexToString(complexDivide(w, z)) + "\n";
	std::cout << "Length/Modulus of your two complex numbers: " + std::to_string(complexModulus(w)) + " | " + std::to_string(complexModulus(z)) + "\n\n";
	//Drill 1.3.1
	std::cout << "Drill 1.3.1\n";
	std::cout << "Polar of your complex numbers: " + complexToPolarString(w) + " | " + complexToPolarString(z) + " Testing of polar to complex(will just repeat the numbers you gave):" + complexToString(polarToCartesian(cartesianToPolar(w))) + " | " + complexToString(polarToCartesian(cartesianToPolar(z))) + "\n";


	//Drill 2.1.1/2.2.1
	std::cout << "Drill 2.1.1/2.2.1\n";
	complexMatrix  mat1 =
	{ {{3,2}, {-1,0}, {2,0}},
	{ {3,3}, {7,-6}, {0,0}},
	{ {4,5}, {5,5}, {-1,-1}} };
	complexMatrix  mat2 =
	{ {{1,0}, {0,0}, {1,0}},
	{ {0,0}, {1,0}, {0,0}},
	{ {1,0}, {0,0}, {2,0} } };
	std::cout << "Matrix 1:\n";
	printComplexMatrix(mat1);
	std::cout << "\n";
	std::cout << "\n";
	std::cout << "Matrix 2:\n";
	printComplexMatrix(mat2);
	std::cout << "\n";
	std::cout << "\n";
	std::cout << "Additive inverse of Matrix 1:\n";
	complexMatrixAdditiveInverse(mat1);
	printComplexMatrix(mat1);
	std::cout << "\n";
	std::cout << "\n";
	std::cout << "Subtraction of Matrix 1 from Matrix 2:\n";
	complexMatrix  addedMat;
	complexMatrixAddition(mat1, mat2, addedMat);
	printComplexMatrix(addedMat);
	std::cout << "\n";
	std::cout << "\n";
	//Drill 2.2.2/2.2.3
	std::cout << "Drill 2.2.2/2.2.3\nMultiplying the additive inverse by Matrix 2\n";
	complexMatrix  multMat;
	complexMatrixMultiplication(mat2, mat1, multMat);
	printComplexMatrix(multMat);
	std::cout << "\n\n";
	//Drill 2.4.1
	complexMatrix vec1 = { {{33,330},{0,34},{324,0}} };
	complexMatrix vec2 = { {{-5,-1},{0,0},{4,0}} };
	complexMatrix abc;
	std::cout << "Drill 2.4.1 \nInner Product of complex matrices\n";
	std::cout << complexToString(complexMatrixInnerProduct(vec1, vec2));
	//Drill 2.4.2
	std::cout << "\nDrill 2.4.2 \nNorm of a complex vector \n";
	std::cout << complexToString(complexMatrixNorm(vec1));
	//Drill 2.4.3
	vec1 = { {{33,330},{0,34},{324,0}} };
	vec2 = { {{-5,-1},{0,0},{4,0}} };

	std::cout << "\nDrill 2.4.3 distance of two complexvectors\n";
	std::cout << "\nFirst Vector:\n";
	printComplexMatrix(vec1);
	std::cout << "\n\n";
	std::cout << "\nSecond Vector:\n";
	printComplexMatrix(vec2);
	std::cout << complexToString(complexMatrixDistance(vec1, vec2));
	//Drill 2.6.1
	vec1 =
	{ {{1,0}, {1,3}, {4,-7}},
	{ {1,-3}, {2,0}, {3,5.000}} ,
	{ {4,7},{3,-5},{1,0}} };
	if (complexMatrixIsHermetian(vec1))
	{
		std::cout << "\nHERMETIAN\n";

	}
	else
	{
		std::cout << "\nNOT HERMETIAN\n";
	}
	//Drill 2.6.2
	std::cout << "\n\nDrill 2.6.2\n\n";
	complexMatrix unitaryTest = { {{0.5,0.5} ,{0,1/pow(3,0.5)}, {1.5/pow(15,0.5),0.5 / pow(15,0.5)}},
								{{-0.5,0.0}, {1/pow(3,0.5),0.0}, {2/pow(15,0.5),1.5/pow(15,0.5)}},
								{ {0.5,0.0},{0.0,-1/pow(3,0.5)},{0.0,2.5/pow(15,0.5)} } };
	if (complexMatrixIsUnitary(unitaryTest))
	{
		std::cout << "\nUNITARY\n";

	}
	else
	{
		std::cout << "\nNOT UNITARY\n";
	}
	//Drill 2.7.1
	std::cout << "\n\nDrill 2.7.1\n\n";
	complexMatrix tens1 =     { {{3,2} ,{5,-1}, {0,2}},
								{{0.0,0.0}, {12,0.0}, {6,-3}},
								{ {2,0.0},{4,4},{9.0,3} } };
	complexMatrix tens2 =       { {{1,0} ,{3,4}, {5,-7}},
								{{10,2.0}, {6,0.0}, {2,5}},
								{ {0.0,0.0},{1.0,0},{2.0,9} } };
	complexMatrix tout;
	complexMatrixTensorProduct(tens1, tens2, tout);
	std::cout << "\n\n";
	printComplexMatrix(tout);
	std::cout << "\n\n";
	//Drill 4.1.1
	complexMatrix ket1 = { { {-3,-1}, {0,-2}, {0,1}, {2,0} } };//row ket 1
	complexMatrix ket2 = { {  {1,0}, {0,1} } };//row ket 2
	complexMatrix ket3 = { {  {0,1}, {-1,0} } };//row ket 3
	complexMatrixScalarMultiplication(pow(2, 0.5) / 2, ket2);
	complexMatrixScalarMultiplication(pow(2, 0.5) / 2, ket3);
	complexMatrix ketout;
	std::cout << "\nDrill 4.1.1\nPosition and Transistion probabilities for kets\n";
	complexMatrixTranspose(ket1);
	complexMatrixTranspose(ket2);
	complexMatrixTranspose(ket3);
	//complexMatrixAdjoint(ket3);
	std::cout << "Probability of finding the particle at position 2 of ket 1: " << ketStateProbability(ket1,2);
	std::cout << "\nProbability of particle state transitioning from ket 2 to ket 3: " << complexToString(ketStateTransitionAmplitude(ket2,ket3)) ;
	//prevent window close
	std::cin >> z.b;
}
	







