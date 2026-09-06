/*!	\file	complex.h
 *
 *	\brief	Définition de la classe complex pour rendre compatible \n
 *			le misérable MS Visual C++ avec la norme ANSI-C99 !
 */

#ifndef COMPLEX_H
#define COMPLEX_H

#include <math.h>



class complex  
{
public:
	double real, imag;
	
	// Constructeur et destructeur
	complex(double a=0, double b=0){
		real = a;
		imag = b;
	}
	virtual ~complex(){};


	// Fonctions membres
	double creal (){
		return real;
	}

	double cimag (){
		return imag;
	}

	complex csqrt (){
		double modd = sqrt(real*real+imag*imag);
		double argg = atan (imag/real);
		return complex(sqrt(modd)*cos(argg/2) , sqrt(modd)*sin(argg/2));
	}

	complex conj (){
		return complex(real , -imag);
	}

	complex cexp (){
		return complex(exp(real)*cos(imag) , exp(real)*sin(imag));
	}
	
	double cabs (){
		return sqrt(real*real+imag*imag);
	}

	// Opérateurs
	complex operator+(complex z){
		return complex(real+z.real,imag+z.imag);
	}

	complex operator-(complex z){
		return complex(real-z.real,imag-z.imag);
	}

	complex operator*(complex z){
		return complex(real*z.real-imag*z.imag , real*z.imag+imag*z.real);
	}

	complex operator*(double r){
		return complex(r*real, r*imag);
	}

	complex operator/(complex z){
		double modd = sqrt(real*real+imag*imag)/sqrt(z.real*z.real+z.imag*z.imag);
		double argg = atan (imag/real) - atan(z.imag/z.real);
		return complex(modd*cos(argg) , modd*sin(argg));
	}

	void operator+=(complex z){
		real += z.real;
		imag += z.imag;
	}

	void operator-=(complex z){
		real -= z.real;
		imag -= z.imag;
	}

	void operator*=(complex z){
		real = z.real*real - z.imag*imag;
		imag = z.real*imag + z.imag*real;
	}
	
	bool operator==(complex z){
		return (real==z.real && imag==z.imag);
	}


	friend complex operator*(double r,complex z);

	friend complex operator/(double r,complex z);

/*------------------*/
};

complex operator*(double r,complex z){
	return complex(r*z.real, r*z.imag);
}

complex operator/(double r,complex z){
	return (complex)r / z;
}

complex operator-(complex z){
	return complex(-(z.real), -(z.imag));
}

#define I				complex(0,1)
#define csqrt(z)		(z).csqrt()
#define cexp(z)			(z).cexp()
#define conj(z) 		(z).conj()
#define cabs(z)			(z).cabs()
#define creal(z)		(z).creal()
#define cimag(z)		(z).cimag()
#define c_omplex(a,b)	complex(a,b) 

#endif // COMPLEX_H
