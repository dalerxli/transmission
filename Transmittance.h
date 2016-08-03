#ifndef TRANSMITTANCE_H_
#define TRANSMITTANCE_H_

#include <algorithm>
#include <iostream>
#include <vector>
#include <time.h>
#include "SubstrateRefractive.h"
#include "const.h"
#include "adept.h"
#include "func.h"


enum SubType { AIR, GLASS, SI, AQUARTZ, CQUARTZ, SILICON, SODA };
using adept::adouble;

class Transmittance {
public:
	SubstrateRefractive *subrefr; //class for calculate substrate refractive index dependence on lambda
	std::vector<double> wavelength; //wavelength mass
	std::vector<double> calc;
	/* Unknows: 
		0 - u, 1 - u1, 2 - W(n-2), 3 - W(n-3),... N - W1
		N+1 - v, N+2 - v1, N+3 - Z(n-2), 2N - Z1 */
	std::vector<double> der; // derivation of func(unknowns)
	std::vector<double> n; // refractive index
	std::vector<double> k; // absorption coefficient
	double wave_infl; // inflation wavelength
	double thickness; // thickness of layer
	std::vector<double> sub; // substrate refractive
	size_t size; // number of given points

	void init_value(double& x, double x0, double x1, double xs);

public:
	Transmittance() : wavelength(0), calc(0), size(0), thickness(0), wave_infl(0) {
		subrefr = new SubstrateRefractive();
	}
	void readfile(const std::string& filename); // read file and get wavelength and cal vectors
	void readtrans(const std::string& filename); // read basefile and get wavelength and cal vectors
	void readthick(const std::string& filename); // read thickness from xml file
	void print_data(std::ostream& file) const;
	void print_der(std::ostream& file) const;
	void print_index(std::ostream& file) const;
	void print_result(std::ostream& file) const;
	void setSubstrate(SubType);
	void setThickness(double thick) { thickness = thick; }
	size_t getSize() const { return size; }
	void setWaveInflation(double wave) { wave_infl = wave;  }
	void set_index(const std::vector<double>& unknowns);
	void setTrueIndex();
	void setTrueUnknowns(std::vector<double>& x);
	void printTrueIndex(std::ostream& file);
	void set_initial(std::vector<double>& x, double nstart, double nend, double kval);

	//compute substrate refractive index only
	void computeCleanSubstrateRefractive();
	//theory formula for given n, k, s, wavelength
	template<class X, class Y>
	X TransTheory(const X& n, const X& k, const Y& s, const Y& wavelength) const;

	template<class X, class Y>
	X TransTheoryInf(const X& n, const X& k, const Y& s, const Y& wavelength) const;

	//functional for minimization
	template<class xdouble>
	void getNK(std::vector<xdouble>& n, std::vector<xdouble>& k, const std::vector<xdouble>& x);
	//adouble MinFunctional(const std::vector<adouble>& x);
	template<class xdouble>
	xdouble MinFunctional(const std::vector<xdouble>& x);

	//gradient of MinFunctional, refresh derivative vector
	//store gradient in der
	double gradient(const std::vector<double>& x_val);

	//optimization procedure
	//find unknown parameters
	
	double optimize(std::vector<double>& x, int max_iter);
	double optimizeCG(std::vector<double>& x, int max_iter);

	//optimization procedure with account of 
	//different initial conditions
	//return best approximation 
	void optimization(
		double ns0, double ns1, double nss,
		double nf0, double nf1, double nfs,
		double k0, double k1, double ks,
		int maxiter);
};

template<class X, class Y>
X Transmittance::TransTheory(const X& n, const X& k, const Y& s, const Y& wavelength) const {
	try {
		X phi = (n*PI4*thickness) / wavelength;
		X A = 16 * s*(n*n + k*k);
		X B = ((n + 1)*(n + 1) + k*k)*((n + 1)*(n + s*s) + k*k);
		X C = ((n*n - 1 + k*k)*(n*n - s*s + k*k) - 2 * k*k*(s*s + 1)) * 2 * cos(phi) - \
			k*(2 * (n*n - s*s + k*k) + (s*s + 1)*(n*n - 1 + k*k)) * 2 * sin(phi);
		X D = ((n - 1)*(n - 1) + k*k)*((n - 1)*(n - s*s) + k*k);
		X alpha = PI4*k / wavelength;
		X x = exp(alpha*(-thickness)); 
		return A*x / (B - C*x + D*x*x);
	}
	catch (std::overflow_error& e) {
		std::cout << e.what() << '\n';
	}
}

template<class X, class Y>
X Transmittance::TransTheoryInf(const X& n, const X& k, const Y& s, const Y& wavelength) const {
	try {
		X phi = (n*PI4*thickness) / wavelength;
		X A = 16 * s*(n*n + k*k);
		X B = ((n + 1)*(n + 1) + k*k)*((n + s)*(n + s) + k*k);
		X C = ((n*n - 1 + k*k)*(n*n - s*s + k*k) + 4 * k*k*s) * 2 * cos(phi) - \
			k*(2 * (n*n - s*s + k*k) + 2*s*(n*n - 1 + k*k)) * 2 * sin(phi);
		X D = ((n - 1)*(n - 1) + k*k)*((n - s)*(n - s) + k*k);
		X alpha = PI4*k / wavelength;
		X x = exp(alpha*(-thickness));
		return A*x / (B - C*x + D*x*x);
	}
	catch (std::overflow_error& e) {
		std::cout << e.what() << '\n';
	}
}

template<class xdouble>
xdouble Transmittance::MinFunctional(const std::vector<xdouble>& x) {
	std::vector<xdouble> n(size);
	std::vector<xdouble> k(size);
	getNK(n, k, x);
	xdouble tmp;
	xdouble y = 0.0;
	for (size_t i = 0; i < size; ++i) {
		tmp = TransTheory(n[i], k[i], sub[i], wavelength[i]);
		y += (tmp - calc[i]) * (tmp - calc[i]);
	}
	return y;
}

template<class xdouble>
void Transmittance::getNK(std::vector<xdouble>& n, std::vector<xdouble>& k,
	const std::vector<xdouble>& x) {
	size_t size = x.size() / 2;
	double h, h1, h2;
	h = wavelength[size - 1] - wavelength[size - 2];
	n[size - 1] = 1.0 + x[size - 1] * x[size - 1];
	k[size - 1] = x[2 * size - 1] * x[2 * size - 1];
	n[size - 2] = n[size - 1] + x[size - 2] * x[size - 2] * h;
	k[size - 2] = k[size - 1] + x[2 * size - 2] * x[2 * size - 2] * h;
	for (size_t i = 3; i <= size; ++i) {
		h1 = wavelength[size - i + 1] - wavelength[size - i];
		h2 = wavelength[size - i + 2] - wavelength[size - i + 1];
		n[size - i] = x[size - i] * x[size - i] * h1*(h1 + h2) * 0.5 + ((h1 + h2) / h2)*n[size - i + 1] - (h1 / h2)*n[size - i + 2];

		if (wavelength[size - i] > wave_infl)
			k[size - i] = x[2 * size - i] * x[2 * size - i] * h1*(h1 + h2) * 0.5 + ((h1 + h2) / h2)*k[size - i + 1] - (h1 / h2)*k[size - i + 2];
		else
			k[size - i] = -1.0*x[2 * size - i] * x[2 * size - i] * h1*(h1 + h2) * 0.5 + ((h1 + h2) / h2)*k[size - i + 1] - (h1 / h2)*k[size - i + 2];
	}
}

#endif // TRANSMITTANCE_H_