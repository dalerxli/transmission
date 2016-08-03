#include <fstream>
#include <string>
#include <sstream>
#include "Transmittance.h"


void Transmittance::readfile(const std::string& filename) {
	std::ifstream file(filename);
	double wave;
	double trans;
	int count = 0;
	while (file >> wave >> trans) {
		if (wave > 549 && wave < 1531) {
			if (count % 10 == 0) {
				wavelength.push_back(wave);
				calc.push_back(trans);
			}
			//wavelength.push_back(wave);
			//calc.push_back(trans);
			count++;
		}
	}
	size = wavelength.size();
	n.assign(size, 0);
	k.assign(size, 0);
	der.assign(2 * size, 0);
	file.close();
}

void Transmittance::readtrans(const std::string& filename) {
	std::ifstream file(filename);
	double wave;
	double trans;
	int count = 0;
	std::string token = "transmittance, T";
	std::string token_data_st = "<DATA>";
	std::string token_data_end = "</DATA>";
	std::string line, tmp;
	size_t pos, pos_start, pos_end;
	getline(file, line);
	while (line.find(token) == std::string::npos)
		getline(file, line);
	pos = line.find(token);
	pos_start = line.find(token_data_st, pos);
	pos_end = line.find(token_data_end, pos_start);
	std::string data = line.substr(pos_start, pos_end - pos_start);
	int comm = 0; // if comm%2 == 0 - wave, else == 1 - transmittance
	for (size_t i = pos_start; i < pos_end; ++i) {
		tmp = "";
		if (line[i] == '\"') {
			i = i + 1; //into a number
			while (line[i] != '\"') {
				tmp += line[i]; // read number to string;
				i++;
			}
			i++; //out from number
			if (comm % 2 == 0) {
				wave = atof(tmp.c_str());
				wavelength.push_back(wave);
			}
			else {
				trans = atof(tmp.c_str());
				calc.push_back(trans * 0.01);
			}
			comm++;
		}
	}
	size = wavelength.size();
	n.assign(size, 0);
	k.assign(size, 0);
	der.assign(2 * size, 0);
	file.close();
}

void Transmittance::readthick(const std::string& filename) {
	std::ifstream file(filename);
	std::string line, tmp;
	size_t pos, pos_start, pos_end, pos_dst, pos_est;
	std::string token = ">thickness<";
	std::string token_data_st = "<AT_VALUE>";
	std::string token_data_end = "</AT_VALUE>";
	std::string dim_st = "<AT_UNIT>";
	std::string dim_end = "</AT_UNIT>";
	getline(file, line);
	while (line.find(token) == std::string::npos)
		getline(file, line);
	pos = line.find(token);
	pos_start = line.find(token_data_st, pos);
	pos_end = line.find(token_data_end, pos_start);
	if (pos_start != std::string::npos) {
		size_t j;
		for (size_t i = pos_start; i < pos_end; ++i) {
			if (line[i] == '>') {
				j = i + 1;
				while (j < pos_end) {
					tmp += line[j];
					j++;
				}
				thickness = atof(tmp.c_str());
			}
		}
		pos_dst = line.find(dim_st, pos_end);
		pos_est = line.find(dim_end, pos_end);
		tmp = "";
		for (size_t i = pos_dst; i < pos_est; ++i) {
			if (line[i] == '>') {
				j = i + 1;
				while (j < pos_est) {
					tmp += line[j];
					j++;
				}
			}
		}
		if (tmp == "mm") thickness *= 1e+06;
	}
	else {
		std::cout << "There is no thickness\n";
	}
	file.close();
}

void Transmittance::setSubstrate(SubType t) {
	if (t == AIR) {
		sub = subrefr->air(wavelength);
	}

	if (t == GLASS) {
		sub = subrefr->glass(wavelength);
	}

	if (t == SI) {
		sub = subrefr->si(wavelength);
	}

	if (t == AQUARTZ) {
		sub = subrefr->AQuartz(wavelength);
	}

	if (t == CQUARTZ) {
		sub = subrefr->CQuartz(wavelength);
	}

	if (t == SILICON) {
		sub = subrefr->silicon(wavelength);
	}

	if (t == SODA) {
		sub = subrefr->sodaglass(wavelength);
	}
}

void Transmittance::print_data(std::ostream& file) const {
	std::cout << "Printing given wavelength and transmittance...\n";
	for (size_t i = 0; i < wavelength.size(); ++i) {
		file << wavelength[i] << ' ' << calc[i] << '\n';
	}
}

void Transmittance::print_der(std::ostream& file) const {
	std::cout << "Printing functional derivation vector...\n";
	for (size_t i = 0; i < size; ++i) {
		file << der[i] << ' '  << der[size + i] << '\n';
	}
}




void Transmittance::print_index(std::ostream& file) const {
	std::cout << "Printing n and k...\n";
	std::vector<double> alpha(k.size());
	for (size_t i = 0; i < k.size(); ++i) {
		alpha[i] = k[i] * PI4 / wavelength[i];
	}
	for (size_t i = 0; i < wavelength.size(); ++i) {
		//std::cout << wavelength[i] << ' ' << n[i] << ' ' << alpha[i]*1e7 << '\n';
		file << wavelength[i] << ' ' << n[i] << ' ' << k[i] << '\n';
	}
}

void Transmittance::print_result(std::ostream& file) const {
	std::cout << "Printing wave and calculated transmittance...\n";
	for (size_t i = 0; i < wavelength.size(); ++i) {
		file << wavelength[i] << ' ' << TransTheory(n[i], k[i], sub[i], wavelength[i]) << '\n';
	}
}

void Transmittance::setTrueIndex() {
	std::cout << "Setting true n and k indexes...\n";
	double E;
	for (size_t i = 0; i < wavelength.size(); ++i) {
		E = 1240.0 / wavelength[i];
		n[i] = std::sqrt(1 + 1 / (0.09195 - 12600 / (wavelength[i] * wavelength[i])));
		if (wavelength[i] > 1000) {
			k[i] = 1.3e-5;
		}
		/*else if (E < 1.40 && E > 0.60) {
			k[i] = exp(6.5944*1e-6*exp(9.0846*E) - 16.102) * wavelength[i] * PI_m4;
		}*/
		else if (wavelength[i] <= 1000.0 && wavelength[i] > 885) {
			k[i] = exp(6.5944*1e-6*exp(9.0846*E) - 16.102) * wavelength[i] * PI_m4;
		}
		else if (E >= 1.40 && E < 1.75) {
			k[i] = exp(20 * E - 41.9) * wavelength[i] / PI4;
		}
		else if (E >= 1.75 && E <= 2.29) {
			k[i] = exp(sqrt(59.56*E - 102.1) - 8.391) * wavelength[i] * PI_m4;
		}
	}
}

void Transmittance::setTrueUnknowns(std::vector<double>& x) {
	double h1, h2;
	size_t size = wavelength.size();
	double h = wavelength[size - 1] - wavelength[size - 2];
	x[size - 1] = sqrt(n[size-1] - 1);
	x[2 * size - 1] = sqrt(k[size-1]);
	x[size - 2] = sqrt((n[size - 2] - n[size - 1]) / h);
	x[2 * size - 2] = sqrt((k[size - 2] - k[size - 1]) / h);
	for (size_t i = 3; i <= size; ++i) {
		h1 = wavelength[size - i + 1] - wavelength[size - i];
		h2 = wavelength[size - i + 2] - wavelength[size - i + 1];
		x[size - i] = sqrt((2 * (h2*n[size - i + 2] - (h1 + h2)*n[size - i + 1] + h1*n[size - i])) / (h1*h2*(h1 + h2)));
		x[2 * size - i] = sqrt((2 * (h2*k[size - i + 2] - (h1 + h2)*k[size - i + 1] + h1*k[size - i])) / (h1*h2*(h1 + h2)));
	}
}


void Transmittance::printTrueIndex(std::ostream& file) {
	for (size_t i = 0; i < wavelength.size(); ++i) {
		double E;
		E = 1240.0 / wavelength[i];
		if (E < 1.40 && E > 0.60) {
			file << wavelength[i] << ' ' << \
				std::sqrt(1 + 1 / (0.09195 - 12600 / (wavelength[i] * wavelength[i]))) \
				<< ' ' << exp(6.5944*1e-6*exp(9.084*E) - 16.102) * wavelength[i] / PI4 << '\n';
		}
		else if (E >= 1.40 && E < 1.75) {
			file << wavelength[i] << ' ' << \
				std::sqrt(1 + 1 / (0.09195 - 12600 / (wavelength[i] * wavelength[i]))) \
				<< ' ' << exp(20 * E - 41.9) * wavelength[i] / PI4 << '\n';
		}
		else if (E >= 1.75 && E <= 2.30) {
			file << wavelength[i] << ' ' << \
				std::sqrt(1 + 1 / (0.09195 - 12600 / (wavelength[i] * wavelength[i]))) \
				<< ' ' << exp(sqrt(59.56*E - 102.1) - 8.391) * wavelength[i] / PI4 << '\n';
		}
	
	}

}

void Transmittance::computeCleanSubstrateRefractive() {
	for (size_t i = 0; i < size; i++) {
		n[i] = 1.0 / calc[i] + sqrt(1.0 / (calc[i] * calc[i]) - 1);
	}
}

void Transmittance::set_index(const std::vector<double>& x) {
	size_t size = x.size() / 2;
	double h1, h2;
	std::cout << "setting index..\n";
	double h = wavelength[size - 1] - wavelength[size - 2];
	n[size - 1] = 1 + x[size - 1] * x[size - 1];
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

void Transmittance::set_initial(std::vector<double>& x, double nstart, double nend, double kval) {
	std::cout << "setting initial...\n";
	double h1, h2;
	size_t size = wavelength.size();
	size_t size2 = size * 0.2;
	double k1 = 0.1 * kval;
	double k2 = 1e-10 * kval;
	double ntga = (nstart - nend) / (wavelength[size - 1] - wavelength[0]);
	double ktga1 = (k1 - k2) / (wavelength[size - 1] - wavelength[size2]);
	double ktga2 = (kval - k1) / (wavelength[size2] - wavelength[0]);
	double h = wavelength[size - 1] - wavelength[size - 2];
	x[size - 1] = sqrt(nend - 1);
	x[2 * size - 1] = sqrt(k2);
	n[size - 1] = nend; k[size - 1] = k2;
	x[size - 2] = sqrt(ntga); x[2 * size - 2] = sqrt(ktga1);
	n[size - 2] = nend + h*ntga;
	k[size - 2] = k2 + h*ktga1;
	for (size_t i = 3; i <= size; ++i) {
		h1 = wavelength[size - i + 1] - wavelength[size - i];
		h2 = wavelength[size - i + 2] - wavelength[size - i + 1];
		n[size - i] = n[size - i + 1] + h1*ntga;
		if (i < size - size2 + 1)
			k[size - i] = k[size - i + 1] + h1*ktga1;
		else 
			k[size - i] = k[size - i + 1] + h1*ktga2;

		x[size - i] = sqrt(std::max(0.0, 2*(h1*n[size - i + 2] - (h1 + h2)*n[size - i + 1] + h2*n[size - i]) / (h1*h2*(h1 + h2))));
		x[2 * size - i] = sqrt(std::max(0.0, 2 * (h1*k[size - i + 2] - (h1 + h2)*k[size - i + 1] + h2*k[size - i]) / (h1*h2*(h1 + h2))));
	}
}

double Transmittance::gradient(const std::vector<double>& x_val) { 
	adept::Stack stack; // Where the derivative information is stored
	// Import adouble from adept
	using adept::adouble;
	size_t size = x_val.size();
	std::vector<adouble> x(size); // Initialize active input variables
	for (size_t i = 0; i < size; i++) {
		x[i] = x_val[i];
	}
	stack.new_recording(); // Start recording
	adouble y = MinFunctional(x); // Call version overloaded for adouble args
	y.set_gradient(1.0); // Defines y as the objective function
	stack.compute_adjoint(); // Run the adjoint algorithm
	for (size_t i = 0; i < size; i++) {
		der[i] = x[i].get_gradient(); // Store the first gradient
	}
	return y.value(); // Return the result of the simple computation
}

//Follow Birgin 1999
double Transmittance::optimize(std::vector<double>& x, int max_iter) {
	std::cout << "starting optimize Birgin 1999..\n";
	size_t size = x.size();
	std::vector<double> x_tmp(size);
	std::vector<double> x_old(size);
	std::vector<double> der_old(size);
	std::vector<double> s(size);
	std::vector<double> y(size);
	std::vector<double> d(size);
	double alpha, lambda, lambda_tmp;
	double func_val, func_val_tmp, func_tmp;
	double delta;
	double beta;
	gradient(x);
	double error = MinFunctional(x);
	double der_value = norm(der);
	//double error = norm(der);
	alpha = ALPHA_MIN + (ALPHA_MAX - ALPHA_MIN) * mRandNorm<double>();
	int k = 0;
	while (error > 4e-06 && der_value > 1e-16 && k < max_iter) {
		// d and x_tmp
		for (size_t i = 0; i < size; i++) {
			d[i] = -alpha*der[i];
			x_tmp[i] = x[i] + d[i];
		}
		//compute step length
		if (k == 0) func_val = MinFunctional(x);
		else if (k < M - 1) {
			func_val_tmp = MinFunctional(x);
			if (func_val < func_val_tmp)
				func_val = func_val_tmp;
		}
		delta = scalar(der, d);
		lambda = 1;
		func_tmp = MinFunctional(x_tmp);
		while ((func_tmp > func_val + GAMMA * lambda * delta)) {
			if (lambda <= SIGMA_MIN)
				lambda = lambda * 0.5;
			else {
				lambda_tmp = (-0.5 * lambda * lambda * delta) / \
					(func_tmp - MinFunctional(x) - alpha*delta);
				if (lambda_tmp >= SIGMA_MIN && lambda_tmp <= SIGMA_MAX * alpha) {
					lambda = lambda_tmp;
				}
				else
					lambda = lambda * 0.5;
			}
			for (size_t i = 0; i < size; i++) {
				x_tmp[i] = x[i] + lambda * d[i];
			}
			func_tmp = MinFunctional(x_tmp);
		}
		for (size_t i = 0; i < size; i++) {
			x_old[i] = x[i];
			x[i] = x[i] + lambda * d[i];
			s[i] = x[i] - x_old[i];
			der_old[i] = der[i];
		}
		gradient(x);
		for (size_t i = 0; i < size; i++) {
			y[i] = der[i] - der_old[i];
		}
		beta = scalar(s, y);
		if (beta < 0)
			alpha = ALPHA_MAX;
		else
			alpha = std::min(ALPHA_MAX, std::max(ALPHA_MIN, scalar(s, s) / beta));
		der_value = norm(der);
		error = MinFunctional(x);
		if (!(k % 1000)) std::cout << k << ' ' << error << ' ' << der_value << '\n';	
		k++;
	}
	std::cout << "final value\n";
	std::cout << MinFunctional(x) << '\n';
	return error;
}

double Transmittance::optimizeCG(std::vector<double>& x, int max_iter) {
	size_t size = x.size();
	std::vector<double> x_tmp(size);
	std::vector<double> x_old(size);
	std::vector<double> der_old(size);
	std::vector<double> y(size);
	std::vector<double> d(size);
	double alpha, lambda, lambda_tmp;
	double func_val;
	double delta;
	double beta;
	gradient(x);
	double error = norm(der);
	int k = 0;
	for (size_t i = 0; i < size; i++) {
		d[i] = -der[i];
	}
	while (error > 1e-5) {
		if (k > max_iter) break;
		//compute step length
		delta = scalar(der, d);
		//lambda = ALPHA_MIN + (ALPHA_MAX - ALPHA_MIN) * mRandNorm<double>();
		lambda = 1;
		for (size_t i = 0; i < size; i++) {
			x_tmp[i] = x[i] + lambda*d[i];
		}
		func_val = MinFunctional(x);
		while (MinFunctional(x_tmp) > func_val + GAMMA * lambda * delta) {
			lambda *= (SIGMA_MIN + (SIGMA_MAX - SIGMA_MIN) * mRandNorm<double>());
			for (size_t i = 0; i < size; i++) {
				x_tmp[i] = x[i] + lambda * d[i];
			}
		}		for (size_t i = 0; i < size; i++) {
			x_old[i] = x[i];
			x[i] = x[i] + lambda * d[i];
			der_old[i] = der[i];
		}
		gradient(x);
		for (size_t i = 0; i < size; i++) {
			y[i] = der[i] - der_old[i];
		}
		beta = scalar(der, y) / scalar(der_old, der_old);
		for (size_t i = 0; i < size; i++) {
			d[i] = -der[i] + beta*d[i];
		}
		error = norm(der);
		if (k%100 == 0)
			std::cout << k << ' ' << error << ' ' << ' ' << func_val << ' ' << lambda << '\n';
		k++;
	}
	std::cout << "final value\n";
	std::cout << MinFunctional(x) << '\n';
	return error;
}


void Transmittance::optimization(
	double ns0, double ns1, double nss,
	double nf0, double nf1, double nfs,
	double k0, double k1, double ks,
	int maxiter) {

	double error = 100, error_tmp;
	double ns = ns0, nf = nf0, kk = k0;
	std::vector<double> x(getSize() * 2);
	std::vector<double> x_best;
	for (ns = ns0; ns <= ns1; ns += nss) {
		for (nf = nf0; nf <= nf1; nf += nfs) {
			for (kk = k0; kk <= k1; kk += ks) {
				if (nf >= ns) {
					set_initial(x, nf, ns, kk);
					error_tmp = optimize(x, maxiter);
					set_index(x);
					std::cout << ns << ' ' << nf << ' ' << kk << ' ' << error_tmp << ' ';
					if (error > error_tmp) {
						error = error_tmp;
						x_best = x;
						std::cout << "!!!\n";
					}
				}
			}
		}
	}
	set_index(x_best);
	print_index(std::cout);
}
