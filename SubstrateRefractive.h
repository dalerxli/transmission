#ifndef SUBSTRATEREFRACTIVE_H_
#define SUBSTRATEREFRACTIVE_H_

#include <vector>

class SubstrateRefractive {
public:
	std::vector<double> air(const std::vector<double>& wavelength) {
		std::vector<double> vec(wavelength.size(), 1);
		return vec;
	}

	std::vector<double> glass(const std::vector<double>& wavelength) {
		std::vector<double> vec(wavelength.size());
		for (size_t i = 0; i < wavelength.size(); ++i) {
			vec[i] = std::sqrt(1.0 + 1.0 / (0.7568 - 7930.0 / (wavelength[i] * wavelength[i])));
		}
		return vec;
	}

	std::vector<double> sodaglass(const std::vector<double>& wavelength) {
		std::vector<double> vec(wavelength.size());
		for (size_t i = 0; i < wavelength.size(); ++i) {
			vec[i] = 1.5130-0.003169*1e-6*wavelength[i]*wavelength[i] + \
				0.003962 / (1e-6*wavelength[i] * wavelength[i]);
		}
		return vec;
	}

	std::vector<double> si(const std::vector<double>& wavelength) {
		std::vector<double> vec(wavelength.size());
		for (size_t i = 0; i < wavelength.size(); ++i) {
			vec[i] = 3.71382 - 8.69123*1e-5*wavelength[i] - 2.47125*1e-8*wavelength[i] * wavelength[i] \
				+ 1.04677*1e-11*wavelength[i] * wavelength[i] * wavelength[i];
		}
		return vec;
	}

	std::vector<double> AQuartz(const std::vector<double>& wavelength) {
		std::vector<double> vec(wavelength.size());
		for (size_t i = 0; i < wavelength.size(); ++i) {
			vec[i] = std::sqrt(1.0 + 0.18394 * wavelength[i] * wavelength[i] /
				(wavelength[i] * wavelength[i] - 18231.83));
		}
		return vec;
	}

	std::vector<double> CQuartz(const std::vector<double>& wavelength) {
		std::vector<double> vec(wavelength.size());
		for (size_t i = 0; i < wavelength.size(); ++i) {
			vec[i] = std::sqrt(1.0 + 1.34157 * wavelength[i] * wavelength[i] /
				(wavelength[i] * wavelength[i] - 8385.79));
		}
		return vec;
	}
	


	std::vector<double> silicon(const std::vector<double>& wavelength) {
		std::vector<double> vec(wavelength.size());
		for (size_t i = 0; i < wavelength.size(); ++i) {
			vec[i] = std::sqrt(11.67316 + 1.0/(1e-6*wavelength[i]*wavelength[i]) + \
				0.00448263 / (1e-6*wavelength[i] * wavelength[i] - 1.228118322025));
		}
		return vec;
	}
};

#endif // SUBSTRATEREFRACTIVE_H_