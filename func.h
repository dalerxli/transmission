#ifndef FUNC_H_
#define FUNC_H_

#include <vector>

template<class X>
X norm(const std::vector<X>& x) {
	X ans = X();
	for (size_t i = 0; i < x.size(); ++i) {
		ans += x[i] * x[i];
	}
	return std::sqrt(ans);
}

template<class X>
X scalar(const std::vector<X>& x, const std::vector<X>& y) {
	X ans = X();
	for (size_t i = 0; i < x.size(); ++i) {
		ans += x[i] * y[i];
	}
	return ans;
}

template<class X>
double mRandNorm() {
	X res = 0;
	while (res == 0 || res == 1) {
		res = ((X)rand()) / (X)(RAND_MAX);
	}
	return res;
}

template<class X>
void print_vec(const std::vector<X>& x) {
	for (size_t i = 0; i < x.size() / 2; ++i) {
		std::cout << x[i] << ' ' << x[i + x.size() / 2] << '\n';
	}

}
#endif //FUNC_H_