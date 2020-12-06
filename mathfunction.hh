#ifndef _mathfunction_
#define _mathfunction_

// --> This library provide mathematic function which are not
// --> included in the standard library <cmath>
// --> Some functions already implemented in c++17 are included.

#include <cmath>
#include <limits>
#include <iostream>
#include <utility>
#include <functional>
#include <vector>

namespace mathconstant {
  const double pi = 3.14159265358979323;
  const double sqrt2 = std::sqrt(2);
  const double sqrt2pi = std::sqrt(2 * pi);
}

namespace mathfunc {
  // --- 数値解析用 --- //
  double power(double x, int N);
  double differential(std::function<double(double)> func, double x, double h = 0.001);
  double differential2(std::function<double(double)> func, double x, double h = 0.001);
  double error_propagation(std::function<double(double*)> func, double* x, double* x_e, const int num_arg, double h = 0.001); // multi variable function
 // ---- 数値計算アルゴリズム ---- //
  double newton_method(std::function<double(double)> func, double init, double epsilon = 1e-12);
  double find_extremum_x(std::function<double(double)> func, double init, double epsilon = 1e-12);
  int64_t gcd(int64_t a, int64_t b);
  int64_t lcm(int64_t a, int64_t b);
  // ---- 特殊関数 ---- //
  double lower_incomp_gamma(double a, double x);
  double normalized_lower_incomp_gamma(double a, double x);
  double incomp_beta(double x, double a, double b);
  double beta(double a, double b);
  // ---- 統計学用の関数 ---- //
  double chisquared_pdf(double x, double deg);
  double chisquared_cdf(double x, double deg);
  double poisson_pdf(unsigned int k, double lambda);
  double chisquared_lower_limit(double alpha, double deg, double init = 0);
  double binomial_pdf(unsigned int n, unsigned int k, double p);
  double binomial_cdf(unsigned int n, unsigned int k, double p);
  double normal_pdf(double x, double mu = 0, double sigma = 1);
  double normal_cdf(double x, double mu = 0, double sigma = 1);
  double geometric_pdf(int k, double p);
  double geometric_cdf(int k, double p);
};


#endif
