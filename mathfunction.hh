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
#include <bitset>

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
  double simpson_rule(std::function<double(double)> func, double a, double b, uint32_t division);
  int64_t gcd(int64_t a, int64_t b);
  int64_t lcm(int64_t a, int64_t b);
  uint64_t llcombination(double n, double r);
  double combination(double n, double r);
  bool is_prime(uint64_t n);
  template <size_t SIEVE_SIZE>
  bool is_prime_soe(uint64_t n);
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

  /**
   * @brief エラトステネスの篩により与えられた整数が素数かどうか判定する
   *
   * エラトステネスの篩により与えられた整数が素数かどうか判定する。
   * 篩のサイズより大きな整数は is_prime により判定する。
   *
   * @tparam SIEVE_SIZE 篩のサイズ
   * @param n 素数かどうか判定する整数
   *
   * @return 与えられた整数が素数ならばtrue、素数でなければfalse
   */
  template <size_t SIEVE_SIZE>
  bool is_prime_soe(uint64_t n)
  {
    static std::bitset<SIEVE_SIZE + 1> prime_sieve;
    static bool first_call = true;

    if (n > SIEVE_SIZE) {
      return is_prime(n);
    }

    if (first_call) {
      first_call = false;
      for (size_t i = 2; i < prime_sieve.size(); ++i) {
        prime_sieve[i] = true;
      }

      for (size_t i = 2; i * i <= SIEVE_SIZE; ++i) {
        if (prime_sieve[i]) {
          for (size_t j = i * 2; j <= SIEVE_SIZE; j += i) {
            prime_sieve[j] = false;
          }
        }
      }
    }

    return prime_sieve[n];
  }
};


#endif
