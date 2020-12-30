#include "mathfunction.hh"

double mathfunc::power(double x, int N) {
  double result = 1;
  for (int i = 0; i < N; ++i) {
    result *= x;
  }
  return result;
}

double mathfunc::differential(std::function<double(double)> func, double x, double h) {
  return (func(x - 2 * h) - 8 * func(x - h) + 8 * func(x + h) - func(x + 2 * h)) / (12 * h);
}

double mathfunc::differential2(std::function<double(double)> func, double x, double h) {
  // return ( func(x - h) - 2 * func(x) + func(x + h) ) / mathfunc::power(h, 2);
  return (16 * (func(x + h) + func(x - h)) - (func(x + 2 * h) + func(x - 2 * h)) - 30 * func(x)) / (12 * h * h);
}

double mathfunc::newton_method(std::function<double(double)> func, double init, double epsilon) {
  double x_old = init;
  for (int count = 0; count < 100; ++count) {
    double y_now = func(x_old);
    double y_dif = mathfunc::differential(func, x_old);
    double x_new = x_old - y_now / y_dif;
    if (std::fabs(x_new - x_old) < epsilon) {
      return x_new;
    } else {
      x_old = x_new;
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

double mathfunc::find_extremum_x(std::function<double(double)> func, double init, double epsilon) {
  double x_old = init;
  for (int count = 0; count < 100; ++count) {
    double y_now = mathfunc::differential(func, x_old);
    double y_dif = mathfunc::differential2(func, x_old);
    double x_new = x_old - y_now / y_dif;
    if (std::fabs(x_new - x_old) < epsilon) {
      return x_new;
    } else {
      x_old = x_new;
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

double mathfunc::simpson_rule(std::function<double(double)> func, double a, double b, uint32_t division) {
  if (division <= 0) return 0;
  if (division % 2 != 0) ++division;
  double h = (b - a) / division;
  double result = 0;
  result += func(a) + func(b) + 4 * func(b - h);
  for (auto i = 1; i <= division - 2; i += 2) {
    result += 4 * func(a + i * h);
    result += 2 * func(a + (i + 1) * h);
  }
  return result * h / 3;
}

int64_t mathfunc::gcd(int64_t a, int64_t b) {
  if (a < b) return gcd(b, a);
  while (b != 0) {
    long int R = a % b;
    a = b; b = R;
  }
  return a;
}

int64_t mathfunc::lcm(int64_t a, int64_t b) {
  int64_t gcd = mathfunc::gcd(a, b);
  return (a / gcd) * b;
}

bool mathfunc::is_prime(uint64_t n) {
  if (n == 0 || n == 1) return false;
  if (n == 2 || n == 3) return true;
  if (n % 2 == 0 || n % 3 == 0) return false;
  for (uint64_t i = 5; i * i <= n; i += 6) {
    if (n % i == 0) return false;
    if (n % (i + 2) == 0) return false;
  }
  return true;
}

double mathfunc::error_propagation(std::function<double(double*)> func, double* x, double* x_e, const int num_arg, double h) {
  std::vector<double> dif;
  dif.reserve(num_arg);
  for (int i = 0; i < num_arg; ++i) {
    x[i] += h;
    double dif_1p = func(x);
    x[i] += h;
    double dif_2p = func(x);
    x[i] -= 3 * h;
    double dif_1n = func(x);
    x[i] -= h;
    double dif_2n = func(x);
    x[i] += 2 * h;
    dif.push_back((dif_2n - 8 * dif_1n + 8 * dif_1p - dif_2p) / (12 * h));
  }
  double result = 0;
  for (int i = 0; i < num_arg; ++i) {
    result += mathfunc::power(dif[i] * x_e[i], 2);
  }
  return std::sqrt(result);
}

double mathfunc::combination(double n, double r) {
  double ln_result = std::lgamma(n + 1) - std::lgamma(r + 1) - std::lgamma(n - r + 1);
  return std::exp(ln_result);
}

uint64_t mathfunc::llcombination(double n, double r) {
  return std::llround(mathfunc::combination(n, r));
}

double mathfunc::chisquared_pdf(double x, double deg) {
  if (deg < 250) {
    return 0.5 * std::pow(x / 2, deg / 2 - 1) / std::tgamma(deg / 2) * std::exp(-x / 2);
  }
  double log_chisquared = (deg / 2 - 1) * std::log(x) - x / 2 - deg / 2 * std::log(2) - std::lgamma(deg / 2);
  return std::exp(log_chisquared);

}

double mathfunc::chisquared_cdf(double x, double deg) {
  return normalized_lower_incomp_gamma(deg / 2, x / 2);
}

double mathfunc::poisson_pdf(unsigned int k, double lambda) {
  if (k > 0 && lambda > 0)
    return std::exp(k * std::log(lambda) - lambda - std::lgamma(k + 1));
  if (lambda >= 0)
    return std::exp(-lambda);

  return std::numeric_limits<double>::quiet_NaN();
}

double mathfunc::binomial_pdf(unsigned int n, unsigned int k, double p) {
  if (n >= k && 0 < p && p < 1) {
    double log_prob = std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1) + k * std::log(p) + (n - k) * std::log(1 - p);
    return std::exp(log_prob);
  }
  if (p == 1 && n == k) return 1;
  if (p == 1 && n != k) return 0;
  if (p == 0 && k == 0) return 1;
  if (p == 0 && k != 0) return 0;
  return std::numeric_limits<double>::quiet_NaN();
}

double mathfunc::binomial_cdf(unsigned int n, unsigned int k, double p) {
  if (n == k) return 1.0;
  if (n < k) return 0.0;
  double result = mathfunc::incomp_beta(1 - p, n - k, k + 1);
  if (result > 1) {
    return 1;
  } else {
    return result;
  }

}

double mathfunc::normal_pdf(double x, double mu, double sigma) {
  return std::exp(-0.5 * std::pow((x - mu) / sigma, 2)) / (mathconstant::sqrt2pi * sigma);
}

double mathfunc::normal_cdf(double x, double mu, double sigma) {
  return 0.5 * (1 + std::erf((x - mu) / (mathconstant::sqrt2 * sigma)));
}

double mathfunc::lower_incomp_gamma(double a, double x) {
  double result = 0;
  double term = 1;
  int k = 0;
  while (1) {
    if (k == 0) {
      term *= 1. / a;
    } else {
      term *= x / (a + k);
    }
    result += term;

    if (term < 1e-15) break;
    k++;
  }// <k>
  return result * std::pow(x, a) * std::exp(-x);
}

double mathfunc::normalized_lower_incomp_gamma(double a, double x) {
  double term = 1;
  double series_sum = 0;
  for (int k = 0; k < 4095; ++k) {
    if (k == 0) {
      term *= 1. / a;
    } else {
      term *= x / (a + k);
    }
    series_sum += term;
    if (series_sum >= std::numeric_limits<double>::infinity()) return 1.0;
    if (term < 1e-15) break;
  }
  double log_result = a * std::log(x) - x - std::lgamma(a) + std::log(series_sum);
  return std::exp(log_result);
}

double mathfunc::chisquared_lower_limit(double alpha, double deg, double init) {
  double x_old;
  double epsilon = 1e-12;
  if (deg < 0.5) {
    x_old = 1e-80;
    epsilon = 1e-12;
  } else if (deg <= 2) {
    x_old = 1e-15;
  } else if (deg > 2 && alpha < 0.1) {
    x_old = deg - 2 - 1e-6;
  } else if (deg > 2 && alpha >= 0.1) {
    x_old = deg - 2 + 1e-6;
  } else if (deg > 250 && alpha < 0.5) {
    x_old = deg - std::sqrt(deg);
  } else if (deg > 250 && alpha >= 0.5) {
    x_old = deg + std::sqrt(deg);
  }

  if (init != 0) x_old = init;

  double count = 0;
  while (1) {
    double x_new = x_old - (mathfunc::chisquared_cdf(x_old, deg) - alpha) / mathfunc::chisquared_pdf(x_old, deg);
    if (std::fabs(1 - x_new / x_old) < epsilon) {
      return x_new;
    } else {
      x_old = x_new;
    }
    count++;
    if (count > 10000) {
      std::cout << "Error : Lower value Not Converged : [alpha, deg] = "
        << alpha << " " << deg << std::endl;
      break;
    }
  }

  return std::numeric_limits<double>::quiet_NaN();
}

double mathfunc::beta(double a, double b) {
  return std::exp(std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b));
}

double mathfunc::incomp_beta(double x, double a, double b) {
  bool flip = false;
  if (x > (a + 1) / (a + b + 1)) {
    x = 1 - x;
    std::swap(a, b);
    flip = true;
  }
  double y = 1 - x;
  auto alpha_ = [&](int m) {
    if (m == 1) return 1.0;
    --m;
    double numerator = (a + m - 1) * (a + b + m - 1) * m * (b - m) * x * x;
    double denominator = mathfunc::power(a + 2 * m - 1, 2);
    return numerator / denominator;
  };
  auto beta_ = [&](int m) {
    if (m == 1) return a * (a - (a + b) * x + 1) / (a + 1);
    --m;
    double term1 = m;
    double term2 = m * (b - m) * x / (a + 2 * m - 1);
    double term3 = (a + m) * (a - (a + b) * x + 1 + m * (2 - x)) / (a + 2 * m + 1);
    return term1 + term2 + term3;
  };

#if __cplusplus >= 201402L
  auto series = [=](double m = 1) {
    auto f = [=](auto& self, double m) {
      if (m == 20)
        return beta_(20);
      return alpha_(m) / (beta_(m) + self(self, m + 1));
    };
    return f(f, m);
  };
#else
  std::function<double(double)> series = [&](double m = 1) {
    if (m == 20) return beta_(20);
    double result = alpha_(m) / (beta_(m) + series(m + 1));
    return result;
  };
#endif

  if (flip == false) {
    return std::pow(x, a) * std::pow(y, b) * series(1) / mathfunc::beta(a, b);
  } else {
    return 1.0 - std::pow(x, a) * std::pow(y, b) * series(1) / mathfunc::beta(a, b);
  }

}

double mathfunc::geometric_pdf(int k, double p) {
  if (k < 1) return 0;
  return std::pow(1 - p, k - 1) * p;
}

double mathfunc::geometric_cdf(int k, double p) {
  if (k < 1) return 0;
  return 1 - std::pow(1 - p, k);
}
