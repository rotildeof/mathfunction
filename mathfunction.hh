#ifndef _mathfunction_
#define _mathfunction_

// --> This libraty provide mathematic function which are not
// --> included in the standard library <cmath>

#include <cmath>
#include <limits>
#include <iostream>
#include <utility>
#include <functional>

namespace mathconstant{
  const double pi = 3.14159265358979323;
  const double sqrt2 = std::sqrt(2);
  const double sqrt2pi = std::sqrt( 2 * pi );
}

namespace mathfunc{
  double power(double x, int N);
  double differential(std::function<double(double)> func, double x, double h = 0.001);
  double differential2(std::function<double(double)> func, double x, double h = 0.001);
  double lower_incomp_gamma(double a, double x);
  double normalized_lower_incomp_gamma(double a, double x);
  double incomp_beta(double x, double a, double b);
  double beta(double a, double b);

  double chisquared_pdf(double x, double deg);
  double chisquared_cdf(double x, double deg);
  double poisson_pdf(unsigned int k, double lambda);
  double chisquared_lower_limit(double alpha, double deg, double init = 0);
  double binomial_pdf(unsigned int n, unsigned int k, double p);
  double binomial_cdf(unsigned int n, unsigned int k, double p);
  double normal_pdf(double x, double mu = 0, double sigma = 1);
  double normal_cdf(double x, double mu = 0, double sigma = 1);

};

double mathfunc::power(double x, int N){
  double result = 1;
  for(int i = 0 ; i < N ; ++i){
    result *= x;
  }
  return result;
}

double mathfunc::differential(std::function<double(double)> func, double x, double h){
  return (func( x -  2 * h ) - 8 * func(x - h) + 8 * func(x + h) - func(x + 2 * h) ) / (12 * h);
}

double mathfunc::differential2(std::function<double(double)> func, double x, double h){
  return ( func(x - h) - 2 * func(x) + func(x + h) ) / mathfunc::power(h, 2);
}

double mathfunc::chisquared_pdf(double x, double deg){
  if(deg < 250){
    return 0.5 * std::pow(x / 2, deg / 2 - 1) / std::tgamma(deg / 2) * std::exp(- x / 2);
  }
  double log_chisquared = (deg / 2 - 1) * std::log(x) - x / 2 - deg / 2 * std::log(2) - std::lgamma(deg / 2);
  return std::exp(log_chisquared);
  
}

double mathfunc::chisquared_cdf(double x, double deg){
  return normalized_lower_incomp_gamma(deg / 2, x / 2);
}

double mathfunc::poisson_pdf(unsigned int k, double lambda){
  if( k > 0 && lambda > 0)
    return std::exp( k * std::log(lambda) - lambda - std::lgamma(k + 1) );
  if( lambda >= 0)
    return std::exp(-lambda);

  return std::numeric_limits<double>::quiet_NaN();
}

double mathfunc::binomial_pdf(unsigned int n, unsigned int k, double p){
  if(n >= k && 0 < p && p < 1){
    double log_prob = std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1) + k * std::log(p) + (n - k) * std::log(1 - p);
    return std::exp(log_prob);
  }
  if( p == 1 && n == k) return 1;
  if( p == 1 && n != k) return 0;
  if( p == 0 && k == 0) return 1;
  if( p == 0 && k != 0) return 0;
  return std::numeric_limits<double>::quiet_NaN();
}

double mathfunc::binomial_cdf(unsigned int n, unsigned int k, double p){
  if(n == k) return 1.0;
  if(n < k) return 0.0;
  double result = mathfunc::incomp_beta(1 - p, n - k, k + 1);
  if(result > 1){
    return 1;
  }else{
    return result;
  }

}

double mathfunc::normal_pdf(double x, double mu, double sigma){
  return std::exp( -0.5 * std::pow( (x - mu) / sigma , 2 ) ) / (mathconstant::sqrt2pi * sigma) ;
}

double mathfunc::normal_cdf(double x, double mu, double sigma){
  return 0.5 * (1 + std::erf( (x - mu) / (mathconstant::sqrt2 * sigma) ) );
}

double mathfunc::lower_incomp_gamma(double a, double x){
  double result = 0;
  double term = 1;
  int k = 0;
  while(1){
    if(k == 0){
      term *= 1./ a;
    }else{
      term *= x / (a + k);
    }
    result += term;

    if(term < 1e-15) break;
    k++;
  }// <k>
  return result * std::pow(x, a) * std::exp( - x );
}

double mathfunc::normalized_lower_incomp_gamma(double a, double x){
  double term = 1;
  double series_sum = 0;
  int k = 0;
  while(1){
    if(k == 0){
      term *= 1. / a;
    }else{
      term *= x / (a + k);
    }
    series_sum += term;
    if(series_sum >= std::numeric_limits<double>::infinity() ) return 1.0;
    if(term < 1e-15) break;

    k++;
    if(k > 4095) break;
  }
  double log_result = a * std::log(x) - x - std::lgamma(a) + std::log(series_sum);
  return std::exp(log_result);
  
}

double mathfunc::chisquared_lower_limit(double alpha, double deg, double init){
  double x_old;
  double epsilon = 1e-12;
  if(deg < 0.5){
    x_old = 1e-80;
    epsilon = 1e-12;
  }else if(deg <= 2){
    x_old = 1e-15;
  }else if(deg > 2 && alpha < 0.1){
    x_old = deg - 2 - 1e-6;
  }else if(deg > 2 && alpha >= 0.1){
    x_old = deg - 2 + 1e-6;
  }else if(deg > 250 && alpha < 0.5){
    x_old = deg - std::sqrt(deg);
  }else if(deg > 250 && alpha >= 0.5){
    x_old = deg + std::sqrt(deg);
  }
  
  if(init != 0) x_old = init;
  
  double count = 0;
  while(1){
    double x_new = x_old - (mathfunc::chisquared_cdf(x_old, deg) - alpha ) / mathfunc::chisquared_pdf(x_old, deg);
    if(std::fabs(1 - x_new / x_old) < epsilon){
      return x_new;
    }else{
      x_old = x_new;
    }
    count++;
    if(count > 10000){
      std::cout << "Error : Lower value Not Converged : [alpha, deg] = "
		<< alpha << " " << deg << std::endl;
      break;
    }
  }

  return std::numeric_limits<double>::quiet_NaN();
}

double mathfunc::beta(double a, double b){
  return std::exp(std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b));
}

double mathfunc::incomp_beta(double x, double a, double b){
  bool flip = false;
  if(x > (a + 1) / (a + b + 1)){
    x = 1 - x;
    std::swap(a, b);
    flip = true;
  }
  double y = 1 - x;
  auto alpha_ = [&](int m){
		  if(m == 1) return 1.0;
		  --m;
		  double numerator = (a + m - 1) * (a + b + m - 1) * m * (b - m) * x * x;
		  double denominator = mathfunc::power(a + 2 * m - 1, 2);
		  return numerator / denominator;
		};
  auto beta_ = [&](int m){
		 if(m == 1) return a * (a - (a + b) * x + 1) / (a + 1);
		 --m;
		 double term1 = m;
		 double term2 = m * (b - m) * x / (a + 2 * m - 1);
		 double term3 = (a + m) * (a - (a + b) * x + 1 + m * (2 - x)) / (a + 2 * m + 1);
		 return term1 + term2 + term3;
	       };
  std::function<double(double)> series = [&](double m = 1){
					   if(m == 20) return beta_(20);
					   double result = alpha_(m) / (beta_(m) + series(m+1) );
					   return result;
					 };
  if(flip == false){
    return std::pow(x, a) * std::pow(y, b) * series(1) / mathfunc::beta(a, b);
  }else{
    return 1.0 - std::pow(x, a) * std::pow(y, b) * series(1) / mathfunc::beta(a, b);
  }
  
}


#endif
