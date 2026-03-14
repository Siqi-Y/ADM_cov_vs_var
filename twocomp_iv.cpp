
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Two-compartment IV bolus: Cp(t) = A*exp(-alpha*t) + B*exp(-beta*t)
// [[Rcpp::export]]
NumericMatrix Cpp_twocomp_iv(NumericVector cl, NumericVector v1,
                             NumericVector v2, NumericVector q,
                             NumericVector ti) {
  int n_ind  = cl.size();
  int n_time = ti.size();
  NumericMatrix out(n_ind, n_time);

  for (int iter = 0; iter < n_ind; ++iter) {
    double k   = cl[iter] / v1[iter];
    double k12 = q[iter]  / v1[iter];
    double k21 = q[iter]  / v2[iter];

    double sum_rates = k + k12 + k21;
    double sqrt_disc = sqrt(sum_rates * sum_rates - 4.0 * k * k21);

    double alpha = 0.5 * (sum_rates + sqrt_disc);
    double beta  = 0.5 * (sum_rates - sqrt_disc);

    // Coefficients for IV bolus (dose/V1 normalized out)
    double A = (alpha - k21) / (alpha - beta) / v1[iter];
    double B = (k21 - beta)  / (alpha - beta) / v1[iter];

    for (int i = 0; i < n_time; ++i) {
      out(iter, i) = A * exp(-alpha * ti[i]) + B * exp(-beta * ti[i]);
    }
  }
  return out;
}

