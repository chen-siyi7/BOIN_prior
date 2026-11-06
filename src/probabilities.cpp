#include <Rcpp.h>
#include <cmath>
#include <limits>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix terminal_mass_cpp(IntegerMatrix n, IntegerMatrix x,
                                NumericVector coefficient, NumericMatrix p) {
  int k = n.nrow(), j = n.ncol(), curves = p.nrow();
  if (x.nrow() != k || x.ncol() != j || p.ncol() != j ||
      coefficient.size() != k) stop("Incompatible terminal likelihood inputs.");
  NumericMatrix mass(k, curves);
  for (int c = 0; c < curves; ++c) {
    if (c % 32 == 0) checkUserInterrupt();
    for (int s = 0; s < k; ++s) {
      double value = std::log(coefficient[s]);
      for (int d = 0; d < j; ++d) {
        double prob = p(c,d);
        if (!std::isfinite(prob) || prob < 0 || prob > 1)
          stop("A toxicity probability is outside [0,1].");
        if (x(s,d)) value += x(s,d) * std::log(prob);
        if (n(s,d) > x(s,d)) value += (n(s,d)-x(s,d)) * std::log1p(-prob);
      }
      mass(s,c) = std::exp(value);
    }
  }
  return mass;
}

// [[Rcpp::export]]
NumericVector adjacent_double(NumericVector value, bool upward = true) {
  NumericVector out(value.size());
  double direction = upward ? R_PosInf : R_NegInf;
  for (R_xlen_t i=0; i<value.size(); ++i)
    out[i] = std::nextafter(value[i], direction);
  return out;
}
