#include <Rcpp.h>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <limits>
using namespace Rcpp;
using u128 = unsigned __int128;
using i128 = __int128;
static constexpr uint64_t scale = uint64_t(1) << 63;
struct Interval { uint64_t low, up; };


static Interval from_double(double value) {
  if (!std::isfinite(value) || value < 0 || value > 1) stop("Probability outside [0,1]");
  if (value == 0) return {0, 0};
  uint64_t bits;
  std::memcpy(&bits, &value, sizeof(bits));
  const int exponent = int((bits >> 52) & 2047);
  if (exponent == 0) return {0, 1};
  const uint64_t mantissa = (uint64_t(1) << 52) | (bits & ((uint64_t(1) << 52) - 1));
  const int shift = exponent - 1012;
  if (shift >= 0) return {mantissa << shift, mantissa << shift};
  if (-shift >= 64) return {0, 1};
  const uint64_t denominator = uint64_t(1) << (-shift);
  return {mantissa / denominator, mantissa / denominator + (mantissa % denominator != 0)};
}
static uint64_t multiply_low(uint64_t a, uint64_t b) {
  return uint64_t((u128(a) * b) >> 63);
}
static uint64_t multiply_up(uint64_t a, uint64_t b) {
  return uint64_t((u128(a) * b + scale - 1) >> 63);
}
static Interval multiply(Interval a, Interval b) {
  return {multiply_low(a.low, b.low), multiply_up(a.up, b.up)};
}
static Interval mass_at(IntegerMatrix n, IntegerMatrix x, NumericVector coefficient,
                        int state, const std::vector<Interval>& probability) {
  Interval mass{scale, scale};
  for (int j = 0; j < n.ncol(); ++j) {
    if (x(state, j) < 0 || x(state, j) > n(state, j)) stop("Invalid counts");
    const Interval p = probability[j];
    const Interval q{scale - p.up, scale - p.low};
    for (int k = 0; k < x(state, j); ++k) mass = multiply(mass, p);
    for (int k = x(state, j); k < n(state, j); ++k) mass = multiply(mass, q);
  }
  const double c = coefficient[state];
  if (!std::isfinite(c) || c < 1 || c != std::floor(c) || c >= 9007199254740992.0)
    stop("Invalid path coefficient");
  const u128 lo = u128(mass.low) * uint64_t(c);
  const u128 up = u128(mass.up) * uint64_t(c);
  if (lo > scale) stop("Lower mass exceeds one");
  return {uint64_t(lo), uint64_t(std::min(up, u128(scale)))};
}
static double lower_double(i128 value, int exponent = 63) {
  return std::nextafter(std::ldexp(double(value), -exponent), -INFINITY);
}
static double upper_double(i128 value, int exponent = 63) {
  return std::nextafter(std::ldexp(double(value), -exponent), INFINITY);
}

// [[Rcpp::export]]
List point_mass_bounds(IntegerMatrix n, IntegerMatrix x, NumericVector coefficient,
                       NumericMatrix points) {
  if (x.nrow() != n.nrow() || x.ncol() != n.ncol() || points.ncol() != n.ncol() ||
      coefficient.size() != n.nrow()) stop("Dimension mismatch");
  NumericMatrix low(n.nrow(), points.nrow()), up(n.nrow(), points.nrow());
  for (int c = 0; c < points.nrow(); ++c) {
    std::vector<Interval> probability;
    for (int j = 0; j < points.ncol(); ++j) probability.push_back(from_double(points(c, j)));
    for (int s = 0; s < n.nrow(); ++s) {
      const Interval mass = mass_at(n, x, coefficient, s, probability);
      low(s, c) = std::max(0.0, lower_double(mass.low));
      up(s, c) = std::min(1.0, upper_double(mass.up));
    }
    if (c % 100 == 0) checkUserInterrupt();
  }
  return List::create(_["low"] = low, _["up"] = up);
}

// [[Rcpp::export]]
List clinical_score_bounds(IntegerMatrix n, IntegerMatrix x, NumericVector coefficient,
  NumericMatrix points, NumericMatrix weight_low, NumericMatrix weight_up,
  IntegerMatrix utility_quarters) {
  const int k = n.nrow(), j = n.ncol(), curves = points.nrow(), priors = weight_low.ncol();
  if (x.nrow() != k || x.ncol() != j || coefficient.size() != k || points.ncol() != j ||
      weight_low.nrow() != curves || weight_up.nrow() != curves ||
      weight_up.ncol() != priors || utility_quarters.nrow() != curves ||
      utility_quarters.ncol() != j + 1) stop("Dimension mismatch");
  const int size = k * (j + 1) * priors;
  std::vector<i128> low(size, 0), up(size, 0);
  for (int c = 0; c < curves; ++c) {
    std::vector<Interval> probability, weights;
    for (int a = 0; a < j; ++a) probability.push_back(from_double(points(c, a)));
    for (int p = 0; p < priors; ++p) {
      if (weight_low(c,p) > weight_up(c,p)) stop("Inverted weight interval");
      weights.push_back({from_double(weight_low(c,p)).low, from_double(weight_up(c,p)).up});
    }
    for (int s = 0; s < k; ++s) {
      const Interval mass = mass_at(n, x, coefficient, s, probability);
      for (int p = 0; p < priors; ++p) {
        const Interval weighted = multiply(mass, weights[p]);
        for (int a = 0; a <= j; ++a) {
          const int u = utility_quarters(c,a);
          if (u < -4 || u > 4) stop("Utility outside checked range");
          const int index = s + k * a + k * (j + 1) * p;
          low[index] += i128(u) * (u >= 0 ? weighted.low : weighted.up);
          up[index] += i128(u) * (u >= 0 ? weighted.up : weighted.low);
        }
      }
    }
    if (c % 100 == 0) checkUserInterrupt();
  }
  NumericVector lower(size), upper(size);
  for (int i = 0; i < size; ++i) {
    lower[i] = lower_double(low[i], 65);
    upper[i] = upper_double(up[i], 65);
  }
  lower.attr("dim") = IntegerVector::create(k, j + 1, priors);
  upper.attr("dim") = IntegerVector::create(k, j + 1, priors);
  return List::create(_["low"] = lower, _["up"] = upper);
}
