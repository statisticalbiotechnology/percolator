#pragma once
#include <vector>
#include <memory>

struct MonotoneParams {
  // Shared
  double clip_lo = 0.0;
  double clip_hi = 1.0;
  // I-spline TRR specific
  double ridge_lambda = 1e-3;
  int ispline_degree = 3;
  bool include_intercept = false;
  int intercept_col = -1; // if include_intercept, usually 0
  std::vector<double> knots; // knot policy fills this
};

class MonotoneRegressor {
public:
  virtual ~MonotoneRegressor() = default;

  // Fit y ~ f(x) with optional weights and output range clamp
  virtual std::vector<double> fit_xy(const std::vector<double>& x,
                                     const std::vector<double>& y,
                                     double clip_lo = 0.0,
                                     double clip_hi = 1.0) = 0;

  // Fit y ~ f(index) (no x); useful for the tdc path when no scores
  virtual std::vector<double> fit_y(const std::vector<double>& y,
                                    double clip_lo = 0.0,
                                    double clip_hi = 1.0) = 0;

  virtual void set_params(const MonotoneParams& p) { params_ = p; }
  virtual MonotoneParams params() const { return params_; }

protected:
  MonotoneParams params_;
};

enum class RegressorType {
  PAVA,          // isotonic/PAVA implementation
  ISPLINE_TRR    // trust-region reflective constrained ridge I-spline
};

std::unique_ptr<MonotoneRegressor> make_monotone_regressor(RegressorType type, const MonotoneParams& params);