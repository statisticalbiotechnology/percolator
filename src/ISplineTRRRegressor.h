#include "MonotoneRegressor.h"
#include <iostream> // for std::cerr, std::endl
#include <Eigen/Dense>
#include <limits>
#include <numeric>
#include <cassert>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Utility: project x to [lb, ub]
inline void projectToBox(VectorXd& x, const VectorXd& lb, const VectorXd& ub) {
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] < lb[i]) x[i] = lb[i];
    else if (x[i] > ub[i]) x[i] = ub[i];
  }
}


// Quadratic objective model:
//   f(x) = x^T G x - 2 h^T x + c
//   g(x) = 2 (G x - h)
struct QuadraticModel {
  MatrixXd G;
  VectorXd h;
  double c;
};

// Solves x'Gx - 2h'x = min subject to lb <= x <= ub for a symmetric
// positive definite G with a primal active-set method. Each iteration
// minimises exactly over the free variables (Cholesky), moves towards that
// minimum until a bound blocks, or, when nothing blocks, releases the bound
// variable whose gradient points furthest into the box. For a strictly
// convex problem this ends at the exact optimum after a finite number of
// iterations.
struct BoxQPResult {
  VectorXd x;
  double fval;
  int iters;
  bool converged;
};

inline BoxQPResult solveBoxQP(const QuadraticModel& q, const VectorXd& lb,
                              const VectorXd& ub, VectorXd x) {
  const int p = static_cast<int>(x.size());
  projectToBox(x, lb, ub);
  std::vector<char> atBound(p);
  for (int i = 0; i < p; ++i) atBound[i] = (x[i] == lb[i] || x[i] == ub[i]);
  const int maxIters = 10 * p + 100;
  int it = 0;
  bool converged = false;
  for (; it < maxIters; ++it) {
    std::vector<int> F;
    for (int i = 0; i < p; ++i) {
      if (!atBound[i]) F.push_back(i);
    }
    const int m = static_cast<int>(F.size());
    // minimum over the free variables, with the bound variables fixed
    MatrixXd GFF(m, m);
    VectorXd rhs(m);
    for (int a = 0; a < m; ++a) {
      double t = q.h[F[a]];
      for (int j = 0; j < p; ++j) {
        if (atBound[j]) t -= q.G(F[a], j) * x[j];
      }
      rhs[a] = t;
      for (int c = 0; c < m; ++c) GFF(a, c) = q.G(F[a], F[c]);
    }
    VectorXd xF(m);
    if (m > 0) {
      Eigen::LLT<MatrixXd> llt(GFF);
      if (llt.info() != Eigen::Success) break;
      xF = llt.solve(rhs);
    }
    // move towards it until the first bound blocks
    double alpha = 1.0;
    int blocking = -1;
    for (int a = 0; a < m; ++a) {
      const int i = F[a];
      const double step = xF[a] - x[i];
      if (step < 0.0 && xF[a] < lb[i]) {
        const double r = (lb[i] - x[i]) / step;
        if (r < alpha) { alpha = r; blocking = i; }
      } else if (step > 0.0 && xF[a] > ub[i]) {
        const double r = (ub[i] - x[i]) / step;
        if (r < alpha) { alpha = r; blocking = i; }
      }
    }
    for (int a = 0; a < m; ++a) {
      const int i = F[a];
      x[i] = (alpha == 1.0) ? xF[a] : x[i] + alpha * (xF[a] - x[i]);
    }
    if (blocking >= 0) {
      x[blocking] = (x[blocking] - lb[blocking] < ub[blocking] - x[blocking])
                        ? lb[blocking] : ub[blocking];
      atBound[blocking] = 1;
      continue;
    }
    // optimal on the free set; check the bound variables
    const VectorXd Gx = q.G * x;
    int release = -1;
    double mostInward = 1e-12;
    for (int i = 0; i < p; ++i) {
      if (!atBound[i]) continue;
      const double g = 2.0 * (Gx[i] - q.h[i]);
      const double inward = (x[i] == lb[i]) ? -g : g;
      if (inward > mostInward) {
        mostInward = inward;
        release = i;
      }
    }
    if (release < 0) {
      converged = true;
      break;
    }
    atBound[release] = 0;
  }
  const double fval = x.dot(q.G * x) - 2.0 * q.h.dot(x) + q.c;
  return {x, fval, it, converged};
}

namespace ispline_detail {

inline std::vector<double> normalize_to_unit_interval(const std::vector<double>& x) {
  if (x.empty()) return {};

  auto mm = std::minmax_element(x.begin(), x.end());
  const double xmin = *mm.first;
  const double xmax = *mm.second;
  const double span = xmax - xmin;
  if (!(span > 0.0)) {
    return std::vector<double>(x.size(), 0.5);
  }

  std::vector<double> out(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    out[i] = (x[i] - xmin) / span;
  }
  return out;
}

} // namespace ispline_detail

// Build I-spline design matrix from proper B-spline basis.
//
// I-spline columns are cumulative tail sums of B-splines:
//   I_j(x) = sum_{m >= j} B_m(x),  j = 1, ..., n_bspline-1
// With clamped (repeated) boundary knots, each I_j is a smooth monotone
// function from 0 to 1 spanning multiple knot intervals.  Non-negative
// coefficients on these columns guarantee a monotone non-decreasing fit.
//
// B-splines are evaluated via the de Boor triangular recursion (O(d^2)
// per data point, independent of the number of basis functions).
//
// The knot vector must be a clamped B-spline knot vector: degree+1
// copies of the left boundary, sorted internal knots, degree+1 copies
// of the right boundary.
inline Eigen::MatrixXd build_ispline_design(
    const std::vector<double>& scores,
    int degree,
    const std::vector<double>& knots,
    bool include_intercept)
{
  const int n = static_cast<int>(scores.size());
  const int n_knots = static_cast<int>(knots.size());
  const int n_bspline = n_knots - degree - 1;
  const int d = degree;

  if (n_bspline <= 1) {
    Eigen::MatrixXd X(n, 1);
    X.setOnes();
    return X;
  }

  const int n_ispline = n_bspline - 1;
  const int p = (include_intercept ? 1 : 0) + n_ispline;
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(n, p);

  int col0 = 0;
  if (include_intercept) {
    X.col(0).setOnes();
    col0 = 1;
  }

  const double* t = knots.data();

  for (int i = 0; i < n; ++i) {
    double x = scores[i];

    // Find knot span k: t[k] <= x < t[k+1], clamped to valid range [d, n_knots-d-2]
    int k;
    if (x >= t[n_knots - d - 1]) {
      k = n_knots - d - 2;
    } else if (x <= t[d]) {
      k = d;
    } else {
      int lo = d, hi = n_knots - d - 1;
      while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (t[mid + 1] <= x) lo = mid + 1;
        else hi = mid;
      }
      k = lo;
    }

    // De Boor triangular recursion: compute B_{k-d}(x), ..., B_k(x)
    static constexpr int kMaxDegree = 31;
    assert(d <= kMaxDegree && "ispline_degree must be <= 31");
    double B[kMaxDegree + 1];
    B[0] = 1.0;
    for (int j = 1; j <= d; ++j) {
      double saved = 0.0;
      for (int r = 0; r < j; ++r) {
        double left_val  = x - t[k + 1 - j + r];
        double right_val = t[k + 1 + r] - x;
        double denom = left_val + right_val;
        double temp = (denom > 0.0) ? B[r] / denom : 0.0;
        B[r] = saved + right_val * temp;
        saved = left_val * temp;
      }
      B[j] = saved;
    }
    // B[r] == B_{k-d+r, d}(x)  for r = 0, ..., d

    // Cumulative sum from right: cum[r] = sum_{s=r}^{d} B[s]
    double cum[kMaxDegree + 1];
    cum[d] = B[d];
    for (int r = d - 1; r >= 0; --r) {
      cum[r] = cum[r + 1] + B[r];
    }

    // Fill I-spline columns.
    // I_j(x) = sum_{m >= j} B_m(x)  for j = 1 .. n_bspline-1
    //   = 1           if j <= k-d    (all non-zero B-splines included)
    //   = cum[j-(k-d)] if k-d < j <= k (partial sum of non-zero B-splines)
    //   = 0           if j > k        (no non-zero B-splines included)
    for (int j = 1; j < n_bspline; ++j) {
      double val;
      if (j <= k - d) {
        val = 1.0;
      } else if (j > k) {
        val = 0.0;
      } else {
        val = cum[j - (k - d)];
      }
      X(i, col0 + j - 1) = val;
    }
  }

  return X;
}


// Build a clamped B-spline knot vector with degree+1 repeated boundary
// knots and internal knots placed at data quantiles.
static std::vector<double> make_default_knots(const std::vector<double>& scores, int degree) {
  const int n = (int)scores.size();
  if (n == 0) return {};

  std::vector<double> sorted = scores;
  std::sort(sorted.begin(), sorted.end());

  double lo = sorted.front();
  double hi = sorted.back();
  if (hi <= lo) hi = lo + 1e-6;

  int num_internal = std::min(200, (int)std::sqrt(n));

  // Place internal knots at quantiles, skipping duplicates and boundaries
  std::vector<double> internal;
  internal.reserve(num_internal);
  double prev = lo;
  for (int k = 1; k <= num_internal; ++k) {
    double q = (double)k / (num_internal + 1);
    size_t idx = std::min((size_t)(q * (n - 1)), (size_t)(n - 1));
    double val = sorted[idx];
    if (val > lo + 1e-12 && val < hi - 1e-12 && val > prev + 1e-12) {
      internal.push_back(val);
      prev = val;
    }
  }

  // Clamped knot vector: (degree+1) copies of lo, internal knots, (degree+1) copies of hi
  const int order = degree + 1;
  std::vector<double> knots;
  knots.reserve(2 * order + (int)internal.size());
  for (int i = 0; i < order; ++i) knots.push_back(lo);
  for (double v : internal) knots.push_back(v);
  for (int i = 0; i < order; ++i) knots.push_back(hi);

  return knots;
}

class ISplineTRRRegressor final : public MonotoneRegressor {
public:
  using MonotoneRegressor::MonotoneRegressor;

  std::vector<double> fit_xy(const std::vector<double>& x,
                             const std::vector<double>& y,
                             double clip_lo = 0.0, double clip_hi = 1.0) override {
    assert(x.size() == y.size());
    const int n = (int)x.size();

    std::vector<double> x_work;
    const std::vector<double>* x_ptr = &x;
    if (params_.y_decreasing_in_x) {
      x_work = x;
      for (auto& v : x_work) v = -v;
      x_ptr = &x_work;
    }
    const std::vector<double> x_norm = ispline_detail::normalize_to_unit_interval(*x_ptr);

    auto knots = params_.knots;
    if (knots.size() <= 0) {
      knots = make_default_knots(x_norm, params_.ispline_degree);
    }
    // Build X
    Eigen::MatrixXd X = build_ispline_design(
        x_norm, params_.ispline_degree, knots, params_.include_intercept);
    const int p = (int)X.cols();

    // Precompute normal-equation terms for objective:
    // (1/n)||X beta - y||^2 + lambda_ridge ||beta||^2 + lambda_smooth ||D beta_spline||^2
    // where D is the second-difference matrix on the I-spline coefficients.
    const Eigen::Map<const Eigen::VectorXd> y_vec(y.data(), n);
    const double inv_n = 1.0 / std::max(1, n);
    QuadraticModel q;
    q.G = inv_n * (X.transpose() * X);
    q.h = inv_n * (X.transpose() * y_vec);
    q.c = inv_n * y_vec.squaredNorm();
    const double lambda = std::max(0.0, params_.ridge_lambda);
    if (lambda > 0.0) {
      q.G.diagonal().array() += lambda;
    }

    // Second-difference smoothing penalty on I-spline coefficients.
    // Penalizes curvature (gamma_{j+2} - 2*gamma_{j+1} + gamma_j)^2,
    // encouraging smooth transitions instead of step functions.
    const int col0_s = params_.include_intercept ? 1 : 0;
    const int n_ispline = p - col0_s;
    const double smooth_lambda = std::max(0.0, params_.smooth_lambda);
    if (smooth_lambda > 0.0 && n_ispline > 2) {
      const int nd = n_ispline - 2;
      Eigen::MatrixXd D = Eigen::MatrixXd::Zero(nd, p);
      for (int di = 0; di < nd; ++di) {
        D(di, col0_s + di)     =  1.0;
        D(di, col0_s + di + 1) = -2.0;
        D(di, col0_s + di + 2) =  1.0;
      }
      q.G += smooth_lambda * (D.transpose() * D);
    }

    // Bounds: non-negative except possibly intercept
    Eigen::VectorXd lb = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd ub = Eigen::VectorXd::Constant(p, std::numeric_limits<double>::infinity());
    if (params_.intercept_col >= 0 && params_.intercept_col < p)
      lb[params_.intercept_col] = -std::numeric_limits<double>::infinity();

    BoxQPResult sol = solveBoxQP(q, lb, ub, Eigen::VectorXd::Zero(p));
    // Predict + clamp
    std::vector<double> yhat(n);
    for (int i = 0; i < n; ++i) {
      double v = X.row(i).dot(sol.x);
      if (v < clip_lo) v = clip_lo;
      if (v > clip_hi) v = clip_hi;
      yhat[i] = v;
    }
    return yhat;
  }
  
  std::vector<double> fit_y(const std::vector<double>& y,
                            double clip_lo = 0.0, double clip_hi = 1.0) override {
    // No x provided: use rank as a x-value.
    std::vector<double> x(y.size());
    std::iota(x.rbegin(), x.rend(), 1.0);  // x[0]=N, ..., x[N-1]=1
    return fit_xy(x, y, clip_lo, clip_hi);
  }

};
