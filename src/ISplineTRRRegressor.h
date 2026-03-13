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

struct CostGrad {
  double f;
  VectorXd g;
};

inline CostGrad costGrad(const QuadraticModel& q, const VectorXd& x) {
  const VectorXd Gx = q.G * x;
  CostGrad cg;
  cg.f = x.dot(Gx) - 2.0 * q.h.dot(x) + q.c;
  cg.g = 2.0 * (Gx - q.h);
  return cg;
}

// Identify the "free" set: variables not at a bound OR where the gradient points interior
inline std::vector<int> freeSet(const VectorXd& x, const VectorXd& g,
                                const VectorXd& lb, const VectorXd& ub, double tol = 1e-12) {
  std::vector<int> F;
  F.reserve(x.size());
  for (int i = 0; i < x.size(); ++i) {
    bool atL = std::abs(x[i] - lb[i]) <= tol, atU = std::abs(x[i] - ub[i]) <= tol;
    bool moveIn = (!atL && !atU) || (atL && g[i] < 0.0) || (atU && g[i] > 0.0);
    if (moveIn) F.push_back(i);
  }
  return F;
}

// Steihaug–Toint truncated CG: solve (H)p = -g on free vars, with ||p|| <= Delta
struct TRStep { VectorXd p; double pred_red; bool hit_boundary; };

inline TRStep steihaugCG(const MatrixXd& H,
                  const VectorXd& x,
                  const VectorXd& g,
                  const std::vector<int>& F,
                  double Delta) {
  // Work only on free coordinates via views.
  auto Hmul = [&](const VectorXd& vF) {
    // Build a full vector v with zeros outside F, apply H, then extract F.
    VectorXd v = VectorXd::Zero(x.size());
    for (int k = 0; k < (int)F.size(); ++k) v[F[k]] = vF[k];
    VectorXd Hv = H * v;
    VectorXd out(vF.size());
    for (int k = 0; k < (int)F.size(); ++k) out[k] = Hv[F[k]];
    return out;
  };
  // Gradient on free set: gF
  VectorXd gF(F.size());
  for (int k = 0; k < (int)F.size(); ++k) gF[k] = g[F[k]];

  // CG init
  VectorXd pF = VectorXd::Zero(F.size());
  VectorXd rCG = -gF;                  // residual
  VectorXd d = rCG;                    // search direction
  double rr = rCG.squaredNorm();
  bool hit_boundary = false;

  const double eps = 1e-12;
  while (rr > eps) {
    VectorXd Hd = Hmul(d);
    double dHd = d.dot(Hd);
    if (dHd <= 0) {
      // Non-positive curvature, step to boundary
      // Solve ||pF + tau d|| = Delta  for tau > 0
      double a = d.squaredNorm();
      double bq = 2.0 * pF.dot(d);
      double c = pF.squaredNorm() - Delta * Delta;
      double tau = (-bq + std::sqrt(std::max(0.0, bq*bq - 4*a*c))) / (2*a);
      pF = pF + tau * d;
      hit_boundary = true;
      break;
    }
    double alpha = rr / dHd;
    VectorXd pTrial = pF + alpha * d;
    if (pTrial.norm() >= Delta) {
      // Intersect with trust-region boundary
      double a = d.squaredNorm();
      double bq = 2.0 * pF.dot(d);
      double c = pF.squaredNorm() - Delta * Delta;
      double tau = (-bq + std::sqrt(std::max(0.0, bq*bq - 4*a*c))) / (2*a);
      pF = pF + tau * d;
      hit_boundary = true;
      break;
    }
    // Regular CG step
    pF = pTrial;
    VectorXd rCG_new = rCG - alpha * Hd;
    double rr_new = rCG_new.squaredNorm();
    double beta = rr_new / rr;
    d = rCG_new + beta * d;
    rCG = rCG_new;
    rr = rr_new;
  }

  // Predicted reduction: -g_F^T pF - 0.5 pF^T H pF
  VectorXd HpF = Hmul(pF);
  double pred = -gF.dot(pF) - 0.5 * pF.dot(HpF);

  // Lift pF back to full space
  VectorXd p = VectorXd::Zero(x.size());
  for (int k = 0; k < (int)F.size(); ++k) p[F[k]] = pF[k];
  return {p, pred, hit_boundary};
}

// Reflect/clamp to bounds along step p; return alpha in (0,1] that hits first bound
inline double stepToBox(const VectorXd& x, const VectorXd& p,
                        const VectorXd& lb, const VectorXd& ub) {
  double alpha = 1.0;
  for (int i = 0; i < x.size(); ++i) {
    if (p[i] > 0.0) alpha = std::min(alpha, (ub[i] - x[i]) / p[i]);
    else if (p[i] < 0.0) alpha = std::min(alpha, (lb[i] - x[i]) / p[i]);
  }
  return std::max(0.0, alpha);
}

// Top-level TRR solver
struct TRROpts {
  double Delta0 = 1.0;
  double DeltaMax = 1e6;
  double eta = 1e-3;     // accept if rho >= eta
  double gtol = 1e-8;    // projected gradient tolerance (inf-norm)
  int max_iters = 200;
};

struct TRRResult {
  VectorXd x;
  double fval;
  int iters;
  bool converged;
};

inline TRRResult trr_boxed_qp(const QuadraticModel& q,
                       const VectorXd& lb,
                       const VectorXd& ub,
                       VectorXd x0,
                       const TRROpts& opts = {}) {
  VectorXd x = x0; projectToBox(x, lb, ub);
  double Delta = opts.Delta0;
  const MatrixXd H = 2.0 * q.G;

  auto proj_grad_inf = [&](const VectorXd& x, const VectorXd& g) {
    double normInf = 0.0;
    for (int i = 0; i < x.size(); ++i) {
      double gi = g[i];
      if      (x[i] <= lb[i] + 1e-12) gi = std::min(0.0, gi); // cannot go below lb
      else if (x[i] >= ub[i] - 1e-12) gi = std::max(0.0, gi); // cannot go above ub
      normInf = std::max(normInf, std::abs(gi));
    }
    return normInf;
  };

  for (int it = 0; it < opts.max_iters; ++it) {
    CostGrad cg = costGrad(q, x);
    if (proj_grad_inf(x, cg.g) < opts.gtol)
      return {x, cg.f, it, true};

    // Free set and TR subproblem
    auto F = freeSet(x, cg.g, lb, ub);
    if (F.empty()) return {x, cg.f, it, true}; // stuck at KKT on bounds

    TRStep st = steihaugCG(H, x, cg.g, F, Delta);
    if (st.p.norm() == 0.0) return {x, cg.f, it, true};

    // Respect bounds by clamping along p
    double alphaBox = stepToBox(x, st.p, lb, ub);
    VectorXd x_trial = x + alphaBox * st.p;

    // Actual & predicted reduction
    double f_old = cg.f;
    double f_new = costGrad(q, x_trial).f;
    double ared = f_old - f_new;
    double pred = st.pred_red * alphaBox; // linear model along clamped step (good enough)
    double rho = (pred > 0.0) ? (ared / pred) : -std::numeric_limits<double>::infinity();

    // Update trust region
    if (rho < 0.25) Delta *= 0.25;
    else if (rho > 0.75 && std::abs(alphaBox - 1.0) < 1e-12 && st.hit_boundary)
      Delta = std::min(2.0 * Delta, opts.DeltaMax);

    // Accept/reject
    if (rho >= opts.eta) x = x_trial; // accept
    // else reject and keep x, with smaller Delta
  }

  CostGrad cg_end = costGrad(q, x);
  return {x, cg_end.f, opts.max_iters, false};
}

namespace ispline_detail {
inline double cubic_ispline_impl(double x, double left, double right) {
  if (x < left) return 0.0;
  if (x >= right) return 1.0;
  double u = (x - left) / (right - left);
  return 3 * u * u - 2 * u * u * u;
}

inline double cubic_ispline_integral_impl(double x, double left, double right) {
  if (x <= left) return 0.0;
  const double width = right - left;
  if (!(width > 0.0)) return 0.0;
  if (x >= right) return 0.5 * width;

  const double u = (x - left) / width;
  return width * (u * u * u - 0.5 * u * u * u * u);
}

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

inline std::vector<double> normalized_rank_axis(size_t n) {
  std::vector<double> u(n);
  if (n == 0) return u;

  const double inv_n = 1.0 / static_cast<double>(n);
  for (size_t i = 0; i < n; ++i) {
    u[i] = static_cast<double>(i + 1) * inv_n;
  }
  return u;
}
}

inline Eigen::MatrixXd build_ispline_design(
    const std::vector<double>& scores,
    int /*degree*/,
    const std::vector<double>& knots_in,
    bool include_intercept)
{
  // Copy + sanitize knots: sort, unique, ensure at least two
  std::vector<double> knots = knots_in;
  if (knots.size() < 2) {
    // Derive a minimal knot span from data
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    for (double s : scores) { xmin = std::min(xmin, s); xmax = std::max(xmax, s); }
    if (!std::isfinite(xmin) || !std::isfinite(xmax)) { xmin = 0.0; xmax = 1.0; }
    if (xmax <= xmin) xmax = xmin + 1e-6; // avoid zero width
    knots = { xmin, xmax };
  }
  std::sort(knots.begin(), knots.end());
  knots.erase(std::unique(knots.begin(), knots.end()), knots.end());
  if (knots.size() < 2) {
    // If still < 2 after unique, expand slightly
    double k0 = knots.empty() ? 0.0 : knots.front();
    knots = { k0, k0 + 1e-6 };
  }

  const int n = static_cast<int>(scores.size());
  const int m = static_cast<int>(knots.size()) - 1; // number of intervals/bases
  const int p = m + (include_intercept ? 1 : 0);

  Eigen::MatrixXd X(n, p);
  int col0 = 0;
  if (include_intercept) {
    X.col(0).setOnes();
    col0 = 1;
  }

  for (int j = 0; j < m; ++j) {
    const double kj   = knots[j];
    const double kj1  = knots[j+1];
    const double wid  = std::max(1e-12, kj1 - kj); // protect against zero width
    for (int i = 0; i < n; ++i) {
      double x = scores[i];
      double val = ispline_detail::cubic_ispline_impl(x, kj, kj1);
      X(i, col0 + j) = val;
    }
  }

  return X;
}

inline Eigen::MatrixXd build_convex_q_design(
    const std::vector<double>& u,
    int /*degree*/,
    const std::vector<double>& knots_in,
    bool include_intercept)
{
  std::vector<double> knots = knots_in;
  if (knots.size() < 2) {
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    for (double s : u) { xmin = std::min(xmin, s); xmax = std::max(xmax, s); }
    if (!std::isfinite(xmin) || !std::isfinite(xmax)) { xmin = 0.0; xmax = 1.0; }
    if (xmax <= xmin) xmax = xmin + 1e-6;
    knots = { xmin, xmax };
  }
  std::sort(knots.begin(), knots.end());
  knots.erase(std::unique(knots.begin(), knots.end()), knots.end());
  if (knots.size() < 2) {
    double k0 = knots.empty() ? 0.0 : knots.front();
    knots = { k0, k0 + 1e-6 };
  }

  const int n = static_cast<int>(u.size());
  const int m = static_cast<int>(knots.size()) - 1;
  const int p = m + 1 + (include_intercept ? 1 : 0); // intercept + linear slope + integrated I-splines

  Eigen::MatrixXd X(n, p);
  int col = 0;
  if (include_intercept) {
    X.col(col++).setOnes();
  }

  for (int i = 0; i < n; ++i) {
    X(i, col) = u[i];
  }
  ++col;

  for (int j = 0; j < m; ++j) {
    const double kj = knots[j];
    const double kj1 = knots[j + 1];
    for (int i = 0; i < n; ++i) {
      X(i, col + j) = ispline_detail::cubic_ispline_integral_impl(u[i], kj, kj1);
    }
  }

  return X;
}

inline Eigen::MatrixXd build_convex_q_derivative_design(
    const std::vector<double>& u,
    int /*degree*/,
    const std::vector<double>& knots_in,
    bool include_intercept)
{
  std::vector<double> knots = knots_in;
  if (knots.size() < 2) {
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    for (double s : u) { xmin = std::min(xmin, s); xmax = std::max(xmax, s); }
    if (!std::isfinite(xmin) || !std::isfinite(xmax)) { xmin = 0.0; xmax = 1.0; }
    if (xmax <= xmin) xmax = xmin + 1e-6;
    knots = { xmin, xmax };
  }
  std::sort(knots.begin(), knots.end());
  knots.erase(std::unique(knots.begin(), knots.end()), knots.end());
  if (knots.size() < 2) {
    double k0 = knots.empty() ? 0.0 : knots.front();
    knots = { k0, k0 + 1e-6 };
  }

  const int n = static_cast<int>(u.size());
  const int m = static_cast<int>(knots.size()) - 1;
  const int p = m + 1 + (include_intercept ? 1 : 0);

  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(n, p);
  int col = 0;
  if (include_intercept) {
    ++col; // derivative of intercept is zero
  }

  X.col(col).setOnes(); // derivative of linear term u
  ++col;

  for (int j = 0; j < m; ++j) {
    const double kj = knots[j];
    const double kj1 = knots[j + 1];
    for (int i = 0; i < n; ++i) {
      X(i, col + j) = ispline_detail::cubic_ispline_impl(u[i], kj, kj1);
    }
  }

  return X;
}


static std::vector<double> make_default_knots(const std::vector<double>& scores, int degree) {
  const int n = (int)scores.size();
  if (n == 0) return {};

  // 1. Sort copy of scores
  std::vector<double> sorted = scores;
  std::sort(sorted.begin(), sorted.end());

  // 2. Boundaries
  double lo = sorted.front();
  double hi = sorted.back();

  // 3. Number of internal knots
  int num_knots = std::min(200, (int)std::sqrt(n));

  // 4. Place them at equally spaced quantiles
  std::vector<double> knots;
  knots.reserve(num_knots + 2* (degree+1));
  knots.push_back(lo);
  for (int k = 1; k <= num_knots; ++k) {
    double q = (double)k / (num_knots+1); // exclude endpoints
    size_t idx = (size_t)(q * (n-1));
    knots.push_back(sorted[idx]);
  }
  knots.push_back(hi);

  // 5. Depending on your basis implementation, you may need to pad
  //    the boundary knots 'degree' times on each side.
  for (int d = 0; d < degree; ++d) {
    knots.insert(knots.begin(), lo);
    knots.push_back(hi);
  }

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
    // ||X beta - y||^2 + lambda ||beta||^2
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

    // Bounds: non-negative except possibly intercept
    Eigen::VectorXd lb = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd ub = Eigen::VectorXd::Constant(p, std::numeric_limits<double>::infinity());
    if (params_.intercept_col >= 0 && params_.intercept_col < p)
      lb[params_.intercept_col] = -std::numeric_limits<double>::infinity();

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(p);
    TRROpts opts; opts.Delta0 = 1.0; opts.gtol = 1e-8; opts.max_iters = 1000;
    TRRResult sol = trr_boxed_qp(q, lb, ub, x0, opts);
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

  std::vector<double> fit_q_xy(const std::vector<double>& x,
                               const std::vector<double>& q_values,
                               double clip_lo = 0.0, double clip_hi = 1.0) {
    if (x.size() != q_values.size()) {
      throw std::invalid_argument("ISplineTRRRegressor::fit_q_xy: x and q_values size mismatch");
    }
    if (q_values.empty()) {
      return {};
    }

    std::vector<size_t> order(q_values.size());
    std::iota(order.begin(), order.end(), 0);
    const bool descending_x = params_.y_decreasing_in_x;
    std::stable_sort(order.begin(), order.end(),
                     [&](size_t lhs, size_t rhs) {
                       if (x[lhs] == x[rhs]) {
                         return lhs < rhs;
                       }
                       return descending_x ? (x[lhs] > x[rhs]) : (x[lhs] < x[rhs]);
                     });

    std::vector<double> q_sorted(q_values.size());
    for (size_t i = 0; i < order.size(); ++i) {
      q_sorted[i] = q_values[order[i]];
    }

    const auto pep_sorted = fit_q_sorted(q_sorted, clip_lo, clip_hi);
    std::vector<double> out(q_values.size());
    for (size_t i = 0; i < order.size(); ++i) {
      out[order[i]] = pep_sorted[i];
    }
    return out;
  }

  std::vector<double> fit_q_y(const std::vector<double>& q_values,
                              double clip_lo = 0.0, double clip_hi = 1.0) {
    return fit_q_sorted(q_values, clip_lo, clip_hi);
  }

  inline double cubic_ispline(double x, double left, double right) const {
    return ispline_detail::cubic_ispline_impl(x, left, right);
  }

private:
  std::vector<double> fit_q_sorted(const std::vector<double>& q_values,
                                   double clip_lo,
                                   double clip_hi) const {
    const int n = static_cast<int>(q_values.size());
    if (n == 0) {
      return {};
    }

    const std::vector<double> u = ispline_detail::normalized_rank_axis(q_values.size());
    auto knots = params_.knots;
    if (knots.empty()) {
      knots = make_default_knots(u, params_.ispline_degree);
    }

    Eigen::MatrixXd Xq = build_convex_q_design(
        u, params_.ispline_degree, knots, params_.include_intercept);
    Eigen::MatrixXd Xdq = build_convex_q_derivative_design(
        u, params_.ispline_degree, knots, params_.include_intercept);
    const int p = static_cast<int>(Xq.cols());

    const Eigen::Map<const Eigen::VectorXd> q_vec(q_values.data(), n);
    const double inv_n = 1.0 / std::max(1, n);
    QuadraticModel qobj;
    qobj.G = inv_n * (Xq.transpose() * Xq);
    qobj.h = inv_n * (Xq.transpose() * q_vec);
    qobj.c = inv_n * q_vec.squaredNorm();

    const double lambda = std::max(1e-2, params_.ridge_lambda);
    if (lambda > 0.0) {
      // Penalize the fitted slope directly so p(u) = q(u) + u q'(u) remains
      // stable and does not spike on narrow spline intervals.
      qobj.G += lambda * inv_n * (Xdq.transpose() * Xdq);

      // Keep a light coefficient ridge as a numerical backstop.
      qobj.G.diagonal().array() += 1e-8;
      if (params_.include_intercept &&
          params_.intercept_col >= 0 &&
          params_.intercept_col < p) {
        qobj.G(params_.intercept_col, params_.intercept_col) -= 1e-8;
      }
    }

    Eigen::VectorXd lb = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd ub = Eigen::VectorXd::Constant(p, std::numeric_limits<double>::infinity());
    if (params_.include_intercept && params_.intercept_col >= 0 && params_.intercept_col < p) {
      lb[params_.intercept_col] = clip_lo;
      ub[params_.intercept_col] = clip_hi;
    }

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(p);
    if (params_.include_intercept && params_.intercept_col >= 0 && params_.intercept_col < p) {
      x0[params_.intercept_col] = std::max(clip_lo, std::min(clip_hi, q_values.front()));
    }

    TRROpts opts;
    opts.Delta0 = 1.0;
    opts.gtol = 1e-8;
    opts.max_iters = 1000;
    TRRResult sol = trr_boxed_qp(qobj, lb, ub, x0, opts);

    std::vector<double> pep(n);
    for (int i = 0; i < n; ++i) {
      const double qhat = Xq.row(i).dot(sol.x);
      const double qprime = Xdq.row(i).dot(sol.x);
      double local_pep = qhat + u[i] * qprime;
      if (local_pep < clip_lo) local_pep = clip_lo;
      if (local_pep > clip_hi) local_pep = clip_hi;
      pep[i] = local_pep;
    }
    return pep;
  }
};
