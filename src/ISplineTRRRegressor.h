#include "MonotoneRegressor.h"
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


// Compute cost and gradient for augmented LS: f = ||A x - b||^2, g = 2 A^T (A x - b)
struct CostGrad {
  double f;
  VectorXd g;
};

inline CostGrad costGrad(const MatrixXd& A, const VectorXd& b, const VectorXd& x) {
  VectorXd r = A * x - b;
  CostGrad cg;
  cg.f = r.squaredNorm();
  cg.g = 2.0 * (A.transpose() * r);
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
// Here H = 2*A_F^T*A_F, but we never form H explicitly; we do Hv = 2*A_F^T*(A_F*v).
struct TRStep { VectorXd p; double pred_red; bool hit_boundary; };

inline TRStep steihaugCG(const MatrixXd& A, const VectorXd& b,
                  const VectorXd& x, const std::vector<int>& F, double Delta) {
  // Work only on free coordinates via views
  // Define operators on free set: apply H_F v = 2*A_F^T*(A_F v)
  auto Hmul = [&](const VectorXd& vF) {
    // Build a full vector v with zeros outside F, apply A, then A^T, then extract F
    VectorXd v = VectorXd::Zero(x.size());
    for (int k = 0; k < (int)F.size(); ++k) v[F[k]] = vF[k];
    VectorXd Av = A * v;                  // size m
    VectorXd AtAv = 2.0 * (A.transpose() * Av); // size n
    VectorXd out(vF.size());
    for (int k = 0; k < (int)F.size(); ++k) out[k] = AtAv[F[k]];
    return out;
  };
  // Gradient on free set: gF
  VectorXd r = A * x - b;
  VectorXd g = 2.0 * (A.transpose() * r);
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

inline TRRResult trr_boxed_ls(const MatrixXd& Atilde, const VectorXd& btilde,
                       const VectorXd& lb, const VectorXd& ub, VectorXd x0,
                       const TRROpts& opts = {}) {
  VectorXd x = x0; projectToBox(x, lb, ub);
  double Delta = opts.Delta0;

  auto proj_grad_inf = [&](const VectorXd& x) {
    CostGrad cg = costGrad(Atilde, btilde, x);
    double normInf = 0.0;
    for (int i = 0; i < x.size(); ++i) {
      double gi = cg.g[i];
      if      (x[i] <= lb[i] + 1e-12) gi = std::min(0.0, gi); // cannot go below lb
      else if (x[i] >= ub[i] - 1e-12) gi = std::max(0.0, gi); // cannot go above ub
      normInf = std::max(normInf, std::abs(gi));
    }
    return normInf;
  };

  for (int it = 0; it < opts.max_iters; ++it) {
    CostGrad cg = costGrad(Atilde, btilde, x);
    if (proj_grad_inf(x) < opts.gtol)
      return {x, cg.f, it, true};

    // Free set and TR subproblem
    auto F = freeSet(x, cg.g, lb, ub);
    if (F.empty()) return {x, cg.f, it, true}; // stuck at KKT on bounds

    TRStep st = steihaugCG(Atilde, btilde, x, F, Delta);
    if (st.p.norm() == 0.0) return {x, cg.f, it, true};

    // Respect bounds by clamping along p
    double alphaBox = stepToBox(x, st.p, lb, ub);
    VectorXd x_trial = x + alphaBox * st.p;

    // Actual & predicted reduction
    double f_old = cg.f;
    double f_new = costGrad(Atilde, btilde, x_trial).f;
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

  CostGrad cg_end = costGrad(Atilde, btilde, x);
  return {x, cg_end.f, opts.max_iters, false};
}

namespace ispline_detail {
inline double cubic_ispline_impl(double x, double left, double right) {
  if (x < left) return 0.0;
  if (x >= right) return 1.0;
  double u = (x - left) / (right - left);
  return 3 * u * u - 2 * u * u * u;
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


class ISplineTRRRegressor final : public MonotoneRegressor {
public:
  using MonotoneRegressor::MonotoneRegressor;

  std::vector<double> fit_xy(const std::vector<double>& x,
                             const std::vector<double>& y,
                             double clip_lo = 0.0, double clip_hi = 1.0) override {
    assert(x.size() == y.size());
    const int n = (int)x.size();

    // Weights: use 1s unless you already carry them
    std::vector<double> w(n, 1.0);

    // Build X
    Eigen::MatrixXd X = build_ispline_design(
        x, params_.ispline_degree, params_.knots, params_.include_intercept);
    const int p = (int)X.cols();

    // Weighted + ridge-augmented system
    Eigen::MatrixXd A(n + p, p);
    Eigen::VectorXd b(n + p);
    for (int i = 0; i < n; ++i) {
      double s = 1.0; // sqrt(w[i]) if you have weights
      A.row(i) = X.row(i) * s;
      b[i]     = y[i] * s;
    }
    A.bottomRows(p).setZero();
    for (int j = 0; j < p; ++j)
      A(n + j, j) = std::sqrt(std::max(0.0, params_.ridge_lambda));
    b.tail(p).setZero();

    // Bounds: non-negative except possibly intercept
    Eigen::VectorXd lb = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd ub = Eigen::VectorXd::Constant(p, std::numeric_limits<double>::infinity());
    if (params_.intercept_col >= 0 && params_.intercept_col < p)
      lb[params_.intercept_col] = -std::numeric_limits<double>::infinity();

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(p);
    TRROpts opts; opts.Delta0 = 1.0; opts.gtol = 1e-8; opts.max_iters = 200;
    TRRResult sol = trr_boxed_ls(A, b, lb, ub, x0, opts);

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
    std::iota(x.begin(), x.end(), 0.0);
    return fit_xy(x, y, clip_lo, clip_hi);
  }

   inline double cubic_ispline(double x, double left, double right) const {
    return ispline_detail::cubic_ispline_impl(x, left, right);
  }
};

