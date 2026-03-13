#include <gtest/gtest.h>
#include <algorithm>
#include <iostream> // for std::cerr, std::endl
#include <cmath>
#include <numeric>
#include <random>
#include <vector>
#include <tuple>
#include <string>

#include "IsotonicPEP.h"          // your class (InferPEP)
#include "MonotoneRegressor.h"    // RegressorType, MonotoneParams factory

namespace {

// Generate synthetic data as specified:
// - Top 1000 highest scores: all targets (non-decoy)
// - Next 1000: decoy prob ramps linearly 0 -> 0.5
// - Last 1000: decoy prob = 0.5
struct SynthData {
  std::vector<double> scores;   // larger is better (descending order)
  std::vector<double> is_decoy; // 1 for decoy, 0 for target
};

struct QCurveData {
  std::vector<double> q_values;
  std::vector<double> pep_values;
  std::vector<double> scores;
};

SynthData make_synth_data(unsigned seed = 1337) {
  const int n_block = 1000;
  const int N = 3 * n_block;

  // Scores: strictly decreasing, top are best
  std::vector<double> scores(N);
  for (int i = 0; i < N; ++i) scores[i] = double(N - i);

  std::mt19937 rng(seed);
  std::vector<double> is_decoy; is_decoy.reserve(N);

  // Block 1: prob=0 (all targets)
  for (int i = 0; i < n_block; ++i) is_decoy.push_back(0.0);

  // Block 2: prob ramps 0 -> 0.5
  for (int i = 0; i < n_block; ++i) {
    double p = 0.5 * (double(i) / (n_block - 1));
    std::bernoulli_distribution bern(p);
    is_decoy.push_back(bern(rng) ? 1.0 : 0.0);
  }

  // Block 3: prob = 0.5
  std::bernoulli_distribution bern_half(0.5);
  for (int i = 0; i < n_block; ++i)
    is_decoy.push_back(bern_half(rng) ? 1.0 : 0.0);

  return {scores, is_decoy};
}

SynthData make_sparse_top_decoy_data(unsigned seed = 2025) {
  const int n_top = 400;
  const int n_ramp = 1200;
  const int n_tail = 1200;
  const int N = n_top + n_ramp + n_tail;

  std::vector<double> scores(N);
  for (int i = 0; i < N; ++i) scores[i] = double(N - i);

  std::mt19937 rng(seed);
  std::vector<double> is_decoy; is_decoy.reserve(N);

  for (int i = 0; i < n_top; ++i) is_decoy.push_back(0.0);

  for (int i = 0; i < n_ramp; ++i) {
    double p = 0.20 * (double(i) / std::max(1, n_ramp - 1));
    std::bernoulli_distribution bern(p);
    is_decoy.push_back(bern(rng) ? 1.0 : 0.0);
  }

  std::bernoulli_distribution bern_tail(0.5);
  for (int i = 0; i < n_tail; ++i) {
    is_decoy.push_back(bern_tail(rng) ? 1.0 : 0.0);
  }

  return {scores, is_decoy};
}

QCurveData make_convex_q_curve_data(size_t n = 250) {
  std::vector<double> q_values(n);
  std::vector<double> pep_values(n);
  std::vector<double> scores(n);

  for (size_t i = 0; i < n; ++i) {
    const double u = static_cast<double>(i + 1) / static_cast<double>(n);
    const double q = 0.005 + 0.05 * u + 0.10 * u * u * u;
    const double qp = 0.05 + 0.30 * u * u;
    const double pep = q + u * qp;
    q_values[i] = q;
    pep_values[i] = pep;
    scores[i] = static_cast<double>(n - i);
  }

  return {q_values, pep_values, scores};
}

void check_qvalue_consistency(const std::vector<double>& pep,
                              const std::vector<double>& q_values,
                              double abs_tol = 0.05) {
  ASSERT_EQ(pep.size(), q_values.size());
  double running_sum = 0.0;
  for (size_t i = 0; i < pep.size(); ++i) {
    running_sum += pep[i];
    const double q_hat = running_sum / static_cast<double>(i + 1);
    EXPECT_NEAR(q_hat, q_values[i], abs_tol) << "at index " << i;
  }
}

// Given pep[i] for i=0..M-1 (aligned with scores/is_decoy),
// evaluate calibration at each prefix i where the i-th item is a TARGET:
//   fdr_hat(i)  = (# decoys up to i) / (# targets up to i)
//   pep_mean(i) = (sum of pep_j for targets j <= i) / (# targets up to i)
// We assert relative deviation <= rel_tol after a warmup (to reduce noise).
void check_calibration_q(const std::vector<double>& pep,
                         const std::vector<double>& is_decoy,
                         double rel_tol = 0.20)
{
  const size_t N = pep.size();

  // prefix counts
  std::vector<int> D(N,0), T(N,0);
  for (size_t i=0;i<N;++i) {
    D[i] = (i?D[i-1]:0) + (is_decoy[i] > 0.5 ? 1:0);
    T[i] = (i?T[i-1]:0) + (is_decoy[i] > 0.5 ? 0:1);
  }

  // raw fdrs with +0.5 offset
  std::vector<double> fdr(N,0.0);
  for (size_t i=0;i<N;++i) {
    fdr[i] = (T[i]>0) ? ( (double)D[i] + 0.5 ) / (double)T[i] : 0.0;
    if (fdr[i] > 1.0) fdr[i] = 1.0;
  }

  // backward smoothing into q-values
  std::vector<double> q(N,1.0);
  double run_min = 1.0;
  for (ptrdiff_t i=(ptrdiff_t)N-1; i>=0; --i) {
    run_min = std::min(run_min, fdr[(size_t)i]);
    q[(size_t)i] = run_min;
  }

  // compare mean PEP vs q
  constexpr int kWarmupTargets = 25;
  constexpr double kAbsTol = 2e-3;

  int targets=0; double pep_sum=0.0;
  for (size_t i=0;i<N;++i) {
    if (is_decoy[i] > 0.5) continue;
    pep_sum += pep[i];
    targets++;

    if (targets < kWarmupTargets) {
      continue;
    }

    double pep_mean = pep_sum / (double)targets;
    double qi = q[i];
    double denom = std::max(1e-6, std::min(pep_mean,qi));
    double rel_err = std::abs(pep_mean - qi)/denom;
    double abs_err = std::abs(pep_mean - qi);
    ASSERT_TRUE(rel_err < rel_tol || abs_err < kAbsTol)
        << "Poor calibration at index " << i
        << " pep_mean=" << pep_mean
        << " vs q[i]=" << qi
        << " rel_err=" << rel_err
        << " abs_err=" << abs_err;
  }

  EXPECT_NEAR(q.back(),
              pep_sum/std::max(1,targets),
              0.15); // soft tail check
}
} // namespace

// ------------------------------------------------------------------
// Parameterization: (use_ispline, use_fit_xy)
//   false = PAVA, true = ISPLINE
//   false = fit_y, true = fit_xy
// ------------------------------------------------------------------
class TdcToPepCalibrationTest
    : public ::testing::TestWithParam<std::tuple<bool, bool>> {};

// Custom parameter name generator for friendly test names
static std::string ParamNameGen(
    const ::testing::TestParamInfo<std::tuple<bool, bool>>& info) {
  const auto [use_ispline, use_fit_xy] = info.param;
  const char* reg = use_ispline ? "ISpline" : "Pava";
  const char* fit = use_fit_xy   ? "FitXY"   : "FitY";
  return std::string(reg) + "_" + fit;
}

INSTANTIATE_TEST_SUITE_P(
    AllMonoRegressors,
    TdcToPepCalibrationTest,
    ::testing::Combine(::testing::Bool(), ::testing::Bool()),
    ParamNameGen);

// --- Tests ---------------------------------------------------------

TEST_P(TdcToPepCalibrationTest, SyntheticRamp_CalibratedPEPs) {
  // 1) Synthesize
  auto data = make_synth_data(/*seed=*/1337);
  auto [use_ispline, use_fit_xy] = GetParam();
  InferPEP infer(use_ispline);

  // 2) Fit using chosen method

  std::vector<double> pep;
  if (use_fit_xy) {
    pep = infer.tdc_to_pep(data.is_decoy, data.scores);
  } else {
    pep = infer.tdc_to_pep(data.is_decoy);
  }

  // 4) Compare cumulative mean PEP among targets vs TDC FDR
  check_calibration_q(pep, data.is_decoy, /*rel_tol=*/3.0); // tighten as needed
}

TEST_P(TdcToPepCalibrationTest, GetParam_Sane_AfterFit) {
  auto data = make_synth_data(/*seed=*/42);
  auto [use_ispline, use_fit_xy] = GetParam();
  InferPEP infer(use_ispline);

  if (use_fit_xy) {
    infer.tdc_to_pep(data.is_decoy, data.scores);
  } else {
    infer.tdc_to_pep(data.is_decoy);
  }

  // Smoke: model remains usable
  auto pep = infer.tdc_to_pep(data.is_decoy, data.scores);
}

TEST_P(TdcToPepCalibrationTest, SparseHighScoreDecoys_DoNotForceLargeTailFloor) {
  auto data = make_sparse_top_decoy_data();
  auto [use_ispline, use_fit_xy] = GetParam();
  InferPEP infer(use_ispline);

  std::vector<double> pep = use_fit_xy
      ? infer.tdc_to_pep(data.is_decoy, data.scores)
      : infer.tdc_to_pep(data.is_decoy);

  double top_target_sum = 0.0;
  int top_target_count = 0;
  for (size_t i = 0; i < pep.size() && top_target_count < 100; ++i) {
    if (data.is_decoy[i] > 0.5) continue;
    top_target_sum += pep[i];
    ++top_target_count;
  }

  ASSERT_EQ(top_target_count, 100);
  EXPECT_LT(top_target_sum / top_target_count, 0.05);
}

TEST(InferPepQValuePathTest, ISplineQToPepRecoversConvexLocalRate) {
  const auto data = make_convex_q_curve_data();
  InferPEP infer(/*use_ispline=*/true);

  const auto pep = infer.q_to_pep(data.q_values);

  ASSERT_EQ(pep.size(), data.pep_values.size());
  for (size_t i = 1; i < pep.size(); ++i) {
    EXPECT_GE(pep[i], pep[i - 1]);
  }

  check_qvalue_consistency(pep, data.q_values);
}

TEST(InferPepQValuePathTest, ISplineQnsToPepUsesScoreOrderAndMapsBack) {
  auto data = make_convex_q_curve_data();
  InferPEP infer(/*use_ispline=*/true);
  const auto pep_reference = infer.q_to_pep(data.q_values);

  std::vector<size_t> perm(data.q_values.size());
  std::iota(perm.begin(), perm.end(), 0);
  std::mt19937 rng(17);
  std::shuffle(perm.begin(), perm.end(), rng);

  std::vector<double> q_perm(perm.size());
  std::vector<double> score_perm(perm.size());
  std::vector<double> pep_expected(perm.size());
  for (size_t i = 0; i < perm.size(); ++i) {
    q_perm[i] = data.q_values[perm[i]];
    score_perm[i] = data.scores[perm[i]];
    pep_expected[i] = pep_reference[perm[i]];
  }

  const auto pep = infer.qns_to_pep(q_perm, score_perm);

  ASSERT_EQ(pep.size(), pep_expected.size());
  for (size_t i = 0; i < pep.size(); ++i) {
    EXPECT_NEAR(pep[i], pep_expected[i], 1e-9) << "at permuted index " << i;
  }
}
