#include "Globals.h"
#include <vector>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <numeric>    // for std::iota
#include <cassert>
#include <iostream>
#include <chrono>
using clock_type = std::chrono::high_resolution_clock;

#include <Eigen/Dense>
#include <Eigen/Sparse> 

#include "MonotoneRegressor.h"
#include "IsotonicPEP.h"

InferPEP::InferPEP(bool use_ispline)
{
    if (use_ispline) {
        regressor_ptr_ = make_monotone_regressor(RegressorType::ISPLINE_TRR, MonotoneParams());
        if (VERB > 1) {
            std::cerr << "Performing isotonic regression using I-Splines" << std::endl;
        }                
    } else {
        regressor_ptr_ = make_monotone_regressor(RegressorType::PAVA, MonotoneParams());
        if (VERB > 1) {
            std::cerr << "Performing isotonic regression using PAVA" << std::endl;
        }
    }
}

std::vector<double> InferPEP::q_to_pep(const std::vector<double>& q_values) {
    qs = q_values;

    std::vector<double> qn(q_values.size());
    std::vector<int> indices(q_values.size());
    std::iota(indices.begin(), indices.end(), 1);

    for (size_t i = 0; i < q_values.size(); ++i) {
        qn[i] = q_values[i] * indices[i];
        assert((i == q_values.size() - 1) || (q_values[i] <= q_values[i + 1]));
    }
    std::vector<double> raw_pep(qn.size());
    raw_pep[0] = qn[0];
    for (size_t i = 1; i < qn.size(); ++i) {
        raw_pep[i] = qn[i] - qn[i - 1];
    }   
    return regressor_ptr_->fit_y(raw_pep);
}

std::vector<double> InferPEP::qns_to_pep(const std::vector<double>& q_values, const std::vector<double>& scores) {
    qs = q_values;

    std::vector<double> qn(q_values.size());
    std::vector<int> indices(q_values.size());
    std::iota(indices.begin(), indices.end(), 1);

    for (size_t i = 0; i < q_values.size(); ++i) {
        qn[i] = q_values[i] * indices[i];
        assert((i == q_values.size() - 1) || (q_values[i] <= q_values[i + 1]));
    }

    std::vector<double> raw_pep(qn.size());
    raw_pep[0] = qn[0];
    for (size_t i = 1; i < qn.size(); ++i) {
        raw_pep[i] = qn[i] - qn[i - 1];
    }

    return regressor_ptr_->fit_xy(scores, raw_pep);
}

std::vector<double>
InferPEP::tdc_to_pep(const std::vector<double>& is_decoy,
                     const std::vector<double>& scores) {
  using namespace std::chrono;
  auto start = clock_type::now();
  if (VERB > 2) std::cerr << "[TIMING] entering tdc_to_pep\n";

  const double epsilon = 1e-20;

  // Build target vector (prepend as your original code does to anchor the fit)
  std::vector<double> is_dec = is_decoy;
  is_dec.insert(is_dec.begin(), 0.5);

  std::vector<double> sc;
  if (!scores.empty()) {
    sc = scores;
    sc.insert(sc.begin(), sc.front());
  }

  // Spline / ridge settings (reuse yours)
  const int degree = ispline_degree_;
  const double lambda = ridge_lambda_tdc_;
  const bool include_intercept = include_intercept_;
  const int intercept_col = include_intercept ? 0 : -1;
  const std::vector<double>& knots = ispline_knots_;

  std::vector<double> w(is_dec.size(), 1.0);

  // Fit decoy-rate in (epsilon, 1-epsilon); clamp outputs
  std::vector<double> decoy_rate;
  if (!scores.empty()) {
    regressor_ptr_->fit_xy(scores, decoy_rate, /*clip_lo=*/epsilon, /*clip_hi=*/1.0 - epsilon);
    // decoy_rate = fit_ispline_trr_predict(
    //    sc, is_dec, w, lambda, degree, knots, include_intercept,
    //    intercept_col, /*clip_lo=*/epsilon, /*clip_hi=*/1.0 - epsilon);
    // drop the prepended element
    // decoy_rate.erase(decoy_rate.begin());
  } else {
    regressor_ptr_->fit_y(decoy_rate, /*clip_lo=*/epsilon, /*clip_hi=*/1.0 - epsilon);
  }

  // Convert to PEP = p(decoy) / (1 - p(decoy))
  std::vector<double> pep_iso; pep_iso.reserve(decoy_rate.size());
  for (double dp : decoy_rate) {
    double pep = dp / (1.0 - dp);
    if (pep > 1.0) pep = 1.0;
    pep_iso.push_back(pep);
  }

  auto end = clock_type::now();
  if (VERB > 2) {
    double duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cerr << "[TIMING] tdc_to_pep duration: " << duration << " seconds\n";
  }
  return pep_iso;
}




