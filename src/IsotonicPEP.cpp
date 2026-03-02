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

namespace {
std::vector<double> q_values_to_raw_pep(const std::vector<double>& q_values) {
  const size_t n = q_values.size();
  std::vector<double> raw_pep(n);
  if (n == 0) {
    return raw_pep;
  }

  double prev_qn = 0.0;
  for (size_t i = 0; i < n; ++i) {
    assert((i + 1 == n) || (q_values[i] <= q_values[i + 1]));
    const double qn = q_values[i] * static_cast<double>(i + 1);
    raw_pep[i] = qn - prev_qn;
    prev_qn = qn;
  }
  return raw_pep;
}
}  // namespace

InferPEP::InferPEP(bool use_ispline)
{
    if (use_ispline) {
        regressor_ptr_ = make_monotone_regressor(RegressorType::ISPLINE_TRR);
        if (VERB > 1) {
            std::cerr << "Performing isotonic regression using I-Splines" << std::endl;
        }                
    } else {
        regressor_ptr_ = make_monotone_regressor(RegressorType::PAVA);
        if (VERB > 1) {
            std::cerr << "Performing isotonic regression using PAVA" << std::endl;
        }
    }
}

std::vector<double> InferPEP::q_to_pep(const std::vector<double>& q_values) {
    qs = q_values;
    const std::vector<double> raw_pep = q_values_to_raw_pep(q_values);
    return regressor_ptr_->fit_y(raw_pep);
}

std::vector<double> InferPEP::qns_to_pep(const std::vector<double>& q_values, const std::vector<double>& scores) {
    qs = q_values;
    const std::vector<double> raw_pep = q_values_to_raw_pep(q_values);
    return regressor_ptr_->fit_xy(scores, raw_pep);
}

std::vector<double>
InferPEP::tdc_to_pep(const std::vector<double>& is_decoy,
                     const std::vector<double>& scores) {
  const bool print_timing = (VERB > 2);
  const auto start = print_timing ? clock_type::now() : clock_type::time_point{};
  if (VERB > 2) std::cerr << "[TIMING] entering tdc_to_pep\n";

  const double epsilon = 1e-20;

  // Decide whether we will add an anchor (only in the fit_xy path)
  const bool will_use_fit_xy = !scores.empty();

  std::vector<double> is_dec;
  std::vector<double> sc;
  if (will_use_fit_xy) {
    is_dec.reserve(is_decoy.size() + 1);
    sc.reserve(scores.size() + 1);
    is_dec.push_back(0.5);     // neutral label
    sc.push_back(scores.front());      // duplicate top score
    is_dec.insert(is_dec.end(), is_decoy.begin(), is_decoy.end());
    sc.insert(sc.end(), scores.begin(), scores.end());
  } else {
    is_dec = is_decoy;
  }

  // Fit decoy rate p(decoy | x)
  const std::vector<double> decoy_rate = will_use_fit_xy
      ? regressor_ptr_->fit_xy(sc, is_dec, /*clip_lo=*/epsilon, /*clip_hi=*/1.0 - epsilon)
      : regressor_ptr_->fit_y(is_dec,      /*clip_lo=*/epsilon, /*clip_hi=*/1.0 - epsilon);

  const size_t offset = will_use_fit_xy ? 1 : 0;
  std::vector<double> pep_iso(decoy_rate.size() - offset);
  for (size_t i = offset; i < decoy_rate.size(); ++i) {
    double p = decoy_rate[i];
    double pep = p / (1.0 - p);  
    pep = std::max(0.0, std::min(1.0, pep));
    pep_iso[i - offset] = pep;
  }

  if (print_timing) {
    const double duration = std::chrono::duration_cast<std::chrono::duration<double>>(clock_type::now() - start).count();
    std::cerr << "[TIMING] tdc_to_pep duration: " << duration << " seconds\n";
  }

  return pep_iso;
}


