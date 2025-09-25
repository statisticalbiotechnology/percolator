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

  // Work copies
  std::vector<double> is_dec = is_decoy;
  std::vector<double> sc = scores;

  // Decide whether we will add an anchor (only in the fit_xy path)
  const bool will_use_fit_xy = !scores.empty();
  bool used_anchor = false;

  if (will_use_fit_xy) {
    // Prepend anchor to BOTH x and y so sizes stay equal
    is_dec.insert(is_dec.begin(), 0.5);     // neutral label
    sc.insert(sc.begin(), sc.front());      // duplicate top score
    used_anchor = true;
  }

  // Fit decoy rate p(decoy | x)
  std::vector<double> decoy_rate;
  if (will_use_fit_xy) {
    decoy_rate = regressor_ptr_->fit_xy(sc, is_dec, /*clip_lo=*/epsilon, /*clip_hi=*/1.0 - epsilon);
  } else {
    decoy_rate = regressor_ptr_->fit_y(is_dec,      /*clip_lo=*/epsilon, /*clip_hi=*/1.0 - epsilon);
  }

  // Drop the anchor we added, if any
  if (used_anchor) {
    // sc.size() == scores.size() + 1 when used_anchor is true
    decoy_rate.erase(decoy_rate.begin());
  }

  std::vector<double> pep_iso;
  pep_iso.reserve(decoy_rate.size());
  for (size_t i = 0; i < decoy_rate.size(); ++i) {
    double p = decoy_rate[i];
    double pep = p / (1.0 - p);  
    pep = std::max(0.0, std::min(1.0, pep));
    pep_iso.push_back(pep);
  }

  if (VERB > 2) {
    double duration = std::chrono::duration_cast<std::chrono::duration<double>>(clock_type::now() - start).count();
    std::cerr << "[TIMING] tdc_to_pep duration: " << duration << " seconds\n";
  }

  return pep_iso;
}




