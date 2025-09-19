#ifndef ISOTONIC_REGRESSION_H_
#define ISOTONIC_REGRESSION_H_

#include <chrono>
#include <iostream> // for std::cerr, std::endl
#include <memory>   // for std::unique_ptr, std::make_unique
#include "MonotoneRegressor.h"

class MonotoneRegressor;

class InferPEP {
    public:
        InferPEP(bool use_ispline);
    
        std::vector<double> q_to_pep(const std::vector<double>& q_values);
        // Calibrate PEPs from q-values (qns) or from target/decoy labels (TDC)
        std::vector<double> qns_to_pep(const std::vector<double>& q_values, const std::vector<double>& scores);
        std::vector<double> tdc_to_pep(const std::vector<double>& is_decoy, const std::vector<double>& scores = {});

        double interpolate(const double q_value, const double q1, const double q2, const double pep1, const double pep2) const {
            double interp_pep = pep1 + (q_value - q1) * (pep2 - pep1) / (q2 - q1);
            return interp_pep;
        }
        
    private:
        std::unique_ptr<MonotoneRegressor> regressor_ptr_;
        std::vector<double> qs;
        // std::vector<double> pep_iso;
    };

    #endif /* ISOTONICPEP_H_ */