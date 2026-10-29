/*******************************************************************************
 * SVMlin
 * Copyright (c) 2006 Vikas Sindhwani at the University of Chicago.
 * Adapted to Percolator by Lukas Käll at the University of Washington
 * Sped up by John Halloran at the University of California, Davis, as detailed in:
 ******************************
 * A Matter of Time: Faster Percolator Analysis via Efficient SVM Learning for 
 * Large-Scale Proteomics
 * John T. Halloran and David M. Rocke
 * Journal of Proteome Research 2018 17 (5), 1978-1982
 *******************************************************************************/
#ifndef _svmlin_H
#define _svmlin_H
#include <vector>
#include <ctime>

using namespace std;

/* OPTIMIZATION CONSTANTS */
#define CGITERMAX 10000 /* maximum number of CGLS iterations */
#define SMALL_CGITERMAX 10 /* for heuristic 1 in reference [2] */
#define EPSILON   1e-7 /* most tolerances are set to this value */
#define BIG_EPSILON 0.01 /* for heuristic 2 in reference [2] */
#define RELATIVE_STOP_EPS 1e-9 /* for L2-SVM-MFN relative stopping criterion */
#define MFNITERMAX 50 /* maximum number of MFN iterations */

#define VERBOSE_CGLS 0

/* Training examples of the SVM. The feature rows are copied into one
   contiguous row-major matrix, with the intercept (1.0) as the last column,
   so that the solver streams through memory instead of following pointers
   to scattered rows. */
class AlgIn {
  public:
    explicit AlgIn(const int numFeat) : n(numFeat) {}
    int m = 0; /* number of examples */
    int n; /* number of features, including the intercept */
    int positives = 0;
    int negatives = 0;
    std::vector<double> rows; /* m x n */
    std::vector<double> Y; /* labels */

    const double* row(int i) const {
      return rows.data() + static_cast<size_t>(i) * static_cast<size_t>(n);
    }
    void clear() {
      m = positives = negatives = 0;
      rows.clear();
      Y.clear();
    }
    /* appends an example with n-1 features */
    void add(const double* features, double label) {
      rows.insert(rows.end(), features, features + (n - 1));
      rows.push_back(1.0);
      Y.push_back(label);
      ++m;
    }
};

struct vector_double { /* defines a vector of doubles */
    int d; /* number of elements */
    double* vec = nullptr; /* ptr to vector elements*/
    ~vector_double(){
      delete[] vec;
    }
};

struct vector_int { /* defines a vector of ints for index subsets */
    int d; /* number of elements */
    int* vec = nullptr; /* ptr to vector elements */
    ~vector_int(){
      delete[] vec;
    }
};

struct options {
    /* user options */
    double lambda; /* regularization parameter */
    double lambda_u; /* regularization parameter over unlabeled examples */
    double epsilon; /* all tolerances */
    int cgitermax; /* max iterations for CGLS */
    int mfnitermax; /* max iterations for L2_SVM_MFN */

};

class Delta { /* used in line search */
  public:
    Delta() {
      delta = 0.0;
      index = 0;
      s = 0;
    }
    ;
    double delta;
    int index;
    int s;
};
inline bool operator<(const Delta& a, const Delta& b) {
  return (a.delta < b.delta);
}

struct SolverWorkspace {
  SolverWorkspace(int m, int n);
  std::vector<double> z, q, r, p;
  std::vector<Delta> deltas;
};

/* svmlin algorithms and their subroutines */

/* Conjugate Gradient for Sparse Linear Least Squares Problems */
/* Solves: min_w 0.5*Options->lamda*w'*w + 0.5*sum_{i in Subset} Data->C[i] (Y[i]- w' x_i)^2 */
/* over a subset of examples x_i specified by vector_int Subset */
int CGLS(const AlgIn& set, const double lambda, const int cgitermax,
         const double epsilon, const vector_int& Subset,
         vector_double& Weights, vector_double& Outputs,
         double cpos, double cneg, SolverWorkspace& workspace);

/* Linear Modified Finite Newton L2-SVM*/
/* Solves: min_w 0.5*Options->lamda*w'*w + 0.5*sum_i Data->C[i] max(0,1 - Y[i] w' x_i)^2 */
int L2_SVM_MFN(const AlgIn& set, options& Options,
               vector_double& Weights,
               vector_double& Outputs, double cpos, double cneg);
double line_search(double* w, double* w_bar, double lambda, double* o,
                         double* o_bar, const double* Y, int d, int l,
                          double cpos, double cneg, Delta* deltas);
#endif
