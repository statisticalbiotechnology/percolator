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
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <set>
#include <vector>
#include <ctype.h>
using namespace std;
#include "Globals.h"
#include "ssl.h"

#include <stdarg.h>
#include <cstring>
#include "Timer.h"

extern "C" {
  extern double ddot_(int *, double *, int *, double *, int *); // compute the dot product of two vectors
  extern int daxpy_(int *, double *, double *, int *, double *, int *); // compute y := alpha * x + y
  extern int dscal_(int *, double *, double *, int *); // Compute y := alpha * y
}

#define VERBOSE 1
#define LOG2(x) 1.4426950408889634*log(x)
// for compatibility issues, not using log2

AlgIn::AlgIn(const unsigned int size, const int numFeat) {
  vals = new double*[size];
  Y = new double[size];
  n = numFeat;
  positives = 0;
  negatives = 0;
}
AlgIn::~AlgIn() {
  delete[] vals;
  delete[] Y;
}

/* Buffers of one L2_SVM_MFN solve, sized for all m examples and reused by
   all its CGLS calls and line searches. */
SolverWorkspace::SolverWorkspace(int m, int n)
    : set2(static_cast<size_t>(m) * static_cast<size_t>(n)),
      z(static_cast<size_t>(m)), q(static_cast<size_t>(m)),
      r(static_cast<size_t>(n)), p(static_cast<size_t>(n)),
      deltas(static_cast<size_t>(m)) {}

/* q = set2 * p for the active rows of set2 (row length n). Gives the same
   result as the reference BLAS dgemv with trans='T', alpha=1 and beta=0,
   which this replaces: every row is summed in the same order, but four rows
   are processed together so that their additions overlap instead of waiting
   on each other. */
static void activeRowsTimesVector(int active, const double* set2, int n,
                                  const double* p, double* q) {
  int i = 0;
  for (; i + 4 <= active; i += 4) {
    const double* x0 = set2 + static_cast<size_t>(i) * static_cast<size_t>(n);
    const double* x1 = x0 + n;
    const double* x2 = x1 + n;
    const double* x3 = x2 + n;
    double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0;
    for (int k = 0; k < n; k++) {
      t0 += x0[k] * p[k];
      t1 += x1[k] * p[k];
      t2 += x2[k] * p[k];
      t3 += x3[k] * p[k];
    }
    // dgemv adds the dot product to the zeroed output vector
    q[i] = 0.0 + t0;
    q[i + 1] = 0.0 + t1;
    q[i + 2] = 0.0 + t2;
    q[i + 3] = 0.0 + t3;
  }
  for (; i < active; i++) {
    const double* x = set2 + static_cast<size_t>(i) * static_cast<size_t>(n);
    double t = 0.0;
    for (int k = 0; k < n; k++) {
      t += x[k] * p[k];
    }
    q[i] = 0.0 + t;
  }
}

double cglsFun1(int active, int* J, const double* Y,
                double* set2, int n, double* q, 
                double* p, double cpos, double cneg){
  double omega_q = 0.0;
  int i = 0;

  activeRowsTimesVector(active, set2, n, p, q);

  for (i = 0; i < active; i++) {
    omega_q += ((Y[J[i]]==1)? cpos : cneg) * (q[i]) * (q[i]);
  }

  return(omega_q);
}

void cglsFun2(int active, int* J, const double* Y,
              double* set2, int n0, int n, double* q, 
              double* o, double* z, double* r, 
              double cpos, double cneg){
  int i;
  int inc = 1;
  
  for (i = 0; i < active; i++) {
    o[J[i]] += q[i];
    z[i] -= ((Y[J[i]]==1)? cpos : cneg) * q[i];
    daxpy_(&n, &(z[i]),
           set2 + static_cast<size_t>(i) * static_cast<size_t>(n), &inc, r,
           &inc);
  }
}

int CGLS(const AlgIn& data, const double lambda, const int cgitermax,
         const double epsilon, const vector_int& Subset,
         vector_double& Weights, vector_double& Outputs,
         double cpos, double cneg, SolverWorkspace& workspace) {
  if (VERBOSE_CGLS) {
    cout << "CGLS starting..." << endl;
  }
  /* Disassemble the structures */
  Timer tictoc;
  int active = Subset.d;
  int* J = Subset.vec;
  double** set = data.vals;
  const double* Y = data.Y;
  int n = data.n;
  double* beta = Weights.vec;
  double* o = Outputs.vec;
  // initialize z
  double* z = workspace.z.data();
  double* q = workspace.q.data();
  int ii = 0;
  int i;
  int n0 = n-1;
  int inc = 1;
  double one = 1;
  double negLambda = -lambda;
  size_t rowStart = 0;
  // contiguous copy of the active rows, with the intercept as last column
  double* set2 = workspace.set2.data();
  double* r = workspace.r.data();
  for (i = n; i--;) {
    r[i] = 0.0;
  }
  for (i = 0; i < active; i++) {
    ii = J[i];
    z[i] = ((Y[ii]==1)? cpos : cneg) * (Y[ii] - o[ii]);
    rowStart = static_cast<size_t>(i) * static_cast<size_t>(n);
    memcpy(set2 + rowStart, set[ii], sizeof(double)*static_cast<std::size_t>(n0));
    set2[rowStart + n0] = 1.0;
    daxpy_(&n, &(z[i]), set2 + rowStart, &inc, r, &inc);
  }
  double* p = workspace.p.data();
  daxpy_(&n, &negLambda, beta, &inc, r, &inc);
  memcpy(p, r, sizeof(double)*static_cast<std::size_t>(n));
  double omega1 = ddot_(&n, r, &inc, r, &inc);
  double omega_p = omega1;
  double omega_q = 0.0;
  double inv_omega2 = 1 / omega1;
  double scale = 0.0;
  double omega_z = 0.0;
  double gamma = 0.0;
  int cgiter = 0;
  int optimality = 0;
  double epsilon2 = epsilon * epsilon;
  // iterate
  while (cgiter < cgitermax) {
    cgiter++;
    omega_q = cglsFun1(active, J, Y, set2, n, q, p, cpos, cneg);
    gamma = omega1 / (lambda * omega_p + omega_q);
    inv_omega2 = 1 / omega1;

    memcpy(r, beta, sizeof(double)*static_cast<std::size_t>(n));
    dscal_(&n, &negLambda, r, &inc);

    daxpy_(&n, &gamma, p, &inc, beta, &inc);
    dscal_(&active, &gamma, q, &inc);

    cglsFun2(active, J, Y, set2,
             n0, n, q, o, z, r, cpos, cneg);

    omega_z = ddot_(&active, z, &inc, z, &inc);
    omega1 = ddot_(&n, r, &inc, r, &inc);

    if (VERBOSE_CGLS) {
      cout << "..." << cgiter << " ( " << omega1 << " )";
    }

    if (omega1 < epsilon2 * omega_z) {
      optimality = 1;
      break;
    }

    scale = omega1 * inv_omega2;
    dscal_(&n, &scale, p, &inc);
    daxpy_(&n, &one, r, &inc, p, &inc);
    omega_p = ddot_(&n, p, &inc, p, &inc);
  }
  if (VERBOSE_CGLS) {
    cout << "...Done." << endl;
  }
  tictoc.stop();
  if (VERB > 4) {
    cerr << "CGLS converged in " << cgiter << " iteration(s) and "
        << tictoc.getCPUTimeStr() << " CPU seconds." << endl;
  }
  return optimality;
}

int L2_SVM_MFN(const AlgIn& data, options& Options,
               vector_double& Weights,
               vector_double& Outputs, double cpos, double cneg) {
  /* Disassemble the structures */
  Timer tictoc;
  double** set = data.vals;
  const double* Y = data.Y;
  int n = Weights.d;
  const int m = data.m;
  double lambda = Options.lambda;
  double epsilon = BIG_EPSILON;
  int cgitermax = SMALL_CGITERMAX;
  double* w = Weights.vec;
  double* o = Outputs.vec;
  double F_old = 0.0;
  double F = 0.0;
  double diff = 0.0;
  int ini = 0;
  int n0 = n-1;
  int inc = 1;
  vector_int ActiveSubset;
  ActiveSubset.vec = new int[m];
  ActiveSubset.d = m;
  // initialize
  F = 0.5 * lambda * ddot_(&n, w, &inc, w, &inc);
  int active = 0;
  int inactive = m - 1; // l-1
  for (int i = 0; i < m; i++) {
    diff = 1 - Y[i] * o[i];
    if (diff > 0) {
      ActiveSubset.vec[active] = i;
      active++;
      // C[i]
      F += 0.5 * ((Y[i]==1)? cpos : cneg) * diff * diff;
    } else {
      ActiveSubset.vec[inactive] = i;
      inactive--;
    }
  }
  ActiveSubset.d = active;
  int iter = 0;
  int opt = 0;
  int opt2 = 0;
  vector_double Weights_bar;
  vector_double Outputs_bar;
  double* w_bar = new double[static_cast<std::size_t>(n)];
  double* o_bar = new double[static_cast<std::size_t>(m)];
  Weights_bar.vec = w_bar;
  Outputs_bar.vec = o_bar;
  Weights_bar.d = n;
  Outputs_bar.d = m;
  double delta = 0.0;
  int ii = 0;
  SolverWorkspace workspace(m, n);
  while (iter < Options.mfnitermax) {
    iter++;
    if (VERB > 4) {
      cerr << "L2_SVM_MFN Iteration# " << iter << " (" << active
          << " active examples, " << " objective_value = " << F << ")"
          << endl;
    }
    memcpy(w_bar, w, sizeof(double)*static_cast<std::size_t>(n));
    memcpy(o_bar, o, sizeof(double)*static_cast<std::size_t>(m));
    opt = CGLS(data,
               lambda,
               cgitermax,
               epsilon,
               ActiveSubset,
               Weights_bar,
               Outputs_bar, cpos, cneg, workspace);
    for (int i = active; i < m; i++) {
      ii = ActiveSubset.vec[i];
      o_bar[ii] = ddot_(&n0, set[ii], &inc, w_bar, &inc) + w_bar[n - 1];
    }
    if (ini == 0) {
      cgitermax = CGITERMAX;
      ini = 1;
    };
    opt2 = 1;
    for (int i = 0; i < m; i++) {
      ii = ActiveSubset.vec[i];
      if (i < active) {
        opt2 = (opt2 && (Y[ii] * o_bar[ii] <= 1 + epsilon));
      } else {
        opt2 = (opt2 && (Y[ii] * o_bar[ii] >= 1 - epsilon));
      }
      if (opt2 == 0) {
        break;
      }
    }
    if (opt && opt2) { // l
      if (epsilon == BIG_EPSILON) {
        epsilon = Options.epsilon;
        if (VERB > 4) {
          cerr << "  epsilon = " << BIG_EPSILON
            << " case converged (speedup heuristic 2). Continuing with epsilon="
            << EPSILON << endl;
        }
        continue;
      } else {
        memcpy(w, w_bar, sizeof(double)*static_cast<std::size_t>(n));
        memcpy(o, o_bar, sizeof(double)*static_cast<std::size_t>(m));
        if (VERB > 3) {
          tictoc.stop();
          cerr << "L2_SVM_MFN converged (optimality) in " << iter
              << " iteration(s) and " << tictoc.getCPUTimeStr() << " CPU seconds. \n"
              << endl;
        }
        return 1;
      }
    }
    delta = line_search(w, w_bar, lambda, o, o_bar, Y, n, m, cpos, cneg,
                        workspace.deltas.data());
    F_old = F;
    double delta2 = 1-delta;
    dscal_(&n, &delta2, w, &inc);
    daxpy_(&n, &delta, w_bar, &inc, w, &inc);
    F = 0.5 * lambda * ddot_(&n, w, &inc, w, &inc);
    active = 0;
    inactive = m - 1;
    for (int i = 0; i < m; i++) {
      o[i] += delta * (o_bar[i] - o[i]);
      diff = 1 - Y[i] * o[i];
      if (diff > 0) {
        ActiveSubset.vec[active] = i;
        active++;
        F += 0.5 * ((Y[i]==1)? cpos : cneg) * diff * diff;
      } else {
        ActiveSubset.vec[inactive] = i;
        inactive--;
      }
    }
    ActiveSubset.d = active;
    if (fabs(F - F_old) < RELATIVE_STOP_EPS * fabs(F_old)) {
      return 2;
    }
  }
  return 0;
}

double line_search(double* w, double* w_bar, double lambda, double* o,
                   double* o_bar, const double* Y, int d, /* data dimensionality -- 'n' */
                   int l, double cpos, double cneg, Delta* deltas){
  int i = 0;
  double omegaL = 0.0;
  double omegaR = 0.0;
  double diff = 0.0;
  for (int i = d; i--;) {
    diff = w_bar[i] - w[i];
    omegaL += w[i] * diff;
    omegaR += w_bar[i] * diff;
  }
  omegaL = lambda * omegaL;
  omegaR = lambda * omegaR;
  double L = omegaL;
  double R = omegaR;
  int ii = 0;
  double d2 = 0.0;

  int p = 0;
  for (i = 0; i < l; i++) {
    diff = Y[i] * (o_bar[i] - o[i]);
    if (Y[i] * o[i] < 1) {
      d2 = ((Y[i]==1)? cpos : cneg) * (o_bar[i] - o[i]);
      L += (o[i] - Y[i]) * d2;
      R += (o_bar[i] - Y[i]) * d2;
      if (diff > 0) {
        deltas[p].delta = (1 - Y[i] * o[i]) / diff;
        deltas[p].index = i;
        deltas[p].s = -1;
        p++;
      }
    } else {
      if (diff < 0) {
        deltas[p].delta = (1 - Y[i] * o[i]) / diff;
        deltas[p].index = i;
        deltas[p].s = 1;
        p++;
      }
    }
  }
  sort(deltas, deltas + p);
  double delta_prime = 0.0;
  for (i = 0; i < p; i++) {
    delta_prime = L + deltas[i].delta * (R - L);
    if (delta_prime >= 0) {
      break;
    }
    ii = deltas[i].index;
    diff = (deltas[i].s) * ((Y[ii]==1)? cpos : cneg) * (o_bar[ii] - o[ii]);
    L += diff * (o[ii] - Y[ii]);
    R += diff * (o_bar[ii] - Y[ii]);
  }
  return (-L / (R - L));
}

