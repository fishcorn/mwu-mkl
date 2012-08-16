#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>

#include "mwu_main.h"

namespace {

  inline double access(double * K, int r, int c, int l, int n) {
    return K[l*n*n + c*n + r];
  }

  inline double * get_column(double * K, int c, int l, int n) {
    return K + l*n*n + c*n;
  }

  struct primal_var 
  {
    primal_var(int m_, int n_) 
      : x11(m_,0.0), x12(m_,0.0), x22(m_,0.0), ea(0.0), m(m_), n(n_)
    {}

    void normalize(double R)
    {
      double trace = std::accumulate(x11.begin(), x11.end(), 0.0, std::plus<double>())
        + std::accumulate(x22.begin(), x22.end(), 0.0, std::plus<double>())
        + m*(n-1)*ea;
      trace /= R;
      for( unsigned int i = 0; i < x11.size(); ++i ) {
        x11[i] /= trace;
        x12[i] /= trace;
        x22[i] /= trace;
      }
      ea /= trace;
    }

    std::vector<double> x11, x12, x22;
    double ea;
    const int m, n;
  };

  bool is_pos(double x) { return (x > 0.0); }

  bool is_neg(double x) { return (x < 0.0); }

  bool oracle(std::vector<double> & yperp0,  // OUTPUT
              const std::vector<double> & g, // INPUT
              int * y,                       // INPUT
              int n,                         // INPUT
              int m)                         // INPUT
  {
    int pidx = std::find_if(y, y + n, &is_pos) - y;
    int nidx = std::find_if(y, y + n, &is_neg) - y;
    double pmax = g[pidx], nmax = g[nidx];

    for( int j = 0; j < n; ++j ) {
      if (y[j] > 0) {
        if (g[j] > pmax) {
          pidx = j;
          pmax = g[j];
        }
      }
      if (y[j] < 0) {
        if (g[j] > nmax) {
          nidx = j;
          nmax = g[j];
        }
      }
    }

    yperp0[pidx] = 0.5;
    yperp0[nidx] = 0.5;
      
    double gyy = (pmax + nmax)/2;

    return (gyy >= -2);
  }

  void exponentiateM(primal_var & X,                     // OUTPUT
                     std::vector<double> & g,            // OUTPUT
                     const std::vector<double> & alGal,  // INPUT
                     const std::vector<double> & Galpha, // INPUT
                     const int n,                        // INPUT
                     const int m,                        // INPUT
                     const double rho,                   // INPUT
                     const double epsp,                  // INPUT
                     const int t,                        // INPUT
                     const double R,                     // INPUT
                     const double cutoff)                // INPUT
  {
    std::vector<double> normu(m);
    std::transform(alGal.begin(), alGal.end(), normu.begin(), (double (*)(double)) &std::sqrt);

    double coeff = -epsp/(2*rho); // always negative

    std::vector<double> ps(m,0.0);

    for( int i = 0; i < m; ++i ) {
      ps[i] = std::abs(coeff)*normu[i];
    }

    // For large x, sinh(x) and cosh(x) are essentially exp(x) 
    // -- this will also overflow for large x, so quash it
    double quash = *std::max_element(ps.begin(), ps.end());
    if (quash > cutoff) {
      for( int i = 0; i < m; ++i ) {
        double epsquash = std::exp(ps[i] - quash)/2;
        X.x11[i] = X.x22[i] = epsquash;
        X.x12[i] = -epsquash;
      }
      X.ea = std::exp(-quash); // probably insignificant
    }
    else {
      // Factoring out exp(ph), as it gets factored out by the normalize anyway
      for( int i = 0; i < m; ++i ) {
        X.x11[i] = X.x22[i] = std::cosh(ps[i]);
        X.x12[i] = -std::sinh(ps[i]);
      }
      X.ea = 1;
    }

    X.normalize(R);
    for( int j = 0; j < n; ++j ) {
      g[j] = 0;
      for( int i = 0; i < m; ++i ) {
        g[j] += 2*Galpha[i*n + j]*X.x12[i]/normu[i];
      }
    }
  }

  bool try_solve(double * alpha,  // OUTPUT: Support vector
                 primal_var & X,  // OUTPUT: primal variable
                 double * K,      // INPUT:  Kernel as one-dim array (columns, then layers)
                 int * y,         // INPUT:  Labels, +/-1
                 double c,        // INPUT:  Desired output trace
                 const std::vector<double> & r,
                                  // INPUT:  Kernel traces
                 int n,           // INPUT:  Number of data points
                 int m,           // INPUT:  Number of kernels
                 double eps,      // INPUT:  Epsilon parameter
                 double ratio,    // INPUT:  Iteration multiplier
                 double cutoff,   // INPUT:  Exponentiation cutoff
                 double C,        // INPUT:  Soft margin parameter
                 double norm1or2, // INPUT:  Is the soft margin 1-norm (1) or 2-norm(2) 
                 //                          or is it a hard margin (0)
                 int verbose      // INPUT:  Be noisy or not (boolean)
                 )
  {
    double R = 2.0;
    double rho = sqrt(2+c)/2;
    double eps0 = eps/(2*rho);
    double epsp = -log1p(-eps0);

    const int tau = std::ceil(ratio*2*std::log(n)/(eps0*eps0));

    std::vector<double> g(n,0.0);
    std::vector<double> Galpha(n*m,0.0);
    std::vector<double> alGal(m,0.0);

    for( int t = 0; t < tau; ++t) {
      std::vector<double> yperp0(n,0.0);
      if (!oracle(yperp0, g, y, n, m)) {
        alpha[0] = t/tau;
        alpha[1] = t;
        return false;
      }

      std::vector<double>::iterator yidx1 = find_if(yperp0.begin(), yperp0.end(), &is_pos);
      std::vector<double>::iterator yidx2 = find_if(yidx1 + 1, yperp0.end(), &is_pos);
      int j1 = yidx1 - yperp0.begin();
      int j2 = yidx2 - yperp0.begin();
      double yp01 = *yidx1, yp02 = *yidx2;

      alpha[j1] += yp01;
      alpha[j2] += yp02;

      for( int i = 0; i < m; ++i ) {
        double * Kij1 = get_column(K, j1, i, n);
        double * Kij2 = get_column(K, j2, i, n);
        for( int k = 0; k < n; ++k ) {
          Galpha[i*n + k] +=  (c * Kij1[k] / r[i]) * yp01 * y[j1] * y[k];
          Galpha[i*n + k] +=  (c * Kij2[k] / r[i]) * yp02 * y[j2] * y[k];
        }
        if (norm1or2 == 2) {
          Galpha[i*n + j1] += yp01 / C;
          Galpha[i*n + j2] += yp02 / C;
        }

        alGal[i] += Galpha[i*n + j1] * yp01;
        alGal[i] += Galpha[i*n + j2] * yp02;
      }

      exponentiateM(X, g, alGal, Galpha, n, m, rho, epsp, t, R, cutoff);
    }

    for( int i = 0; i < n; ++i ) {
      alpha[i] /= tau;
    }

    return true;
  }
}

void run_mwu_cpp(double * Sigma,  /* OUTPUT: The kernel weights                             */
                 double * alpha,  /* OUTPUT: Support vector                                 */
                 double * bsvm,   /* OUTPUT: Bias                                           */
                 int * posw,      /* OUTPUT: Support indicators                             */
                 double * K,      /* INPUT:  Kernel as one-dim array (columns, then layers) */
                 int * y,         /* INPUT:  Labels, +/-1                                   */
                 int n,           /* INPUT:  Number of data points                          */
                 int m,           /* INPUT:  Number of kernels                              */
                 double eps,      /* INPUT:  Epsilon parameter                              */
                 double ratio,    /* INPUT:  Iteration multiplier                           */
                 double cutoff,   /* INPUT:  Exponentiation cutoff                          */
                 double C,        /* INPUT:  Margin parameter                               */
                 int norm1or2,    /* INPUT:  Is the soft margin 1-norm (1) or 2-norm (2)    */
                 /*                          or is it a hard margin (0)                     */
                 int verbose      /* INPUT:  Be noisy or not (boolean)                      */
                 )
{
  // compute traces of all the kernels
  std::vector<double> r(m,0);
  for( int i = 0; i < m; ++i ) {
    for( int j = 0; j < n; ++j ) {
      r[i] += access(K, j, j, i, n);
    }
  }

  primal_var X(m, n);
  double c = m;

  if (!try_solve(alpha, X, K, y, c, r, n, m, eps, ratio, cutoff, C, norm1or2, verbose)) {
    Sigma[0] = -1e-300;
    return;
  }

  // compute posw
  for( int j = 0; j < n; ++j ) {
    if (alpha[j] != 0)
      posw[j] = 1;
  }

  // compute Sigma
  double Sigma_sum = std::accumulate(X.x12.begin(), X.x12.end(), 0.0, std::plus<double>());
  for( int i = 0; i < m; ++i ) {
    Sigma[i] = c * X.x12[i] / (Sigma_sum * r[i]);
  }

  // compute the values of the support points
  std::vector<double> suppval(n,0);
  for( int i = 0; i < m; ++i ) {
    for( int j = 0; j < n; ++j ) {
      if (posw[j] == 0)
        continue;
      double * Kij = get_column(K, j, i, n);
      for( int k = 0; k < n; ++k ) {
        if (posw[k] == 0)
          continue;
        suppval[k] += Sigma[i] * Kij[k] * alpha[j] * y[j];
      }
    }
  }

  // find the averages over the positive and negative supports
  double pavg = 0.0, navg = 0;
  int pcount = 0, ncount = 0;
  for (int j = 0; j < n; ++j) {
    if (posw[j] == 0)
      continue;
    if (y[j] > 0) {
      ++pcount;
      pavg += suppval[j];
    }
    if (y[j] < 0) {
      ++ncount;
      navg += suppval[j];
    }
  }
  pavg /= pcount;
  navg /= ncount;

  // compute final alpha and bsvm
  double scale = std::abs(pavg - navg)/2;
  *bsvm = -(pavg + navg)/(2 * scale);
  for (int j = 0; j < n; ++j) {
    alpha[j] /= scale;
  }

  return;
}
