#include "FGRadPot.hpp"
#include <cassert>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <iostream>

namespace FGRP {

//******************************************************************************
// Some helper functions:
//------------------------------------------------------------------------------
// E1(x) := -std::expint(-x)
static inline double E1(double x) {
  // return -std::expint(-x);
  // c++-17 provides expint [note definition] - but not widely available
  return gsl_sf_expint_E1(x);
}

//------------------------------------------------------------------------------
// Calculates x * std::cosh(x) - sinh(x) accounting for low-x instability
static inline double xCoshxMinusSinhx(double x) {
  // x * std::cosh(x) - sinh(x) - large cancellation for low x
  // Fix numerical instability by using expansion at low x
  // Good to at few parts in 10^10
  if ((x < 0.3) && (x > -0.3)) {
    const auto x3 = x * x * x;
    return x3 / 3.0 + (x3 * x * x) / 30.0 + (x3 * x3 * x) / 840.0;
  }
  return x * std::cosh(x) - sinh(x);
}

//******************************************************************************
double V_Uehling(double Z, double r, double rN) {
  const auto f = (Z * alpha) / (3.0 * M_PI * r);
  const auto phi = t_integral(Uehling::J_Ueh_gsl, {r, rN});
  return f * phi;
}

//------------------------------------------------------------------------------
double H_Magnetic(double Z, double r, double rN) {
  const auto f = (Z * alpha * alpha) / (4.0 * M_PI * r * r);
  const auto phi = t_integral(Magnetic::J_mag_gsl, {r, rN});
  return f * phi;
}
//------------------------------------------------------------------------------
double V_SEl(double Z, double r, double rN) {
  // Note: Doesn't include B_l(Z)
  const auto za = Z * alpha;
  const auto f = -Z * za * za * za;
  return f * SE::F_SEl(Z, r, rN);
}

//------------------------------------------------------------------------------
double V_SEh(double Z, double r, double rN) {
  // Note: Doesn't include A_l(Z)
  const auto f = -(Z * alpha) / (M_PI * r);
  const auto phi = t_integral(SE::J_SE_gsl, {r, rN, Z}, 1.0e-3);
  return f * phi;
}

//------------------------------------------------------------------------------
double V_WK(double Z, double r, double) {
  // V. V. Flambaum and J. S. M. Ginges, Phys. Rev. A 72, 052115 (2005).
  // final is rN - finite size NOT included
  const double f_eff = 2.0 / 3.0; // effective scaling
  const double za = Z * alpha;
  const double f = (2.0 * za) / (3.0 * M_PI * r);
  const double x = 0.092 * za * za;
  const double y = 1.62 * r / lam_c;
  return f_eff * f * x / (1.0 + y * y * y * y);
}

//******************************************************************************
// Function that performs the t integral
double t_integral(double (*f)(double, void *), std::vector<double> params,
                  double eps) {
  const double abs_err_lim = 1.0e-6;
  const double rel_err_lim = eps;
  const unsigned long max_num_subintvls = 1000;
  gsl_set_error_handler_off();

  gsl_function f_gsl;
  f_gsl.function = f;
  f_gsl.params = params.data();

  double result{0.0};
  double abs_err{0.0};
  gsl_integration_workspace *gsl_int_wrk =
      gsl_integration_workspace_alloc(max_num_subintvls + 1);

  double rel_err_targ = rel_err_lim;
  while (rel_err_targ < 1.0) {
    // If gsl_integration_qagiu fails to converge, answer is typically rubbish
    // So, try again, with lower threshold.
    // If never converges, set to zero
    // Not really sure why need to do this; documentation says should be OK..

    gsl_integration_qagiu(&f_gsl, 1.0, abs_err_lim, rel_err_targ,
                          max_num_subintvls, gsl_int_wrk, &result, &abs_err);

    if (std::abs(abs_err / result) < rel_err_targ ||
        std::abs(abs_err) <= abs_err_lim) {
      break;
    } else {
      result = 0.0;
      abs_err = 0.0;
      rel_err_targ *= 5.0;
    }
  }
  gsl_integration_workspace_free(gsl_int_wrk);

  return result;
}

//******************************************************************************
// Scaling factors:
double Fit::A(double Z, int l) {
  const double x = (Z - 80.0) * alpha;
  if (l > 1)
    return 0.0;
  double f = 0.0;
  for (auto i = 0ul; i < a01.size(); ++i) {
    f += a01[i] * std::pow(x, i);
  }
  return f;
}

double Fit::B(double Z, int l) {
  const double za = Z * alpha;
  if (l > 2)
    return 0.0;
  double f = 0.0;
  const auto &b = l == 2 ? b2 : b01;
  for (auto i = 0ul; i < b.size(); ++i) {
    f += b[i] * std::pow(za, i);
  }
  return f;
}

//******************************************************************************
// Helper functions:
//******************************************************************************
// Uehling
//------------------------------------------------------------------------------
double Uehling::G_Ueh(double xn, double x) {
  if (xn <= 0.0)
    return 1.0;
  const auto f1 = 3.0 / (xn * xn * xn);
  if (x >= xn)
    return f1 * xCoshxMinusSinhx(xn);
  else
    return f1 * std::exp(x) * (x - std::exp(-xn) * (1.0 + xn) * std::sinh(x));
}

//------------------------------------------------------------------------------
// Form of function required by GSL integrations; p is {R, RN}
double Uehling::J_Ueh_gsl(double t, void *p) {
  const auto params = static_cast<double *>(p);
  const auto r = params[0];
  const auto rN = params[1];
  assert(t >= 1.0);
  const double xn = 2.0 * t * rN / lam_c;
  const double x = 2.0 * t * r / lam_c;
  const auto t2 = t * t;
  const auto FUeh = (2.0 * t2 + 1.0) * std::sqrt(t2 - 1.0) / (t2 * t2);
  return std::exp(-x) * FUeh * G_Ueh(xn, x);
}

//******************************************************************************
// Magnetic self-energy:
//------------------------------------------------------------------------------
double Magnetic::G_mag(double xn, double x) {
  if (xn <= 0.0)
    return 1.0;
  const double chi = std::min(xn, x);
  const double eta = std::max(xn, x);
  const double f = 3.0 / (xn * xn * xn);
  return f * std::exp(x - eta) * xCoshxMinusSinhx(chi);
}
//------------------------------------------------------------------------------
double Magnetic::J_mag_gsl(double t, void *p) {
  assert(t > 1.0);
  const auto params = static_cast<double *>(p);
  const auto r = params[0];
  const auto rN = params[1];
  const double xn = 2.0 * t * rN / lam_c;
  const double x = 2.0 * t * r / lam_c;
  const auto t2 = t * t;

  const auto Fmag = 1.0 / (t2 * std::sqrt(t2 - 1.0));
  const double chi = std::min(xn, x);
  const double eta = std::max(xn, x);
  const auto rat = xn == 0.0 ? 1.0 : (chi / xn);
  return Fmag * ((1.0 + eta) * std::exp(-x) * G_mag(xn, x) - (rat * rat * rat));
}

//******************************************************************************
// Electric self-energy:
//------------------------------------------------------------------------------
double SE::xi(double x) { return (1.0 + x) * std::exp(-x); }

//------------------------------------------------------------------------------
double SE::F_SEl(double Z, double r, double rN) {
  if (rN <= 0.0)
    return std::exp(-Z * r);

  const auto fact = 1.5 / (r * rN * rN * rN * Z * Z);
  auto func = [=](auto x) {
    const double a1 = Z * std::abs(r - x);
    const double a2 = Z * (r + x);
    return x * (xi(a1) - xi(a2));
  };

  return fact * r_integral(func, 0.0, rN, 500);
}

//------------------------------------------------------------------------------
double SE::I1(double Z, double t) {
  const double t2 = t * t;
  const double f = 1.0 / std::sqrt(t2 - 1.0);
  const double a = 1.0 / t2 - 1.5;
  const double b = 1.0 - 0.5 / t2;
  const double arg = 1.0 / (Z * alpha) + 0.5;
  const double c = std::log(t2 - 1.0) + 4.0 * std::log(arg);
  return f * (a + b * c);
}

//------------------------------------------------------------------------------
double SE::I2(double Z, double r, double rn, double t) {

  const double rA = 0.07 * Z * Z * alpha * alpha * alpha;

  if (rn <= 0.0) {
    // NB: This is here for completeness, but
    // Rn=0.0 case actually taken care of in J_SE_gsl - much more numerically
    // stable that way (no large subtraction)
    const double x = 2.0 * t * r / alpha;
    return std::exp(-x) * (rA / (r + rA));
  }

  const auto expr = std::exp(2.0 * rA * t / alpha);
  const double fact = 1.5 * rA / (rn * rn * rn) * expr;
  const auto ttoa = 2.0 * t / alpha;

  auto func = [=](auto x) {
    const double a1 = (std::abs(r - x) + rA) * ttoa;
    const double a2 = (r + x + rA) * ttoa;
    return x * (E1(a1) - E1(a2));
  };

  return fact * r_integral(func, 0.0, rn, 500);
}

//------------------------------------------------------------------------------
double SE::J_SE_gsl(double t, void *p) {
  const auto params = static_cast<double *>(p);
  const auto r = params[0];
  const auto rN = params[1];
  const auto Z = params[2];
  const double x = 2.0 * t * r / alpha;
  const double xn = 2.0 * t * rN / alpha;

  // Note: this is the numerically unstable part:
  // Both sides cancel down to roughly std::exp(-x)
  const auto eGminusI =
      std::exp(-x) * Uehling::G_Ueh(xn, x) - SE::I2(Z, r, rN, t);

  return SE::I1(Z, t) * eGminusI;
}

} // namespace FGRP
