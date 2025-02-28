//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NaClUCl3FluidPropertiesTest.h"
#include "SinglePhaseFluidPropertiesTestUtils.h"
#include <cmath>

/**
 * Verify calculation of the NaCl–UCl₃ fluid properties.
 *
 * The fluid property correlations used are:
 *
 * - Density:
 *   \[
 *   \rho(T)=4212.6-1.0686\,T
 *   \]
 *
 * - Specific heat capacity:
 *   \[
 *   c_p(T)=8.900439\times10^{3}-1.377936\times10^{1}\,T+6.400369\times10^{-3}\,T^2-\frac{8.443758\times10^{8}}{T^2}
 *   \]
 *
 * - Thermal conductivity:
 *   \[
 *   k(T)=5.6820-8.7832\times10^{-3}\,T+4.0967\times10^{-6}\,T^2-\frac{5.7642\times10^{5}}{T^2}
 *   \]
 *
 * - Dynamic viscosity:
 *   \[
 *   \mu(T)=1.505\times10^{-4}\,\exp\!\left(\frac{2.666\times10^{4}}{8.314\,T}\right)
 *   \]
 *
 * - Specific enthalpy is computed as:
 *   \[
 *   h(T)=A\,(T-T_{mo})+\frac{B}{2}\,(T^2-T_{mo}^2)+\frac{C}{3}\,(T^3-T_{mo}^3)-D\left(\frac{1}{T}-\frac{1}{T_{mo}}\right)
 *   \]
 *   with \(A=8.900439\times10^{3}\), \(B=-1.377936\times10^{1}\), \(C=6.400369\times10^{-3}\),
 *   \(D=-8.443758\times10^{8}\) and \(T_{mo}=796.15\,\mathrm{K}\).
 *
 * - The internal energy is defined as \(e=h-\frac{p}{\rho}\).
 *
 * The test loops over pressures \(\{1\times10^{5}, 1\times10^{6}, 5\times10^{6}\}\) Pa and
 * temperatures \(\{850,900,950\}\) K.
 */
TEST_F(NaClUCl3FluidPropertiesTest, properties)
{
  const Real tol = REL_TOL_SAVED_VALUE;
  const std::vector<Real> pressures = {1e5, 1e6, 5e6};
  const std::vector<Real> temperatures = {850., 900., 950.};

  // Constants for cp integration
  const Real T_mo = 796.15;
  const Real A_cp = 8.900439e+03;
  const Real B_cp = -1.377936e+01;
  const Real C_cp = 6.400369e-03;
  const Real D_cp = -8.443758e+08;

  for (size_t ip = 0; ip < pressures.size(); ip++)
  {
    const Real p = pressures[ip];
    for (size_t iT = 0; iT < temperatures.size(); iT++)
    {
      const Real T = temperatures[iT];

      // Compute reference density:
      const Real rho_ref = 4212.6 - 1.0686 * T;
      // Compute reference specific heat:
      const Real cp_ref =
          8.900439e+03 + (-1.377936e+01) * T + 6.400369e-03 * T * T + (-8.443758e+08) / (T * T);
      // Compute reference thermal conductivity:
      const Real k_ref = 5.6820 + (-8.7832e-03) * T + 4.0967e-06 * T * T + (-5.7642e+05) / (T * T);
      // Compute reference dynamic viscosity:
      const Real mu_ref = 1.505e-04 * std::exp(2.666e+04 / (8.314 * T));
      // Compute reference specific enthalpy by integrating cp:
      const Real h_ref = A_cp * (T - T_mo) + 0.5 * B_cp * (T * T - T_mo * T_mo) +
                         (C_cp / 3.0) * (T * T * T - T_mo * T_mo * T_mo) -
                         D_cp * (1.0 / T - 1.0 / T_mo);
      // Compute reference specific internal energy:
      const Real v_ref = 1.0 / rho_ref;
      const Real e_ref = h_ref - p * v_ref;

      // Needed to compute other properties
      const Real v = _fp->v_from_p_T(p, T);
      const Real rho = _fp->rho_from_p_T(p, T);
      const Real e = _fp->e_from_p_T(p, T);
      const Real h = _fp->h_from_p_T(p, T);

      // Density and specific volume:
      REL_TEST(_fp->rho_from_p_T(p, T), rho_ref, tol);
      REL_TEST(_fp->v_from_p_T(p, T), 1.0 / rho_ref, tol);

      // Specific enthalpy:
      REL_TEST(_fp->h_from_p_T(p, T), h_ref, tol);

      // Specific internal energy:
      REL_TEST(_fp->e_from_p_T(p, T), e_ref, tol);
      REL_TEST(_fp->e_from_p_rho(p, rho), e_ref, tol);

      // Temperature (from v,e, from p,ρ, and from p,h):
      REL_TEST(_fp->T_from_v_e(v, e), T, tol);
      REL_TEST(_fp->T_from_p_rho(p, rho), T, tol);
      REL_TEST(_fp->T_from_p_h(p, h), T, tol);

      // Thermal conductivity:
      REL_TEST(_fp->k_from_p_T(p, T), k_ref, tol);

      // Isobaric specific heat:
      REL_TEST(_fp->cp_from_p_T(p, T), cp_ref, tol);

      // Dynamic viscosity:
      REL_TEST(_fp->mu_from_p_T(p, T), mu_ref, tol);
    }
  }
}

/**
 * Verify the derivatives of the NaCl–UCl₃ properties via finite differences.
 * The test temperature is chosen within the valid range (e.g., 900 K).
 */
TEST_F(NaClUCl3FluidPropertiesTest, derivatives)
{
  const Real tol = REL_TOL_DERIVATIVE;
  const Real p = 3e6; // 3 MPa
  const Real T = 900; // 900 K
  const Real rho = _fp->rho_from_p_T(p, T);
  const Real v = 1.0 / rho;
  const Real h = _fp->h_from_p_T(p, T);
  const Real e = _fp->e_from_p_T(p, T);

  DERIV_TEST(_fp->rho_from_p_T, p, T, tol);
  DERIV_TEST(_fp->e_from_p_T, p, T, tol);
  DERIV_TEST(_fp->v_from_p_T, p, T, tol);
  DERIV_TEST(_fp->h_from_p_T, p, T, tol);
  DERIV_TEST(_fp->k_from_p_T, p, T, tol);
  DERIV_TEST(_fp->cp_from_p_T, p, T, tol);
  DERIV_TEST(_fp->mu_from_p_T, p, T, tol);

  DERIV_TEST(_fp->T_from_v_e, v, e, tol);

  DERIV_TEST(_fp->T_from_p_rho, p, rho, tol);
  DERIV_TEST(_fp->e_from_p_rho, p, rho, tol);
  DERIV_TEST(_fp->T_from_p_h, p, h, tol);
}
