//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NaClUCl3FluidProperties.h"

registerMooseObject("FluidPropertiesApp", NaClUCl3FluidProperties);

InputParameters
NaClUCl3FluidProperties::validParams()
{
  InputParameters params = SinglePhaseFluidProperties::validParams();
  params.addClassDescription(
      "Fluid properties for NaCl-UCl3 salt eutectic (66.7 mol% NaCl, 33.3 mol% UCl3)");
  return params;
}

NaClUCl3FluidProperties::NaClUCl3FluidProperties(const InputParameters & parameters)
  : SinglePhaseFluidProperties(parameters)
{
}

std::string
NaClUCl3FluidProperties::fluidName() const
{
  return "NaClUCl3";
}

Real
NaClUCl3FluidProperties::molarMass() const
{
  // Average molar mass: 0.333*MW(UCl3) + 0.667*MW(NaCl)
  // MW(UCl3) ≈ 344.35 g/mol, MW(NaCl) ≈ 58.44 g/mol
  // Average ≈ 0.333×344.35 + 0.667×58.44 = 153.72 g/mol = 0.15372 kg/mol
  return 0.15372;
}

Real
NaClUCl3FluidProperties::rho_from_p_T(Real /*p*/, Real T) const
{
  // Density: ρ(T) = 4212.6 - 1.0686*T
  return 4212.6 - 1.0686 * T;
}

void
NaClUCl3FluidProperties::rho_from_p_T(
    Real p, Real T, Real & rho, Real & drho_dp, Real & drho_dT) const
{
  rho = rho_from_p_T(p, T);
  drho_dp = 0;
  drho_dT = -1.0686;
}

Real
NaClUCl3FluidProperties::v_from_p_T(Real p, Real T) const
{
  return 1.0 / rho_from_p_T(p, T);
}

void
NaClUCl3FluidProperties::v_from_p_T(Real p, Real T, Real & v, Real & dv_dp, Real & dv_dT) const
{
  Real rho = rho_from_p_T(p, T);
  v = 1.0 / rho;
  dv_dp = 0;
  dv_dT = 1.0686 / (rho * rho);
}

Real
NaClUCl3FluidProperties::cp_from_p_T(Real /*p*/, Real T) const
{
  // Specific heat capacity: c_p(T) = A + B*T + C*T^2 + D/T^2
  const Real A = 8.900439e+03;
  const Real B = -1.377936e+01;
  const Real C = 6.400369e-03;
  const Real D = -8.443758e+08;
  return A + B * T + C * T * T + D / (T * T);
}

void
NaClUCl3FluidProperties::cp_from_p_T(Real p, Real T, Real & cp, Real & dcp_dp, Real & dcp_dT) const
{
  cp = cp_from_p_T(p, T);
  dcp_dp = 0;
  const Real D = -8.443758e+08;
  // d(cp)/dT = B + 2*C*T - 2*D/T^3
  dcp_dT = -1.377936e+01 + 2.0 * 6.400369e-03 * T - 2.0 * D / (T * T * T);
}

Real
NaClUCl3FluidProperties::k_from_p_T(Real /*p*/, Real T) const
{
  // Thermal conductivity: k(T) = A + B*T + C*T^2 + D/T^2
  const Real A = 5.6820;
  const Real B = -8.7832e-03;
  const Real C = 4.0967e-06;
  const Real D = -5.7642e+05;
  return A + B * T + C * T * T + D / (T * T);
}

void
NaClUCl3FluidProperties::k_from_p_T(Real p, Real T, Real & k, Real & dk_dp, Real & dk_dT) const
{
  k = k_from_p_T(p, T);
  dk_dp = 0;
  dk_dT = -8.7832e-03 + 2.0 * 4.0967e-06 * T - 2.0 * (-5.7642e+05) / (T * T * T);
}

Real
NaClUCl3FluidProperties::mu_from_p_T(Real /*p*/, Real T) const
{
  // Dynamic viscosity: μ(T) = A * exp(Ea/(R*T))
  const Real A = 1.505e-04;
  const Real Ea = 2.666e+04;
  const Real R = 8.314;
  return A * std::exp(Ea / (R * T));
}

void
NaClUCl3FluidProperties::mu_from_p_T(Real p, Real T, Real & mu, Real & dmu_dp, Real & dmu_dT) const
{
  mu = mu_from_p_T(p, T);
  dmu_dp = 0;
  dmu_dT = -mu * (2.666e+04 / (8.314 * T * T));
}

Real
NaClUCl3FluidProperties::h_from_p_T(Real p, Real T) const
{
  // Specific enthalpy as the integral of c_p from T_mo to T:
  // h(T) = A*(T-T_{mo}) + (B/2)*(T^2-T_{mo}^2) + (C/3)*(T^3-T_{mo}^3) - D*(1/T - 1/T_{mo})
  const Real A_cp = 8.900439e+03;
  const Real B_cp = -1.377936e+01;
  const Real C_cp = 6.400369e-03;
  const Real D_cp = -8.443758e+08;
  return A_cp * (T - _T_mo) + 0.5 * B_cp * (T * T - _T_mo * _T_mo) +
         (C_cp / 3.0) * (T * T * T - _T_mo * _T_mo * _T_mo) - D_cp * (1.0 / T - 1.0 / _T_mo);
}

void
NaClUCl3FluidProperties::h_from_p_T(Real p, Real T, Real & h, Real & dh_dp, Real & dh_dT) const
{
  h = h_from_p_T(p, T);
  dh_dp = 0;
  dh_dT = cp_from_p_T(p, T);
}

Real
NaClUCl3FluidProperties::e_from_p_T(Real p, Real T) const
{
  // e = h - p*v, where v = 1/ρ
  Real h = h_from_p_T(p, T);
  Real v = v_from_p_T(p, T);
  return h - p * v;
}

void
NaClUCl3FluidProperties::e_from_p_T(Real p, Real T, Real & e, Real & de_dp, Real & de_dT) const
{
  Real h, dh_dp, dh_dT, v, dv_dp, dv_dT;
  h_from_p_T(p, T, h, dh_dp, dh_dT);
  v_from_p_T(p, T, v, dv_dp, dv_dT);
  e = h - p * v;
  de_dp = dh_dp - v - p * dv_dp;
  de_dT = dh_dT - p * dv_dT;
}

Real
NaClUCl3FluidProperties::e_from_p_rho(Real p, Real rho) const
{
  return e_from_p_T(p, T_from_p_rho(p, rho));
}

void
NaClUCl3FluidProperties::e_from_p_rho(
    Real p, Real rho, Real & e, Real & de_dp, Real & de_drho) const
{
  Real T, dT_dp, dT_drho;
  T_from_p_rho(p, rho, T, dT_dp, dT_drho);
  Real de_dp_T, de_dT;
  e_from_p_T(p, T, e, de_dp_T, de_dT);
  de_dp = de_dp_T + de_dT * dT_dp;
  de_drho = de_dT * dT_drho;
}

Real
NaClUCl3FluidProperties::T_from_v_e(Real v, Real /*e*/) const
{
  // Invert ρ = 4212.6 - 1.0686*T  with ρ = 1/v:
  return (4212.6 - 1.0 / v) / 1.0686;
}

void
NaClUCl3FluidProperties::T_from_v_e(Real v, Real e, Real & T, Real & dT_dv, Real & dT_de) const
{
  T = T_from_v_e(v, e);
  dT_de = 0;
  dT_dv = 1.0 / (1.0686 * v * v);
}

Real
NaClUCl3FluidProperties::T_from_p_rho(Real /*p*/, Real rho) const
{
  return (4212.6 - rho) / 1.0686;
}

void
NaClUCl3FluidProperties::T_from_p_rho(
    Real p, Real rho, Real & T, Real & dT_dp, Real & dT_drho) const
{
  T = T_from_p_rho(p, rho);
  dT_dp = 0;
  dT_drho = -1.0 / 1.0686;
}

Real
NaClUCl3FluidProperties::T_from_p_h(Real p, Real h) const
{
  auto lambda = [&](Real p, Real current_T, Real & new_h, Real & dh_dp, Real & dh_dT)
  { h_from_p_T(p, current_T, new_h, dh_dp, dh_dT); };
  Real T = FluidPropertiesUtils::NewtonSolve(
               p, h, _T_initial_guess, _tolerance, lambda, name() + "::T_from_p_h")
               .first;
  if (std::isnan(T))
    mooseError("Conversion from pressure and enthalpy to temperature failed to converge.");
  return T;
}

void
NaClUCl3FluidProperties::T_from_p_h(Real p, Real h, Real & T, Real & dT_dp, Real & dT_dh) const
{
  T = T_from_p_h(p, h);
  dT_dp = 0;
  Real h1, dh_dp, dh_dT;
  h_from_p_T(p, T, h1, dh_dp, dh_dT);
  dT_dh = 1.0 / dh_dT;
}

Real
NaClUCl3FluidProperties::bulk_modulus_from_p_T(Real /*p*/, Real /*T*/) const
{
  // Isentropic bulk modulus
  return 5e9;
}

Real
NaClUCl3FluidProperties::c_from_v_e(Real v, Real /*e*/) const
{
  return std::sqrt(bulk_modulus_from_p_T(0, 0) * v);
}

ADReal
NaClUCl3FluidProperties::c_from_v_e(const ADReal & v, const ADReal & /*e*/) const
{
  return std::sqrt(bulk_modulus_from_p_T(0, 0) * v);
}
