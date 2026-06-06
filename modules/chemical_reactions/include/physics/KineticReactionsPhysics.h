//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ReactionNetworkPhysicsBase.h"

/**
 * Creates the equations for a kinetic (mass-action) diffusion-reaction network.
 *
 * For each species s a continuous Galerkin equation is formed:
 *   dC_s/dt = div(D_s grad C_s) + S_s(C)
 * where the net production S_s is the sum over every reaction j of
 *   nu_{s,j} * k_j * prod_r C_r^{a_{r,j}}
 * with nu_{s,j} the net stoichiometry of s in reaction j (products minus
 * reactants), k_j the forward rate constant, and a_{r,j} the stoichiometry of
 * reactant r. Because the reactant stoichiometry appears both as the rate-law
 * exponent and (for a reactant) as the net-stoichiometry multiplier, the
 * factor-of-2 convention for like-species reactions (e.g. 2A -> P giving
 * dA/dt = -2 k [A]^2) is handled automatically.
 *
 * The equations are assembled entirely from existing framework objects:
 *   - ADTimeDerivative           (time derivative)
 *   - ADMatDiffusion             (Fickian diffusion, obeys the mesh coord_type)
 *   - ADParsedMaterial           (the net reaction rate S_s as a functor)
 *   - ADMatBodyForce             (applies -S_s to the residual)
 *   - ADGenericConstantMaterial  (the per-species diffusivities)
 */
class KineticReactionsPhysics : public ReactionNetworkPhysicsBase
{
public:
  static InputParameters validParams();

  KineticReactionsPhysics(const InputParameters & parameters);

private:
  virtual void addFEKernels() override;
  virtual void addMaterials() override;

  /// Per-species diffusion coefficients, ordered as the 'species' parameter
  const std::vector<Real> _diffusivities;
  /// Forward rate constant of each reaction (from the [k=...] metadata)
  std::vector<Real> _rate_constants;
  /// Net stoichiometry of each species in each reaction (products - reactants)
  std::vector<std::vector<Real>> _net_stoichiometry;
  /// FParser expression for the net reaction rate S_s of each species
  std::vector<std::string> _rate_expression;
  /// Species coupled into the rate expression S_s of each species
  std::vector<std::vector<VariableName>> _rate_coupled_species;

  /// Material property name holding the diffusivity of a species
  std::string diffusivityName(const VariableName & species) const { return "D_" + species; }
  /// Material property name holding the net reaction rate of a species
  std::string rateName(const VariableName & species) const { return "S_" + species; }
};
