//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KineticReactionsPhysics.h"

#include <sstream>
#include <iomanip>
#include <set>

registerMooseAction("ChemicalReactionsApp", KineticReactionsPhysics, "add_kernel");
registerMooseAction("ChemicalReactionsApp", KineticReactionsPhysics, "add_material");
registerReactionNetworkPhysicsBaseTasks("ChemicalReactionsApp", KineticReactionsPhysics);

namespace
{
/// Format a Real with enough precision for an FParser expression
std::string
fp(Real v)
{
  std::ostringstream os;
  os << std::setprecision(12) << v;
  return os.str();
}
}

InputParameters
KineticReactionsPhysics::validParams()
{
  InputParameters params = ReactionNetworkPhysicsBase::validParams();
  params.addClassDescription(
      "Forms the equations for a kinetic (mass-action) diffusion-reaction network using a "
      "continuous Galerkin finite element discretization. Each reaction proceeds at a rate "
      "k * product(C_reactant^stoichiometry), and the per-species net production is assembled "
      "from existing framework kernels and materials.");

  params.renameParam("solver_variables", "species", "The list of species to solve for");
  // A kinetic network solves for every species; auxiliary output species are optional
  params.makeParamNotRequired("auxiliary_variables");
  params.set<std::vector<AuxVariableName>>("auxiliary_variables") = {};
  params.renameParam("variable_order", "order", "Order of the species variables");
  params.renameParam("equation_scaling", "scaling", "");
  params.setDocString("reactions",
                      "The list of kinetic reactions, one per line. Each reaction must carry its "
                      "forward rate constant in its metadata, e.g. 'e_aq + Hp -> H [k=1.7e10]'. "
                      "Use species names matching 'species' (no +/- charge notation).");

  params.addRequiredParam<std::vector<Real>>(
      "diffusivities", "Diffusion coefficient of each species, in the same order as 'species'");

  return params;
}

KineticReactionsPhysics::KineticReactionsPhysics(const InputParameters & parameters)
  : ReactionNetworkPhysicsBase(parameters),
    _diffusivities(getParam<std::vector<Real>>("diffusivities")),
    _rate_constants(_num_reactions, 0.0),
    _rate_expression(_num_solver_species),
    _rate_coupled_species(_num_solver_species)
{
  if (_diffusivities.size() != _num_solver_species)
    paramError("diffusivities",
               "Expected one diffusivity per species (" + std::to_string(_num_solver_species) +
                   "), got " + std::to_string(_diffusivities.size()));

  // Check that every species named in a reaction is a known species
  std::set<VariableName> known(_solver_species.begin(), _solver_species.end());
  known.insert(_aux_species.begin(), _aux_species.end());

  // Build, per reaction, the mass-action rate expression and its reactant list
  std::vector<std::string> reaction_rate_expr(_num_reactions);
  std::vector<std::vector<VariableName>> reaction_reactants(_num_reactions);
  std::vector<std::map<VariableName, Real>> net_stoich(_num_reactions);

  for (const auto j : index_range(_reactions))
  {
    const auto & reaction = _reactions[j];

    if (!reaction.hasMetaData<Real>("k"))
      paramError("reactions",
                 "Reaction: '" + _reactions_input[j] +
                     "'\n is missing a forward rate constant [k=...] in its metadata");
    _rate_constants[j] = std::stod(reaction.getMetaData("k"));

    // rate = k * product over reactant terms of (species ^ stoichiometry)
    std::string rexpr = fp(_rate_constants[j]);
    for (const auto & term : reaction.reactants)
    {
      if (!known.count(term.species))
        paramError("reactions",
                   "Reaction: '" + _reactions_input[j] + "'\n references unknown species '" +
                       term.species +
                       "'. It must be listed in 'species' or 'auxiliary_variables'.");
      reaction_reactants[j].push_back(term.species);
      rexpr += "*" + term.species;
      if (term.coefficient != 1.0)
        rexpr += "^" + fp(term.coefficient);
    }
    reaction_rate_expr[j] = rexpr;
    net_stoich[j] = reaction.getUniqueStoichiometricCoefficients();
    // Untracked products (e.g. the solvent H2O) are allowed: they simply get
    // no equation. Only reactants must be variables, since they drive the rate
    // and are coupled into the parsed material (checked above).
  }

  // Assemble, per species, the net production expression S_s = sum_j nu_{s,j} * rate_j
  for (const auto i : index_range(_solver_species))
  {
    const auto & species = _solver_species[i];
    std::set<VariableName> coupled;

    for (const auto j : make_range(_num_reactions))
    {
      const auto it = net_stoich[j].find(species);
      if (it == net_stoich[j].end() || it->second == 0.0)
        continue;
      const Real nu = it->second;

      // nu * rate_j ; fp() carries the sign of nu
      const std::string term = fp(nu) + "*" + reaction_rate_expr[j];
      if (_rate_expression[i].empty())
        _rate_expression[i] = term;
      else
        _rate_expression[i] += (nu > 0.0 ? "+" : "") + term;

      for (const auto & r : reaction_reactants[j])
        coupled.insert(r);
    }
    _rate_coupled_species[i].assign(coupled.begin(), coupled.end());
  }
}

void
KineticReactionsPhysics::addMaterials()
{
  // Constant per-species diffusivities
  {
    InputParameters params = getFactory().getValidParams("ADGenericConstantMaterial");
    assignBlocks(params, _blocks);
    std::vector<std::string> names;
    for (const auto & species : _solver_species)
      names.push_back(diffusivityName(species));
    params.set<std::vector<std::string>>("prop_names") = names;
    params.set<std::vector<Real>>("prop_values") = _diffusivities;
    getProblem().addMaterial("ADGenericConstantMaterial", prefix() + "diffusivities", params);
  }

  // One parsed material per species holding its net reaction rate S_s
  for (const auto i : index_range(_solver_species))
  {
    if (_rate_expression[i].empty())
      continue;
    InputParameters params = getFactory().getValidParams("ADParsedMaterial");
    assignBlocks(params, _blocks);
    params.set<std::string>("property_name") = rateName(_solver_species[i]);
    params.set<std::vector<VariableName>>("coupled_variables") = _rate_coupled_species[i];
    params.set<std::string>("expression") = _rate_expression[i];
    getProblem().addMaterial("ADParsedMaterial", prefix() + rateName(_solver_species[i]), params);
  }
}

void
KineticReactionsPhysics::addFEKernels()
{
  for (const auto i : index_range(_solver_species))
  {
    const auto & species = _solver_species[i];

    // Time derivative
    {
      InputParameters params = _factory.getValidParams("ADTimeDerivative");
      assignBlocks(params, _blocks);
      params.set<NonlinearVariableName>("variable") = species;
      _problem->addKernel("ADTimeDerivative", prefix() + species + "_time", params);
    }

    // Fickian diffusion (skip immobile species). ADMatDiffusion obeys coord_type.
    if (_diffusivities[i] != 0.0)
    {
      InputParameters params = _factory.getValidParams("ADMatDiffusion");
      assignBlocks(params, _blocks);
      params.set<NonlinearVariableName>("variable") = species;
      params.set<MaterialPropertyName>("diffusivity") = diffusivityName(species);
      _problem->addKernel("ADMatDiffusion", prefix() + species + "_diff", params);
    }

    // Net reaction source/sink: residual -S_s, applied via the parsed material
    if (!_rate_expression[i].empty())
    {
      InputParameters params = _factory.getValidParams("ADMatBodyForce");
      assignBlocks(params, _blocks);
      params.set<NonlinearVariableName>("variable") = species;
      params.set<MaterialPropertyName>("material_property") = rateName(species);
      _problem->addKernel("ADMatBodyForce", prefix() + species + "_rxn", params);
    }
  }
}
