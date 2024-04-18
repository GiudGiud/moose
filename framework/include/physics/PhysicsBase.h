//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Action.h"
#include "InputParametersChecksUtils.h"

// We include these headers for all the derived classes that will be building objects
#include "FEProblemBase.h"
#include "Factory.h"
#include "MultiMooseEnum.h"

#define registerPhysicsBaseTasks(app_name, derived_name)                                           \
  registerMooseAction(app_name, derived_name, "init_physics");                                     \
  registerMooseAction(app_name, derived_name, "copy_vars_physics")

/**
 * Base class to help creating an entire physics
 */
class PhysicsBase : public Action, public InputParametersChecksUtils<PhysicsBase>
{
public:
  static InputParameters validParams();

  PhysicsBase(const InputParameters & parameters);

  /// Forwards from the action tasks to the implemented addXYZ() in the derived classes
  /// If you need more than these:
  /// - register your action to the new task using
  ///   registerMooseAction("AppName", ActionClass, "task_name");
  /// - override actOnAdditionalTasks and add your additional work there
  virtual void act() override final;

  /// Routine to add additional setup work on additional registered tasks to a Physics
  virtual void actOnAdditionalTasks() {}

  /// Add new blocks to the Physics
  void addBlocks(const std::vector<SubdomainName> & blocks);

  /// Return the blocks this physics is defined on
  const std::vector<SubdomainName> & blocks() const { return _blocks; }

  /// Check if an external object has the same block restriction
  bool checkBlockRestrictionIdentical(const std::string & object_name,
                                      const std::vector<SubdomainName> & blocks,
                                      const bool error_if_not_identical = true) const;

  /// Provide additional parameters for the relationship managers
  virtual InputParameters getAdditionalRMParams() const { return emptyInputParameters(); };

  /// Get a Physics from the ActionWarehouse with the requested type and name
  template <typename T>
  const T * getCoupledPhysics(const PhysicsName & phys_name, const bool allow_fail = false) const;
  /// Get all Physics from the ActionWarehouse with the requested type
  template <typename T>
  const std::vector<T *> getCoupledPhysics(const bool allow_fail = false) const;

  /// Utilities to merge two Physics of the same type together
  /// Check that parameters are compatible for a merge with another Physics
  virtual bool checkParametersMergeable(const InputParameters & /*param*/, bool /*warn*/) const
  {
    mooseError("Not implemented");
  }
  /// Merge these parameters into existing parameters of this Physics
  virtual void mergeParameters(const InputParameters & /*params*/)
  {
    mooseError("Not implemented");
  }

protected:
  /// Return whether the Physics is solved using a transient
  bool isTransient() const;
  /// Return the maximum dimension of the blocks the Physics is active on
  unsigned int dimension() const;

  /// Get the factory for this physics
  /// The factory lets you get the parameters for objects
  virtual Factory & getFactory() { return _factory; }
  virtual Factory & getFactory() const { return _factory; }
  /// Get the problem for this physics
  /// Useful to add objects to the simulation
  virtual FEProblemBase & getProblem()
  {
    mooseAssert(_problem, "Requesting the problem too early");
    return *_problem;
  }
  virtual const FEProblemBase & getProblem() const
  {
    mooseAssert(_problem, "Requesting the problem too early");
    return *_problem;
  }

  /// Tell the app if we want to use Exodus restart
  void prepareCopyVariablesFromMesh() const;
  /// Copy variables from the mesh file
  void copyVariablesFromMesh(const std::vector<VariableName> & variables_to_copy);

  /// Use prefix() to disambiguate names
  std::string prefix() const { return name() + "_"; }

  /// Return the list of nonlinear variables in this physics
  const std::vector<VariableName> & nonlinearVariableNames() const { return _nl_var_names; };
  /// Keep track of the name of a nonlinear variable defined in the Physics
  void saveNonlinearVariableName(const VariableName & var_name)
  {
    _nl_var_names.push_back(var_name);
  }

  /// Check whether a nonlinear variable already exists
  bool nonlinearVariableExists(const VariableName & var_name, bool error_if_aux) const;

  /// Add a new required task for all physics deriving from this class
  /// NOTE: This does not register the task, you still need to call registerMooseAction
  void addRequiredPhysicsTask(const std::string & task) { _required_tasks.insert(task); }

  /// System number for the system owning the variables
  const unsigned int _sys_number;

  /// Whether to output additional information
  const bool _verbose;

  /// Whether to add a default preconditioning.
  /// The implementation of the default is defined by the derived class
  const MooseEnum & _preconditioning;

  /// Keep track of the subdomains the Physics is defined on
  std::vector<SubdomainName> _blocks;

  /// Utilities to process and forward parameters
  void assignBlocks(InputParameters & params, const std::vector<SubdomainName> & blocks) const;
  /// Check if a vector contains all the mesh blocks
  bool allMeshBlocks(const std::vector<SubdomainName> & blocks) const;

  /// Routine to help create maps
  template <typename T, typename C>
  std::map<T, C> createMapFromVectors(std::vector<T> keys, std::vector<C> values) const;
  template <typename T>
  std::map<T, MooseEnum> createMapFromVectorAndMultiMooseEnum(std::vector<T> keys,
                                                              MultiMooseEnum values) const;

private:
  /// Gathers additional parameters for the relationship managers from the Physics
  /// then calls the parent Action::addRelationshipManagers with those parameters
  using Action::addRelationshipManagers;
  virtual void addRelationshipManagers(Moose::RelationshipManagerType input_rm_type) override;

  /// Process some parameters that require the problem to be created. Executed on init_physics
  void initializePhysics();
  /// Additional initialization work that should happen very early, as soon as the problem is created
  virtual void initializePhysicsAdditional() {}

  /// The default implementation of these routines will do nothing as we do not expect all Physics
  /// to be defining an object of every type
  virtual void addNonlinearVariables() {}
  virtual void addAuxiliaryVariables() {}
  virtual void addInitialConditions() {}
  virtual void addFEKernels() {}
  virtual void addFVKernels() {}
  virtual void addNodalKernels() {}
  virtual void addDiracKernels() {}
  virtual void addDGKernels() {}
  virtual void addScalarKernels() {}
  virtual void addInterfaceKernels() {}
  virtual void addFVInterfaceKernels() {}
  virtual void addFEBCs() {}
  virtual void addFVBCs() {}
  virtual void addNodalBCs() {}
  virtual void addPeriodicBCs() {}
  virtual void addFunctions() {}
  virtual void addAuxiliaryKernels() {}
  virtual void addMaterials() {}
  virtual void addFunctorMaterials() {}
  virtual void addUserObjects() {}
  virtual void addPostprocessors() {}
  virtual void addVectorPostprocessors() {}
  virtual void addReporters() {}
  virtual void addOutputs() {}
  virtual void addPreconditioning() {}
  virtual void addExecutioner() {}
  virtual void addExecutors() {}

  /// Check the list of required tasks for missing tasks
  void checkRequiredTasks() const;

  /// Whether the physics is to be solved as a transient. It can be advantageous to solve
  /// some physics directly to steady state
  MooseEnum _is_transient;

  /// Vector of the nonlinear variables in the Physics
  std::vector<VariableName> _nl_var_names;

  /// Dimension of the physics, which we expect for now to be the dimension of the mesh
  /// NOTE: this is not known at construction time, only after initializePhysics which is a huge bummer
  unsigned int _dim = libMesh::invalid_uint;

  /// Manually keeps track of the tasks required by each physics as tasks cannot be inherited
  std::set<std::string> _required_tasks;
};

template <typename T>
const T *
PhysicsBase::getCoupledPhysics(const PhysicsName & phys_name, const bool allow_fail) const
{
  const auto all_T_physics = _awh.getActions<T>();
  for (const auto * const physics : all_T_physics)
  {
    if (physics->name() == phys_name)
      return physics;
  }
  if (!allow_fail)
    mooseError("Requested Physics '",
               phys_name,
               "' does not exist or is not of type '",
               MooseUtils::prettyCppType<T>(),
               "'");
  else
    return nullptr;
}

template <typename T>
const std::vector<T *>
PhysicsBase::getCoupledPhysics(const bool allow_fail) const
{
  const auto all_T_physics = _awh.getActions<T>();
  if (!allow_fail && all_T_physics.empty())
    mooseError("No Physics of requested type '", MooseUtils::prettyCppType<T>(), "'");
  else
    return all_T_physics;
}

template <typename T, typename C>
std::map<T, C>
PhysicsBase::createMapFromVectors(std::vector<T> keys, std::vector<C> values) const
{
  std::map<T, C> map;
  // No values have been specified.
  if (!values.size())
  {
    return map;
    // If we cant return a map of default C, dont try it
    // if constexpr (std::is_same_v<MooseEnum, T> || std::is_same_v<MultiMooseEnum, T>)
    //   return map;

    // C def;
    // for (const auto & k : keys)
    //   map[k] = def;
    // return map;
  }
  std::transform(keys.begin(),
                 keys.end(),
                 values.begin(),
                 std::inserter(map, map.end()),
                 [](T a, C b) { return std::make_pair(a, b); });
  return map;
}

template <typename T>
std::map<T, MooseEnum>
PhysicsBase::createMapFromVectorAndMultiMooseEnum(std::vector<T> keys, MultiMooseEnum values) const
{
  std::map<T, MooseEnum> map;
  // No values have been specified. We cant form a map of empty MooseEnum
  if (!values.size())
    return map;
  std::transform(keys.begin(),
                 keys.end(),
                 values.begin(),
                 std::inserter(map, map.end()),
                 [values](T a, MooseEnumItem b)
                 {
                   // Create a MooseEnum from the available values in the MultiMooseEnum and an
                   // actual current active item from that same MultiMooseEnum
                   MooseEnum single_value(values.getRawNames(), b.name());
                   return std::make_pair(a, single_value);
                 });
  return map;
}
