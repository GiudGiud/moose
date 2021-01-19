//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVPenaltyFreeSlipBC.h"

registerMooseObject("NavierStokesApp", INSFVPenaltyFreeSlipBC);

InputParameters
INSFVPenaltyFreeSlipBC::validParams()
{
  InputParameters params = INSFVSlipWallBC::validParams();
  params.addClassDescription(
      "Implements a free slip boundary condition using a penalty formulation.");
  params.addRequiredCoupledVar("u", "The velocity in the x direction.");
  params.addCoupledVar("v", 0, "The velocity in the y direction.");
  params.addCoupledVar("w", 0, "The velocity in the z direction.");
  MooseEnum momentum_component("x=0 y=1 z=2", "x");
  params.addParam<MooseEnum>("momentum_component",
                             momentum_component,
                             "The component of the momentum equation that this BC applies to.");
  params.addParam<Real>("penalty", 1e6, "The penalty factor");
  return params;
}

INSFVPenaltyFreeSlipBC::INSFVPenaltyFreeSlipBC(const InputParameters & params)
  : INSFVSlipWallBC(params),
    _u(adCoupledValue("u")),
    _v(adCoupledValue("v")),
    _w(adCoupledValue("w")),
    _comp(getParam<MooseEnum>("momentum_component")),
    _p(getParam<Real>("penalty"))

{
}

ADReal
INSFVPenaltyFreeSlipBC::computeQpResidual()
{
  const FaceInfo * const fi = _face_info;

  /// Obtain the variable names from the parameters
  const auto & velx_name = parameters().getParamHelper("u", parameters(),
      static_cast<std::vector<VariableName, std::allocator<VariableName> > *>(0));
  const auto & vely_name = parameters().getParamHelper("v", parameters(),
      static_cast<std::vector<VariableName, std::allocator<VariableName> > *>(0));
  const auto & velz_name = parameters().getParamHelper("w", parameters(),
      static_cast<std::vector<VariableName, std::allocator<VariableName> > *>(0));

  /// Get face value for velocity, using the variable's interpolation method
  const auto & vx_face = velx_name.empty() ? _u[_qp] :
      dynamic_cast<const MooseVariableFV<Real> *>(
      &_subproblem.getVariable(_tid, velx_name[0]))->getBoundaryFaceValue(*fi);

  const auto & vy_face = vely_name.empty() ? _v[_qp] :
      dynamic_cast<const MooseVariableFV<Real> *>(
      &_subproblem.getVariable(_tid, vely_name[0]))->getBoundaryFaceValue(*fi);

  const auto & vz_face = velz_name.empty() ? _w[_qp] :
      MetaPhysicL::raw_value(dynamic_cast<const MooseVariableFV<Real> *>(
      &_subproblem.getVariable(_tid, velz_name[0]))->getBoundaryFaceValue(*fi));

  return _p * _normal(_comp) * (_normal * ADRealVectorValue(vx_face, vy_face, vz_face));
}
