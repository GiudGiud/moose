//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PINSFVPorosityJumpv2.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", PINSFVPorosityJumpv2);

InputParameters
PINSFVPorosityJumpv2::validParams()
{
  InputParameters params = FVInterfaceKernel::validParams();
  params.addClassDescription("Creates a pressure drop across from a porosity interface based on "
      "the Bernouilli equation.");
  params.addRequiredCoupledVar(NS::pressure, "Pressure variable");
  params.addRequiredParam<MooseFunctorName>(NS::porosity, "Porosity functor");
  params.addRequiredParam<MooseFunctorName>(NS::density, "Density functor");
  params.addRequiredParam<MooseFunctorName>(NS::superficial_velocity_x,
                                            "Superficial x-velocity functor");
  params.addParam<MooseFunctorName>(NS::superficial_velocity_y, "Superficial y-velocity functor");
  params.addParam<MooseFunctorName>(NS::superficial_velocity_z, "Superficial z-velocity functor");
  MooseEnum momentum_component("x=0 y=1 z=2");
  params.addRequiredParam<MooseEnum>("momentum_component",
                                     momentum_component,
                                     "The component of the momentum equation that this kernel "
                                     "applies to.");
  return params;
}

PINSFVPorosityJumpv2::PINSFVPorosityJumpv2(const InputParameters & params)
  : FVInterfaceKernel(params),
  _dim(_subproblem.mesh().dimension()),
  _p_elem(adCoupledValue(NS::pressure)),
  _p_neighbor(adCoupledNeighborValue(NS::pressure)),
  _eps(getFunctor<ADReal>(NS::porosity)),
  _rho(getFunctor<ADReal>(NS::density)),
  _u(getFunctor<ADReal>(NS::superficial_velocity_x)),
  _v(isParamValid(NS::superficial_velocity_y) ? &getFunctor<ADReal>(NS::superficial_velocity_y)
                                              : nullptr),
  _w(isParamValid(NS::superficial_velocity_z) ? &getFunctor<ADReal>(NS::superficial_velocity_z)
                                              : nullptr),
  _index(getParam<MooseEnum>("momentum_component"))
{
#ifndef MOOSE_GLOBAL_AD_INDEXING
  mooseError("INSFV is not supported by local AD indexing. In order to use INSFV, please run the "
             "configure script in the root MOOSE directory with the configure option "
             "'--with-ad-indexing-type=global'");
#endif

  if (_dim >= 2 && !isParamValid(NS::superficial_velocity_y))
    paramError("v",
               "In two or more dimensions, the v superficial velocity must be supplied.");

  if (_dim >= 3 && !isParamValid(NS::superficial_velocity_z))
    paramError("w",
               "In three-dimensions, the w superficial velocity must be supplied.");
}

void
PINSFVPorosityJumpv2::computeResidual(const FaceInfo & fi)
{
  setupData(fi);

  const auto var_elem_num = _elem_is_one ? _var1.number() : _var2.number();
  const auto var_neigh_num = _elem_is_one ? _var2.number() : _var1.number();

  // COMPUTE RESIDUALS for each side
  const auto elem_arg = elemFromFace();
  const auto neighbor_arg = neighborFromFace();

  Real eps_1 = _eps(elem_arg).value();
  Real eps_2 = _eps(neighbor_arg).value();
  Real k = eps_1/eps_2 * (eps_2 - eps_1);

  // Superficial velocity is continuous, no need to chose a side
  // though can we avoid oscillations by picking a side?
  ADReal v_square = _u(elem_arg) * _u(elem_arg);
  if (_dim >= 2)
    v_square += (*_v)(elem_arg) * (*_v)(elem_arg);
  if (_dim == 3)
    v_square += (*_w)(elem_arg) * (*_w)(elem_arg);

  ADReal pressure_drop = -0.5 * _rho(elem_arg) * k * v_square;
  ADReal p_average;
  Moose::FV::interpolate(Moose::FV::InterpMethod::Average,
                         p_average,
                         _p_elem[_qp],
                         _p_neighbor[_qp],
                         *_face_info,
                         true);

  const auto residual1 = eps_1 * _normal(_index) * (p_average + pressure_drop/2);
  const auto residual2 = eps_2 * _normal(_index) * (p_average - pressure_drop/2);
  _console << _p_elem[_qp] << " " << _p_neighbor[_qp] << " " << pressure_drop << std::endl;
  _console << residual1 << " " << residual2 << std::endl;

  const auto r1 = MetaPhysicL::raw_value(fi.faceArea() * fi.faceCoord() * residual1);
  const auto r2 = MetaPhysicL::raw_value(fi.faceArea() * fi.faceCoord() * residual2);

  processResidual(r1, var_elem_num, false);
  processResidual(-r2, var_neigh_num, true);
}

void
PINSFVPorosityJumpv2::computeJacobian(const FaceInfo & fi)
{
  setupData(fi);

  const auto & elem_dof_indices = _elem_is_one ? _var1.dofIndices() : _var2.dofIndices();
  const auto & neigh_dof_indices =
      _elem_is_one ? _var2.dofIndicesNeighbor() : _var1.dofIndicesNeighbor();
  mooseAssert((elem_dof_indices.size() == 1) && (neigh_dof_indices.size() == 1),
              "We're currently built to use CONSTANT MONOMIALS");

  // COMPUTE RESIDUALS for each side
  const auto elem_arg = elemFromFace();
  const auto neighbor_arg = neighborFromFace();

  Real eps_1 = _eps(elem_arg).value();
  Real eps_2 = _eps(neighbor_arg).value();
  Real k = eps_1/eps_2 * (eps_2 - eps_1);

  // Superficial velocity is continuous, no need to chose a side
  // though can we avoid oscillations by picking a side?
  ADReal v_square = _u(elem_arg) * _u(elem_arg);
  if (_dim >= 2)
    v_square += (*_v)(elem_arg) * (*_v)(elem_arg);
  if (_dim == 3)
    v_square += (*_w)(elem_arg) * (*_w)(elem_arg);

  ADReal pressure_drop = -0.5 * _rho(elem_arg) * k * v_square;
  ADReal p_average;
  Moose::FV::interpolate(Moose::FV::InterpMethod::Average,
                         p_average,
                         _p_elem[_qp],
                         _p_neighbor[_qp],
                         *_face_info,
                         true);

  // eps_1 = 1;
  // eps_2 = 1;
  const auto residual1 = eps_1 * _normal(_index) * (p_average + pressure_drop/2);
  const auto residual2 = eps_2 * _normal(_index) * (p_average - pressure_drop/2);
  _console << _p_elem[_qp] << " " << _p_neighbor[_qp] << " " << pressure_drop << std::endl;
  _console << residual1 << " " << residual2 << std::endl;

  const auto r1 = fi.faceArea() * fi.faceCoord() * residual1;
  const auto r2 = fi.faceArea() * fi.faceCoord() * residual2;

  processDerivatives(r1, elem_dof_indices[0]);
  processDerivatives(-r2, neigh_dof_indices[0]);
}
