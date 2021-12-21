//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PINSFVPorosityJump.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", PINSFVPorosityJump);

InputParameters
PINSFVPorosityJump::validParams()
{
  InputParameters params = FVScalarLagrangeMultiplierInterface::validParams();
  params.addClassDescription("Creates a pressure drop across from a porosity interface based on "
      "the Bernouilli equation.");
  params.addRequiredCoupledVar(NS::pressure, "Pressure variable");
  params.addRequiredParam<MooseFunctorName>(NS::porosity, "Porosity functor");
  params.addRequiredParam<MooseFunctorName>(NS::density, "Density functor");
  params.addRequiredParam<MooseFunctorName>(NS::superficial_velocity_x,
                                            "Superficial x-velocity functor");
  params.addParam<MooseFunctorName>(NS::superficial_velocity_y, "Superficial y-velocity functor");
  params.addParam<MooseFunctorName>(NS::superficial_velocity_z, "Superficial z-velocity functor");

  return params;
}

PINSFVPorosityJump::PINSFVPorosityJump(const InputParameters & params)
  : FVScalarLagrangeMultiplierInterface(params),
  _dim(_subproblem.mesh().dimension()),
  _p_elem(adCoupledValue(NS::pressure)),
  _p_neighbor(adCoupledNeighborValue(NS::pressure)),
  _eps(getFunctor<ADReal>(NS::porosity)),
  _rho(getFunctor<ADReal>(NS::density)),
  _u(getFunctor<ADReal>(NS::superficial_velocity_x)),
  _v(isParamValid(NS::superficial_velocity_y) ? &getFunctor<ADReal>(NS::superficial_velocity_y)
                                              : nullptr),
  _w(isParamValid(NS::superficial_velocity_z) ? &getFunctor<ADReal>(NS::superficial_velocity_z)
                                              : nullptr)
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

ADReal
PINSFVPorosityJump::computeQpResidual()
{
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

  ADReal pressure_drop = 0.5 * _rho(elem_arg) * k * v_square;

  return _p_elem[_qp] - _p_neighbor[_qp] + pressure_drop;
}
