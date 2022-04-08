//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVMassWallFlowDiode.h"

registerMooseObject("NavierStokesApp", INSFVMassWallFlowDiode);

InputParameters
INSFVMassWallFlowDiode::validParams()
{
  InputParameters params = INSFVAdvectionInterfaceKernel::validParams();
  params.addClassDescription("If the flow is in the direction of the sideset, matches an mass "
                             "advection kernel. Else, becomes a wall boundary condition");
  params.addRequiredParam<MooseFunctorName>(NS::density, "Density functor");
  params.addParam<bool>(
      "invert_diode",
      false,
      "Whether to use the direction opposite the sideset normal for allowing flow");

  return params;
}

INSFVMassWallFlowDiode::INSFVMassWallFlowDiode(const InputParameters & params)
  : INSFVAdvectionInterfaceKernel(params),
    _rho(getFunctor<ADReal>(NS::density)),
    _inverted(1 - 2 * getParam<bool>("invert_diode"))
{
}

ADReal
INSFVMassWallFlowDiode::computeQpResidual()
{
  const auto elem_face = elemFromFace();
  const auto neighbor_face = neighborFromFace();

  auto v = _rc_vel_provider.getVelocity(_velocity_interp_method, *_face_info, _tid);

  if (_inverted * v * _normal < 0)
    v = 0;

  ADReal rho_interface;
  Moose::FV::interpolate(_advected_interp_method,
                         rho_interface,
                         _rho(elem_face),
                         _rho(neighbor_face),
                         v,
                         *_face_info,
                         true);
  return _normal * v * rho_interface;
}
