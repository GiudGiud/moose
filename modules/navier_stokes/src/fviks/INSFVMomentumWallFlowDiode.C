//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVMomentumWallFlowDiode.h"
#include "MooseVariableFV.h"

registerMooseObject("NavierStokesApp", INSFVMomentumWallFlowDiode);

InputParameters
INSFVMomentumWallFlowDiode::validParams()
{
  InputParameters params = INSFVAdvectionInterfaceKernel::validParams();
  params += INSFVMomentumResidualObject::validParams();

  params.addClassDescription("If the flow is in the direction of the sideset, matches an advection "
                             "kernel. Else, becomes a free-slip wall boundary condition");
  params.addRequiredParam<MooseFunctorName>(NS::density, "Density functor");
  params.addParam<bool>(
      "invert_diode",
      false,
      "Whether to use the direction opposite the sideset normal for allowing flow");
  return params;
}

INSFVMomentumWallFlowDiode::INSFVMomentumWallFlowDiode(const InputParameters & params)
  : INSFVAdvectionInterfaceKernel(params),
    INSFVMomentumResidualObject(*this),
    _rho(getFunctor<ADReal>(NS::density)),
    _inverted(1 - 2 * getParam<bool>("invert_diode"))
{
}

ADReal
INSFVMomentumWallFlowDiode::computeQpResidual()
{
  const auto elem_face = elemFromFace();
  const auto neighbor_face = neighborFromFace();

  auto v = _rc_vel_provider.getVelocity(_velocity_interp_method, *_face_info, _tid);

  if (_inverted * v * _normal < 0)
    v = 0;

  const auto interp_coeffs = Moose::FV::interpCoeffs(_advected_interp_method, *_face_info, true, v);
  _ae = _normal * v * _rho(elem_face) * interp_coeffs.first;
  // Minus sign because we apply a minus sign to the residual in computeResidual
  _an = -_normal * v * _rho(neighbor_face) * interp_coeffs.second;

  return _ae * var1()(elem_face) - _an * var2()(neighbor_face);
}

void
INSFVMomentumWallFlowDiode::gatherRCData(const FaceInfo & fi)
{
  _face_info = &fi;
  _normal = fi.normal();

  const auto saved_velocity_interp_method = _velocity_interp_method;
  _velocity_interp_method = Moose::FV::InterpMethod::Average;
  // Fill-in the coefficients _ae and _an (but without multiplication by A)
  computeQpResidual();
  _velocity_interp_method = saved_velocity_interp_method;

  _rc_uo.addToA(&fi.elem(), _index, _ae * (fi.faceArea() * fi.faceCoord()));
  _rc_uo.addToA(fi.neighborPtr(), _index, _an * (fi.faceArea() * fi.faceCoord()));
}
