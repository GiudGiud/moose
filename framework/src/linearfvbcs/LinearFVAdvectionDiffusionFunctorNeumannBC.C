//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVAdvectionDiffusionFunctorNeumannBC.h"

registerMooseObject("MooseApp", LinearFVAdvectionDiffusionFunctorNeumannBC);

InputParameters
LinearFVAdvectionDiffusionFunctorNeumannBC::validParams()
{
  InputParameters params = LinearFVAdvectionDiffusionBC::validParams();
  params.addClassDescription("Adds a Neumann BC which can be used for the assembly of linear "
                             "finite volume system and whose face values are determined using a "
                             "functor. This boundary condition is "
                             "only designed to work with advection-diffusion problems.");
  params.addRequiredParam<MooseFunctorName>("functor", "The functor for this boundary condition.");
  return params;
}

LinearFVAdvectionDiffusionFunctorNeumannBC::LinearFVAdvectionDiffusionFunctorNeumannBC(
    const InputParameters & parameters)
  : LinearFVAdvectionDiffusionBC(parameters), _functor(getFunctor<Real>("functor"))
{
}

Real
LinearFVAdvectionDiffusionFunctorNeumannBC::computeBoundaryValue() const
{
  mooseError("Not implemented for NeumannBC");
}

Real
LinearFVAdvectionDiffusionFunctorNeumannBC::computeBoundaryNormalGradient() const
{
  mooseError("Not implemented for NeumannBC");
}

Real
LinearFVAdvectionDiffusionFunctorNeumannBC::computeBoundaryValueMatrixContribution() const
{
  // Ths will not contribute to the matrix from the value considering that
  // the value is independent of the solution.
  return 0.0;
}

Real
LinearFVAdvectionDiffusionFunctorNeumannBC::computeBoundaryValueRHSContribution() const
{
  // Fetch the boundary value from the provided functor.
  return 0;
}

Real
LinearFVAdvectionDiffusionFunctorNeumannBC::computeBoundaryGradientMatrixContribution() const
{
  return 0.0;
}

Real
LinearFVAdvectionDiffusionFunctorNeumannBC::computeBoundaryGradientRHSContribution() const
{
  // The boundary term from the central difference approximation of the
  // normal gradient.
  return _functor(singleSidedFaceArg(_current_face_info), determineState());
}
