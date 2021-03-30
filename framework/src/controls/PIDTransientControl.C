//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "PIDTransientControl.h"
#include "Function.h"
#include "Transient.h"

registerMooseObject("MooseApp", PIDTransientControl);

defineLegacyParams(PIDTransientControl);

InputParameters
PIDTransientControl::validParams()
{
  InputParameters params = Control::validParams();
  params.addClassDescription(
      "Sets the value of a 'Real' input parameters based on a Proportional Integral "
      "Derivative control to make a postprocessor match a target value.");
  params.addRequiredParam<PostprocessorName>(
      "postprocessor", "The postprocessor to use for controlling the specified parameter.");
  params.addRequiredParam<FunctionName>("target", "The target value 1D time function for the postprocessor");
  params.addParam<Real>("K_integral", 1, "The coefficient multiplying the integral term");
  params.addParam<Real>("K_proportional", 1, "The coefficient multiplying the difference term");
  params.addParam<Real>("K_derivative", 0.1, "The coefficient multiplying the derivative term");
  params.addRequiredParam<std::string>(
      "parameter",
      "The input parameter(s) to control. Specify a single parameter name and all "
      "parameters in all objects matching the name will be updated");
  params.addParam<Real>("start_time", -std::numeric_limits<Real>::max(), "The time to start the PID controller at");
  params.addParam<Real>("stop_time", std::numeric_limits<Real>::max(), "The time to stop the PID controller at");
  params.addParam<bool>("reset_every_timestep",
      false,
      "Reset the PID integral when changing timestep, to keep the PID within Picard iterations");

  return params;
}

PIDTransientControl::PIDTransientControl(const InputParameters & parameters)
  : Control(parameters),
  _current(getPostprocessorValueByName(getParam<PostprocessorName>("postprocessor"))),
  _target(getFunction("target")),
  _Kint(getParam<Real>("K_integral")),
  _Kpro(getParam<Real>("K_proportional")),
  _Kder(getParam<Real>("K_derivative")),
  _start_time(getParam<Real>("start_time")),
  _stop_time(getParam<Real>("stop_time")),
  _reset_every_timestep(getParam<bool>("reset_every_timestep"))
{
  _integral = 0;

  // For Picard iteration resets
  _value_old = 0;
  _integral_old = 0;
}

void
PIDTransientControl::execute()
{
  Point dummy;

  if (_t >= _start_time && _t < _stop_time)
  {
    // Compute the new parameter based on the PID control of the postprocessor
    Real value = getControllableValue<Real>("parameter");

    // Save integral and controlled value for the first Picard iteration
    if (dynamic_cast<Transient *>(_app.getExecutioner())->picardSolve().numPicardIts() == 0)
    {
      _integral_old = _integral;
      _value_old = value;
    }

    // If there were Picard iterations during the transient and they failed,
    // need to reset the controlled value and the error integral
    if (dynamic_cast<Transient *>(_app.getExecutioner())->picardSolve().numPicardIts() == 0)
    {
      _integral = _integral_old;
      value = _value_old;
    }

    // Reset the error integral if PID is used only within each timestep
    if (_reset_every_timestep)
      _integral = 0;

    // Compute the three error terms and add them to the controlled value
    _integral += (_current - _target.value(_t, dummy)) * _dt;
    value += _Kint * _integral + _Kpro * (_current- _target.value(_t, dummy));
    if (_dt > 0)
      value += _Kder * (_current- _target.value(_t, dummy)) / _dt;

    std::cout << "PID "<<getControllableValue<Real>("parameter")<< " (" << _Kint * _integral<<" " << _Kpro * (_current- _target.value(_t, dummy)) << " " <<
       _Kder * (_current- _target.value(_t, dummy)) / _dt << ") -> " << value << std::endl;

    // Set the new value using the Controllable system
    setControllableValue<Real>("parameter", value);
  }
}
