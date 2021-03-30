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
  params.addParam<std::string>(
      "parameter",
      "The input parameter(s) to control. Specify a single parameter name and all "
      "parameters in all objects matching the name will be updated");
  params.addParam<std::string>("parameter_pp", "The postprocessor parameter(s) to control.");
  params.addParam<Real>("start_time", -std::numeric_limits<Real>::max(), "The time to start the PID controller at");
  params.addParam<Real>("stop_time", std::numeric_limits<Real>::max(), "The time to stop the PID controller at");
  params.addParam<bool>("reset_every_timestep",
      false,
      "Reset the PID integral when changing timestep, to keep the PID within Picard iterations");
  params.addParam<bool>("reset_integral_windup",
      true,
      "Reset the PID integral when the error crosses zero and the integral is larger than the error.");

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
  _reset_every_timestep(getParam<bool>("reset_every_timestep")),
  _reset_integral_windup(getParam<bool>("reset_integral_windup"))
{
  _integral = 0;
  _previous_value = 0;

  // For Picard iteration resets
  _value_old = 0;
  _integral_old = 0;

  if (isParamValid("parameter") && isParamValid("parameter_pp"))
    paramError("parameter_pp",
               "Either a controllable parameter or a postprocessor to control should be specified, not both.");
  if (!isParamValid("parameter") && !isParamValid("parameter_pp"))
    mooseError("A parameter or a postprocessor to control should be specified.");
}

void
PIDTransientControl::execute()
{
  Point dummy;

  if (_t >= _start_time && _t < _stop_time)
  {
    std::cout << _t << " " << _t_step << " " << dynamic_cast<Transient *>(_app.getExecutioner())->picardSolve().numPicardIts() << std::endl;
    // Get current value of the controllable parameter
    Real value;
    if (isParamValid("parameter"))
      value = getControllableValue<Real>("parameter");
    else
      value = getPostprocessorValueByName(getParam<std::string>("parameter_pp"));

    // Save integral and controlled value for the first Picard iteration
    // if the Picard iteration fails, a smaller time step will be used but _t_step is unchanged
    if (_t_step != _t_step_old)
    {
      // Reset the error integral if PID is used only within each timestep
      if (_reset_every_timestep)
        _integral = 0;

      _integral_old = _integral;
      _value_old = value;
      _t_step_old = _t_step;
      _old_delta = 0;
    }

    // If there were Picard iterations during the transient and they failed,
    // need to reset the controlled value and the error integral
    if (dynamic_cast<Transient *>(_app.getExecutioner())->picardSolve().numPicardIts() == 1)
    {
      _integral = _integral_old;
      std::cout << "Resetting integral" << std::endl;
      value = _value_old;
    }

    // If there are Picard iterations, dont use the current received postprocessor value,
    // use the value from the previous iteration
    if (dynamic_cast<Transient *>(_app.getExecutioner())->picardSolve().numPicardIts() > 1)
    {
      value = _previous_value;
    }

    // Compute the delta between the current value of the postprocessor and the desired value
    Real delta = _current - _target.value(_t, dummy);

    // Reset integral of error if error crosses zero
    if (_reset_integral_windup && delta * _old_delta < 0) // && std::abs(_Kint *_integral) > std::abs(_Kpro * delta))
    { //std::cout << "Integral reset " <<  _Kint *_integral << " " << _Kpro * delta << std::endl;
      _integral = 0;
    }

    // Compute the three error terms and add them to the controlled value
    _integral += delta * _dt;
    value += _Kint * _integral + _Kpro * delta;
    if (_dt > 0)
      value += _Kder * delta / _dt;

    std::cout << "PID "<< value - _Kder * (_current- _target.value(_t, dummy)) / _dt-(_Kint * _integral + _Kpro * (_current- _target.value(_t, dummy)))
       << " (" << _Kint * _integral<<" " << _Kpro * delta << " " <<
       _Kder * delta / _dt << ") -> " << value << std::endl;

    if (isParamValid("parameter"))
      // Set the new value using the Controllable system
      setControllableValue<Real>("parameter", value);
    else
    {
      // Set the new postprocessor value
      _fe_problem.setPostprocessorValueByName(getParam<std::string>("parameter_pp"), value);

      if (dynamic_cast<Transient *>(_app.getExecutioner())->picardSolve().numPicardIts() > 1)
        // Set the previous postprocessor value, to have ChangeOverPicardPostprocessor correct
        _fe_problem.setPostprocessorValueByName(getParam<std::string>("parameter_pp"), value, 1);
    }

    // Save value for Picard iteration purposes
    _previous_value = value;

    // Keep track of previous delta to avoid integral windup
    _old_delta = delta;
  }
}
