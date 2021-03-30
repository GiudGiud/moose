//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "Control.h"

// Forward declarations
class PIDTransientControl;

template <>
InputParameters validParams<PIDTransientControl>();

/**
 * An time-dependent control for changing an input parameter to make a target
 * postprocessor match a desired value.
 */
class PIDTransientControl : public Control
{
public:
  /**
   * Class constructor
   * @param parameters Input parameters for this Control object
   */
  static InputParameters validParams();

  PIDTransientControl(const InputParameters & parameters);

  virtual void execute() override;

private:
  /// The current value of the target postprocessor
  const PostprocessorValue & _current;
  /// The target 1D time-dependent function for the postprocessor
  const Function & _target;
  /// Integral of the error
  Real _integral;
  /// The coefficient multiplying the integral of the error
  const Real _Kint;
  /// The coefficient multiplying the error
  const Real _Kpro;
  /// The coefficient multiplying the derivative of the error
  const Real _Kder;
  /// The time to start the PID controller on
  const Real _start_time;
  /// The time to stop using the PID controller on
  const Real _stop_time;
  /// Whether to reset reset the PID error when changing timestep, for a 'Picard' PID
  const bool _reset_every_timestep;
  /// Data to save to be able to recover from failed Picard solves
  Real _integral_old;
  Real _value_old;
};
