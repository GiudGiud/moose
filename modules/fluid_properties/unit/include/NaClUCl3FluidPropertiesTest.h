//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseObjectUnitTest.h"
#include "NaClUCl3FluidProperties.h"

class NaClUCl3FluidPropertiesTest : public MooseObjectUnitTest
{
public:
  NaClUCl3FluidPropertiesTest() : MooseObjectUnitTest("FluidPropertiesApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("NaClUCl3FluidProperties");
    _fe_problem->addUserObject("NaClUCl3FluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<NaClUCl3FluidProperties>("fp");
  }

  const NaClUCl3FluidProperties * _fp;
};
