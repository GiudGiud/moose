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
#include "UserObject.h"
#include "BoundaryRestrictableRequired.h"
#include "MaterialPropertyInterface.h"
#include "Coupleable.h"
#include "MooseVariableDependencyInterface.h"
#include "TransientInterface.h"
#include "ElementIDInterface.h"

class SideUserObject : public UserObject,
                       public BoundaryRestrictableRequired,
                       public MaterialPropertyInterface,
                       public Coupleable,
                       public MooseVariableDependencyInterface,
                       public TransientInterface,
                       public ElementIDInterface
{
public:
  static InputParameters validParams();

  SideUserObject(const InputParameters & parameters);

protected:
  MooseMesh & _mesh;

  const MooseArray<Point> & _q_point;
  const QBase * const & _qrule;
  const MooseArray<Real> & _JxW;
  const MooseArray<Real> & _coord;
  const MooseArray<Point> & _normals;

  const Elem * const & _current_elem;
  /// current side of the current element
  const unsigned int & _current_side;

  const Elem * const & _current_side_elem;
  const Real & _current_side_volume;

  const BoundaryID & _current_boundary_id;

  /**
   * Determine the single sided face argument when evaluating a functor on a face.
   * This is used to perform evaluations of material properties with the actual face values of
   * their dependences, rather than interpolate the material property to the boundary.
   * @param fi the FaceInfo for this face
   * @param limiter_type the limiter type, to be specified if more than the default average
   *        interpolation is required for the parameters of the functor
   * @param correct_skewness whether to perform skew correction at the face
   */
  Moose::SingleSidedFaceArg singleSidedFaceArg(
      const MooseFunctorName & functor_name,
      const FaceInfo & fi,
      Moose::FV::LimiterType limiter_type = Moose::FV::LimiterType::CentralDifference,
      bool correct_skewness = false) const;
};
