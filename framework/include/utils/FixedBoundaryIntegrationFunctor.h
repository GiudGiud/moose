//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseMesh.h"
#include "MooseTypes.h"
#include "MooseError.h"
#include "MooseFunctor.h"
#include "Moose.h"
#include "Limiter.h"
#include "MathFVUtils.h"

#include "libmesh/elem.h"
#include "libmesh/remote_elem.h"

#include <unordered_map>
#include <functional>

/**
 * A functor that computes the integral of a functor on a fixed boundary.
 */
template <typename T>
class FixedBoundaryIntegralFunctor : public Moose::FunctorBase<T>
{
public:
  FixedBoundaryIntegralFunctor(const std::string & name,
                               const Moose::FunctorEnvelope<T> & functor,
                               const std::set<ExecFlagType> & clearance_schedule,
                               const MooseMesh & mesh,
                               const BoundaryID boundary_id);

  virtual ~FixedBoundaryIntegralFunctor() = default;

  bool isExtrapolatedBoundaryFace(const FaceInfo & fi,
                                  const Elem * elem,
                                  const Moose::StateArg & time) const override;

  bool hasBlocks(SubdomainID /* id */) const override { return true; }
  bool hasFaceSide(const FaceInfo & /*fi*/, const bool /*fi_elem_side*/) const override
  {
    return true;
  }

  using typename Moose::FunctorBase<T>::FunctorType;
  using typename Moose::FunctorBase<T>::ValueType;
  using typename Moose::FunctorBase<T>::FunctorReturnType;

protected:
  virtual ValueType globalValue(const Moose::StateArg & time) const;
  ValueType evaluate(const Moose::ElemArg & elem_arg, const Moose::StateArg & time) const override;
  ValueType evaluate(const Moose::FaceArg & face, const Moose::StateArg & time) const override;
  ValueType evaluate(const Moose::ElemQpArg & elem_qp, const Moose::StateArg & time) const override;
  ValueType evaluate(const Moose::ElemSideQpArg & elem_side_qp,
                     const Moose::StateArg & time) const override;
  ValueType evaluate(const Moose::ElemPointArg & elem_point,
                     const Moose::StateArg & time) const override;
  ValueType evaluate(const Moose::NodeArg & elem_point,
                     const Moose::StateArg & time) const override;

  /// Functor to evaluate to compute the boundary integral
  const Moose::Functor<T> & _functor;

  /// Boundary id of the boundary to integrate on
  const BoundaryID _bid;

  /// The mesh that this functor operates on
  const MooseMesh & _mesh;
};

template <typename T>
FixedBoundaryIntegralFunctor<T>::FixedBoundaryIntegralFunctor(
    const std::string & name,
    const Moose::FunctorEnvelope<T> & functor,
    const std::set<ExecFlagType> & clearance_schedule,
    const MooseMesh & mesh,
    const BoundaryID boundary_id)
  : Moose::FunctorBase<T>(name, clearance_schedule),
    _functor(functor),
    _bid(boundary_id),
    _mesh(mesh)
{
}

template <typename T>
bool
FixedBoundaryIntegralFunctor<T>::isExtrapolatedBoundaryFace(const FaceInfo & /*fi*/,
                                                            const Elem *,
                                                            const Moose::StateArg &) const
{
  return false;
}

template <typename T>
typename FixedBoundaryIntegralFunctor<T>::ValueType
FixedBoundaryIntegralFunctor<T>::globalValue(const Moose::StateArg & time) const
{
  using namespace Moose::FV;
  auto & binfo = _mesh.getMesh().get_boundary_info();

  // Compute the surface integral
  T sum = 0;
  for (const auto & [elem_id, side_id, bc_id] :
       binfo.build_side_list(libMesh::BoundaryInfo::BCTupleSortBy::BOUNDARY_ID))
  {
    if (bc_id != _bid)
      continue;
    const auto elem = _mesh.elemPtr(elem_id);
    auto fi = _mesh.faceInfo(elem, side_id);
    if (!fi)
    {
      const auto neighbor = elem->neighbor_ptr(side_id);
      mooseAssert(neighbor, "We should have a neighbor");
      fi = _mesh.faceInfo(neighbor, neighbor->which_neighbor_am_i(elem));
    }
    mooseAssert(fi, "We should have a face info");
    Moose::FaceArg face_arg = {fi, Moose::FV::LimiterType::CentralDifference, true, true, nullptr};
    sum += _functor(face_arg, time) * fi->faceArea() * fi->faceCoord();
  }
  return sum;
}

template <typename T>
typename FixedBoundaryIntegralFunctor<T>::ValueType
FixedBoundaryIntegralFunctor<T>::evaluate(const Moose::ElemArg & /*elem_arg*/,
                                          const Moose::StateArg & time) const
{
  return globalValue(time);
}

template <typename T>
typename FixedBoundaryIntegralFunctor<T>::ValueType
FixedBoundaryIntegralFunctor<T>::evaluate(const Moose::FaceArg & /*face*/,
                                          const Moose::StateArg & time) const
{
  return globalValue(time);
}

template <typename T>
typename FixedBoundaryIntegralFunctor<T>::ValueType
FixedBoundaryIntegralFunctor<T>::evaluate(const Moose::ElemQpArg & /*elem_qp*/
                                          ,
                                          const Moose::StateArg & time) const
{
  return globalValue(time);
}

template <typename T>
typename FixedBoundaryIntegralFunctor<T>::ValueType
FixedBoundaryIntegralFunctor<T>::evaluate(const Moose::ElemSideQpArg & /*elem_side_qp*/,
                                          const Moose::StateArg & time) const
{
  return globalValue(time);
}

template <typename T>
typename FixedBoundaryIntegralFunctor<T>::ValueType
FixedBoundaryIntegralFunctor<T>::evaluate(const Moose::ElemPointArg & /*elem_point_arg*/,
                                          const Moose::StateArg & time) const
{
  return globalValue(time);
}

template <typename T>
typename FixedBoundaryIntegralFunctor<T>::ValueType
FixedBoundaryIntegralFunctor<T>::evaluate(const Moose::NodeArg & /*node_arg*/,
                                          const Moose::StateArg & time) const
{
  return globalValue(time);
}

/**
 * A functor that computes the average of another functor on a boundary
 */
template <typename T>
class BoundaryAverageFunctor : public FixedBoundaryIntegralFunctor<T>
{
public:
  BoundaryAverageFunctor(const std::string & name,
                         const Moose::FunctorEnvelope<T> & functor,
                         const std::set<ExecFlagType> & clearance_schedule,
                         const MooseMesh & mesh,
                         const BoundaryID boundary_id);

  virtual ~BoundaryAverageFunctor() = default;

  using typename Moose::FunctorBase<T>::FunctorType;
  using typename Moose::FunctorBase<T>::ValueType;
  using typename Moose::FunctorBase<T>::FunctorReturnType;
  using FixedBoundaryIntegralFunctor<T>::_mesh;
  using FixedBoundaryIntegralFunctor<T>::_bid;
  using FixedBoundaryIntegralFunctor<T>::_functor;

protected:
  ValueType globalValue(const Moose::StateArg & time) const override;
};

template <typename T>
BoundaryAverageFunctor<T>::BoundaryAverageFunctor(const std::string & name,
                                                  const Moose::FunctorEnvelope<T> & functor,
                                                  const std::set<ExecFlagType> & clearance_schedule,
                                                  const MooseMesh & mesh,
                                                  const BoundaryID boundary_id)
  : FixedBoundaryIntegralFunctor<T>(name, functor, clearance_schedule, mesh, boundary_id)
{
}

template <typename T>
typename BoundaryAverageFunctor<T>::ValueType
BoundaryAverageFunctor<T>::globalValue(const Moose::StateArg & time) const
{
  using namespace Moose::FV;
  auto & binfo = _mesh.getMesh().get_boundary_info();

  // Compute the surface integral
  T sum = 0;
  T area = 0;
  for (const auto & [elem_id, side_id, bc_id] :
       binfo.build_side_list(libMesh::BoundaryInfo::BCTupleSortBy::BOUNDARY_ID))
  {
    if (bc_id != _bid)
      continue;
    const auto elem = _mesh.elemPtr(elem_id);
    auto fi = _mesh.faceInfo(elem, side_id);
    if (!fi)
    {
      const auto neighbor = elem->neighbor_ptr(side_id);
      mooseAssert(neighbor, "We should have a neighbor");
      fi = _mesh.faceInfo(neighbor, neighbor->which_neighbor_am_i(elem));
    }
    mooseAssert(fi, "We should have a face info");
    mooseAssert(fi, "We should have a face info");
    Moose::FaceArg face_arg = {fi, Moose::FV::LimiterType::CentralDifference, true, true, nullptr};
    sum += _functor(face_arg, time) * fi->faceArea() * fi->faceCoord();
    area += fi->faceArea() * fi->faceCoord();
  }
  return sum / area;
}
