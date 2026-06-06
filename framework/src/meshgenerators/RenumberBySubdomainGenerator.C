//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RenumberBySubdomainGenerator.h"
#include "CastUniquePointer.h"
#include "MooseMeshUtils.h"

#include "libmesh/boundary_info.h"

registerMooseObject("MooseApp", RenumberBySubdomainGenerator);

InputParameters
RenumberBySubdomainGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addParam<std::vector<SubdomainName>>(
      "blocks_to_renumber",
      "Elements and nodes within these blocks will be renumbered. If none are specified, all "
      "subdomains are affected by the renumbering");
  params.addParam<std::vector<BoundaryName>>(
      "boundary_ordering",
      "If specified, within each subdomain the elements and nodes that lie on these boundaries are "
      "numbered first, in the order the boundaries are listed, with interior elements and nodes "
      "numbered afterwards");

  params.addClassDescription("Changes the element and node IDs so that elements and nodes are "
                             "contiguous within a subdomain. Note that DoF ordering may be "
                             "affected as well, and that the mesh renumbering will be turned off.");

  return params;
}

RenumberBySubdomainGenerator::RenumberBySubdomainGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters), _input(getMesh("input"))
{
}

std::unique_ptr<MeshBase>
RenumberBySubdomainGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  // Not impossible to do.
  if (!mesh->is_serial())
    mooseError("Not implemented for non-serialized distributed meshes");

  // Orphaned nodes would cause problems on renumbering, we are looping on the nodes attached
  // to elements
  mesh->remove_orphaned_nodes();
  mesh->renumber_nodes_and_elements();

  // Get the blocks provided by the user
  std::optional<std::vector<SubdomainName>> blocks =
      isParamValid("blocks_to_renumber")
          ? getParam<std::vector<SubdomainName>>("blocks_to_renumber")
          : std::vector<SubdomainName>();
  std::vector<SubdomainID> block_ids(blocks->size());
  std::stringstream missing_block;

  for (const auto i : index_range(block_ids))
  {
    const SubdomainName & name = blocks.value()[i];

    // Convert the SubdomainName to an id and store
    const auto id = MooseMeshUtils::getSubdomainID(name, *mesh);
    block_ids[i] = id;

    // Block does not exist - store for a future error
    if (id == Moose::INVALID_BLOCK_ID)
      missing_block << name << " ";
  }
  if (missing_block.str().size())
    paramError("blocks_to_renumber",
               "The following blocks were requested to be renumbered, but do not exist: ",
               missing_block.str());

  // User did not specify blocks, just get them all
  if (blocks->empty())
  {
    std::set<subdomain_id_type> block_ids_set;
    mesh->subdomain_ids(block_ids_set);
    block_ids.reserve(block_ids_set.size());
    block_ids.assign(block_ids_set.begin(), block_ids_set.end());

    // Does not have useful data
    blocks.reset();
  }

  // Resolve the boundaries used for the inner ordering, preserving the user-specified order
  std::vector<boundary_id_type> boundary_order;
  if (isParamValid("boundary_ordering"))
  {
    const auto & boundary_names = getParam<std::vector<BoundaryName>>("boundary_ordering");
    boundary_order.resize(boundary_names.size());
    std::stringstream missing_boundary;
    for (const auto i : index_range(boundary_names))
    {
      const auto id = MooseMeshUtils::getBoundaryID(boundary_names[i], *mesh);
      boundary_order[i] = id;
      if (id == Moose::INVALID_BOUNDARY_ID)
        missing_boundary << boundary_names[i] << " ";
    }
    if (missing_boundary.str().size())
      paramError("boundary_ordering",
                 "The following boundaries were requested for the inner ordering, but do not "
                 "exist: ",
                 missing_boundary.str());
  }

  // Build, for each requested boundary, the set of elements and nodes that lie on it. Nodes are
  // tracked by pointer rather than ID because the node IDs change as we renumber below
  std::map<boundary_id_type, std::unordered_set<dof_id_type>> boundary_elems;
  std::map<boundary_id_type, std::unordered_set<const Node *>> boundary_nodes;
  if (!boundary_order.empty())
  {
    const std::set<boundary_id_type> requested(boundary_order.begin(), boundary_order.end());
    const auto & boundary_info = mesh->get_boundary_info();

    // Elements (and the nodes of their boundary sides) from the sidesets
    for (const auto & [elem_id, side, bid] : boundary_info.build_active_side_list())
      if (requested.count(bid))
      {
        boundary_elems[bid].insert(elem_id);
        const auto * const elem = mesh->elem_ptr(elem_id);
        for (const auto local_node : elem->nodes_on_side(side))
          boundary_nodes[bid].insert(elem->node_ptr(local_node));
      }

    // Nodes from the nodesets
    for (const auto & [node_id, bid] : boundary_info.build_node_list())
      if (requested.count(bid))
        boundary_nodes[bid].insert(mesh->node_ptr(node_id));
  }

  // Renumber all elements with an ID we can recognize (so we can tell an already renumbered elem)
  const auto max_elem_id = mesh->max_elem_id();
  const auto max_node_id = mesh->max_node_id();
  std::unordered_map<dof_id_type, dof_id_type> new_elem_ids;
  new_elem_ids.reserve(max_elem_id);
  std::unordered_set<const Node *> renumbered_nodes;
  renumbered_nodes.reserve(max_node_id);

  // Renumber block IDs one at a time
  // We have to move them out of range, then back to range
  dof_id_type elem_count = 0;
  dof_id_type node_count = 0;

  // Assign the next available element ID to an element, skipping IDs that are still taken (which
  // can happen when not all the subdomains are being renumbered). The actual renumbering happens
  // later so we keep the mapping here.
  auto renumber_elem = [&](const Elem * elem)
  {
    // prevent re-renumbering
    if (new_elem_ids.count(elem->id()))
      return;
    while (mesh->query_elem_ptr(elem_count))
      elem_count++;
    // We can't mess with the range while looping in them
    new_elem_ids[elem->id()] = elem_count++;
  };

  // Assign the next available node ID to a node, skipping IDs that are still taken
  auto renumber_node = [&](Node & node)
  {
    // prevent re-renumbering
    if (renumbered_nodes.count(&node))
      return;
    while (mesh->query_node_ptr(node_count))
      node_count++;
    renumbered_nodes.insert(&node);
    mesh->renumber_node(node.id(), node_count++);
  };

  for (const auto i : index_range(block_ids))
  {
    const auto subdomain = block_ids[i];

    // Inner ordering: number the elements and nodes that lie on the requested boundaries first,
    // in the order the boundaries were listed
    for (const auto bid : boundary_order)
    {
      const auto & elems_on_bid = boundary_elems[bid];
      const auto & nodes_on_bid = boundary_nodes[bid];
      for (auto elem : mesh->active_subdomain_elements_ptr_range(subdomain))
      {
        if (elems_on_bid.count(elem->id()))
          renumber_elem(elem);
        for (auto & node : elem->node_ref_range())
          if (nodes_on_bid.count(&node))
            renumber_node(node);
      }
    }

    // Number the remaining (interior) elements and nodes of the subdomain
    for (auto elem : mesh->active_subdomain_elements_ptr_range(subdomain))
    {
      renumber_elem(elem);
      for (auto & node : elem->node_ref_range())
        renumber_node(node);
    }
  }

  // Now change the IDs
  for (const auto [key, value] : new_elem_ids)
    mesh->renumber_elem(key, value);

  // Update the max ids
  mesh->contract();

  // Try to disallow renumbering from now on since we just renumbered
  mesh->allow_renumbering(false);
  // No gain if exodus just renumbers this. Better to tell the user
  if (mesh->n_nodes() != mesh->max_node_id() || mesh->n_elem() != mesh->max_elem_id())
    mooseWarning("Mesh is not contiguously numbered after renumbering. The numbering may be erased "
                 "by outputs that require contiguous numbering such as Exodus.\nNumber of nodes: " +
                 std::to_string(mesh->n_nodes()) +
                 "\nMax node ID: " + std::to_string(mesh->max_node_id() - 1) +
                 "\nNumber of elements: " + std::to_string(mesh->n_elem()) +
                 "\nMax elem ID: " + std::to_string(mesh->max_elem_id() - 1));

  mesh->unset_is_prepared();
  return dynamic_pointer_cast<MeshBase>(mesh);
}
