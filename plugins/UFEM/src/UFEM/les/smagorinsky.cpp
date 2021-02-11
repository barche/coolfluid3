// Copyright (C) 2010-2011 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Based on WALE files

#include <cmath>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include "common/Component.hpp"
#include "common/Builder.hpp"
#include "common/FindComponents.hpp"
#include "common/OptionT.hpp"
#include "common/OptionArray.hpp"
#include "common/PropertyList.hpp"

#include "math/LSS/SolveLSS.hpp"
#include "math/LSS/ZeroLSS.hpp"
#include "mesh/LagrangeP1/ElementTypes.hpp"

#include "solver/actions/Iterate.hpp"
#include "solver/actions/NodeValence.hpp"
#include "solver/CriterionTime.hpp"
#include "solver/actions/AdvanceTime.hpp"
#include "solver/Time.hpp"
#include "solver/Tags.hpp"

#include "smagorinsky.hpp"
#include "../Tags.hpp"

namespace cf3 {
namespace UFEM {
namespace les {

using namespace solver::actions::Proto;
using boost::proto::lit;

////////////////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < Smagorinsky, common::Action, LibUFEMLES > Smagorinsky_builder;

////////////////////////////////////////////////////////////////////////////////////////////

Smagorinsky::Smagorinsky(const std::string& name) :
  ProtoAction(name),
  smagorinsky(boost::proto::as_child(m_smagorinsky_op))
{
  options().add("cs", m_smagorinsky_op.op.m_cs)
    .pretty_name("Cs")
    .description("Smagorinsky constant")
    .link_to(&m_smagorinsky_op.op.m_cs)
    .mark_basic();
    
  options().add("use_dynamic_smagorinsky", m_smagorinsky_op.op.use_dynamic_smagorinsky)
    .pretty_name("Activate dynamic smagorinsky")
    .description("Activate the dynamic smagorinsky parameter")
    .link_to(&m_smagorinsky_op.op.use_dynamic_smagorinsky)
    .mark_basic();

  options().add("use_anisotropic_correction", m_smagorinsky_op.op.use_anisotropic_correction)
    .pretty_name("Use anisotropic correction")
    .description("Correct the LES length scale for grid anisotropy")
    .link_to(&m_smagorinsky_op.op.use_anisotropic_correction)
    .mark_basic();

  options().add("use_mason_wallDamping", m_smagorinsky_op.op.use_mason_wallDamping)
    .pretty_name("Use Mason Wall Damping")
    .description("Reduce the LES length scale near the wall according to Mason")
    .link_to(&m_smagorinsky_op.op.use_mason_wallDamping)
    .mark_basic();

  options().add("masonCoef", m_smagorinsky_op.op.m_masonCoef)
    .pretty_name("n-damping function shape")
    .description("Integer that provides the coefficient or damping function shape for the Mason wall damping theory")
    .link_to(&m_smagorinsky_op.op.m_masonCoef)
    .mark_basic();

  options().add("initial_conditions", m_initial_conditions)
    .pretty_name("Initial Conditions")
    .description("The component that is used to manage the initial conditions in the solver this action belongs to")
    .link_to(&m_initial_conditions)
    .attach_trigger(boost::bind(&Smagorinsky::trigger_initial_conditions, this));
    
  options().add("velocity_tag", "navier_stokes_u_solution")
    .pretty_name("Velocity Tag")
    .description("Tag for the field containing the velocity")
    .attach_trigger(boost::bind(&Smagorinsky::trigger_set_expression, this));
    
  link_physics_constant("kappa", m_smagorinsky_op.op.m_kappa);
  link_physics_constant("z0", m_smagorinsky_op.op.m_z0);

  trigger_set_expression();
}

void Smagorinsky::trigger_set_expression()
{
  FieldVariable<0, VectorField> u("Velocity", options().value<std::string>("velocity_tag"));
  FieldVariable<1, ScalarField> nu_eff("EffectiveViscosity", "navier_stokes_viscosity");
  FieldVariable<2, ScalarField> valence("Valence", "node_valence");
  FieldVariable<3, ScalarField> cs_elts("cs", "cs"); // cs value for each element of the field
  FieldVariable<4, ScalarField> sgsLambda_elts("sgsLambda", "sgsLambda"); // sgsLambda value for each element of the field
  
  PhysicsConstant nu_visc("kinematic_viscosity");

  // List of applicable elements
  typedef boost::mpl::vector3<
    mesh::LagrangeP1::Hexa3D,
    mesh::LagrangeP1::Tetra3D,
    mesh::LagrangeP1::Prism3D
  > AllowedElementTypesT;

  m_reset_viscosity = create_component<ProtoAction>("ResetViscosity");
  m_reset_viscosity->set_expression(nodes_expression(
    group(
      nu_eff = 0.,
      cs_elts = 0.,
      sgsLambda_elts = 0.
      )
    ));

  m_node_connectivity = create_static_component<mesh::NodeConnectivity>("NodeConnectivity");
  m_smagorinsky_op.op.m_node_connectivity = m_node_connectivity.get();

  // m_set_node_velocity = create_component<ProtoAction>("NodeVelocity");
  // m_set_node_velocity->set_expression(nodes_expression(
  //     u = coordinates
  // ));

  set_expression(elements_expression(AllowedElementTypesT(), smagorinsky(u, nu_eff, valence, nu_visc, lit(m_smagorinsky_op.op.m_cs), cs_elts, sgsLambda_elts)));
}

void Smagorinsky::execute()
{
  m_reset_viscosity->execute();
  // m_set_node_velocity->execute();
  ProtoAction::execute();
}


void Smagorinsky::on_regions_set()
{
  if(is_not_null(m_node_valence))
  {
    m_node_valence->options().set("regions", options().option("regions").value());
  }
  m_reset_viscosity->options().set("regions", options().option("regions").value());
  // m_set_node_velocity->options().set("regions", options().option("regions").value());

  if(!m_loop_regions.empty())
  {
    mesh::Mesh& mesh = common::find_parent_component<mesh::Mesh>(*m_loop_regions.front());
    m_node_connectivity->initialize(common::find_components_recursively_with_filter<mesh::Elements>(mesh, mesh::IsElementsVolume()));
  }
}


void Smagorinsky::trigger_initial_conditions()
{
  if(is_null(m_initial_conditions))
    return;

  if(is_null(m_node_valence))
  {
    m_node_valence = m_initial_conditions->create_initial_condition("node_valence", "cf3.solver.actions.NodeValence");
    on_regions_set();
  }
}

} // les
} // UFEM
} // cf3
