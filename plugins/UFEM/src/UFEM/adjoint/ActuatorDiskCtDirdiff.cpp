// Copyright (C) 2010-2011 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <cmath>
#include "common/Core.hpp"
#include "common/FindComponents.hpp"
#include "common/Foreach.hpp"
#include "common/Log.hpp"
#include "common/OptionList.hpp"
#include "common/Signal.hpp"
#include "common/Builder.hpp"
#include "common/OptionT.hpp"
#include "common/EventHandler.hpp"

#include "math/LSS/System.hpp"

#include "mesh/Region.hpp"
#include "mesh/LagrangeP0/LibLagrangeP0.hpp"
#include "mesh/LagrangeP0/Quad.hpp"
#include "mesh/LagrangeP0/Line.hpp"

#include "ActuatorDiskCtDirdiff.hpp"

#include "solver/actions/Proto/SurfaceIntegration.hpp"
#include "solver/actions/Proto/VolumeIntegration.hpp"
#include "solver/actions/Proto/ProtoAction.hpp"
#include "solver/actions/Proto/Expression.hpp"

#include <iomanip>

namespace cf3
{

namespace UFEM
{

namespace adjoint
{

using namespace solver::actions::Proto;

////////////////////////////////////////////////////////////////////////////////////////////

common::ComponentBuilder < ActuatorDiskCtDirdiff, common::Action, LibUFEMAdjoint > ActuatorDiskCtDirdiff_Builder;

////////////////////////////////////////////////////////////////////////////////////////////

using boost::proto::lit;

ActuatorDiskCtDirdiff::ActuatorDiskCtDirdiff(const std::string& name) :
  UnsteadyAction(name),
  rhs(options().add("lss", Handle<math::LSS::System>())
    .pretty_name("LSS")
    .description("The linear system for which the boundary condition is applied")),
  system_matrix(options().option("lss"))
{
  options().add("area", m_area)
    .pretty_name("Area")
    .description("Area of the disk")
    .link_to(&m_area)
    .mark_basic(); // is this is enabled, the option can be accessed directly from Python, otherwise .options is needed
  
  options().add("ct", m_ct)
    .pretty_name("Ct")
    .description("Trust coefficient of the disk")
    .link_to(&m_ct)
    .mark_basic(); // if ct < 0 --> a in function of u_mean_disk

  options().add("th", m_th)
    .pretty_name("Thickness")
    .description("Thickness of the disk")
    .link_to(&m_th)
    .mark_basic();
	
  options().add("result", Real())
      .pretty_name("Result")
      .description("Result of the integration (read-only)")
      .mark_basic();

  
  // The component that  will set the force
  create_static_component<ProtoAction>("SetDirdiffForce")->options().option("regions").add_tag("norecurse");
  create_static_component<ProtoAction>("SetSensUDisk")->options().option("regions").add_tag("norecurse");

  // Initialize the expression
  trigger_setup();
}

ActuatorDiskCtDirdiff::~ActuatorDiskCtDirdiff()
{
}

void ActuatorDiskCtDirdiff::on_regions_set()
{
  auto regions = options()["regions"].value<std::vector<common::URI>>();

  // Set the regions when the option is set
  get_child("SetDirdiffForce")->options().set("regions", std::vector<common::URI>({regions[0]}));
  get_child("SetSensUDisk")->options().set("regions", std::vector<common::URI>({regions[0]}));
}

void ActuatorDiskCtDirdiff::trigger_setup()
{
  FieldVariable<0, VectorField> SensF("SensForce", "sensitivity_body_force");
  FieldVariable<1, VectorField> SensU("SensU", "sensitivity_solution");
  FieldVariable<2, ScalarField> Ct("ThrustCoefficient", "actuator_disk");
  FieldVariable<3, VectorField> uDisk("MeanDiskSpeed", "actuator_disk");
  FieldVariable<4, VectorField> SensUDisk("MeanDirdiffDiskSpeed","dirdiff_actuator_disk");

  Handle<ProtoAction> set_SensUDisk(get_child("SetSensUDisk"));
  set_SensUDisk->set_expression(nodes_expression
  (
    group
    (
      SensUDisk[0] = lit(m_SensU_mean_disk)
    )
  ));

  Handle<ProtoAction> set_dirdiff_force(get_child("SetDirdiffForce"));
  set_dirdiff_force->set_expression(nodes_expression
  (
    group
    (
      SensF[0] = 0.5 * uDisk[0] * uDisk[0] / lit(m_th)
     )
  ));

}

void ActuatorDiskCtDirdiff::execute()
{
  FieldVariable<0, VectorField> SensU("SensU", "sensitivity_solution");
  m_SensU_mean_disk = 0;

  // surface_integral(m_SensU_mean_disk, std::vector<Handle<mesh::Region>>({m_loop_regions[1]}), ((SensU*normal)[0]));

  boost::mpl::vector<mesh::LagrangeP1::Triag2D, mesh::LagrangeP1::Tetra3D, mesh::LagrangeP1::Quad2D> etypes;
  volume_integral(m_SensU_mean_disk, std::vector<Handle<mesh::Region>>({m_loop_regions[0]}), U[0], etypes);
  m_SensU_mean_disk /= (m_area * m_th);
  
  // m_SensU_mean_disk /= m_area;

  // m_F = -m_U_mean_disk * 
  // CFinfo << std::setprecision(20) <<"force set to " << m_f << ", a: " << m_a << "m_U_mean_disk :" << m_U_mean_disk <<  " pow2 " << m_u_mean_disk2 << " pow3 " << m_u_mean_disk3 << CFendl;
  options().set("result", m_SensU_mean_disk);
  Handle<ProtoAction> set_SensUDisk(get_child("SetSensUDisk"));
  Handle<ProtoAction> set_dirdiff_force(get_child("SetDirdiffForce"));
  set_dirdiff_force->execute();
  set_SensUDisk->execute();

}

} // namespace adjoint

} // namespace UFEM

} // namespace cf3
