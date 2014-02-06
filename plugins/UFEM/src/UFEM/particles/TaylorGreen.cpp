// Copyright (C) 2010-2011 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/array.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>

#include "common/Core.hpp"
#include "common/FindComponents.hpp"
#include "common/Foreach.hpp"
#include "common/Log.hpp"
#include "common/OptionList.hpp"
#include "common/PE/Comm.hpp"
#include "common/Signal.hpp"
#include "common/Builder.hpp"
#include "common/OptionT.hpp"
#include "common/EventHandler.hpp"

#include "math/Consts.hpp"

#include "TaylorGreen.hpp"

#include "solver/actions/Proto/ProtoAction.hpp"
#include "solver/actions/Proto/Expression.hpp"
#include "solver/Tags.hpp"

namespace cf3
{

namespace UFEM
{

namespace particles
{

using namespace solver::actions::Proto;
using boost::proto::lit;

namespace detail
{

using math::Consts::pi;

inline Real square(const Real x)
{
  return x*x;
}

Real TaylorGreenModel::ux() const
{
  return Ua - (Vs*::cos((pi()*(-(t*Ua) + x))/D)*::sin((pi()*(-(t*Va) + y))/D))/::exp((2*nu*pi()*pi()*t)/D*D);
}

Real TaylorGreenModel::uy() const
{
  return Va + (Vs*::cos((pi()*(-(t*Va) + y))/D)*::sin((pi()*(-(t*Ua) + x))/D))/::exp((2*nu*pi()*pi()*t)/D*D);
}

Real TaylorGreenModel::vx() const
{
  return (32*(D*D*D)*::exp((6*nu*(pi()*pi())*t)/(D*D))*(rhop*rhop)*Ua + 4*::exp((2*nu*(pi()*pi())*t)/(D*D))*pi()*taup*(Vs*Vs)*(D*pi()*square(rhof + 2*rhop)*taup*Ua*(::cos((2*pi()*(t*Ua - x))/D)
    + ::cos((2*pi()*(t*Va - y))/D)) + 4*(rhof - rhop)*((D*D)*rhop + nu*(pi()*pi())*(rhof + 2*rhop)*taup)*::sin((2*pi()*(t*Ua - x))/D))
    + 16*D*::exp((4*nu*(pi()*pi())*t)/(D*D))*rhop*((D*D)*rhop + 2*nu*(pi()*pi())*(-rhof + rhop)*taup)*Vs*(::sin((pi()*(t*(Ua + Va) - x - y))/D)
    - ::sin((pi()*(t*(Ua - Va) - x + y))/D)) + 3*D*(pi()*pi())*rhof*(rhof + 2*rhop)*(taup*taup)*(Vs*Vs*Vs)*(::sin((pi()*(t*(Ua + 3*Va) - x - 3*y))/D)
    + ::sin((pi()*(3*t*Ua + t*Va - 3*x - y))/D) - ::sin((pi()*(3*t*Ua - t*Va - 3*x + y))/D) - ::sin((pi()*(t*Ua - 3*t*Va - x + 3*y))/D)))
    / (32*(D*D*D)*::exp((6*nu*(pi()*pi())*t)/(D*D))*(rhop*rhop) + 4*D*::exp((2*nu*(pi()*pi())*t)/(D*D))*(pi()*pi())*square(rhof + 2*rhop)*(taup*taup)*(Vs*Vs)*(::cos((2*pi()*(-(t*Ua) + x))/D)
    + ::cos((2*pi()*(-(t*Va) + y))/D)));
}

Real TaylorGreenModel::vy() const
{
  return (32*(D*D*D)*::exp((6*nu*(pi()*pi())*t)/(D*D))*(rhop*rhop)*Va + 4*::exp((2*nu*(pi()*pi())*t)/(D*D))*pi()*taup*(Vs*Vs)*(D*pi()*square(rhof + 2*rhop)*taup*Va*(::cos((2*pi()*(t*Ua - x))/D)
    + ::cos((2*pi()*(t*Va - y))/D)) + 4*(rhof - rhop)*((D*D)*rhop + nu*(pi()*pi())*(rhof + 2*rhop)*taup)*::sin((2*pi()*(t*Va - y))/D))
    - 16*D*::exp((4*nu*(pi()*pi())*t)/(D*D))*rhop*((D*D)*rhop + 2*nu*(pi()*pi())*(-rhof + rhop)*taup)*Vs*(::sin((pi()*(t*(Ua + Va) - x - y))/D)
    + ::sin((pi()*(t*(Ua - Va) - x + y))/D)) - 3*D*(pi()*pi())*rhof*(rhof + 2*rhop)*(taup*taup)*(Vs*Vs*Vs)*(::sin((pi()*(t*(Ua + 3*Va) - x - 3*y))/D)
    + ::sin((pi()*(3*t*Ua + t*Va - 3*x - y))/D) + ::sin((pi()*(3*t*Ua - t*Va - 3*x + y))/D) + ::sin((pi()*(t*Ua - 3*t*Va - x + 3*y))/D)))
    / (32*(D*D*D)*::exp((6*nu*(pi()*pi())*t)/(D*D))*(rhop*rhop) + 4*D*::exp((2*nu*(pi()*pi())*t)/(D*D))*(pi()*pi())*square(rhof + 2*rhop)*(taup*taup)*(Vs*Vs)*(::cos((2*pi()*(-(t*Ua) + x))/D)
    + ::cos((2*pi()*(-(t*Va) + y))/D)));
}

void update_tg_values_func(const TaylorGreenModel& d, RealVector& out)
{
  out[0] = d.ux();
  out[1] = d.uy();
  out[2] = d.vx();
  out[3] = d.vy();
}
static boost::proto::terminal< void(*)(const TaylorGreenModel&, RealVector&) >::type const update_tg_values = {&update_tg_values_func};

}

common::ComponentBuilder < TaylorGreen, common::Action, LibUFEMParticles > TaylorGreen_Builder;

TaylorGreen::TaylorGreen(const std::string& name) :
  ProtoAction(name),
  m_tg_values(4)
{
  options().add(solver::Tags::time(), m_time)
    .pretty_name("Time")
    .description("Component that keeps track of time for this simulation")
    .attach_trigger(boost::bind(&TaylorGreen::trigger_time, this))
    .link_to(&m_time);
    
  options().add("ua", m_tg_model.Ua)
    .pretty_name("Ua")
    .description("Constant X velocity")
    .link_to(&(m_tg_model.Ua));
      
  options().add("va", m_tg_model.Va)
    .pretty_name("Va")
    .description("Constant Y velocity")
    .link_to(&(m_tg_model.Va));

  options().add("vs", m_tg_model.Vs)
    .pretty_name("Vs")
    .description("Vortex velocity")
    .link_to(&(m_tg_model.Vs));

  options().add("d", m_tg_model.D)
    .pretty_name("D")
    .description("Vortex diameter")
    .link_to(&(m_tg_model.D));

  options().add("tau_p", m_tg_model.taup)
    .pretty_name("Tau_p")
    .description("Particle relaxation time")
    .link_to(&(m_tg_model.taup));

  options().add("rho_p", m_tg_model.rhop)
    .pretty_name("Rho_p")
    .description("Particle density")
    .link_to(&(m_tg_model.rhop));
    
  FieldVariable<0, VectorField> u("FluidVelocityTG", "taylor_green");
  FieldVariable<1, VectorField> v("ParticleVelocityTG", "taylor_green");
  PhysicsConstant nu("kinematic_viscosity");
  PhysicsConstant rho_f("density");
  
  set_expression(nodes_expression_2d
  (
    group
    (
      lit(m_tg_model.nu) = nu,
      lit(m_tg_model.rhof) = rho_f,
      lit(m_tg_model.x) = coordinates[0],
      lit(m_tg_model.y) = coordinates[1],
      detail::update_tg_values(lit(m_tg_model), lit(m_tg_values)),
      u[0] = lit(m_tg_values)[0],
      u[1] = lit(m_tg_values)[1],
      v[0] = lit(m_tg_values)[2],
      v[1] = lit(m_tg_values)[3]
    )
  ));
}

TaylorGreen::~TaylorGreen()
{
}


void TaylorGreen::trigger_time()
{
  if(is_null(m_time))
      return;

  m_time->options().option("current_time").attach_trigger(boost::bind(&TaylorGreen::trigger_current_time, this));
  trigger_current_time();
}

void TaylorGreen::trigger_current_time()
{
  m_tg_model.t = m_time->current_time();
}

} // namespace particles

} // namespace UFEM

} // namespace cf3
