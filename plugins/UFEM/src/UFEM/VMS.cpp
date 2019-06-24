/// Copyright (C) 2010-2011 von Karman Institute for Fluid Dynamics, Belgium
///
/// This software is distributed under the terms of the
/// GNU Lesser General Public License version 3 (LGPLv3).
/// See doc/lgpl.txt and doc/gpl.txt for the license text.

/// derived from ScalarAdvection.cpp, SUPG.hpp, NavierStokesAssembly.hpp, NavierStokes.cpp

/// Source: Variational Multiscale (VMS) Residual-based turbulence modeling 
/// for LES of incompressible flows (Y. Bazilevs, Elsevier, (2007) 173-201)

#include "VMS.hpp"

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include "common/Component.hpp"
#include "common/Builder.hpp"
#include "common/OptionT.hpp"
#include "common/OptionArray.hpp"
#include "common/PropertyList.hpp"

#include "math/LSS/SolveLSS.hpp"
#include "math/LSS/ZeroLSS.hpp"

#include "solver/actions/Proto/ProtoAction.hpp"
#include "solver/actions/Proto/Expression.hpp"
#include "solver/actions/Iterate.hpp"
#include "solver/CriterionTime.hpp"
#include "solver/actions/AdvanceTime.hpp"
#include "solver/Time.hpp"
#include "solver/Tags.hpp"

#include "Tags.hpp"

namespace cf3 {
namespace UFEM {

using namespace common;
using namespace solver;
using namespace solver::actions;
using namespace solver::actions::Proto;

using boost::proto::lit;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ComponentBuilder < VMS, LSSActionUnsteady, LibUFEM > VMS_builder;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

VMS::VMS(const std::string& name) :
  LSSActionUnsteady(name),
  m_alphaM(0.5), // Values from (43,84) p. 182 from Bazilevs 2007
  m_alphaF(0.5),
  m_gamma(0.5)
{
  options().add("alphaM", m_alphaM)
    .description("alphaM")
    .pretty_name("alphaM coefficient for the alpha-method.")
    .link_to(&m_alphaM)
    .mark_basic();

  options().add("alphaF", m_alphaF)
    .description("alphaF")
    .pretty_name("alphaF coefficient for the alpha-method.")
    .link_to(&m_alphaF)
    .mark_basic();

  options().add("gamma", m_gamma)
    .pretty_name("gamma")
    .description("Gamma coefficient for the alpha-method.")
    .link_to(&m_gamma);

  // options().add("scalar_name", "Temperature")
  //   .pretty_name("Scalar Name")
  //   .description("Internal (and default visible) name to use for the scalar")
  //   .attach_trigger(boost::bind(&VMS::trigger_scalar_name, this));

  // options().add("velocity_tag", "navier_stokes_solution")
  //   .pretty_name("Velocity Tag")
  //   .description("Tag for the velocity field")
  //   .attach_trigger(boost::bind(&VMS::trigger_scalar_name, this));

  set_solution_tag("variational_multiscale_solution");

  create_component<ProtoAction>("Predictor");
  create_component<ProtoAction>("CorrectorInitialiser");
  create_component<math::LSS::ZeroLSS>("ZeroLSS");
  create_component<ProtoAction>("Assembly");
  create_component<BoundaryConditions>("BoundaryConditions")->set_solution_tag(solution_tag());
  create_component<math::LSS::SolveLSS>("SolveLSS");
  create_component<ProtoAction>("Update");

  get_child("BoundaryConditions")->mark_basic();

  /// Set the default scalar name
  set_expression();
}

struct ComputeTauVMSImpl
{
  typedef void result_type;

  /// Compute the coefficients for the full Navier-Stokes equations
  template<typename UT, typename NUT>
  void operator()(const UT& u, const NUT& nu_eff, const Real& dt, Real& tau_ps, Real& tau_su, Real& tau_bulk) const
  {
    // Average viscosity
    const Real element_nu = fabs(detail::mean(nu_eff.value()));
    compute_coefficients(u, element_nu, dt, tau_ps, tau_su, tau_bulk);
  }
};

MakeSFOp<ComputeTauVMSImpl>::type compute_tau = {};


/// Convenience type for a compute_tau operation, grouping the stored operator and its proto counterpart
struct ComputeTauVMS
{
  ComputeTauVMS() :
    apply(boost::proto::as_child(data))
  {
  }
  
  // Stores the operator
  solver::actions::Proto::MakeSFOp<ComputeTauVMSImpl>::stored_type data;
  
  // Use as apply(velocity_field, nu_eff_field, dt, tau_ps, tau_su, tau_bulk)
  solver::actions::Proto::MakeSFOp<ComputeTauVMSImpl>::reference_type apply;
};

void VMS::set_expression()
{
  /// Make sure variables aren't registered multiple times
  if(is_not_null(m_physical_model))
  {
    if(is_not_null(m_physical_model->variable_manager().get_child(solution_tag())))
      m_physical_model->variable_manager().remove_component(solution_tag());
  }

  /// List of applicable elements
  typedef boost::mpl::vector<
    mesh::LagrangeP1::Quad2D,
    mesh::LagrangeP1::Triag2D
    // mesh::LagrangeP1::Triag2D,
    // mesh::LagrangeP1::Hexa3D,
    // mesh::LagrangeP1::Tetra3D,
    // mesh::LagrangeP1::Prism3D
  > AllowedElementTypesT;

  /// Scalar name is obtained from an option
  FieldVariable<1, VectorField> u("Velocity",options().value<std::string>("velocity_tag"));
  FieldVariable<2, VectorField> u_dot("Velocity_dot",options().value<std::string>("velocity_dot_tag"));
  FieldVariable<3, ScalarField> p("Pressure", "navier_stokes_pressure");
  FieldVariable<4, ScalarField> nu_eff("EffectiveViscosity", "navier_stokes_viscosity");
  // FieldVariable<0, ScalarField> phi(options().value<std::string>("scalar_name"), solution_tag());
  // PhysicsConstant nu_lam("kinematic_viscosity");
  // ConfigurableConstant<Real> relaxation_factor_scalar("relaxation_factor_scalar", "factor for relaxation in case of coupling", 1.);


  Handle<ProtoAction>(get_child("Predictor"))->set_expression( nodes_expression(
    group
    (
      u = u,
      u_dot = (m_gamma-1)/m_gamma * u_dot,
      p = p
    )
  ));

  Handle<ProtoAction>(get_child("CorrectorInitialiser"))->set_expression( nodes_expression(
    group
    (
      u_dot += m_alphaM * (u_dot - u_dot),
      u += m_alphaF * (u - u),
      p = p
    )
  ));

  /// Set the proto expression that handles the assembly
  Handle<ProtoAction>(get_child("Assembly"))->set_expression(elements_expression(
    AllowedElementTypesT(),
    group
    (
      _A = _0, _a = _0,
      // compute_tau(),
      // compute_tau(u, nu_eff, lit(dt()), lit(tau_su)),
      element_quadrature
      (
        /// K (p.183 - eq.102)
        _A(u[_i],u[_j]) += \
          m_alphaF * m_gamma * lit(dt()) * transpose(nabla(u)[_i]) * nu_eff * nabla(u)[_j] \
          + m_alphaF * m_gamma * lit(dt()) * transpose(nabla(u)[_i]) * tau_c * nabla(u)[_j],
        
        _A(u[_i],u[_i]) += m_alphaM * transpose(N(u)) * N(u) \
          + m_alphaM * (transpose(u * tau_m * nabla(u)) * N(u)) \
          + m_alphaF * m_gamma * lit(dt()) * transpose(N(u)) * u * nabla(u) \
          + m_alphaF * m_gamma * lit(dt()) * transpose(nabla(u) * nu_eff) * nabla(u) \
          + m_alphaF * m_gamma * lit(dt()) * transpose(u * nabla(u) * tau_m) * (u * nabla(u)),
        
        /// G (eq.104)
        _A(u[_i],p) +=  - transpose(nabla(u)[_i]) *  N(p) \
          + transpose(u[_i] * tau_m * nabla(u)) * nabla(p),
                
        /// D (eq.106)
        _A(p,u[_i]) += m_alphaF * m_gamma * lit(dt()) * transpose(N(p)) * nabla(u)[_i] \
          + m_alphaF * m_gamma * lit(dt()) * transpose(nabla(p) * tau_m) * u[_i] * nabla(u) \
          + m_alphaM * transpose(nabla(p)[_i] * tau_m) * N(u),
        
        /// L (eq.108)
        _A(p,p) += tau_m * transpose(nabla(p)) * nabla(p) // Continuity, PSPG
        
        /// RHS (Rm=_a[u],Rc=_a[p]; p.182 - eq.92-93)
        //_a[u[_i]] += -,
        //_a[p] += -,
      ),
      // system_matrix += invdt() * _T + m_theta * _A,
      system_matrix += _A,
      system_rhs += -_A * _x + _a
    )
  ));

  /// Set the proto expression for the update step
  Handle<ProtoAction>(get_child("Update"))->set_expression( nodes_expression(
    group(
      u_dot += (u_dot - u_dot1),
      u += m_gamma * lit(dt()) * u_dot,
      p += (p - p1),
      u_dot += solution(u_dot),
      u += solution(u),
      p += solution(p)
    )
  ));
}

void VMS::execute() /// derived from NavierStokesSemiImplicit.cpp
{
  typedef std::pair<Uint,Uint> BlockrowIdxT;
  
  // a->reset(0.);
  // delta_p_sum->reset(0.);
  // u->assign(*u_lss->solution());
  // p->assign(*p_lss->solution());    
  // for(Uint i = 0; i != nb_iterations; ++i)
  // {
  //   u_lss->rhs()->reset(0.);

  //   /// Computing m_velocity_assembly in the loop. Modified for the ABL implementation.
  //   u_lss->matrix()->reset(0.);
  //   m_velocity_assembly->execute();

  //   /// Velocity system: compute delta_a_star
  //   m_u_rhs_assembly->execute();
  //   if(i == 0) /// Apply velocity BC the first inner iteration
  //   {
  //     u_lss->rhs()->scale(m_time->dt());
  //     velocity_bc->execute();
  //     /// The velocity BC deals with velocity, so we need to write this in terms of acceleration
  //     u_lss->rhs()->scale(1./m_time->dt());
  //   }
  //   else /// Zero boundary condition after first pass. Modified for the ABL implementation.
  //   {
  //     velocity_bc->execute();
  //     u_lss->dirichlet_apply(true, true);
  //   }

  //   inner_bc->execute(); /// Inner Boundary Condition for the RHS implementation of ABL.
  //   // u_lss->matrix()->print(std::cout); /// Print m_velocity_assembly

  //   u_lss->solution()->reset(0.);
  //   CFdebug << "Solving velocity LSS..." << CFendl;
  //   solve_u_lss->execute();
  //   u_lss->solution()->sync();

  //   /// Pressure system: compute delta_p
  //   p_lss->rhs()->reset(0.);
  //   m_p_rhs_assembly->execute();
  //   p_lss->solution()->reset(0.);
  //   /// Apply BC if the first iteration, set RHS to 0 otherwise
  //   if(i ==0)
  //   {
  //     pressure_bc->execute();
  //   }
  //   else
  //   {
  //     BOOST_FOREACH(const BlockrowIdxT& diri_idx, Handle<math::LSS::TrilinosCrsMatrix>(p_lss->matrix())->get_dirichlet_nodes())
  //     {
  //       p_lss->rhs()->set_value(diri_idx.first, diri_idx.second, 0.);
  //     }
  //   }
  //   CFdebug << "Solving pressure LSS..." << CFendl;
  //   if(i == 0 && is_not_null(m_p_strategy_first))
  //   {
  //     p_lss->options().set("solution_strategy_component", m_p_strategy_first);
  //   }
  //   else if (i > 0 && is_not_null(m_p_strategy_second))
  //   {
  //     p_lss->options().set("solution_strategy_component", m_p_strategy_second);
  //   }
  //   solve_p_lss->execute();
  //   p_lss->solution()->sync();
    
  //   /// Compute delta_a
  //   u_lss->rhs()->reset(0.);
  //   m_apply_aup->execute(); /// Compute Aup*delta_p (stored in u_lss RHS)
  //   /// delta_a is delta_a_star for the dirichlet nodes
  //   BOOST_FOREACH(const BlockrowIdxT& diri_idx, Handle<math::LSS::TrilinosCrsMatrix>(u_lss->matrix())->get_dirichlet_nodes())
  //   {
  //     u_lss->rhs()->set_value(diri_idx.first, diri_idx.second, 0.);
  //   }
  //   Thyra::apply(*lumped_m_op, Thyra::NOTRANS, *aup_delta_p, delta_a.ptr(), -1., 1.); /// delta_a = delta_a_star - Ml_inv*Aup*delta_p
  //   u_lss->solution()->sync(); /// delta_a is a link to u_lss->solution(), so it needs a sync after matrix apply

  //   const math::LSS::Vector& da = *u_lss->solution();
  //   const math::LSS::Vector& dp = *p_lss->solution();
    
  //   a->update(da);
  //   u->update(da, m_time->dt());
  //   p->update(dp);
  //   delta_p_sum->update(dp);
  // }
    
  // u_lss->solution()->assign(*u);
  // p_lss->solution()->assign(*p);
}

} /// UFEM
} /// cf3
