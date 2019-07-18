/// Copyright (C) 2010-2011 von Karman Institute for Fluid Dynamics, Belgium
///
/// This software is distributed under the terms of the
/// GNU Lesser General Public License version 3 (LGPLv3).
/// See doc/lgpl.txt and doc/gpl.txt for the license text.

/// derived from ScalarAdvection.cpp, SUPG.hpp, NavierStokesAssembly.hpp, NavierStokes.cpp

/// Source: Variational Multiscale (VMS) Residual-based turbulence modeling 
/// for LES of incompressible flows (Y. Bazilevs, Elsevier, (2007) 173-201)

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include "common/Component.hpp"
#include "common/Builder.hpp"
#include "common/OptionT.hpp"
#include "common/OptionArray.hpp"
#include "common/PropertyList.hpp"

#include "math/LSS/SolveLSS.hpp"
#include "math/LSS/ZeroLSS.hpp"
#include "math/LSS/Trilinos/TrilinosCrsMatrix.hpp"
#include "math/LSS/Trilinos/TekoBlockedOperator.hpp"
// #include "math/LSS/Trilinos/TrilinosVector.hpp"

#include "solver/actions/Proto/ProtoAction.hpp"
#include "solver/actions/Proto/Expression.hpp"
#include "solver/actions/Iterate.hpp"
#include "solver/CriterionTime.hpp"
#include "solver/actions/AdvanceTime.hpp"
#include "solver/Time.hpp"
#include "solver/Tags.hpp"

#include "Tags.hpp"

#include "VMS.hpp"

namespace cf3 {
namespace UFEM {

using namespace common;
using namespace solver;
using namespace solver::actions;
using namespace solver::actions::Proto;

using boost::proto::lit;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// needed for each class being defined
ComponentBuilder < VMS, LSSActionUnsteady, LibUFEM > VMS_builder;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

VMS::VMS(const std::string& name) :
  LSSActionUnsteady(name),
  m_alphaM(0.5), // Values from (43,84) p. 182 from Bazilevs 2007
  m_alphaF(0.5),
  m_gamma(0.5),
  nb_iterations(2),
  u1("u1", "velocityNew"), // computed
  u1Dot("u1Dot", "accelNew"),
  p1("p1", "pressureNew"), // computed
  uaF("uaF", "velocityAlphaF"), // computed
  uaMDot("uaMDot", "accelAlphaM"), // computed
  f("f","source_term")
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

  options().add("nb_iterations", nb_iterations)
    .pretty_name("number of iterations")
    .description("Number of inner loop iterations");
  
  set_solution_tag("vms_solution");

  predictor = create_component<ProtoAction>("Predictor");
  correctorInitialiser = create_component<ProtoAction>("CorrectorInitialiser");
  zeroLSS = create_component<math::LSS::ZeroLSS>("ZeroLSS");
  assembly = create_component<ProtoAction>("Assembly");
  /// Handle to convert ProtoAction to BoundaryCondition
  bc = create_component<BoundaryConditions>("BoundaryConditions");
  bc->mark_basic();
  bc->set_solution_tag(solution_tag());
  solveLSS = create_component<math::LSS::SolveLSS>("SolveLSS");
  update = create_component<ProtoAction>("Update");

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
  FieldVariable<1, VectorField> u("u", solution_tag());
  FieldVariable<2, VectorField> uDot("uDot","vms_uDot");
  FieldVariable<3, ScalarField> p("p", solution_tag());
  FieldVariable<4, ScalarField> nu_eff("EffectiveViscosity","navier_stokes_viscosity");
  
  /// Solution fields computed from assembly
  FieldVariable<11, VectorField> Du1Dot("Du1Dot","vms_solution1");
  FieldVariable<12, ScalarField> Dp1("Dp1","vms_solution1");

  predictor->set_expression( nodes_expression(
    group(
      u1 = u,
      u1Dot = (m_gamma-1)/m_gamma * uDot,
      p1 = p,
      _cout << "yop: u1: " << u1 << ", u1Dot: " << u1Dot << "\n"
    )
  ));

  correctorInitialiser->set_expression( nodes_expression(
    group(
      uaMDot = uDot + m_alphaM * (u1Dot - uDot),
      uaF = u + m_alphaF * (u1 - u),
      p1 = p1,
      _cout << "yop: uaMDot: " << uaMDot << ", uaF: " << uaF << "\n"
    )
  ));

  /// Set the proto expression that handles the assembly
  assembly->set_expression( elements_expression(
    AllowedElementTypesT(), 
    group(
      _A = _0, _a = _0,
      // compute_tau(),
      // compute_tau(u, nu_eff, lit(dt()), lit(tau_su)),

/*       element_quadrature
      (
        /// K (p.183 - eq.102)
        _A(Du1Dot[_i],Du1Dot[_j]) += \
          m_alphaF * m_gamma * lit(dt()) * transpose(nabla(Du1Dot)[_i]) * nu_eff * nabla(Du1Dot)[_j] \
          + m_alphaF * m_gamma * lit(dt()) * transpose(nabla(Du1Dot)[_i]) * tau_c * nabla(Du1Dot)[_j],
        
        _A(Du1Dot[_i],Du1Dot[_i]) += m_alphaM * transpose(N(Du1Dot)) * N(Du1Dot) \
          + m_alphaM * (transpose(Du1Dot * tau_m * nabla(Du1Dot)) * N(Du1Dot)) \
          + m_alphaF * m_gamma * lit(dt()) * transpose(N(Du1Dot)) * Du1Dot * nabla(Du1Dot) \
          + m_alphaF * m_gamma * lit(dt()) * transpose(nabla(Du1Dot) * nu_eff) * nabla(Du1Dot) \
          + m_alphaF * m_gamma * lit(dt()) * transpose(Du1Dot * nabla(Du1Dot) * tau_m) * (Du1Dot * nabla(Du1Dot)),
        
        /// G (eq.104)
        _A(Du1Dot[_i],Dp1) +=  - transpose(nabla(Du1Dot)[_i]) *  N(Dp1) \
          + transpose(Du1Dot[_i] * tau_m * nabla(Du1Dot)) * nabla(Dp1),
                
        /// D (eq.106)
        _A(Dp1,Du1Dot[_i]) += m_alphaF * m_gamma * lit(dt()) * transpose(N(Dp1)) * nabla(Du1Dot)[_i] \
          + m_alphaF * m_gamma * lit(dt()) * transpose(nabla(Dp1) * tau_m) * Du1Dot[_i] * nabla(Du1Dot) \
          + m_alphaM * transpose(nabla(Dp1)[_i] * tau_m) * N(Du1Dot),
        
        /// L (eq.108) // Continuity, PSPG
        _A(Dp1,Dp1) += transpose(nabla(Dp1)) * tau_m * nabla(Dp1)
        
        /// RHS (Rm=_a[u],Rc=_a[p]; p.182 - eq.92-93) equals 0 => _a = _0
      ),
      // system_matrix += invdt() * _T + m_theta * _A,
      system_matrix += _A,
      system_rhs += -_A * _x + _a
 */

      element_quadrature // Integration over the element
      (
        _A(u[_i],u[_i]) += transpose(nabla(u)) * nabla(u),
        _a[u[_i]] += transpose(N(u)) * f,
        _A(p,p) += transpose(nabla(p)) * nabla(p),
        _a[p] += transpose(N(p)) * f

      ),
      system_matrix +=  _A, // Assemble into the global linear system
      system_rhs += _a,
      _cout << "_A= " << _A << "\n",
      _cout << "_a= " << _a << "\n"
    )
  ));
 
  /// Set the proto expression for the update step
  update->set_expression( nodes_expression(
    group(
      u = solution(u),
      p = solution(p),
      _cout << "yop: u: " << u << ", p: " << p << "\n"
    )
  ));

/*
  /// Set the proto expression for the update step
  update->set_expression( nodes_expression(
    group(
      u1Dot += solution(Du1Dot),
      u1 += m_gamma * lit(dt()) * solution(Du1Dot),
      p1 += solution(Dp1)
    )
  ));
*/
}

void VMS::execute() /// derived from NavierStokesSemiImplicit.cpp
{
  // typedef std::pair<Uint,Uint> BlockrowIdxT;
  
  // predictor->execute();
  for(Uint i = 0; i != nb_iterations; ++i)
  {
    // correctorInitialiser->execute();
    zeroLSS->execute();
    assembly->execute();
    Handle<math::LSS::System> lss(get_child("LSS"));
    if(i == 0) /// Apply velocity BC the first inner iteration
    {
      lss->rhs()->scale(dt());
      bc->execute();
      /// The velocity BC deals with velocity, so we need to write this in terms of acceleration
      lss->rhs()->scale(1./dt());
    }
    else /// Zero boundary condition after first pass
    {
      bc->execute();
      lss->dirichlet_apply(true, true);
    }
    solveLSS->execute();
    update->execute();
  }
}

} /// UFEM
} /// cf3