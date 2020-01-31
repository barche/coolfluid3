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

////////////////////////////////////////////////////////////////////////////////
/// needed for each class being defined
ComponentBuilder < VMS, LSSActionUnsteady, LibUFEM > VMS_builder;
////////////////////////////////////////////////////////////////////////////////

VMS::VMS(const std::string& name) :
  LSSActionUnsteady(name),
  m_alphaM(0.5), // Values from (43,84) p. 182 from Bazilevs 2007
  m_alphaF(0.5),
  m_gamma(0.5),
  m_tau_m(0.),
  m_tau_c(0.),
  m_dynTau(1),
  m_nb_iterations(1),
  m_c1(1.),
  m_uDot("uDot","accel"),
  m_u1("u1", "velocityNew"), // computed
  m_uaF("uaF", "velocityAlphaF"), // computed
  m_uaMDot("uaMDot", "accelAlphaM"), // computed
  m_f("f","source_term")
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

  options().add("dynTau", m_dynTau)
    .pretty_name("dynTau")
    .description("Compute stabilization coefficients tauM & tauF.")
    .link_to(&m_dynTau);

  options().add("nb_iterations", m_nb_iterations)
    .pretty_name("number of iterations")
    .description("Number of inner loop iterations")
    .link_to(&m_nb_iterations);
  
  set_solution_tag("vms_solution");

  m_predictor = create_component<ProtoAction>("Predictor");
  m_correctorInitialiser =create_component<ProtoAction>("CorrectorInitialiser");
  m_assembly = create_component<ProtoAction>("Assembly");
  /// Handle to convert ProtoAction to BoundaryCondition
  m_bc = create_component<BoundaryConditions>("BoundaryConditions");
  m_bc->mark_basic();
  m_bc->set_solution_tag(solution_tag());
  m_solveLSS = create_component<math::LSS::SolveLSS>("SolveLSS");
  m_update = create_component<ProtoAction>("Update");
  m_renew = create_component<ProtoAction>("Renew");

  /// Set the default scalar name
  set_expression();
}

////////////////////////////////////////////////////////////////////////////////
struct ComputeTauVMSImpl : boost::noncopyable
{
  typedef void result_type;

  /// Compute the coefficients for the full Navier-Stokes equations
  template<typename UT, typename NUT>
  void operator()(const UT& u, const NUT& nu_eff, const Real& dt, Real& m_tau_c, Real& m_tau_m, Real& m_c1, const bool m_dynTau) const
  {
    // Average viscosity
    const Real element_nu = fabs(detail::mean(nu_eff.value()));
    compute_coefficients(u, element_nu, dt, m_tau_c, m_tau_m, m_c1, m_dynTau);
  }

  // Compute the VMS stabilization coefficients
  template<typename UT>
  void compute_coefficients(const UT& u, const Real element_nu, const Real& dt, Real& m_tau_c, Real& m_tau_m, Real& m_c1, const bool m_dynTau) const
  {
    typedef typename UT::EtypeT ElementT;
    static const Uint dim = ElementT::dimension;
    typedef mesh::Integrators::GaussMappedCoords<1, ElementT::shape> GaussT;

    // Compute u in the cell's center (~Gauss node average)
    u.compute_values(GaussT::instance().coords.col(0));
    std::cout << "yop: u_cellCenter: " << std::endl << u.eval() << std::endl;
    u.support().compute_jacobian(GaussT::instance().coords.col(0));
    
    const Eigen::Matrix<Real, ElementT::dimensionality, ElementT::dimensionality> Jinv = u.support().jacobian_inverse();
    const Eigen::Matrix<Real, 1, ElementT::dimensionality> gi = Jinv.colwise().sum();
    const Real gg = gi*gi.transpose();
    
    const Eigen::Matrix<Real, ElementT::dimensionality, ElementT::dimensionality> Gij = (Jinv * Jinv.transpose());
    const Real GG = Gij.cwiseProduct(Gij).sum();
    const Real uGu = u.eval() * Gij * detail::transpose(u.eval());

    if(m_dynTau) {
      m_tau_m = 1. / sqrt((4./(dt*dt)) + uGu + m_c1 * element_nu*element_nu * GG);
      m_tau_c = 1. / (m_tau_m * gg);
    }
    else
    {
      m_tau_m = 0.;
      m_tau_c = 0.;
    }

    // std::cout << "yop: **********: " << std::endl;
    // std::cout << "yop: Jinv: " << std::endl << Jinv << std::endl;
    // std::cout << "yop: gi: " << std::endl << gi << std::endl;
    // std::cout << "yop: gg: " << std::endl << gg << std::endl;
    // std::cout << "yop: Gij: " << std::endl << Gij << std::endl;
    // std::cout << "yop: GG: " << std::endl << GG << std::endl;
    // std::cout << "yop: u.eval(): " << std::endl << u.eval() << std::endl;
    // // std::cout << "yop: typeid(u.eval()).name(): " << std::endl << typeid(u.eval()).name() << std::endl;
    // std::cout << "yop: uGu: " << std::endl << uGu << std::endl;
    std::cout << "yop: m_tau_m: " << m_tau_m << std::endl;
    std::cout << "yop: m_tau_c: " << m_tau_c << std::endl;
  }
};

////////////////////////////////////////////////////////////////////////////////
/// Convenience type for a compute_tau operation, grouping the stored operator and its proto counterpart
struct ComputeTauVMS
{
  ComputeTauVMS() :
    apply(boost::proto::as_child(data))
  {
  }
  
  // Stores the operator
  solver::actions::Proto::MakeSFOp<ComputeTauVMSImpl>::stored_type data;
  
  // Use as apply(velocity_field, nu_eff_field, dt, m_tau_c, m_tau_m, m_c1)
  solver::actions::Proto::MakeSFOp<ComputeTauVMSImpl>::reference_type apply;
};

////////////////////////////////////////////////////////////////////////////////
void VMS::on_initial_conditions_set(InitialConditions& initial_conditions)
{
  m_initial_conditions = initial_conditions.create_initial_condition(solution_tag());

  // Use a proto action to set the linearized_velocity easily
  Handle<ProtoAction> accel_ic (initial_conditions.create_initial_condition("uDot", "cf3.solver.ProtoAction"));
  accel_ic->set_expression(nodes_expression(group(
    // initialize uDot=0
    m_uDot[_i] = 1e-3
    //_cout << "yop: m_uDot= " << transpose(m_uDot) << "\n"
    )));
}

////////////////////////////////////////////////////////////////////////////////
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
  FieldVariable<1, VectorField> u("velocity", "vms_u");
  FieldVariable<3, ScalarField> p("pressure", "vms_p");
  FieldVariable<4, ScalarField> nu_eff("effViscosity","vms_nu");
  
  /// Solution fields computed from assembly
  FieldVariable<5, VectorField> m_u1Dot("m_u1Dot", solution_tag());
  FieldVariable<6, ScalarField> m_p1("m_p1", solution_tag());

  if(is_not_null(m_initial_conditions))
  {
    Handle<InitialConditions> solver_ic(m_initial_conditions->parent());
    cf3_assert(is_not_null(solver_ic));
    solver_ic->remove_component(*m_initial_conditions);
    m_initial_conditions = solver_ic->create_initial_condition(solution_tag());
  }

  m_predictor->set_expression( nodes_expression(
    group(
      m_u1 = u,
      m_u1Dot = (lit(m_gamma)-1)/lit(m_gamma) * m_uDot,
      m_p1 = p,
      _cout << "yop - predictor: m_u1: " << transpose(m_u1) << ", m_u1Dot: " << transpose(m_u1Dot) << ", m_uDot: " << transpose(m_uDot) << ", m_p1: " << transpose(m_p1) << "\n"
    )
  ));
  // std::cout << predictor->tree() << std::endl;

  m_correctorInitialiser->set_expression( nodes_expression(
    group(
      m_uaMDot = m_uDot + lit(m_alphaM) * (m_u1Dot - m_uDot),
      m_uaF = u + lit(m_alphaF) * (m_u1 - u),
      m_p1 = m_p1,
      _cout << "yop - corrector: lit(m_alphaM):" << lit(m_alphaM) << ", lit(m_alphaF):" << lit(m_alphaF) << ",  m_uaMDot: " << transpose(m_uaMDot) << ", m_uaF: " << transpose(m_uaF) << ", m_u1Dot: " << transpose(m_u1Dot) <<  ", m_p1: " << transpose(m_p1) << "\n",
      _cout << "nu_eff= " << nu_eff << " ; u= " << transpose(u) << "\n"

    )
  ));

  ComputeTauVMS compute_tau;

  /// Set the proto expression that handles the assembly
  m_assembly->set_expression( elements_expression(
    AllowedElementTypesT(), 
    group(
      _A = _0, _a = _0,
      compute_tau.apply(u, nu_eff, lit(dt()), lit(m_tau_c), lit(m_tau_m), lit(m_c1), lit(m_dynTau)),

      element_quadrature
      (
        /// K (p.183 - eq.102)
        _A(m_u1Dot[_i],m_u1Dot[_i]) += lit(m_alphaM) * transpose(N(m_u1Dot)) * N(m_u1Dot) * lit(0) \
          + lit(m_alphaM) * (transpose(u * m_tau_m * nabla(m_u1Dot)) * N(m_u1Dot)) \
          + lit(m_alphaF) * lit(m_gamma) * lit(dt()) * transpose(N(m_u1Dot)) * u * nabla(m_u1Dot) \
          + lit(m_alphaF) * lit(m_gamma) * lit(dt()) * transpose(nabla(m_u1Dot) * nu_eff) * nabla(m_u1Dot) \
          + lit(m_alphaF) * lit(m_gamma) * lit(dt()) * transpose(u * nabla(m_u1Dot) * m_tau_m) * (u * nabla(m_u1Dot)),

        // _A(m_u1Dot[_i],m_u1Dot[_j]) += \
        //   lit(m_alphaF) * lit(m_gamma) * lit(dt()) * transpose(nabla(m_u1Dot)[_j]) * nu_eff * nabla(m_u1Dot)[_i] \
        //   + lit(m_alphaF) * lit(m_gamma) * lit(dt()) * transpose(nabla(m_u1Dot)[_i]) * m_tau_c * nabla(m_u1Dot)[_j],
        
        /// G (eq.104)
        // _A(m_u1Dot[_i],m_p1) +=  - transpose(nabla(m_u1Dot)[_i]) *  N(m_p1) \
        //   + transpose(transpose(u) * m_tau_m * nabla(m_u1Dot)[_i]) * nabla(m_p1),
        //? Transposed the complete matrix to have 0 on the last column !?
        //? minus added in front of the whole expr to correspond to SUPG !?
        _A(m_u1Dot[_i],m_p1) +=  - transpose( \
          -transpose(nabla(m_u1Dot)[_i]) *  N(m_p1) \
          + transpose(transpose(u) * m_tau_m * nabla(m_u1Dot)[_i]) * nabla(m_p1) \
          ),

        /// D (eq.106)
        _A(m_p1,m_u1Dot[_i]) += lit(m_alphaF) * lit(m_gamma) * lit(dt()) * transpose(N(m_p1)) * nabla(m_u1Dot)[_i] \
          + lit(m_alphaF) * lit(m_gamma) * lit(dt()) * transpose(nabla(m_p1)[_i] * m_tau_m) * (u * nabla(m_u1Dot)) \
          + lit(m_alphaM) * transpose(nabla(m_p1)[_i] * m_tau_m) * N(m_u1Dot),


          // + lit(m_alphaF) * lit(m_gamma) * lit(dt()) * transpose(nabla(m_p1) * m_tau_m) * (transpose(u) * nabla(m_u1Dot)[_i]) \



        /// L (eq.108) // Continuity, PSPG
        _A(m_p1,m_p1) += transpose(nabla(m_p1)) * m_tau_m * nabla(m_p1)

        /// RHS doesn't need to be expressed because equals zero
        /// RHS (Rm=_a[u],Rc=_a[p]; p.182 - eq.92-93) equals 0 => _a = _0
      ),
      /// system_matrix and _rhs expressed properly
      // system_matrix += invdt() * _T + m_theta * _A,
      system_matrix +=  _A, // Assemble into the global linear system
      system_rhs += -_A * _x + _a,
      // system_rhs += _a,
      _cout << "_A= \n" << _A << "\n\n" \
            << " ; lit(m_alphaF)= " << lit(m_alphaF) \
            << " ; lit(m_gamma)= " << lit(m_gamma) \
            << " ; lit(dt())= " << lit(dt()) << "\n" \
            << "_a= " << transpose(_a) << "\n"
            << "visited element " << element_index << "\n"
    )
  ));
 
  /// Set the proto expression for the update step
  m_update->set_expression( nodes_expression(
    group(
      m_u1Dot += solution(m_u1Dot),
      m_u1 += lit(m_gamma) * lit(dt()) * solution(m_u1Dot),
      m_p1 += solution(m_p1),
      _cout << "yop - update: m_u1Dot: " << transpose(m_u1Dot) << ", m_u1: " << transpose(m_u1) << ", m_p1: " << transpose(m_p1) << "\n"
    )
  ));

  /// Renew the proto expression
  m_renew->set_expression( nodes_expression(
    group(
      u = m_u1,
      p = m_p1,
      m_uDot = m_u1Dot,
      _cout << "yop - Renew: p: " << p << ", u: " << transpose(u) << "\n"
    )
  ));
}

////////////////////////////////////////////////////////////////////////////////
void VMS::execute() /// derived from NavierStokesSemiImplicit.cpp
{
  // typedef std::pair<Uint,Uint> BlockrowIdxT;

  m_predictor->execute();
  // std::cout << "yop: nb_iterations: " << m_nb_iterations << std::endl;

  for(Uint i = 0; i != m_nb_iterations; ++i)
  {
    m_correctorInitialiser->execute();
    //! zeroLSS not yet activated
    // This ensures that the linear system matrix is reset to zero each timestep
    m_assembly->create_component<math::LSS::ZeroLSS>("ZeroLSS")->options().set("reset_solution", true);
    // m_zeroLSS->execute();
    m_assembly->execute();
    Handle<math::LSS::System> m_lss(get_child("LSS"));
    if(i == 0) /// Apply velocity BC the first inner iteration
    {
      m_lss->rhs()->scale(dt());
      m_bc->execute();
      /// The velocity BC deals with velocity, so we need to write this in terms of acceleration
      m_lss->rhs()->scale(1./dt());
    }
    else /// Zero boundary condition after first pass
    {
      m_bc->execute();
      m_lss->dirichlet_apply(true, true);
    }
    m_solveLSS->execute();
    m_update->execute();
  }
  m_renew->execute();
}

} /// UFEM
} /// cf3