/// Copyright (C) 2010-2011 von Karman Institute for Fluid Dynamics, Belgium
///
/// This software is distributed under the terms of the
/// GNU Lesser General Public License version 3 (LGPLv3).
/// See doc/lgpl.txt and doc/gpl.txt for the license text.

/// derived from ScalarAdvection.hpp

/// Source: Variational Multiscale (VMS) Residual-based turbulence modeling 
/// for LES of incompressible flows (Y. Bazilevs, Elsevier, (2007) 173-201)

#ifndef cf3_UFEM_VMS_hpp
#define cf3_UFEM_VMS_hpp

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/max.hpp>

#define BOOST_PROTO_MAX_ARITY 10
#ifdef BOOST_MPL_LIMIT_METAFUNCTION_ARITY
  #undef BOOST_MPL_LIMIT_METAFUNCTION_ARITY
#endif
#define BOOST_MPL_LIMIT_METAFUNCTION_ARITY 10

#include <boost/scoped_ptr.hpp>

#include "LibUFEM.hpp"
#include "LSSActionUnsteady.hpp"

#include "SUPG.hpp"

#include "ns_semi_implicit/LSSVectorOps.hpp"

namespace cf3 {

namespace UFEM {

using solver::actions::Proto::SFOp;
using solver::actions::Proto::CustomSFOp;  

/// solver for scalar transport
class UFEM_API VMS : public LSSActionUnsteady
{
public: // functions

  /// Constructor
  /// @param name of the component
  VMS ( const std::string& name );
  
  /// Modification of the NS execute function
  void execute() override;
  
  /// Get the class name
  static std::string type_name () { return "VMS"; }

private:

  /// Create the solver
  void trigger_assembly();

  /// 
  void set_expression();

  /// AlphaM, AlphaF, Gamma coefficients for the alpha-method
  Real m_alphaM, m_alphaF, m_gamma;

  /// Storage of the stabilization coefficients
  Real tau_su, tau_c, tau_m;

  /// Variables
  /// 
  FieldVariable<5, VectorField> u_dot1;
  /// The pressure solution field at time n-1
  FieldVariable<6, ScalarField> p1;
  

  SFOp< CustomSFOp<VectorLSSVector> > a;
  SFOp< CustomSFOp<ScalarLSSVector> > delta_p_sum;
  /// LSS for the velocity & pressure
  Handle<LSSAction> u_lss;
  Handle<LSSActionUnsteady> p_lss;
  // Actions that handle different stages of assembly, used by the set_elements_expressions function
  Handle<common::Component> m_velocity_assembly;
  Handle<common::Component> m_u_rhs_assembly;Handle<common::Component> m_p_rhs_assembly;
  Handle<common::Component> m_p_rhs;
  Handle<common::Component> m_pressure_assembly;
  Handle<common::Component> m_p_strategy_first;
  Handle<common::Component> m_p_strategy_second;
  Handle<common::Component> solve_p_lss;
  Handle<common::Component> m_apply_aup;

  // Time component
  Handle<solver::Time> m_time;



  // /// Time step
  // Real dt;
    
  /// LSS for the pressure
  Handle<InitialConditions> m_initial_conditions;

  // Actions that handle different stages of assembly, used by the set_elements_expressions function
  Handle<common::Component> m_mass_matrix_assembly;
  Handle<common::Component> m_inner_loop;

};


} // UFEM
} // cf3


#endif // cf3_UFEM_VMS_hpp