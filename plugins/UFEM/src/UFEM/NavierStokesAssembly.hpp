// Copyright (C) 2010-2011 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NavierStokes.hpp"

#include <boost/mpl/back_inserter.hpp>
#include <boost/mpl/copy.hpp>

#include "solver/actions/Proto/ProtoAction.hpp"
#include "solver/actions/Proto/Expression.hpp"

#include "NavierStokesSpecializations.hpp"

namespace cf3 {
namespace UFEM {

using namespace common;
using namespace solver;
using namespace solver::actions;
using namespace solver::actions::Proto;

using boost::proto::lit;

template<typename GenericElementsT, typename SpecializedElementsT>
void NavierStokes::set_assembly_expression(const std::string& action_name)
{
  // Get all the relevant types as the concatenation of the generic and specialized element types:
  typedef typename boost::mpl::copy< SpecializedElementsT, boost::mpl::back_inserter<GenericElementsT> >::type AllElementsT;

  // Proto function that applies expressions only to GenericElementsT
  static const typename boost::proto::terminal< RestrictToElementTypeTag<GenericElementsT> >::type for_generic_elements = {};
  // Proto function that applies expressions only to SpecializedElementsT
  static const typename boost::proto::terminal< RestrictToElementTypeTag<SpecializedElementsT> >::type for_specialized_elements = {};

  const Real theta = options().option("theta").value<Real>();
  if(theta < 0. || theta > 1.)
    throw SetupError(FromHere(), "Value " + to_str(theta) + " for theta option of " + uri().path() + " is outside of the valid range from 0 to 1.");

  static const boost::proto::terminal<ElementSystemMatrix<boost::mpl::int_<2>>>::type _B = {}; // SUPG non-time terms
  static const boost::proto::terminal<ElementSystemMatrix<boost::mpl::int_<3>>>::type _U = {}; // SUPG time terms

  // The actual matrix assembly
  m_assembly->add_component(create_proto_action
  (
    action_name,
    elements_expression
    (
      AllElementsT(),
      group
      (
        _A = _0, _T = _0, _B = _0, _U = _0, _a = _0,
        compute_tau.apply(u_adv, nu_eff, lit(dt()), lit(tau_ps), lit(tau_su), lit(tau_bulk)),
        element_quadrature
        (
          _A(p    , u[_i]) += transpose(N(p)) * nabla(u)[_i],
          _A(u[_i], u[_i]) += nu_eff * transpose(nabla(u)) * nabla(u) + transpose(N(u)) * u_adv*nabla(u),
          _A(u[_i], u[_j]) += nu_eff * transpose(nabla(u)[_j]) * nabla(u)[_i],
          _A(u[_i], p)     += transpose(N(u)) * nabla(p)[_i],
          _T(u[_i], u[_i]) += transpose(N(u)) * N(u),
          _B(p    , u[_i]) += tau_ps * transpose(nabla(p)[_i]) * u_adv*nabla(u),
          _B(p    , p)     += tau_ps * transpose(nabla(p)) * nabla(p),
          _B(u[_i], u[_i]) += transpose(tau_su*u_adv*nabla(u)) * u_adv*nabla(u),
          _B(u[_i], p)     += transpose(tau_su*u_adv*nabla(u)) * nabla(p)[_i],
          _B(u[_i], u[_j]) += transpose(tau_bulk*nabla(u)[_i]) * nabla(u)[_j],
          _U(p    , u[_i]) += tau_ps * transpose(nabla(p)[_i]) * N(u),
          _U(u[_i], u[_i]) += transpose(tau_su*u_adv*nabla(u)) * N(u),
          _a[u[_i]]        += transpose(N(u) + tau_su*u_adv*nabla(u)) * g[_i] * density_ratio
        ),
        system_rhs += -(_A + _B) * _x,
        _A(p) = _A(p) / theta,
        _B(p) = _B(p) / theta,
        system_matrix += invdt() * (_T + _U) + theta * (_A + _B)
      )
    )
  ));
}

} // UFEM
} // cf3
