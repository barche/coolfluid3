// Copyright (C) 2010-2011 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Based on WALE files

#ifndef cf3_UFEM_EulerDNS_hpp
#define cf3_UFEM_EulerDNS_hpp

#include "../InitialConditions.hpp"

#include "LibUFEMLES.hpp"
#include "solver/actions/Proto/ProtoAction.hpp"
#include "solver/actions/Proto/Expression.hpp"

namespace cf3 {
namespace UFEM {
namespace les {

namespace detail
{

inline void cell_sizes(const mesh::LagrangeP1::Hexa3D::NodesT& nodes, Real& d1, Real& d2, Real& d3)
{
  d1 = (nodes.row(1) - nodes.row(0)).norm();
  d2 = (nodes.row(3) - nodes.row(0)).norm();
  d3 = (nodes.row(4) - nodes.row(0)).norm();
}

inline void cell_sizes(const mesh::LagrangeP1::Prism3D::NodesT& nodes, Real& d1, Real& d2, Real& d3)
{
  d1 = (nodes.row(1) - nodes.row(0)).norm();
  d2 = (nodes.row(2) - nodes.row(0)).norm();
  d3 = (nodes.row(3) - nodes.row(0)).norm();
}

inline void cell_sizes(const mesh::LagrangeP1::Tetra3D::NodesT& nodes, Real& d1, Real& d2, Real& d3)
{
  d1 = (nodes.row(1) - nodes.row(0)).norm();
  d2 = (nodes.row(2) - nodes.row(0)).norm();
  d3 = (nodes.row(3) - nodes.row(0)).norm();
}

inline void aspect_ratios(const Real d1, const Real d2, const Real d3, Real& a1, Real& a2)
{
  if (d1 >= d2 && d1 >= d3)
  {
    a1 = d2/d1;
    a2 = d3/d1;
  }
  else if (d2 >= d1 && d2 >= d3)
  {
    a1 = d1/d2;
    a2 = d3/d2;
  }
  else
  {
    a1 = d1/d3;
    a2 = d2/d3;
  }
}
  
/// Proto functor to compute the turbulent viscosity
struct ComputeNuDynSmagorinsky
{
  typedef void result_type;
  
  ComputeNuDynSmagorinsky() :
    cs(0.148),
    use_anisotropic_correction(false)
  {
  }

  template<typename UT, typename NUT, typename ValenceT>
  void operator()(const UT& u, NUT& nu, const ValenceT& valence, const Real nu_visc) const
  {
    typedef typename UT::EtypeT ElementT;
    static const Uint dim = ElementT::dimension;
    
    typedef mesh::Integrators::GaussMappedCoords<1, ElementT::shape> GaussT;
    typedef Eigen::Matrix<Real, dim, dim> SMatT;
        
    const SMatT grad_u = u.nabla(GaussT::instance().coords.col(0))*u.value();
    const SMatT S = 0.5*(grad_u + grad_u.transpose());
    const Real S_norm2 = S.squaredNorm();
    const Real S_norm4 = 2 * S_norm2;
    // const SMatT grad_u2 = grad_u*grad_u.transpose();

    // SMatT Sd = 0.5*(grad_u2 + grad_u2.transpose());
    // Sd.diagonal().array() -= grad_u2.trace()/3.;
    // const Real Sd_norm2 = Sd.squaredNorm();

    // Compute the anisotropic cell size adjustment using the method of Scotti et al.
    Real f = 1.;
    if(use_anisotropic_correction)
    {
      Real d1, d2, d3, a1, a2;
      cell_sizes(u.support().nodes(), d1, d2, d3);
      aspect_ratios(d1, d2, d3, a1, a2);
      const Real log_a1 = ::log(a1);
      const Real log_a2 = ::log(a2);
      f = ::cosh(::sqrt(4./27.*(log_a1*log_a1 - log_a1*log_a2 + log_a2*log_a2)));
    }
    
    // Get the isotropic length scale
    const Real delta_iso = ::pow(u.support().volume(), 2./3.);
    const Real delta2_iso = 2 * ::pow(u.support().volume(), 2./3.);

    // Find the dynamic Smagorinsky parameter (depending on space and time)
    const SMatT M = delta2_iso * ::pow(2*S_norm4, 0.5) * S - delta_iso * ::pow(2*S_norm2, 0.5) * S;
    const SMatT L = -2 * cs * M; // Germano identity
    const SMatT Num = - L * M.transpose();
    const SMatT Den = 2 * M * M.transpose();
    const SMatT Frac = Num * Den.inverse().transpose();
    const Real cs_xt = Frac.squaredNorm();

    // Compute the viscosity
    Real nu_t = cs_xt*cs_xt * f*f * delta_iso * ::pow(2*S_norm2, 0.5);
    if(nu_t < 0. || !std::isfinite(nu_t))
      nu_t = 0.;

    const Eigen::Matrix<Real, ElementT::nb_nodes, 1> nodal_vals = (nu_t + nu_visc)*valence.value().array().inverse();
    nu.add_nodal_values(nodal_vals);
  }
  
  // Model constant
  Real cs;
  bool use_anisotropic_correction;
};

}
  
/// DynSmagorinsky LES model component
class DynSmagorinsky : public solver::actions::Proto::ProtoAction
{
public: // functions
  /// Contructor
  /// @param name of the component
  DynSmagorinsky( const std::string& name );

  /// Get the class name
  static std::string type_name () { return "DynSmagorinsky"; }
  
  virtual void execute();

private:
  void trigger_set_expression();
  void trigger_initial_conditions();
  virtual void on_regions_set();
  Handle<InitialConditions> m_initial_conditions;
  Handle<common::Component> m_node_valence;
  Handle<solver::actions::Proto::ProtoAction> m_reset_viscosity;
  
  /// The data stored by the DynSmagorinsky op terminal
  solver::actions::Proto::MakeSFOp<detail::ComputeNuDynSmagorinsky>::stored_type m_dynSmagorinsky_op;
  
  /// Terminal with a reference to the DynSmagorinsky op data
  solver::actions::Proto::MakeSFOp<detail::ComputeNuDynSmagorinsky>::reference_type dynSmagorinsky;
};

} // les
} // UFEM
} // cf3


#endif // cf3_UFEM_EulerDNS_hpp
