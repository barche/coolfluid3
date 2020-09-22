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
#include "mesh/ConnectivityData.hpp"
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

template<typename UT>
inline std::set<Uint> get_neighbourElts(std::set<Uint>& neighbourElts, const UT& u, const mesh::NodeConnectivity* m_node_connectivity)
{
  const auto element_nodes = u.support().element_connectivity();
  // UT neighbour_u = u;

  // CFdebug << "looking up adjacent elements for element " << u.support().element_idx() << CFendl;
  for(Uint node_idx : element_nodes)
  {
    // CFdebug << "Elements next to node " << node_idx << ": ";
    const Uint elements_begin = m_node_connectivity->node_first_elements()[node_idx];
    const Uint elements_end = elements_begin + m_node_connectivity->node_element_counts()[node_idx];
    for(Uint j = elements_begin; j != elements_end; ++j)
    {
      // CFdebug << " (" << m_node_connectivity->node_elements()[j].first << "," << m_node_connectivity->node_elements()[j].second <<")";
      // neighbour_u.set_element(m_node_connectivity->node_elements()[j].second); // only one type of elt considered
      neighbourElts.insert(m_node_connectivity->node_elements()[j].second);
    }
    // CFdebug << CFendl;
  }
  // CFdebug << CFendl;
  
  /// Check the number of neighbouring elts:
  // CFdebug << "Confirmation: Elements are: ";
  // for (auto i = neighbourElts.begin(); i != neighbourElts.end(); i++)
  //   CFdebug << *i << " ";
  // CFdebug << CFendl;

  return neighbourElts;
}
  
// Helper to get the transpose of either a vector or a scalar
template<typename T>
inline Eigen::Transpose<T const> transpose(const T& mat)
{
  return mat.transpose();
}


/// Proto functor to compute the turbulent viscosity
struct ComputeNuSmagorinsky
{
  typedef void result_type;
  
  ComputeNuSmagorinsky() :
    cs(0.148),
    use_anisotropic_correction(false),
    use_dynamic_smagorinsky(false)
  {
  }

  template<typename UT, typename NUT, typename ValenceT>
  void operator()(const UT& u, NUT& nu, const ValenceT& valence, const Real nu_visc, Real& cs) const
  {
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

    // List all neighbouring elements
    std::set<Uint> neighbourElts {};
    get_neighbourElts(neighbourElts, u, m_node_connectivity);

    // Compute the dyn Smag parameters
    typedef typename UT::EtypeT ElementT;
    typedef typename UT::ValueT UValT;
    typedef typename UT::EvalT UEvalT;
    static const Uint dim = ElementT::dimension;
    
    typedef mesh::Integrators::GaussMappedCoords<1, ElementT::shape> GaussT;
    typedef Eigen::Matrix<Real, dim, dim> SMatT;

    const SMatT grad_u = u.nabla(GaussT::instance().coords.col(0))*u.value();
    const SMatT S = 0.5*(grad_u + grad_u.transpose());
    const Real S_norm = S.squaredNorm(); // SquaredNorm is the square of the square root of the sum of the abs squared elts ; in other words, the sum of the abs squared elts.

    // Get the isotropic length scale
    const Real delta_iso = ::pow(u.support().volume(), 2./3.);

    UT neighbour_u = u;
    Real neighbourElts_volume {0.};

    for (auto i = neighbourElts.begin(); i != neighbourElts.end(); i++)
    {
      neighbour_u.set_element(*i);
      neighbourElts_volume += neighbour_u.support().volume();

      if(use_dynamic_smagorinsky) {
        // test filter suffix: X2 (the grid filter has been multiplied by 2)
        const int tf = 2;
        const SMatT S2 = tf * S;
        const Real S2_norm = S2.squaredNorm();
        const Real delta2_iso = ::pow(tf * u.support().volume(), 2./3.);

        // to get a computed value for u.eval():
        u.compute_values(GaussT::instance().coords.col(0));

        // Find the dynamic Smagorinsky parameter (depending on space and time)
        const SMatT M = delta2_iso * ::pow(2*S2_norm, 0.5) * S2 - delta_iso * ::pow(2*S_norm, 0.5) * S;
        // const UValT u2 = u.value() * tf * u.support().volume();
        // const SMatT uu2 = u.value().transpose() * u.value() * tf * u.support().volume(); 
        // const SMatT L = uu2 - u2.transpose() * u2;
        const UEvalT u2 = u.eval() * (tf * u.support().volume() / 26. * tf * u.support().volume()); // 26 elements surrounding the hexa3D elt
        const SMatT uu2 = detail::transpose(u.eval()) * u.eval() * (tf * u.support().volume() / 26. * tf * u.support().volume()); 
        const SMatT L = uu2 - detail::transpose(u2) * u2;
        // const SMatT L = M;

        const SMatT Num = L * M.transpose();
        const SMatT Den = 2 * M * M.transpose();
        const SMatT Frac = - Num * Den.inverse().transpose();
        cs = Frac.squaredNorm();

        // Limiter (cf. Tamas les_dynamicsmagorinsky.c #377-378)
        cs = std::max(cs,0.);
        cs = std::min(std::sqrt(cs), 0.23);

        // CFdebug << "yop: ************************* " << CFendl;
        // CFdebug << "yop: u = " << u.value() << CFendl;
        // CFdebug << "yop: GaussT::instance().coords = " << GaussT::instance().coords << CFendl;
        // CFdebug << "yop: GaussT::instance().coords.col(0) = " << GaussT::instance().coords.col(0) << CFendl;
        // CFdebug << "yop: u.nabla(GaussT::instance().coords.col(0)) = " << u.nabla(GaussT::instance().coords.col(0)) << CFendl;
        // CFdebug << "yop: grad_u = " << grad_u << CFendl;
        
        // CFdebug << "yop: u[0] = " << u.element_vector()[0] << CFendl;
        // CFdebug << "yop: u.value = " << u.value() << CFendl;
        // CFdebug << "yop: u.eval = " << u.eval() << CFendl;
        // CFdebug << "yop: u.nodes = " << u.support().nodes() << CFendl;
        // CFdebug << "yop: vol = " << u.support().volume() << CFendl;
        // CFdebug << "yop: u * tf * vol = " << u2 << CFendl;
        // CFdebug << "yop: u' * u * tf * vol = " << uu2 << CFendl;
        // CFdebug << "yop: u'*u*tf*vol - u*tf*vol * u*tf*vol = " << L << CFendl;
        // CFdebug << "yop: M = " << M << CFendl;
        // CFdebug << "yop: S = " << S << CFendl;
        // CFdebug << "yop: S2 = " << S2 << CFendl;
        // CFdebug << "yop: S_norm = " << S_norm << CFendl;
        // CFdebug << "yop: S2_norm = " << S2_norm << CFendl;
        // CFdebug << "yop: delta = " << delta_iso << CFendl;
        // CFdebug << "yop: delta2 = " << delta2_iso << CFendl;
      }
    }

    // Compute the viscosity
    Real nu_t = cs*cs*f*f*delta_iso * ::pow(2*S_norm, 0.5);
    CFdebug << "yop: cs = " << cs << CFendl;

    if(nu_t < 0. || !std::isfinite(nu_t))
      nu_t = 0.;

//    Morpheus version, giving the same result
//    const Eigen::Matrix<Real, dim, 1> du = grad_u.col(0);
//    const Eigen::Matrix<Real, dim, 1> dv = grad_u.col(1);
//    const Eigen::Matrix<Real, dim, 1> dw = grad_u.col(2);

//    const Real S2= (du[0]*du[0]+dv[1]*dv[1]+dw[2]*dw[2])
//          +1./2.*( (du[1]+dv[0])*(du[1]+dv[0])
//                  +(du[2]+dw[0])*(du[2]+dw[0])
//                  +(dv[2]+dw[1])*(dv[2]+dw[1]));

//      Real Sd2=0.;
//      Real Sdtmp = du[2]*dw[0]+dv[2]*dw[1]+2.*dw[2]*dw[2]-du[0]*du[0]-2.*du[1]*dv[0]-dv[1]*dv[1];
//      Sd2+= 1./9.*Sdtmp*Sdtmp;
//      Sdtmp = du[1]*dv[0]+2.*dv[1]*dv[1]+dv[2]*dw[1]-du[0]*du[0]-2.*du[2]*dw[0]-dw[2]*dw[2];
//      Sd2+= 1./9.*Sdtmp*Sdtmp;
//      Sdtmp = 2.*du[0]*du[0]+du[1]*dv[0]+du[2]*dw[0]-dv[1]*dv[1]-2.*dv[2]*dw[1]-dw[2]*dw[2];
//      Sd2+= 1./9.*Sdtmp*Sdtmp;
//      Sdtmp = du[1]*dw[0]+du[2]*dv[0]+dv[1]*dv[2]+dv[1]*dw[1]+dv[2]*dw[2]+dw[1]*dw[2];
//      Sd2+= 1./2.*Sdtmp*Sdtmp;
//      Sdtmp = du[0]*du[2]+du[0]*dw[0]+du[1]*dv[2]+du[2]*dw[2]+dv[0]*dw[1]+dw[0]*dw[2];
//      Sd2+= 1./2.*Sdtmp*Sdtmp;
//      Sdtmp = du[0]*du[1]+du[0]*dv[0]+du[1]*dv[1]+du[2]*dw[1]+dv[0]*dv[1]+dv[2]*dw[0];
//      Sd2+= 1./2.*Sdtmp*Sdtmp;

//      if(common::PE::Comm::instance().rank() == 0)
//        std::cout << "Tamas S2: " << S2 << ", Sd2: " << Sd2 << ", Bart S2: " << S_norm2 << ", Sd2: " << Sd_norm2 << std::endl;

//      const Real Ls = ::pow(u.support().volume(), 1./3.)*cw;
//      const Real Sd2_32=pow(Sd2,3./2.);
//      const Real Sd2_54_plus_S2_52=pow(Sd2,5./4.)+pow(S2,5./2.);
      //const Real nu_t = (Sd2_54_plus_S2_52!=0.) ? Ls*Ls*Sd2_32/Sd2_54_plus_S2_52 : 0.;
    
    // CFdebug << "yop: nu_t (smag) = " << nu_t << ", nu_visc = " << nu_visc << CFendl;
    const Eigen::Matrix<Real, ElementT::nb_nodes, 1> nodal_vals = (nu_t + nu_visc)*valence.value().array().inverse();
    nu.add_nodal_values(nodal_vals);
  }
  
  // Model constant
  Real cs;
  bool use_anisotropic_correction;
  bool use_dynamic_smagorinsky;
  mesh::NodeConnectivity* m_node_connectivity = nullptr;
};

}
  
/// Smagorinsky LES model component
class Smagorinsky : public solver::actions::Proto::ProtoAction
{
public: // functions
  /// Contructor
  /// @param name of the component
  Smagorinsky( const std::string& name );

  /// Get the class name
  static std::string type_name () { return "Smagorinsky"; }
  
  virtual void execute();

private:
  void trigger_set_expression();
  void trigger_initial_conditions();
  virtual void on_regions_set();
  Handle<InitialConditions> m_initial_conditions;
  Handle<common::Component> m_node_valence;
  Handle<solver::actions::Proto::ProtoAction> m_reset_viscosity;
  Handle<mesh::NodeConnectivity> m_node_connectivity;
  
  /// The data stored by the Smagorinsky op terminal
  solver::actions::Proto::MakeSFOp<detail::ComputeNuSmagorinsky>::stored_type m_smagorinsky_op;
  
  /// Terminal with a reference to the Smagorinsky op data
  solver::actions::Proto::MakeSFOp<detail::ComputeNuSmagorinsky>::reference_type smagorinsky;
};

} // les
} // UFEM
} // cf3


#endif // cf3_UFEM_EulerDNS_hpp
