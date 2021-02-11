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

using namespace solver::actions::Proto;

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

/// Get neighbouring elements and center element's data
template<typename UT>
inline void get_neighbElts(std::set<Uint>& neighbElts, const UT& u, const mesh::NodeConnectivity* node_connectivity)
{
  const auto element_nodes = u.support().element_connectivity();
  // UT neighbElt_u = u;

  // CFdebug << "xxxxxxxxxxxxxxxxxxxxxxxxxx" << CFendl;
  // CFdebug << "looking up adjacent elements for element " << u.support().element_idx() << CFendl;
  for(Uint node_idx : element_nodes)
  {
    // CFdebug << "Elements next to node " << node_idx << ": ";
    const Uint elements_begin = node_connectivity->node_first_elements()[node_idx];
    const Uint elements_end = elements_begin + node_connectivity->node_element_counts()[node_idx];
    for(Uint j = elements_begin; j != elements_end; ++j)
    {
      // CFdebug << " (" << node_connectivity->node_elements()[j].first << "," << node_connectivity->node_elements()[j].second <<")";
      // neighbElt_u.set_element(node_connectivity->node_elements()[j].second); // only one type of elt considered
      neighbElts.insert(node_connectivity->node_elements()[j].second);
      
      /// Reading rank of mpi node
      // int rank;
      // MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
      // std::cout << "mpi rank from process #" << rank << std::endl;
    }
    // CFdebug << "" << CFendl;
  }
  
  /// Check the number of neighbouring elts:
  // CFdebug << "" << CFendl;
  // CFdebug << "Confirmation: Elements are: ";
  // for (auto i = neighbElts.begin(); i != neighbElts.end(); i++)
  //   CFdebug << *i << " ";
  // CFdebug << CFendl;
}
  
// Helper to get the transpose of either a vector or a scalar
template<typename T>
inline Eigen::Transpose<T const> transpose(const T& mat)
{
  return mat.transpose();
}

// Function to remove temporarly the const on a referenced value
template<class T, class U=typename std::remove_const<typename std::remove_reference<T>::type>::type> 
inline U& removeConstRef(T& support)
{
  return const_cast<U &> (support);
}

// Function to print out the nodes/elts coordinates and their respective u values
template<typename UT>
inline void check_nodeVelocity(const UT& u, bool uVal=true)
{
  const auto element_nodes = u.support().element_connectivity();
  CFdebug << "nodes id:";
  for(Uint node_idx : element_nodes)
  {
    CFdebug << " " << node_idx;
  }
  CFdebug << CFendl;
  CFdebug << "nodes coords = " << u.support().nodes() << CFendl;
  if(uVal)
  {
    CFdebug << "the values of u should be equal to the nodes coords." << CFendl;
    CFdebug << "u = " << u.value() << CFendl;
  }
}


/// Proto functor to compute the turbulent viscosity
struct ComputeNuSmagorinsky
{
  typedef void result_type;
  
  ComputeNuSmagorinsky() :
    use_anisotropic_correction(false),
    use_dynamic_smagorinsky(false),
    use_mason_wallDamping(false),
    m_cs(0.148),
    m_masonCoef(3)
  {
  }

  template<typename UT, typename NUT, typename ValenceT, typename CsT, typename SgsT>
  void operator()(UT& u, NUT& nu, const ValenceT& valence, const Real nu_visc, Real& cs, CsT& cs_elts, SgsT& sgsLambda_elts) const
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

    // List all neighbouring elements (center elt incl.)
    std::set<Uint> neighbElts {};
    get_neighbElts(neighbElts, u, m_node_connectivity);

    // CFdebug << "xxxxxxxxxxxxxxxxxxxxxxxxxx" << CFendl;
    // CFdebug << "xxxxxxxxx center Elt:" << u.support().element_idx() << CFendl;
    /// Check that the node velocity equals its coordinates
    // check_nodeVelocity(u);

    // Compute the dyn Smag parameters
    typedef typename UT::EtypeT ElementT;
    static const Uint dim = ElementT::dimension;
    
    typedef mesh::Integrators::GaussMappedCoords<1, ElementT::shape> GaussT;
    typedef Eigen::Matrix<Real, dim, dim> SMatT;
    typedef Eigen::Matrix<Real, 1, dim> SVecT;

    u.compute_values(GaussT::instance().coords.col(0));  // Compute u.eval()

    // Initialisation for static smagorinsky:
    const Real vol_gF = u.support().volume(); // "g"rid "F"ilter = volume of center elt
    const Real delta2_gF = ::pow(vol_gF, 2./3.); // square of isotropic length scale for grid filter
    SMatT grad_u = u.nabla(GaussT::instance().coords.col(0))*u.value();
    SMatT S = 0.5*(grad_u + grad_u.transpose());
    Real S_norm = S.squaredNorm(); // SquaredNorm is the sum of the square of all the matrix entries (Frobenius norm).
  
    if(use_dynamic_smagorinsky) {
      // Initialisation for dynamic smagorinsky (TEST filter > GRID filter):
      Real smagCoefC = 0; // "C" cte equal to cs^2 ; visible in Germano paper eq(5)
      Real vol_tF = 0.; // "t"est "F"ilter = sum of volume on all neighbouring elts + center elt (e.g. 27 for a unit hexa3D in the center of a 3x3x3 cubus)
      Uint storedElt_idx = u.support().element_idx(); // store the id of the centered element
      SVecT u_gtF {}; // grid + test filtered velocity
      u_gtF.setZero();
      SMatT S_gF {}; // grid filtered strain-tensor
      S_gF.setZero();
      SMatT uu_gtF {}; // grid + test filtered velocities product
      uu_gtF.setZero();
      Real S_gF_norm = 0.; // grid filtered strain-tensor magnitude
      Real SS_norm = 0.; // strain tensor product magnitude
      SMatT SSxS {}; // grid + test filtered strain tensor product magnitude times grid filtered strain tensor
      SSxS.setZero();

      // CFdebug << "neighbElts: " << neighbElts.size() << CFendl;
      for (const auto i : neighbElts)
      {
        // set the neighbouring element's data to object u:
        u.set_element(i);
        removeConstRef(u.support()).set_element(i); // needed to switch to other elt

        // Check the that the node velocity equals its coordinates
        // CFdebug << "xxxxxxxxx neighb Elt " << i << ":" << CFendl;
        // check_nodeVelocity(u,false);

        // u.value() gives velocity components for all nodes, u.eval() computed for the elt
        u.compute_values(GaussT::instance().coords.col(0)); // Compute u.eval()

        vol_tF += u.support().volume(); // data[1]
        u_gtF += u.eval() * u.support().volume(); // data[2-4]
        
        grad_u = u.nabla(GaussT::instance().coords.col(0))*u.value();
        S = 0.5*(grad_u + grad_u.transpose());

        S_gF += S * u.support().volume(); // data[5-13]
        uu_gtF += u.eval().transpose() * u.eval() * u.support().volume(); // data[15-23]
        SS_norm = std::sqrt(2 * S.squaredNorm()); // SSij
        SSxS += SS_norm * S * u.support().volume(); // data[24-32]
      }

      // Restore and re-initialise the center element's data to object u and recompute its values:
      u.set_element(storedElt_idx);
      removeConstRef(u.support()).set_element(storedElt_idx);
      u.compute_values(GaussT::instance().coords.col(0));
      grad_u = u.nabla(GaussT::instance().coords.col(0))*u.value();
      S = 0.5*(grad_u + grad_u.transpose());
      S_norm = S.squaredNorm();

      // Normalize the neighbElts stencil sums
      const Real invVol_tF = 1. / vol_tF;
      u_gtF *= invVol_tF;
      S_gF *= invVol_tF;
      S_gF_norm = std::sqrt(2 * S_gF.squaredNorm()); //data[14]
      uu_gtF *= invVol_tF;
      SSxS *= invVol_tF;

      // Find the dynamic Smagorinsky parameter (depending on space and time)
      const Real delta2_tF = ::pow(vol_tF, 2./3.); // square of isotropic length scale for TEST filter  
      const SMatT M = delta2_gF * SSxS - delta2_tF * S_gF_norm * S_gF;
      const SMatT L = uu_gtF - u_gtF.transpose() * u_gtF;
      smagCoefC = L.cwiseProduct(M).sum() / (2. * M.cwiseProduct(M).sum());

      // CFdebug << "xxxxxxxxxxxxxxxxxxxxxxxxxx" << CFendl;
      // CFdebug << "vol_tF = " << vol_tF << CFendl;
      // CFdebug << "u.value() = " << u.value() << CFendl;
      // CFdebug << "vol_gF = " << vol_gF << CFendl;
      // CFdebug << "grad_u = " << grad_u << CFendl;
      // CFdebug << "S = " << S << CFendl;
      // CFdebug << "S_norm = " << S_norm << CFendl;
      // CFdebug << "" << CFendl;

      // CFdebug << "delta2_gF = " << delta2_gF << CFendl;
      // CFdebug << "SSxS = " << SSxS << CFendl;
      // CFdebug << "delta2_tF = " << delta2_tF << CFendl;
      // CFdebug << "S_gF_norm = " << S_gF_norm << CFendl;
      // CFdebug << "S_gF = " << S_gF << CFendl;
      // CFdebug << "M = " << M << CFendl;
      // CFdebug << "" << CFendl;

      // CFdebug << "u_gtF = " << u_gtF << CFendl;
      // CFdebug << "uu_gtF = " << uu_gtF << CFendl;
      // CFdebug << "L = " << L << CFendl;
      // CFdebug << "" << CFendl;

      // Limiter (cf. Tamas les_dynamicsmagorinsky.c #377-378)
      smagCoefC = std::max(smagCoefC,0.); // limiter to avoid NaN
      cs = std::sqrt(smagCoefC); // Come back to the expression of cs:
      cs = std::min(cs, 0.23); // limiter cf. Lilly and Germano's paper
    }

    // CFdebug << "cs = " << cs << CFendl;

    Real sgsLambda = cs * std::sqrt(delta2_gF); // Subgrid Scale length scale

    if(use_mason_wallDamping) {
      // CFdebug << "sgsLambda-orig: " << sgsLambda << CFendl;
      Real zcoord = u.support().coordinates(GaussT::instance().coords.col(0))[1]; // In fact y-coordinate (...[1]) but it is seen as the z-coordinate for us...
      sgsLambda = std::pow(std::pow(sgsLambda, (-m_masonCoef)) + std::pow(m_kappa*(zcoord + m_z0), (-m_masonCoef)), (-1./m_masonCoef));
      // CFdebug << "sgsLambda-Mason: " << sgsLambda << CFendl;
    }

    // Compute the viscosity
    // Real nu_t = cs*cs*delta2_gF *f*f* std::sqrt(2*S_norm);
    Real nu_t = std::pow(sgsLambda, 2) * f*f * std::sqrt(2*S_norm);

    if(nu_t < 0. || !std::isfinite(nu_t))
      nu_t = 0.;

    // CFdebug << "nu_t (smag) = " << nu_t << ", nu_visc = " << nu_visc << CFendl;
    const Eigen::Matrix<Real, ElementT::nb_nodes, 1> nodal_vals = (nu_t + nu_visc)*valence.value().array().inverse();
    const Eigen::Matrix<Real, ElementT::nb_nodes, 1> cs_vals = cs*valence.value().array().inverse();
    const Eigen::Matrix<Real, ElementT::nb_nodes, 1> sgsLambda_vals = sgsLambda*valence.value().array().inverse();
    nu.add_nodal_values(nodal_vals);
    cs_elts.add_nodal_values(cs_vals);
    sgsLambda_elts.add_nodal_values(sgsLambda_vals);
  }
  
  // Model constant
  Real m_cs;
  int m_masonCoef; // Mason coefficient. Set to 3 according to Goit.
  bool use_anisotropic_correction;
  bool use_dynamic_smagorinsky;
  bool use_mason_wallDamping;
  mesh::NodeConnectivity* m_node_connectivity = nullptr;

  /// Default from KEpsilonPhysics.cpp. Else, from the model
  Real m_kappa = 0;
  Real m_z0 = 0;
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
  // Handle<solver::actions::Proto::ProtoAction> m_set_node_velocity;
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
