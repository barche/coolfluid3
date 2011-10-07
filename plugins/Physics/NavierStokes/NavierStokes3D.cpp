// Copyright (C) 2010-2011 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/algorithm/string/predicate.hpp>

#include "Common/CBuilder.hpp"

#include "Physics/Variables.hpp"

#include "NavierStokes3D.hpp"

namespace CF {
namespace Physics {
namespace NavierStokes {

using namespace Common;

////////////////////////////////////////////////////////////////////////////////////////////

Common::ComponentBuilder < NavierStokes::NavierStokes3D,
                           Physics::PhysModel,
                           LibNavierStokes >
                           Builder_NavierStokes3D;

////////////////////////////////////////////////////////////////////////////////////////////

const Real NavierStokes3D::Properties::gamma = 1.4;
const Real NavierStokes3D::Properties::R = 287.05;
const Real NavierStokes3D::Properties::gamma_minus_1 = NavierStokes3D::Properties::gamma-1.;

////////////////////////////////////////////////////////////////////////////////////////////

NavierStokes3D::NavierStokes3D( const std::string& name ) : Physics::PhysModel(name)
{
}

NavierStokes3D::~NavierStokes3D()
{
}

boost::shared_ptr< Physics::Variables > NavierStokes3D::create_variables( const std::string type, const std::string name )
{
  Physics::Variables::Ptr vars = boost::algorithm::contains( type, "." ) ?
        build_component_abstract_type< Physics::Variables >( type, name ) :
        build_component_abstract_type< Physics::Variables >( LibNavierStokes::library_namespace() + "." + type, name );

  add_component( vars );

  return vars;
}

////////////////////////////////////////////////////////////////////////////////////////////

} // NavierStokes
} // Physics
} // CF
