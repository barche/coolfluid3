// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "common/RegistLibrary.hpp"

#include "mesh/tecplot/LibTecplot.hpp"

namespace cf3 {
namespace mesh {
namespace tecplot {

cf3::common::RegistLibrary<LibTecplot> libtecplot;

} // tecplot
} // mesh
} // cf3
