// Copyright (C) 2010 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

////////////////////////////////////////////////////////////////////////////////

#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/regex.hpp>

#include "Common/Foreach.hpp"
#include "Common/Log.hpp"
#include "Common/OptionT.hpp"

#include "Physics/VariablesDescriptor.hpp"

namespace CF {
namespace Physics {

using namespace Common;

////////////////////////////////////////////////////////////////////////////////


struct VariablesDescriptor::Implementation
{
  Implementation(Component& component) :
    m_component(component),
    m_dim(0u),
    m_dim_configured(false)
  {
    m_component.options().add_option< OptionT<Uint> >("dimensions", 0)
      ->pretty_name("Dimensions")
      ->description("Dimensionality of the problem, i.e. the number of components for the spatial coordinates")
      ->mark_basic()
      ->link_to(&m_dim)
      ->attach_trigger(boost::bind(&Implementation::trigger_dimensions, this));
      
    m_component.options().add_option< OptionT<std::string> >("field_name", m_component.name())
      ->pretty_name("Field Name")
      ->description("Name for the field that will store the variables described here")
      ->mark_basic()
      ->link_to(&m_field_name);
  }

  //////////////// Interface implementation /////////////////////
  
  void push_back(const std::string& name, const VariablesDescriptor::Dimensionalities::Type type)
  {
    // Only proceed if the variable did not exist already
    if(!m_indices.insert(std::make_pair(name, m_types.size())).second)
    {
      CFdebug << "Ignoring double registration of variable " << name << CFendl;
      return;
    }
    
    m_types.push_back(type);
    m_offsets.push_back(m_size);
    m_user_names.push_back(name);
    
    m_component.options().add_option< OptionT<std::string> >(variable_property_name(name), name)
        ->pretty_name(name + std::string(" Variable Name"))
        ->description("Variable name for variable " + name)
        ->link_to(&m_user_names.back());
    
    m_size += to_size(type);
  }
  
  Uint size() const
  {
    if(!m_dim_configured)
      throw SetupError(FromHere(), "Attempt to get total size for " + m_component.uri().string() + " before dimension is configured");
    
    return m_size;
  }
  
  Uint size(const std::string& name) const
  {
    if(!m_dim_configured)
      throw SetupError(FromHere(), "Attempt to get dimension for variable " + name + " in " + m_component.uri().string() + " before dimension is configured");
    
    return to_size(m_types[checked_index(name)]);
  }
  
  Uint offset(const std::string& name)
  {
    if(!m_dim_configured)
      throw SetupError(FromHere(), "Attempt to get offset for variable " + name + " in " + m_component.uri().string() + " before dimension is configured");
    
    return m_offsets[checked_index(name)];
  }
  
  const std::string& user_variable_name(const std::string& name) const
  {
    return m_user_names[checked_index(name)];
  }
  
  /// Implementation based on Willem Deconincks code for CField
  void set_variables(const std::string& description)
  {
    const boost::regex e_variable("([[:word:]]+)[[:space:]]*(\\[[[:space:]]*[[:word:]]+[[:space:]]*\\])?");
    const boost::regex e_scalar  ("((s(cal(ar)?)?)?)|1"     ,boost::regex::perl|boost::regex::icase);
    const boost::regex e_vector("v(ec(tor)?)?",boost::regex::perl|boost::regex::icase);
    const boost::regex e_tensor("t(ens(or)?)?",boost::regex::perl|boost::regex::icase);
    
    std::vector<std::string> tokenized_variables;
    boost::split(tokenized_variables, description, boost::is_any_of(","));
    const Uint nb_vars = tokenized_variables.size();
    
    std::vector<std::string> names_to_add; names_to_add.reserve(nb_vars);
    std::vector<Dimensionalities::Type> types_to_add; types_to_add.reserve(nb_vars);
    
    BOOST_FOREACH(std::string var, tokenized_variables)
    {
      boost::match_results<std::string::const_iterator> what;
      if (boost::regex_search(var,what,e_variable))
      { 
        names_to_add.push_back(what[1]); boost::trim(names_to_add.back());
        std::string var_type = what[2];
        boost::trim_if(var_type, boost::is_any_of(" []"));
        
        if(boost::regex_match(var_type,e_scalar))
        {
          types_to_add.push_back(Dimensionalities::SCALAR);
        }
        else if(boost::regex_match(var_type,e_vector))
        {
          types_to_add.push_back(Dimensionalities::VECTOR);
        }
        else if(boost::regex_match(var_type,e_tensor))
        {
          types_to_add.push_back(Dimensionalities::TENSOR);
        }
        else
        {
          throw ParsingFailed(FromHere(), "Type " + var_type + " deduced from " + var + " is not recognized");
        }
      }
      else
      {
        throw ParsingFailed(FromHere(), "Invalid variable: " + var);
      }
    }
    
    for(Uint i = 0; i != nb_vars; ++i)
    {
      push_back(names_to_add[i], types_to_add[i]);
    }
  }
  
  std::string description() const
  {
    if(!m_dim_configured)
      throw SetupError(FromHere(), "Attempt to get field description in " + m_component.uri().string() + " before dimension is configured");
    
    const Uint nb_vars = m_indices.size();
    // Build an ordered list of the variables
    std::vector<std::string> ordered_vars(nb_vars);
    for(IndexMapT::const_iterator it = m_indices.begin(); it != m_indices.end(); ++it)
    {
      ordered_vars[it->second] = it->first;
    }
    
    // Create a string with the description of the variables
    std::stringstream result_str;
    for(Uint i = 0; i != nb_vars; ++i)
    {
      result_str << (i == 0 ? "" : ",") << ordered_vars[i] << "[" << to_size(m_types[i]) << "]";
    }
    
    return result_str.str();
  }

  ///////////// Helper functions ////////
  
  /// Convert the dimensionality type to a real size
  Uint to_size(Dimensionalities::Type dim_type) const
  {
    // special cases
    if(dim_type == Dimensionalities::VECTOR)
      return m_dim;
    
    if(dim_type == Dimensionalities::TENSOR)
      return m_dim*m_dim;
    
    // all others can be converted directly
    int converted = static_cast<int>(dim_type);
    cf_assert(converted > 0);
    
    return static_cast<Uint>(converted);
  }
  
  void trigger_dimensions()
  {
    // Recalculate size and offsets
    m_size = 0;
    const Uint nb_vars = m_indices.size();
    for(Uint i = 0; i != nb_vars; ++i)
    {
      m_offsets[i] = m_size;
      m_size += to_size(m_types[i]);
    }
    
    m_dim_configured = true;
  }
  
  std::string variable_property_name(std::string var_name)
  {
    boost::to_lower(var_name);
    return var_name + "_variable_name";
  }
  
  /// Index for the given internal name, throwing a nice error if it's not found
  Uint checked_index(const std::string& name) const
  {
    IndexMapT::const_iterator found = m_indices.find(name);
    if(found != m_indices.end())
      return found->second;
    
    throw ValueNotFound(FromHere(), "Variable with internal name " + name + " was not found in descriptor " + m_component.uri().string());
  }

  /////////////// data //////////////

  Component& m_component;

  /// Name of the field
  std::string m_field_name;
  
  /// dimensionality of physics
  Uint m_dim;
  /// True if the dimensions have been configured
  bool m_dim_configured;

  /// Total size
  Uint m_size;
  
  /// Mapping from variable internal name to index in the vectors
  typedef std::map<std::string, Uint> IndexMapT;
  IndexMapT m_indices;
  
  /// Type of each variable
  typedef std::vector<Dimensionalities::Type> VarTypesT;
  VarTypesT m_types;
  
  /// Offsets for the variables
  std::vector<Uint> m_offsets;
  
  /// User defined variable names
  std::vector<std::string> m_user_names;
  
};

////////////////////////////////////////////////////////////////////////////////

VariablesDescriptor::Dimensionalities::Convert::Convert()
{
  all_fwd = boost::assign::map_list_of
      ( VariablesDescriptor::Dimensionalities::SCALAR,    "scalar" )
      ( VariablesDescriptor::Dimensionalities::VECTOR,    "vector" )
      ( VariablesDescriptor::Dimensionalities::TENSOR,    "tensor" );

  all_rev = boost::assign::map_list_of
      ("scalar", VariablesDescriptor::Dimensionalities::SCALAR)
      ("vector", VariablesDescriptor::Dimensionalities::VECTOR)
      ("tensor", VariablesDescriptor::Dimensionalities::TENSOR);
}

////////////////////////////////////////////////////////////////////////////////

VariablesDescriptor::VariablesDescriptor( const std::string& name ) :
  Component(name),
  m_implementation(new Implementation(*this))
{
}

VariablesDescriptor::~VariablesDescriptor()
{
}

////////////////////////////////////////////////////////////////////////////////

void VariablesDescriptor::push_back(const std::string& name, const VariablesDescriptor::Dimensionalities::Type type)
{
  m_implementation->push_back(name, type);
}

////////////////////////////////////////////////////////////////////////////////

Uint VariablesDescriptor::size() const
{
  return m_implementation->size();
}

////////////////////////////////////////////////////////////////////////////////

Uint VariablesDescriptor::size(const std::string& name)
{
  return m_implementation->size(name);
}

////////////////////////////////////////////////////////////////////////////////

Uint VariablesDescriptor::offset(const std::string& name)
{
  return m_implementation->offset(name);
}

////////////////////////////////////////////////////////////////////////////////

const std::string& VariablesDescriptor::field_name() const
{
  return m_implementation->m_field_name;
}

////////////////////////////////////////////////////////////////////////////////

void VariablesDescriptor::set_field_name(const std::string& name)
{
  configure_option("field_name", name);
}

////////////////////////////////////////////////////////////////////////////////

const std::string& VariablesDescriptor::user_variable_name(const std::string& name) const
{
  return m_implementation->user_variable_name(name);
}

////////////////////////////////////////////////////////////////////////////////

std::string VariablesDescriptor::description() const
{
  return m_implementation->description();
}

////////////////////////////////////////////////////////////////////////////////

void VariablesDescriptor::set_variables(const std::string& description)
{
  m_implementation->set_variables(description);
}

////////////////////////////////////////////////////////////////////////////////

} // Physics
} // CF
