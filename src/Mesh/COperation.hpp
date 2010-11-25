// Copyright (C) 2010 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef CF_Mesh_COperation_hpp
#define CF_Mesh_COperation_hpp

#include "Common/CBuilder.hpp"
#include "Common/OptionT.hpp"
#include "Common/ComponentPredicates.hpp"

#include "Mesh/CRegion.hpp"
#include "Mesh/CTable.hpp"
#include "Mesh/CField.hpp"
#include "Mesh/CFieldElements.hpp"
#include "Mesh/ElementData.hpp"
#include "Mesh/CElements.hpp"
#include "Mesh/ConnectivityData.hpp"

/////////////////////////////////////////////////////////////////////////////////////

namespace CF {
namespace Mesh {

///////////////////////////////////////////////////////////////////////////////////////

class Mesh_API COperation : public Common::Component
{
public: // typedefs

  /// pointers
  typedef boost::shared_ptr<COperation> Ptr;
  typedef boost::shared_ptr<COperation const> ConstPtr;

public: // functions
  /// Contructor
  /// @param name of the component
  COperation ( const std::string& name );

  /// Virtual destructor
  virtual ~COperation() {};

  /// Get the class name
  static std::string type_name () { return "COperation"; }

  virtual void set_loophelper (CElements& geometry_elements );

  virtual void set_loophelper (CTable<Real>& coordinates );

  virtual void execute (Uint index = 0 )
  {
    throw Common::NotImplemented(FromHere(), "Must create child that overloads this function");
  }

  /// Templated version for high efficiency
  template < typename EType >
  void executeT ( Uint index=0 )
  {
    throw Common::NotImplemented(FromHere(), "Must create child that overloads this function");
  }

  virtual COperation& operation();

  /// Only for use in non-templatized version
  COperation& create_operation(const std::string operation_type);


  ///@todo this function is temporary until statically linked childs are available.
  virtual void bind()
  {
    add_component(operation().get());
  }

private: // data

  Uint m_counter;

};

///////////////////////////////////////////////////////////////////////////////////////

template < typename OP1, typename OP2 >
class COperationMergeT : public COperation
{
public: // typedefs

  typedef boost::shared_ptr<COperationMergeT> Ptr;
  typedef boost::shared_ptr<COperationMergeT const> ConstPtr;

public: // functions
  /// Contructor
  /// @param name of the component
  COperationMergeT ( const std::string& name ) :
    COperation(name),
    m_op1(Common::allocate_component_type<OP1>("operation_1") ),
    m_op2(Common::allocate_component_type<OP2>("operation_2") )
  {
    BuildComponent<none>().build(this);
  }

  /// Virtual destructor
  virtual ~COperationMergeT() {};

  /// Get the class name
  static std::string type_name () { return "COperationMergeT"; }

  void set_loophelper (CElements& geometry_elements )
  {
    m_op1->set_loophelper(geometry_elements);
    m_op2->set_loophelper(geometry_elements);
  }

  const OP1& operation1() const
  {
    return *m_op1;
  }

  OP1& operation1()
  {
    return *m_op1;
  }


  const OP2& operation2() const
  {
    return *m_op2;
  }

  OP2& operation2()
  {
    return *m_op2;
  }

  template < typename EType >
  void executeT (  Uint elem )
  {
    m_op1->template executeT<EType>( elem );
    m_op2->template executeT<EType>( elem );
  }

  void execute (  Uint elem )
  {
    m_op1->execute( elem );
    m_op2->execute( elem );
  }


  ///@todo this function is temporary until statically linked childs are available.
  virtual void bind()
  {
    add_component(operation1().get());
    add_component(operation2().get());
  }

private: // data

  typename OP1::Ptr m_op1;
  typename OP2::Ptr m_op2;
};

typedef COperationMergeT<COperation,COperation> COperationMerge;

///////////////////////////////////////////////////////////////////////////////////////

class Mesh_API COutputField : public COperation
{
public: // typedefs

  typedef boost::shared_ptr<COutputField> Ptr;
  typedef boost::shared_ptr<COutputField const> ConstPtr;

public: // functions
  /// Contructor
  /// @param name of the component
  COutputField ( const std::string& name ) : COperation(name)
  {
    BuildComponent<nosignals>().build(this);
    m_properties["Field"].as_option().attach_trigger ( boost::bind ( &COutputField::trigger_Field,   this ) );
  }

  void trigger_Field()
  {
    Common::CPath field_path (property("Field").value<Common::URI>());
    CFdebug << "field_path = " << field_path.string() << CFendl;
    scalar_field = look_component_type<CField>(field_path);
    scalar_name = scalar_field->field_name();
  }

  /// Virtual destructor
  virtual ~COutputField() {};

  /// Get the class name
  static std::string type_name () { return "COutputField"; }

  /// Configuration Options
  virtual void define_config_properties ()
  {
    m_properties.add_option< Common::OptionT<Common::URI> > ("Field","Field URI to output", Common::URI("cpath://"))->mark_basic();
  }

  void set_loophelper (CElements& geometry_elements )
  {
    data = boost::shared_ptr<LoopHelper> ( new LoopHelper(*scalar_field, geometry_elements ) );
    CFinfo << data->field_elements.full_path().string() << CFendl;
  }

  template < typename SFType >
  void executeT ( Uint elem )
  {
    execute(elem);
  }

  void execute ( Uint elem )
  {
    CFinfo << "   " << scalar_name << "["<<elem<<"] = " << data->scalars[elem][0] << CFendl;
  }

private: // data

  struct LoopHelper
  {
    LoopHelper(CField& scalar_field, CElements& geometry_elements) :
      field_elements(geometry_elements.get_field_elements(scalar_field.field_name())),
      scalars(field_elements.data())
    { }
    CFieldElements& field_elements;
    CTable<Real>& scalars;
  };

  boost::shared_ptr<LoopHelper> data;

  CField::Ptr scalar_field;
  std::string scalar_name;
};

/////////////////////////////////////////////////////////////////////////////////////

class Mesh_API CComputeVolumes : public COperation
{
public: // typedefs

  typedef boost::shared_ptr<CComputeVolumes> Ptr;
  typedef boost::shared_ptr<CComputeVolumes const> ConstPtr;

public: // functions
  /// Contructor
  /// @param name of the component
  CComputeVolumes ( const std::string& name ) : COperation(name)
  {
    BuildComponent<nosignals>().build(this);
    m_properties["Field"].as_option().attach_trigger ( boost::bind ( &CComputeVolumes::trigger_Field,   this ) );
  }

  void trigger_Field()
  {
    Common::CPath field_path (property("Field").value<Common::URI>());
    CFdebug << "field_path = " << field_path.string() << CFendl;
    volume_field = look_component_type<CField>(field_path);
  }

  /// Virtual destructor
  virtual ~CComputeVolumes() {};

  /// Get the class name
  static std::string type_name () { return "CComputeVolumes"; }

  /// Configuration Options
  virtual void define_config_properties ()
  {
    m_properties.add_option< Common::OptionT<Common::URI> > ("Field","Field URI to output", Common::URI("cpath://"))->mark_basic();
  }

  void set_loophelper (CElements& geometry_elements )
  {
    data = boost::shared_ptr<LoopHelper> ( new LoopHelper(*volume_field, geometry_elements ) );
  }

  template < typename SFType >
  void executeT ( Uint elem )
  {
    typename SFType::NodeMatrixT nodes;
    fill( nodes, data->coordinates, data->connectivity_table[elem] );
    data->volumes[elem][0] = SFType::volume( nodes );
  }

  void execute ( Uint elem )
  {
    ElementType::NodesT nodes(data->connectivity_table.row_size(), data->coordinates.row_size());
    fill(nodes, data->coordinates, data->connectivity_table[elem]);
    data->volumes[elem][0] = data->field_elements.element_type().computeVolume( nodes );
  }

private: // data

  struct LoopHelper
  {
    LoopHelper(CField& volume_field, CElements& geometry_elements) :
      field_elements(geometry_elements.get_field_elements(volume_field.field_name())),
      volumes(field_elements.data()),
      coordinates(field_elements.coordinates()),
      connectivity_table(field_elements.connectivity_table())
    { }
    CFieldElements& field_elements;
    CTable<Real>& volumes;
    CTable<Real>& coordinates;
    CTable<Uint>& connectivity_table;
  };

  boost::shared_ptr<LoopHelper> data;

  CField::Ptr volume_field;
};

/////////////////////////////////////////////////////////////////////////////////////

class Mesh_API CSetValue : public COperation
{
public: // typedefs

  typedef boost::shared_ptr<CSetValue> Ptr;
  typedef boost::shared_ptr<CSetValue const> ConstPtr;

public: // functions
  /// Contructor
  /// @param name of the component
  CSetValue ( const std::string& name ) : COperation(name)
  {
    BuildComponent<nosignals>().build(this);
    m_properties["Field"].as_option().attach_trigger ( boost::bind ( &CSetValue::trigger_Field,   this ) );
  }

  void trigger_Field()
  {
    Common::CPath field_path (property("Field").value<Common::URI>());
    CFdebug << "field_path = " << field_path.string() << CFendl;
    field = look_component_type<CField>(field_path);
  }

  /// Virtual destructor
  virtual ~CSetValue() {};

  /// Get the class name
  static std::string type_name () { return "CSetValue"; }

  /// Configuration Options
  virtual void define_config_properties ()
  {
    m_properties.add_option< Common::OptionT<Common::URI> > ("Field","Field URI to output", Common::URI("cpath://"))->mark_basic();
  }

  virtual void set_loophelper (CTable<Real>& coordinates )
  {
    data = boost::shared_ptr<LoopHelper> ( new LoopHelper(*field, coordinates ) );
  }

  template < typename SFType >
  void executeT ( Uint node )
  {
    execute(node);
  }

  void execute ( Uint node )
  {

    CTable<Real>& field_data = data->field_data;
    field_data[node][0] = data->coordinates[node][XX]*data->coordinates[node][XX];
  }

private: // data

  struct LoopHelper
  {
    LoopHelper(CField& field, CTable<Real>& coords) :
      coordinates(coords),
      node_connectivity(*coords.look_component_type<CNodeConnectivity>("../node_connectivity")),
      local_field(coords.get_parent()->get_type<CRegion>()->get_field(field.name())),
      field_data(Common::get_tagged_component_typed<CTable<Real> >(local_field, "field_data"))
    { }
    const CTable<Real>& coordinates;
    const CNodeConnectivity& node_connectivity;
    CField& local_field;
    CTable<Real>& field_data;
  };

  boost::shared_ptr<LoopHelper> data;

  CField::Ptr field;
};

/////////////////////////////////////////////////////////////////////////////////////

} // Mesh
} // CF

/////////////////////////////////////////////////////////////////////////////////////

#endif // CF_Mesh_COperation_hpp
