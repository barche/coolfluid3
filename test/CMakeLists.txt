include_directories( ${CMAKE_SOURCE_DIR} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

add_subdirectory( rapidxml )
add_subdirectory( Eigen )

add_subdirectory( common )
add_subdirectory( math )
add_subdirectory( mesh )
add_subdirectory( physics )
add_subdirectory( solver )
add_subdirectory( Tools )
if(CF3_ENABLE_GUI)
  add_subdirectory( ui )
endif()
add_subdirectory( python )
