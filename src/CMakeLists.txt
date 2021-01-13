
# Find the VTK package and the required components.
#find_package(VTK COMPONENTS
#  
#  vtkCommonColor
#  vtkCommonCore
#  vtkFiltersCore
#  vtkFiltersSources
#  vtkInfovisCore
#  vtkInteractionStyle
#  vtkRenderingContextOpenGL2
#  vtkRenderingCore
#  vtkRenderingFreeType
#  vtkRenderingGL2PSOpenGL2
#  vtkRenderingOpenGL2
#  vtkViewsQt
#  QUIET
#)
find_package(VTK REQUIRED QUIET)

# Stop processing the file if VTK was not found.
if (NOT VTK_FOUND)

  message("Skipping TIGER: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
  
endif()

# Show the VTK version.
message (STATUS "VTK_VERSION: ${VTK_VERSION}")

# Check if VTK has Qt support.
if(vtkGUISupportQt_LOADED OR TARGET VTK::GUISupportQt)

  # Enable moc and autouic for handling Qt projects.
  if (VTK_VERSION VERSION_LESS "8.90.0")
  
    # Instruct CMake to run moc automatically when needed.
    set(CMAKE_AUTOMOC ON)
	
  else()
  
     # Instruct CMake to run moc automatically when needed.
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
	
  endif()
  
  # We have ui files, this will bring in the macro: qt5_wrap_ui
  find_package(Qt5Widgets REQUIRED QUIET)
  
else()

  # Stop processing the file if VTK doesn't have Qt support.
  message(STATUS "GuiSupportQt not found.")
  return()
  
endif()

# 
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
include_directories($ENV{GUROBI_HOME}/include)

# 
file(GLOB UI_FILES *.ui)
file(GLOB QT_WRAP *.h)
file(GLOB CXX_FILES *.cxx)

qt5_wrap_ui(UISrcs ${UI_FILES})

if (VTK_VERSION VERSION_LESS "8.90.0")

  # old system
  include(${VTK_USE_FILE})
  
  # CMAKE_AUTOMOC in ON so the MOC headers will be automatically wrapped.
  add_executable(TIGER MACOSX_BUNDLE ${CXX_FILES} ${UISrcs} ${QT_WRAP})
  qt5_use_modules(TIGER Core Gui)
  
  # Define a list with the requires Gurobi library files
  #list(APPEND GUROBI_LIBRARIES "")
  #list(APPEND GUROBI_LIBRARIES "$ENV{GUROBI_HOME}/lib/gurobi_c++mdd2017.lib")
  #list(APPEND GUROBI_LIBRARIES "$ENV{GUROBI_HOME}/lib/gurobi81.lib")
  set(GUROBI_LIBRARIES "$ENV{GUROBI_HOME}/lib/gurobi_c++mdd2017.lib" "$ENV{GUROBI_HOME}/lib/gurobi81.lib")
  
  # Define a list with all required libraries
  #list(APPEND TIGER_LIBRARIES "")
  #list(APPEND TIGER_LIBRARIES ${VTK_LIBRARIES})
  #list(APPEND TIGER_LIBRARIES ${GUROBI_LIBRARIES})
  set(TIGER_LIBRARIES ${VTK_LIBRARIES} ${GUROBI_LIBRARIES})
  
  target_link_libraries(TIGER ${TIGER_LIBRARIES})
  
  
  #target_link_libraries(TIGER ${VTK_LIBRARIES})
  
else()

  # CMAKE_AUTOMOC in ON so the MOC headers will be automatically wrapped.
  add_executable(TIGER MACOSX_BUNDLE ${CXX_FILES} ${UISrcs} ${QT_WRAP})
  qt5_use_modules(TIGER Core Gui)
  
  
  # Define a list with the requires Gurobi library files
  #list(APPEND GUROBI_LIBRARIES "")
  #list(APPEND GUROBI_LIBRARIES "$ENV{GUROBI_HOME}/lib/gurobi_c++mdd2017.lib")
  #list(APPEND GUROBI_LIBRARIES "$ENV{GUROBI_HOME}/lib/gurobi81.lib")
  set(GUROBI_LIBRARIES "$ENV{GUROBI_HOME}/lib/gurobi_c++mdd2017.lib" "$ENV{GUROBI_HOME}/lib/gurobi81.lib")
  
  # Define a list with all required libraries
  #list(APPEND TIGER_LIBRARIES "")
  #list(APPEND TIGER_LIBRARIES ${VTK_LIBRARIES})
  #list(APPEND TIGER_LIBRARIES ${GUROBI_LIBRARIES})
  set(TIGER_LIBRARIES ${VTK_LIBRARIES} ${GUROBI_LIBRARIES})
  
  target_link_libraries(TIGER ${TIGER_LIBRARIES})
  
  #target_link_libraries(TIGER ${VTK_LIBRARIES})
  
  # vtk_module_autoinit is needed
  vtk_module_autoinit( TARGETS TIGER MODULES ${VTK_LIBRARIES} )
  
endif()
