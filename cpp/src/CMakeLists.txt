#
include_directories(include)

# Build the basic libraries
add_subdirectory(libs)

# 
set(TIGER_BASIC_LIBRARIES 
	TIGER_UTILS_LIB 
	TIGER_TOOLKIT_LIB 
	TIGER_DS_LIB 
	TIGER_ALGORITHMS_LIB 
	TIGER_GEOMETRIES_LIB 
	TIGER_TIC_LIB 
	TIGER_WORKSPACE_LIB 
	TIGER_ANALYSIS_LIB)

# Build the applications
add_subdirectory(apps/qt_vtk)
