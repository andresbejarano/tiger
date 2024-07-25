# Define the utils library
add_library(TIGER_UTILS_LIB STATIC utils.cxx)
target_compile_definitions(TIGER_UTILS_LIB PUBLIC)

# Define the geometric toolkit library
add_library(TIGER_TOOLKIT_LIB STATIC ray.cxx linesegment.cxx plane.cxx)
target_compile_definitions(TIGER_TOOLKIT_LIB PUBLIC)

# Define the algorithms library
add_library(TIGER_ALGORITHMS_LIB STATIC algorithms.cxx)
target_compile_definitions(TIGER_ALGORITHMS_LIB PUBLIC)

# Define the data structures library
add_library(TIGER_DS_LIB STATIC vf.cxx dcel.cxx)
target_compile_definitions(TIGER_DS_LIB PUBLIC)

# Define the geometries library
add_library(TIGER_GEOMETRIES_LIB STATIC geometries.cxx)
target_compile_definitions(TIGER_GEOMETRIES_LIB PUBLIC)

# Define the topological interlocking configuration elements library
add_library(TIGER_TIC_LIB STATIC tessellation.cxx geometricclass.cxx block.cxx assembly.cxx interfaces.cxx)
target_compile_definitions(TIGER_TIC_LIB PUBLIC)

# Define the tic workspace library
add_library(TIGER_WORKSPACE_LIB STATIC message.cxx workspace.cxx)
target_compile_definitions(TIGER_WORKSPACE_LIB PUBLIC)

# Define the analysis library. For now it requires gurobi. I plan to make it optional
# and have a default simulated anealing solution in case gurobi is not available
add_library(TIGER_ANALYSIS_LIB STATIC equilibriumanalysis.cxx)
target_include_directories(TIGER_ANALYSIS_LIB PUBLIC "$ENV{GUROBI_HOME}/include")
set(GUROBI_LIBRARIES "$ENV{GUROBI_HOME}/lib/gurobi_c++mdd2019.lib" "$ENV{GUROBI_HOME}/lib/gurobi91.lib")
target_link_libraries(TIGER_ANALYSIS_LIB ${GUROBI_LIBRARIES})
target_compile_definitions(TIGER_ANALYSIS_LIB PUBLIC)