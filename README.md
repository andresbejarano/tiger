<div style="text-align:center;">
  <img src="https://github.com/andresbejarano/tiger/blob/main/img/tiger_pipeline_example.png" />
</div>

# TIGER - Topological Interlocking GEneratoR

TIGER is a desktop application to generate Topological Interlocking Configurations (TICs) based on free-form 3D geometric domains. It offers tools to generate some 3D surfaces and works with meshes in Wavefront (.obj) format.

TIGER has the implementation of two published TIC generation method:
- **Traditional Method:** A. J. Kanel-Belov, A. V. Dyskin, Y. Estrin, E. Pasternak, and I. A. IvanovPogodaev. Interlocking of convex polyhedra: towards a geometric theory of fragmented solids. arXiv:0812.5089 [math], December 2008. arXiv: 0812.5089.
- **Height-Bisection Method:** Andres Bejarano and Christoph Hoffmann. A Generalized Framework for Designing Topological Interlocking Configurations. International Journal of Architectural Computing, (February 2019). doi:10.1177/1478077119827187. 

To calculate the minimum tension forces required for the resultant TICs to reach static equilibrium, TIGER has an implementation of the Structure Feasibility Analysis method proposed by Whiting *et al.*: Emily Whiting, John Ochsendorf, Fredo Durand, Emily Whiting, John Ochsendorf, and Fredo Durand. Procedural modeling of structurally-sound masonry buildings. In ACM Transactions on Graphics (TOG), volume 28, page 112. ACM, December 2009.

TIGER is an ongoing project built in C++ using Eigen, Qt, VTK, and Gurobi. Before building and compiling TIGER:
- Install Qt. To make this step faster, select only your compiler of preference in the selection components tree.
- Download VTK's source code and compile it. While generating the project using CMake, select to use Qt to generate the QVTKOpenGLNativeWidget.
