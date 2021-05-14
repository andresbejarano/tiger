#ifndef _WORKSPACE_H_
#define _WORKSPACE_H_

#pragma once

#include <tiger/tic/tessellation.h>
#include <tiger/tic/assembly.h>
#include <tiger/tic/interfaces.h>
#include <tiger/workspace/message.h>
#include <tiger/toolkit/plane.h>
#include <map>
#include <tiger/equilibriumanalysis.h>

//
// The class representing a workspace in TIGER. A workspace
//
class Workspace
{

private:

    // The pointer to the assembly with the blocks
    std::shared_ptr<Assembly> m_assembly;

    // The path to the directory to save the content of the workspace
    std::string m_directory;

    // The pointer to the result of the latest Static Equilibrium Analysis
    std::shared_ptr<EquilibriumAnalysis::Result> m_equilibriumResults;

    // The pointer to the interface polygons
    std::shared_ptr<InterfacePolygons> m_interfaces;

    // The pointer to the tessellation
    std::shared_ptr<Tessellation> m_tessellation;

public:

    //
    // Constructor of the class.
    //
    Workspace();

    //
    // Constructor of the class.
    // @param const std::string & directory
    //
    Workspace(const std::string & directory);

    //
    // Destructor of the class.
    //
    ~Workspace();

    //
    // Clips the blocks of an assembly using the adaptive tile offset method.
    // @param const std::string & topFunction The top function for the adaptive clipping.
    // @param const std::string & bottomFunction The bottom function for the adaptive clipping.
    // @param std::vector<VF> & blocks The vector with the pointers to the blocks to be clipped.
    // @param double threshold The threshold for values close to zero.
    // @return bool
    //
    bool AdaptiveTileOffsetClipping(
        const std::string & topFunction, 
        const std::string & bottomFunction, 
        std::vector<VF> & blocks,
        double threshold = 1e-8) const;

    //
    // Checks the requirements to generate an Abaqus script code are defined: Requirements are:
    // - There are blocks in the assembly.
    // - There are physic attributes.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckAbaqusCodeRequirements(Message& msg) const;

    //
    // Checks the requirements to generate a code with the geometry of the blocks are defined: 
    // Requirements are:
    // - There are blocks in the assembly.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckAssemblyCodeRequirements(Message& msg) const;

    //
    // Checks the requirements to locate the centroid of the tessellation to the origin. Requirements 
    // are:
    // - There is a tessellation.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckCentroidToOriginRequirements(Message& msg) const;

    //
    //
    bool CheckClippingRequirements(Message& msg) const;

    //
    // Checks the requirements to calculate the dual of the tessellation. Requirements are:
    // - There is a tessellation.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckDualTessellationRequirements(Message& msg) const;

    //
    // Checks the requirements to run the equilibrium analysis of the configuration. Requirements are:
    // - There are blocks.
    // - There are interface polygons.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckEquilibriumAnalysisRequirements(Message& msg) const;

    //
    //
    bool CheckFlipTessellationRequirements(Message& msg) const;

    //
    // Checks the requirements to fix the blocks at the bottom of the assembly. Requierements are:
    // - There is a tessellation.
    // - There are blocks.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckFixBottomBlocksRequirements(Message& msg) const;

    //
    // Checks the requirements to run the Height-Bisection Method are defined. Requirements are:
    // - There is a tessellation to work with.
    // - All tiles have an even number of sides.
    // - All tiles have a selected center point indicated by the ATTRIB_CENTER attribute.
    // - All half edges have their direction value indicated by the ATTRIB_DIRECTION_VALUE attribute.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckHeightBisectionRequirements(Message& msg) const;

    //
    // Checks the requirements to run the calculation of the interface polygons between the blocks in 
    // an assembly. Requirements are:
    // - There is a tessellation
    // - There is an assembly
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckInterfacePolygonsRequirements(Message& msg) const;

    //
    // Checks the requirements to rotate the tessellation. Requirements are:
    // - There is a tessellation.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckRotateTessellationRequirements(Message& msg) const;

    //
    // Checks the requirements to scale the blocks are defined. Requirements are:
    // - There are blocks.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckScaleBlocksRequirements(Message& msg) const;

    //
    // Checks the requirements to scale the tessellation. Requirements are:
    // - There is a tessellation.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckScaleTessellationRequirements(Message& msg) const;

    //
    // Checks the requirements to set the edge directions at the edges of the tessellation.
    // Requirements are:
    // - There is a tessellation to work with.
    // - Tiles have the center point.
    // - Tiles have an even number of sides.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckSetEdgeDirectionsRequirements(Message& msg) const;

    //
    // Checks the requirements to set the center points at the tiles of the tessellation. 
    // Requirements are:
    // - There is a tessellation to work with.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckSetFaceCentersRequirements(Message& msg) const;

    //
    // Checks the requirements to normalize the vertices of the tessellation. Requirements are:
    // - There is a tessellation to work with.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckTessellationVertexNormalizationRequirements(Message& msg) const;

    //
    // Checks the requirements to run a face subdivision on the faces of the tessellation.
    // Requirements are:
    // - There is a tessellation.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckTileSubdivisionRequirements(Message& msg) const;

    //
    // Checks the requirements to run the Tilting Angle Method are defined. Requirements are:
    // - There is a tessellation to work with.
    // - Tiles have the center point.
    // - Tiles have an even number of sides.
    // - Edges have the direction value.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckTiltingAngleRequirements(Message& msg) const;

    //
    // Checks the requirements to truncate the blocks are defined. Requirements are:
    // - There is a tessellation.
    // - There are blocks.
    // - blocks have not been modified.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckTruncateBlocksRequirements(Message& msg) const;

    //
    // Checks the requirements to normalize the vertices of the tessellation are defined. 
    // Requirements are:
    // - There is a tessellation.
    // @param Message& msg The reference to the message generated while checking the requirements.
    // @return bool Indicates if the requirements are satisfied.
    //
    bool CheckVertexNormalizationRequirements(Message& msg) const;

    //
    // Clears the content of the workspace. This includes elements in the scene and internal objects.
    //
    void Clear();

    //
    // Erases the assembly from the workspace. This also erases the interface polygons and equilirium 
    // results from the workspace since they won't be compatible anymore.
    //
    void EraseAssembly();

    //
    // Erases the edge directions from the workspace. This also erases the assembly, interface 
    // polygons, and equilirbium results from the workspace since they won't be compatible anymore.
    //
    void EraseEdgeDirections();

    //
    // Erases the Static Equilibrium Analysis results from the workspace.
    //
    void EraseEquilibriumResults();

    //
    // Erases the interface polygons from the workspace. This also erases the equilribrium results 
    // from the workspace since they won't be compatible anymore.
    //
    void EraseInterfaces();

    //
    // Erases the tessellation from the workspace. This also erases the tile centers, edge directions,
    // assembly, interface polygons, and equilibrium results from the workspace since they won't be 
    // compatible anymore.
    //
    void EraseTessellation();

    //
    // Erases the tile centers from the workspace. This also erases the edge directions, assembly, 
    // interface polygons, and equilibrium results from the workspace since they won't be compatible 
    // anymore.
    //
    void EraseTileCenters();

    //
    // Populates the containers with the new indices of the elements as they are prepared for the 
    // equilibrium analysis. Not all blocks and interface polygons in the workspace are considered for
    // equilibrium analysis (e.g., blocks at the boundary and interface polygons between blocks at the
    // boundary). Also, it prepares the centroids and loads vectors with the respective information.
    // @param const std::list<std::tuple<size_t, size_t>> & R The list with the indices tuples (as 
    // they are stored in the respective vectors in the workspace) representing the relationships 
    // between blocks and interface polygons.
    // @param std::vector<Eigen::Vector3d> & C The reference to the vector for storing the centroids 
    // of the blocks.
    // @param std::vector<std::vector<double>> & W The reference to the vector for storing the loads 
    // that affect the blocks.
    // @param double density The density of the blocks.
    // @param double gravity The gravity value.
    // @param std::vector<dcel::DCEL> & I The reference to the vector for storing the pointers to the 
    // interface polygons with at least one block enabled.
    // @param std::list<std::tuple<size_t, size_t>> & nR
    //
    void GenerateEquilibriumSystem(
        const std::list<std::tuple<size_t, size_t>> & R,
        std::vector<Eigen::Vector3d> & C,
        std::vector<std::vector<double>> & W, 
        double density, 
        double gravity, 
        std::vector<std::shared_ptr<VF>> & I, 
        std::list<std::tuple<size_t, size_t>> & nR) const;

public:

    //
    // Returns the pointer to the assembly of the workspace.
    // @return std::shared_ptr<Assembly> The pointer to the assemby.
    //
    std::shared_ptr<Assembly> GetAssembly() const;

    //
    // @param const std::shared_ptr<dcel::Face> face
    // @return VF The vertices and faces of the respective block.
    //
    VF GetBlock(const std::shared_ptr<dcel::Face> face) const;

    //
    // @param const std::shared_ptr<dcel::Face> face
    // @param std::list<Eigen::Vector3d> & blockVertices
    // @param std::list<std::vector<size_t>> & blockFaces
    //
    void GetBlockElements(
        const std::shared_ptr<dcel::Face> face,
        std::list<Eigen::Vector3d> & blockVertices,
        std::list<std::vector<size_t>> & blockFaces) const;

private:

    //
    // Populates a list with the tuples representing the relationships between enabled blocks and
    // interface polygons. An enabled block is one that lies within the assembly (not located at the
    // boundary). Each tuple in the list stores the index of the enabled block and the index of the
    // interface polygon (always in that order).
    // @param std::list<std::tuple<size_t, size_t>> & PI The reference to the list for storing the
    // tuples representing the relationships between enabled blocks and interface polygons.
    // @return bool Indicates if at least one tuple was generated.
    //
    bool GetBlockInterfaceRelations(std::list<std::tuple<size_t, size_t>> & PI) const;

    //
    // Generates the blocks based on the intersection of the planes defined by the rotated vectors in
    // the half edges of the tessellation. We assume each face has the ATTRIB_CENTER dynamic
    // attribute and each half edge has the ATTRIB_MIDPOINT, ATTRIB_DIRECTION_VALUE and
    // ATTRIB_DIRECTION_VECTOR dynamic attributes. This method generates and antiprism from each face.
    // Faces and blocks are associated to each other through the ATTRIB_FACE_INDEX (on the blocks) and
    // ATTRIB_block_INDEX (on the faces of the tessellation).
    // @param std::vector<VF> & blocks The reference to the vector to store the geometry of the
    // blocks.
    // @return bool Indicates if at least one block was generated.
    //
    bool GetBlocksFromRotatedVectors(std::vector<VF> & blocks) const;

public:

    //
    // Returns the pointer to the object with the results with the latest Static Equilibrium Analysis.
    // @return std::shared_ptr<EquilibriumAnalysis::Result> The pointer to the object with the results
    // with the latest Static Equilibrium Analysis.
    //
    std::shared_ptr<EquilibriumAnalysis::Result> GetEquilibriumResults() const;

private:

    //
    // Returns the planes of the rotated planes associated to a face in a tessellation. The 
    // rotated vector of the half edge is considered the normal vector of the rotated plane. The 
    // midpoint of the half edge is used as a point in the plane. If a half edge in the face doesn't 
    // have the ATTRIB_ROTATED_VECTOR attribute then the function stops and returns false. It is 
    // assumed all half edges in the face havethe ATTRIB_MIDPOINT attribute. Each half edge in the 
    // face produces one plane.
    // @param std::shared_ptr<dcel::Face> face The pointer to the face.
    // @param std::vector<toolkit::Plane> & planes The reference to the vector where the planes will be stored.
    // @return bool Indicates if the planes could be calculated.
    //
    bool GetFaceRotatedPlanes(std::shared_ptr<dcel::Face> face, std::vector<toolkit::Plane> & planes) const;

public:

    //
    // Returns the pointer to the interface polygons.
    // @return std::shared_ptr<InterfacePolygons> The pointer to the interface polygons.
    //
    std::shared_ptr<InterfacePolygons> GetInterfacePolygons() const;

    //
    // Returns the pointer to the tessellation of the workspace.
    // @return std::shared_ptr<Tessellation> The pointer to the tessellation.
    //
    std::shared_ptr<Tessellation> GetTessellation() const;

    //
    // Indicates whether there are Static Equilibrium Analysis results.
    // @return bool
    //
    bool HasEquilibriumResults() const;

    //
    // Indicates whether there are interface polygons.
    // @return bool 
    //
    bool HasInterfacePolygons() const;

    //
    // Indicates if there is a tessellation loaded in the workspace.
    // @return bool Indicates if there is a tessellation in the workspace.
    //
    bool HasTessellation() const;

    //
    // Indicates if the half edges of the tessellation have the rotated vectors.
    // @return bool
    //
    bool HasRotatedVectors() const;

    //
    // Generates a Topological Interlocking Configuration using the Height-Bisection method. We assume
    // each face in the tessellation has the ATTRIB_CENTER dynamic attribute and each half edge 
    // has the ATTRIB_MIDPOINT, ATTRIB_DIRECTION_VALUE and ATTRIB_DIRECTION_VECTOR dynamic attributes.
    // This method generates an antiprism from each face.
    // @param double height The height parameter value.
    // @param bool boundary Indicates whether to calculate the blocks at the boundary of the geometric
    // Tessellation.
    // @param std::vector<VF> & blocks The reference to the vector to store the geometry to the 
    // blocks.
    // @return bool Indicates if the blocks were generated.
    //
    bool HeightBisectionMethod(
        double height, 
        bool boundary, 
        std::vector<VF> & blocks) const;

    //
    // Runs the equilibrium analysis of the blocks from the assembly.
    // @param double density The density value for the blocks of the assembly.
    // @param double friction Thge friction value for surface of the the blocks.
    // @param const std::string lengthUnit The length unit. Expected values are "m" and "cm".
    // @param bool verbose Indicates whether to print the results on the terminal.
    // @param bool files Indicates whether to save the results on files.
    // @return bool Indicates if there are results from the equilibrium analysis.
    //
    bool RunEquilibriumAnalysis(
        double density, 
        double friction, 
        const std::string lengthUnit, 
        bool verbose = false, 
        bool files = false, 
        double cWeight = 1.0e5,
        double tWeight = 1.0e5,
        double uWeight = 1.0e3,
        double vWeight = 1.0e3);

    //
    //
    bool SaveAsTicFile(const std::string filename, const std::string directory = "") const;

    //
    // Scales the tessellation. This function removes 
    // @param double scale The scale factor.
    //
    void ScaleTessellation(double scale);

private:

    //
    // Sets the blocks of the assembly. Nothing else happens here.
    // @param std::vector<VF> & blocks The reference to the vector with the geometries to the blocks.
    // @param bool disableBoundary Indicates whether to disable the blocks at the boundary or not.
    //
    void SetAssembly(std::vector<VF> & blocks, bool disableBoundary = true);

public:

    // 
    // Returns a vector with the indices of the selected blocks according to 
    // the values in the given vector of strings. Accepted string values are:
    // - "all"
    // - A single index (e.g., "4") in range [0, n), where n is the number of 
    //   blocks.
    // - A range of values: (e.g., "15-20"). Each value must be in range 
    //   [0, n), where n is the number of blocks, and the first value must be
    //   less than the second.
    // @param const std::vector<std::string>& strings
    // @return std::vector<size_t>
    // 
    std::vector<size_t> SelectedBlocks(const std::vector<std::string>& strings) const;

    //
    // @param const std::string & topFunction
    // @param const std::string & bottomFunction
    // @param double threshold
    //
    void SetAssemblyUsingAdaptiveTileOffsetClippedBlocks(
        const std::string & topFunction, 
        const std::string & bottomFunction, 
        double threshold = 1e-8);

    //
    // Generates the blocks using the Height-Bisection method and Set them to the assembly.
    // @param double height The height parameter.
    // @param bool boundary Indicates whether to generate the blocks at the boundary of the geometric 
    // Tessellation.
    //
    void SetAssemblyUsingHeightBisectionMethod(double height, bool boundary);

    //
    // Clips the blocks of the assembly using the Tile Offset Clipping method and set them to the 
    // assembly.
    // @param double extrados The extrados clipping factor.
    // @param double intrados The intrados clipping factor.
    // @param double threshold The threshold for values close to zero.
    //
    void SetAssemblyUsingTileOffsetClippedBlocks(
        double extrados, 
        double intrados, 
        double threshold = 1e-8);

    //
    // Generates the blocks using the Tilting Angle method and Set them to the assembly.
    // @param double angle The angle value (in degrees).
    //
    void SetAssemblyUsingTiltingAngleMethod(double angle);

    //
    // Sets the vertex coordinates and vertex indices of the tessellation. Setting a geometric 
    // Tessellation resets any content currently loaded in the workspace.
    // @param const VF & vf The vertex coordinates and vertex indices of the geometry.
    //
    void SetTessellation(const VF & vf);

    //
    // Sets the direction values, direction vectors and midpoints of the edges of the geometric 
    // Tessellation. Setting the direction values resets any content currently loaded in the workspace.
    // @param int dir The initial direction value.
    //
    void SetTessellationEdgeDirections(int dir);

    //
    // Sets the indicated center on the tiles of the tessellation. Center point are stored in the
    // ATTRIB_CENTER dynamic attribute of each tile. We assume there is a tessellation loaded in
    // workspace.
    //
    void SetTessellationTileCenters();

    //
    // Scales the tessellation. Setting a tessellation scale resets any content currently 
    // loaded in the workspace.
    // @param double scale The scale factor.
    //
    void SetTessellationScale(double scale);

    //
    // Sets the interface polygons to the workspace.
    // @param std::vector<InterfacePolygon> & interfaces The reference to the vector with the pointers
    // to the geometries of the interface polygons.
    //
    void SetInterfacePolygons(std::vector<std::shared_ptr<VF>> & interfaces);

    //
    // @param Tessellation::TILE_SUBDIVISION type The type of face subdivision.
    //
    void SubdivideTessellationTiles(Tessellation::TILE_SUBDIVISION_TYPE type);

    //
    //
    bool TileOffsetClipping(
        double extrados,
        double intrados,
        std::vector<VF> & blocks, 
        double threshold = 1e-8) const;

    //
    // Generates a topological interlocking configuration using the Tilting Angle method. It is 
    // assumed each face has the ATTRIB_CENTER attribute and each half edge has the ATTRIB_MIDPOINT,
    // ATTRIB_DIRECTION_VALUE and ATTRIB_DIRECTION_VECTOR attributes. This method generates an 
    // antiprism from each face.
    // @param double angle The rotation angle.
    // @param std::vector<VF> & blocks The reference to the vector to store the geometry of the 
    // blocks.
    // @return bool
    //
    bool TiltingAngleMethod(double angle, std::vector<VF> & blocks) const;

    //
    // Writes the content of the workspace.
    //
    void Write() const;

    //
    // Writes the content of the workspace.
    // @param std::stringstream & ss The reference to the string stream for writing the content.
    //
    void Write(std::stringstream & ss) const;

    //
    // Writes the file with the JavaScript code to generate the assembly in GeoGebra.
    // @param const std::string filename The name of the file.
    // @param int r
    // @param int g
    // @param int b
    // @param double threshold The threshold for values close to zero.
    //
    void WriteAllAssemblyGeogebraJsFile(
        const std::string filename = "geogebra.js",
        int r = 255,
        int g = 255,
        int b = 255,
        double threshold = 1e-8) const;

    //
    // Writes the file with the centroids of the blocks.
    // @param const std::vector<Eigen::Vector3d> & C The reference to the vector with the centroids of
    // the blocks.
    // @param const std::string filename The name of the file.
    //
    void WriteCentroidsFile(
        const std::vector<Eigen::Vector3d> & C, 
        const std::string filename = "centroids.txt") const;

    //
    // Writes a JavaScript file with the geometric content of the workspace (tessellation, assembly, 
    // and interface polygons). The file content can be used in Geogebra for reproducing the 
    // goemetries.
    // @param const std::string & filename The name of the file
    // @param double threshold The threshold for values close to zero.
    //
    void WriteGeogebraJSGeometryFile(
        const std::string & filename = "geogebra.js", 
        double threshold = 1e-8) const;

    //
    // Write the file with the JavaScript code for GeoGebra with the content of the elements.
    // @param const std::vector<Eigen::Vector3d> & C
    // @param const std::vector<std::shared_ptr<VF>> & I
    // @param const std::string filename
    //
    void WriteGeogebraJsFile(
        const std::vector<Eigen::Vector3d> & C, 
        const std::vector<std::shared_ptr<VF>> & I, 
        const std::string filename = "geogebra.js") const;

    //
    // Writes the file with the JavaScript code to generate the geometry of the tessellation in
    // GeoGebra.
    // @param const std::string & filename The name of the file.
    // @param int r
    // @param int g
    // @param int b
    // @param double threshold The threshold for values close to zero.
    //
    void WriteTessellationGeogebraJsFile(
        const std::string & filename = "geogebra.js", 
        int r = 128, 
        int g = 128, 
        int b = 128, 
        double threshold = 1e-8) const;

    //
    // @param const std::string & filename
    // @param double threshold
    //
    void WriteTessellationObjFile(
        const std::string & filename = "geogebra.obj", 
        double threshold = 1e-8) const;

    //
    // Writes the file with the geometry of the interface polygons.
    // @param const std::vector<std::shared_ptr<VF>> & I The reference to the vector with the pointers
    // to the geometries of the interface polygons.
    // @param const std::string filename The name of the file.
    //
    void WriteInterfacesFile(
        const std::vector<std::shared_ptr<VF>> & I, 
        const std::string filename = "interfaces.txt") const;

    //
    // Writes the file with the geometry of the blocks.
    // @param const std::string filename The name of the file.
    //
    void WriteBlocksFile(const std::string filename = "blocks.txt") const;

    //
    // Writes the file with the relations between blocks and interface polygons.
    // @param const std::list<std::tuple<size_t, size_t>> & R The reference to the vector with the 
    // relations between blocks and interface polygons.
    // @param const std::string filename The name of the file.
    //
    void WriteRelationsFile(
        const std::list<std::tuple<size_t, size_t>> & R, 
        const std::string filename = "relations.txt") const;

};

#endif
