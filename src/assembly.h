#ifndef _ASSEMBLY_H_
#define _ASSEMBLY_H_

#pragma once

#include "tessellation.h"
#include "block.h"

/*
The class representing an assembly of interlocking blocks. Tile tessellations and blocks are a 
1-to-1 match. That is, the index of the tile in the respective data structure is the same index of
the respective block in the assembly.
*/
class Assembly 
{

private:

    // The vector for storing the pointers to the blocks of the assembly
    std::vector<std::shared_ptr<Block>> m_blocks;

public:

    /*
    Constructor of the class.
    */
    Assembly();

    /*
    Constructor of the class.
    @param const std::vector<VF> & blocks The reference to the vector with the geometry of the 
    blocks.
    */
    Assembly(const std::vector<VF> & blocks);

    /*
    Destructor of the class.
    */
    ~Assembly();

    /*
    Clips the blocks using the adaptive tile offset method.
    @param const std::shared_ptr<Tessellation> tessellation The pointer to the tessellation with 
    the information of the faces associated to the blocks.
    @param std::string & topFunction The top adaptive function.
    @param const std::string & bottomFunction The bottom adaptive function.
    @param std::vector<VF> & clipped The vector to store the pointers to the clipped blocks.
    @param double threshold The threshold for values close to zero.
    */
    void AdaptiveTileOffsetClipping(
        const std::shared_ptr<Tessellation> tessellation,
        const std::string & topFunction, 
        const std::string & bottomFunction, 
        std::vector<VF> & clipped,
        double threshold = 1e-8) const;

    /*
    Indicates if there is at least one block that has been modified. A modified block will have the
    ATTRIB_MODIFIED dynamic attribute and its value is true.
    @return bool Indicates if the blocks have been modified.
    */
    bool AreBlocksModified() const;
    
    /*
    Calculates the interface polyons between the blocks of an assembly.
    @param std::shared_ptr<Tessellation> tessellation The pointer to the tessellation.
    @param std::vector<std::shared_ptr<VF>> & interfaces The reference to the vector for storing 
    the pointers to the geometries of the interface polygons.
    @param double threshold The threshold for values close to zero.
    @return size_t The number of interfaces between the blocks.
    */
    size_t CalculateInterfacePolygons(
        std::shared_ptr<Tessellation> tessellation,
        std::vector<std::shared_ptr<VF>> & interfaces,
        double threshold = 1e-8) const;

    /*
    Calculates the overlapping between the blocks of the assembly.
    @param std::shared_ptr<Tessellation> tessellation The pointer to the tessellation.
    @param std::vector<std::shared_ptr<VF>> & overlapping The reference to the vector to store the 
    pointers to the geometries of the overlapping.
    @param double threshold The threshold for values close to zero.
    @return bool Indicates if there is overlapping.
    */
    bool CalculateOverlapping(
        std::shared_ptr<Tessellation> tessellation,
        std::vector<std::shared_ptr<VF>> & overlapping,
        double threshold = 1e-8) const;
    
    /*
    Populates a vector with the centroid of each block in the assembly.
    @param std::vector<Eigen::Vector3d> & C The reference to the vector for storing the centroids.
    @param bool fixZeros Indicates whether to fix values close to zero.
    @param double threshold The threshold for values close to zero.
    */
    void Centroids(
        std::vector<Eigen::Vector3d> & C, 
        bool fixZeros = false, 
        double threshold = 1e-8) const;

    /*
    Clears the content of the assembly.
    */
    void Clear();

    /*
    Returns the number of blocks in the assembly.
    @return size_t The number of blocks in the assembly.
    */
    size_t CountBlocks() const;

    /*
    Returns the number of edges of the blocks of the assembly.
    @return size_t The number of edges of the blocks of the assembly.
    */
    size_t CountEdges() const;

    /*
    Returns the number of triangles required to represent the blocks of the assembly.
    @return size_t The number of triangles required to represent the blocks in the assembly.
    */
    size_t CountTriangles() const;

    /*
    Returns the number of vertices of the blocks of the assembly.
    @return size_t The number of vertices of the blocks of the assembly.
    */
    size_t CountVertices() const;

    /*
    Disable the blocks that intersect with a given plane.
    @param const toolkit::Plane & plane The reference to a plane.
    */
    void DisableIntersectedBlocks(const toolkit::Plane & plane);

    /*
    Enable all blocks of the assembly. blocks store their respective enabled value in the 
    ATTRIB_block_ENABLED dynamic attribute.
    */
    void EnableAll();

    /*
    @param const std::shared_ptr<Tessellation> tessellation
    @param const double extrados
    @param const double intrados
    @param std::vector<VF> & clipped
    @param double threshold
    */
    void TileOffsetClipping(
        const std::shared_ptr<Tessellation> tessellation,
        double extrados,
        double intrados,
        std::vector<VF> & clipped, 
        double threshold = 1e-8) const;

    /*
    Returns the reference to the vector with the pointers to the blocks of the assembly.
    @return const std::vector<block>std::vector<std::shared_ptr<block>> & The reference to the 
    vector with the pointers to the blocks of the assembly.
    */
    const std::vector<std::shared_ptr<Block>> & GetBlocks() const;

    /*
    Indicates if there are no blocks in the assembly.
    @return bool Indicates if there are no blocks in the assembly.
    */
    bool IsEmpty() const;

    /*
    Links the blocks to their respective tiles of the tessellation.
    @param std::shared_ptr<Tessellation> tessellation The pointer to the tessellation.
    */
    //void LinkBlocksToTiles(std::shared_ptr<Tessellation> tessellation);

    /*
    Returns a vector with the loads of the blocks. The are six values in per block vector: three 
    force loads along the X, Y and Z axes, and three momentum loads along the X, Y and Z axes. All 
    but the Y force load are set to zero. Y force load is the volume of the block times the density
    value times -9.8.
    @param double density The density of the blocks.
    @return std::vector<std::vector<double>> The vector with the loads that apply to the blocks.
    */
    std::vector<std::vector<double>> Loads(double density) const;

    /*
    Calculates the central lengths from the neighboring blocks of a given block.
    @param const std::shared_ptr<Tessellation> tessellation The pointer to the tessellation that 
    contains the face associated to the block.
    @param size_t index The index of the block.
    @param std::vector<double> & lenghts The vector to store the central length values from the 
    neighboring blocks of the given block.
    @param int direction The direction of the central length. The function only considers the sign 
    of the given value.
    @param double threshold The threshold for values close to zero.
    @return size_t
    */
    size_t NeighborhoodCentralLengths(
        const std::shared_ptr<Tessellation> tessellation,
        size_t index, 
        std::vector<double> & lenghts, 
        int direction = 1, 
        double threshold = 1e-8) const;

    /*
    Removes an attribute from all blocks of the assembly.
    @param std::string & name The name of the attribute.
    */
    void RemoveBlocksAttribute(const std::string & name);

    /*
    Sets the blocks. This method clears any previous content and copies the pointers given in the
    vector, no new objects are generated.
    @param std::vector<std::shared_ptr<block>> & blocks The reference to the vector with the
    pointers to the blocks.
    */
    //void Set(std::vector<std::shared_ptr<block>> & blocks);

    /*
    Sets the blocks. This method clears any previous content. Dynamic attributes from the VF 
    objects are not copied.
    @param std::vector<std::shared_ptr<block>> & blocks The reference to the vector with the
    pointers to the blocks.
    */
    //void Set(const std::vector<VF> & vfs);

    /*
    Sets the indices of the blocks. Each interface keeps its index using the ATTRIB_INDEX dynamic 
    attribute.
    */
    //void SetBlockIndex();

    /*
    Sets an attribute on all blocks of the assembly.
    @param std::string & name The name of the attribute.
    @param T value The value of the attribute.
    */
    template <typename T>
    void SetBlocksAttribute(const std::string & name, T value);

    /*
    Enables or disables the blocks according to their location. A block is enabled if its
    associated face in the tessellation has all of its neighbors (aka not at the boundary); 
    otherwise, it is disabled (since it is at the boundary).
    @param std::shared_ptr<Tessellation> tessellation The pointer to the tessellation.
    */
    void ToggleBlocksByBoundary(std::shared_ptr<Tessellation> tessellation);

    /*
    Writes the information of the blocks.
    */
    void Write() const;

    /*
    Writes the information of the blocks.
    @param std::stringstream & ss The reference to a string stream.
    */
    void Write(std::stringstream & ss) const;

    /*
    @param std::stringstream && ss
    @param const std::string & prefix
    @param int r
    @param int g
    @param int b
    @param int a
    @param double threshold
    */
    void WriteGeogebraJs(
        std::stringstream && ss, 
        const std::string & prefix,
        int r = 128,
        int g = 128,
        int b = 128,
        int a = 1,
        double threshold = 1e-8) const;

    /*
    Writes the geometry of the blocks in an OBJ file.
    @param const std::string & filepath The path to the file.
    @param double threshold The threshold for values close to zero.
    */
    void WriteObjFile(const std::string & filepath, double threshold = 1e-8) const;

};

#endif
