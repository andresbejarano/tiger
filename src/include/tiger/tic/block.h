#ifndef _BLOCK_H_
#define _BLOCK_H_

#pragma once

#include "tessellation.h"

/*
The class to represent a block in an assembly. A block contains its own 
geometry and attributes related to the block itself.
*/
class Block : public AttributedClass, public GeometricClass
{

private:

    // Indicates whether the block is enabled or disabled. An enabled block 
    // participates in the static equilibrium analysis
    bool m_enabled;

    // The index of the face associated to the block
    size_t m_faceIndex;

public:

    /*
    Constructor of the class.
    */
    Block();

    /*
    Constructor of the class. Copies the given vertices and faces.
    @param VF & vf The reference to the vertices and faces of the geometry.
    @param size_t faceIndex The index of the face associated to the block.
    */
    Block(const VF & vf, size_t faceIndex);

    /*
    Destructor of the class.
    */
    ~Block();

    /*
    Calculates the central length of the block. The central length is the 
    distance from the centroid of the associated face with the intersection 
    point with the geometry along the normal vector of the face.
    @param const std::shared_ptr<Tessellation> tessellation The pointer to the 
    tessellation that contains the tile associated to the block.
    @param int direction The direction of the central length. The function only
    considers the sign of the given value.
    @param double threshold The threshold for values close to zero.
    @return double The central length of the block.
    */
    double CentralLength(
        const std::shared_ptr<Tessellation> tessellation, 
        int direction = 1,
        double threshold = 1e-8) const;

    /*
    Clears the content of the block.
    */
    void Clear();

    /*
    Labels the vertices of the block that are coplanar with the plane defined 
    by the normal vector N and point P. Vertex coplanarity label is stored in 
    the ATTRIB_COPLANAR dynamic attribute of each vertex.
    @param const Eigen::Vector3d & N The reference to the normal vector of the 
    plane.
    @param const Eigen::Vector3d & P The reference to the point in the plane.
    @param double threshold The threshold for values close to zero.
    @return std::vector<bool> The vector with the coplanar indicator for each 
    vertex.
    */
    std::vector<bool> CoplanarVertices(
        const Eigen::Vector3d & N,
        const Eigen::Vector3d & P,
        double threshold = 1e-8) const;

    /*
    Labels the vertices of the block that are coplanar with a plane. Vertex 
    coplanarity label is stored in the ATTRIB_COPLANAR dynamic attribute of 
    each vertex.
    @param toolkit::Plane & plane The reference to a plane.
    @param double threshold The threshold for values close to zero.
    @return std::vector<bool> The vector with the coplanar indicator for each 
    vertex.
    */
    std::vector<bool> CoplanarVertices(const toolkit::Plane & plane, double threshold = 1e-8) const;

    /*
    Disables the block. A disabled block does not participate in the Static 
    Equilibrium Analysis.
    */
    void Disable();

    /*
    Disables the block if its edges intersect with a given plane.
    @param const toolkit::Plane & plane The reference to a plane.
    */
    void DisableIfIntersects(const toolkit::Plane & plane);

    /*
    Enabled the block. An enabled block participates in the Static Equilibrium 
    Analysis.
    */
    void Enable();

    /*
    Returns the index of the face associated to the block.
    @return size_t The index of the face associated to the blocl.
    */
    size_t FaceIndex() const;

    /*
    Returns the index of the first face found to be coplanar with a given plane.
    @param const toolkit::Plane & plane The reference to a plane.
    @param double threshold The threshold for values close to zero.
    @return size_t The index of the face coplanar with the plane.
    */
    size_t GetCoplanarFace(const toolkit::Plane & plane, double threshold = 1e-8) const;

    /*
    Calculates the interface polygon between this block and another block. Keep
    in mind an interface polygon between both blocks exists only if there 
    exists a plane such that each block has a face coplanar to it, and there 
    exists intersection points between such coplanar faces. Otherwise, the 
    interface polygon doesn't exist. (MUST REVISE THIS DEFINITION)
    @param std::shared_ptr<block> block The pointer to a block.
    @param InterfacePolygon & interfacePolygon The reference to the interface 
    polygon between the blocks.
    @param double threshold The threshold for values close to zero.
    @return VF The vertices and face of the interface polygon.
    */
    VF InterfacePolygon(const std::shared_ptr<Block> block, double threshold = 1e-8) const;

private:

    /*
    Calculates the vertices and faces of the interface polygon between this 
    block and a given one.
    @param const std::shared_ptr<block> block The pointer to a block.
    @param std::list<Eigen::Vector3d> & interfaceVertices The reference to the 
    list to store the vertices of the interface.
    @param std::list<std::vector<size_t>> & interfaceFace The reference to the 
    list to store the face of the interface. (Although we assume there is only 
    one interface between two blocks, I keep this just in case we support 
    concave blocks, which could have more than one interface).
    @param const double threshold The threshold for values close to zero.
    */
    void InterfacePolygonElements(
        const std::shared_ptr<Block> block,
        std::list<Eigen::Vector3d> & interfaceVertices,
        std::list<std::vector<size_t>> & interfaceFace,
        double threshold = 1e-8) const;

public:

    /*
    Returns the geometry of the intersection between this block and a given one.
    @param const std::shared_ptr<block> block The pointer to a block.
    @param double threshold The threshold for values close to zero.
    @return VF The vertices and faces of the intersection between the blocks.
    */
    VF Intersect(const std::shared_ptr<Block> block, double threshold = 1e-8) const;

    /*
    Returns the intersection point between the block and the ray defined by the
    centroid and normal vector from the associated face of the block. By 
    default it finds the intersection point at the top of the block (along the 
    direction of the normal vector of the face); giving a negative direction 
    value makes the function to find the intersection point at the bottom of 
    the block (along the opposite direction of the normal vector of the face). 
    Since the ray origins within the block then we assume there is always an 
    intersection point.
    @param const std::shared_ptr<Tessellation> tessellation The pointer to the 
    tessellation.
    @param double & t The parameter value for the intersection point along the 
    ray.
    @param int direction The direction of the central length. The function only
    considers the sign of the given value.
    @param double threshold The threshold for values close to zero.
    @return Eigen::Vector3d The intersection point between the ray and the 
    block.
    */
    Eigen::Vector3d Intersect(
        const std::shared_ptr<Tessellation> tessellation, 
        double & t, 
        int direction = 1, 
        double threshold = 1e-8) const;

    /*
    @param const toolkit::Ray & ray
    @param double & t
    @param double threshold
    */
    Eigen::Vector3d Intersect(const toolkit::Ray & ray, double & t, double threshold = 1e-8) const;

    /*
    Indicates if the block is enabled.
    @return bool Indicates if the block is enabled.
    */
    bool IsEnabled() const;

    /*
    Returns a vector with the loads of the block. The are six values in the 
    vector: three force loads along the X, Y and Z axes, and three momentum 
    loads along the X, Y and Z axes. All but the Y force load are set to zero. 
    Y force load is the volume of the block times the density value times the 
    gravity value. Adjust the gravity value for specific units. it is usually 
    -9.8 m/s^2 if working with meters, or -980 cm/s^2 if working with 
    centimeters.
    @param double density The density of the block.
    @param double gravity The gravity value.
    @return std::vector<double> The vector with the loads of the block.
    */
    std::vector<double> Loads(double density, double gravity = -9.8) const;

    /*
    @param const toolkit::Plane & plane
    @param double extrados
    @param double intrados
    @param double threshold
    @return VF
    */
    VF PlaneOffsetClipping(
        const toolkit::Plane & plane, 
        double extrados, 
        double intrados, 
        double threshold = 1e-8) const;

    /*
    Sets the enabled status of the block.
    @para bool enabled The enabled status for the block.
    */
    void SetEnabled(bool enabled);

    /*
    Clips the geometry of the block using the Tile Offset Clipping method. It 
    is assumed the block has the index of the associated tile in the 
    tessellation stored in the ATTRIB_FACE_INDEX dynamic attribute.
    @param const std::shared_ptr<Tessellation> domain The pointer to the 
    tessellation that contains the associated tile of the block.
    @param const double extrados The extrados factor.
    @param const double intrados The intrados factor.
    @param const double threshold The threshold for values close to zero.
    @return std::shared_ptr<block> The pointer to the clipped block.
    */
    /*std::shared_ptr<block> TileOffsetClip(
        const std::shared_ptr<Tessellation> domain, 
        const double extrados, 
        const double intrados, 
        const double threshold = 1e-8) const;*/

    //Eigen::Vector3d TopActivePoint(const )

    /*
    Labels the vertices of the block as located at the extrados (above) or the 
    intrados (below) with respect of a given plane. That is, the vertices of 
    the block located above a plane in the direction of its normal vector are 
    labeled as extrados; otherwise, they are labeled as intrados. Vertex 
    location labels are stored in the ATTRIB_LOCATION dynamic attribute of the
    vertex.
    @param const toolkit::Plane & plane The reference to a plane.
    @param double threshold The threshold for values close to zero.
    @return std::vector<int> The vector with the location of the vertices.
    */
    std::vector<int> VertexLocation(const toolkit::Plane & plane, double threshold = 1e-8);

    /*
    Truncates the geometry of the given block. It is assumed the block has the 
    ATTRIB_FACE dynamic attribute. The ATTRIB_FACE attribute of the block and 
    the ATTRIB_COLOR of its faces are copied to the truncated block. The block 
    is set with the ATTRIB_MODIFIED dynamic attribute.
    @param const toolkit::Plane & plane The reference to a plane.
    @param double intrados The truncation parameter for the intrados of the 
    block. It must be a value between 0 and 1.
    @param double extrados The truncation parameter for the extrados of the 
    block. It must be a value between 0 and 1.
    @return std::shared_ptr<block> The pointer to the truncated block.
    */
    //std::shared_ptr<block> Truncate(const toolkit::Plane & plane, double intrados, double extrados) const;

    /*
    Write the information of the block.
    */
    void Write() const;

    /*
    Write the information of the block.
    @param std::stringstream & The reference to the string stream where content
    will be written.
    */
    void Write(std::stringstream & ss) const;

};

#endif
