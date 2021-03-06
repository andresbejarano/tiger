#ifndef _BLOCK_H_
#define _BLOCK_H_

#pragma once

#include "tessellation.h"

//
// The class to represent a block in an assembly. A block contains its own 
// geometry and attributes related to the block itself.
//
class Block : public AttributedClass, public GeometricClass
{

private:

    // The density of the block. We require this value for the static 
    // equilibrium analysis and normalizing the resultant forces (if 
    // indicated). By default, its value is 1.0.
    double m_density;

    // Indicates whether the block is enabled or disabled. An enabled block 
    // participates in the static equilibrium analysis
    bool m_enabled;

    // The index of the face associated to the block
    size_t m_faceIndex;

    // The array to store the loads on the block. The first 3 values correspond
    // to the force loads (x,y,z). The following 3 values correspond to the 
    // torque loads (x,y,z)
    std::array<double, 6> m_loads;

public:

    //
    // Constructor of the class.
    //
    Block();

    //
    // Constructor of the class. Copies the given vertices and faces.
    // @param VF & vf The reference to the vertices and faces of the geometry.
    // @param size_t faceIndex The index of the face associated to the block.
    //
    Block(const VF & vf, size_t faceIndex);

    //
    // Destructor of the class.
    //
    ~Block();

    //
    // Adds a force to the current force load values.
    // @param double x The X component of the added force.
    // @param double y The Y component of the added force.
    // @param double z The Z component of the added force.
    //
    void AddForceLoad(double x, double y, double z);

    //
    // Adds a torque to the current torque load values.
    // @param double x The X component of the added force.
    // @param double y The Y component of the added force.
    // @param double z The Z component of the added force.
    //
    void AddTorqueLoad(double x, double y, double z);

    //
    // Calculates the central length of the block. The central length is the 
    // distance from the centroid of the associated face with the intersection 
    // point with the geometry along the normal vector of the face.
    // @param const std::shared_ptr<Tessellation> tessellation The pointer to the 
    // tessellation that contains the tile associated to the block.
    // @param int direction The direction of the central length. The function only
    // considers the sign of the given value.
    // @param double threshold The threshold for values close to zero.
    // @return double The central length of the block.
    //
    double CentralLength(
        const std::shared_ptr<Tessellation> tessellation, 
        int direction = 1,
        double threshold = 1e-8) const;

    //
    // Clears the content of the block.
    //
    void Clear();

    //
    // Labels the vertices of the block that are coplanar with the plane defined 
    // by the normal vector N and point P. Vertex coplanarity label is stored in 
    // the ATTRIB_COPLANAR dynamic attribute of each vertex.
    // @param const Eigen::Vector3d & N The reference to the normal vector of the 
    // plane.
    // @param const Eigen::Vector3d & P The reference to the point in the plane.
    // @param double threshold The threshold for values close to zero.
    // @return std::vector<bool> The vector with the coplanar indicator for each 
    // vertex.
    //
    std::vector<bool> CoplanarVertices(
        const Eigen::Vector3d & N,
        const Eigen::Vector3d & P,
        double threshold = 1e-8) const;

    //
    // Labels the vertices of the block that are coplanar with a plane. Vertex 
    // coplanarity label is stored in the ATTRIB_COPLANAR dynamic attribute of 
    // each vertex.
    // @param toolkit::Plane & plane The reference to a plane.
    // @param double threshold The threshold for values close to zero.
    // @return std::vector<bool> The vector with the coplanar indicator for each 
    // vertex.
    //
    std::vector<bool> CoplanarVertices(const toolkit::Plane & plane, double threshold = 1e-8) const;

    //
    // Disables the block. A disabled block does not participate in the Static 
    // Equilibrium Analysis.
    //
    void Disable();

    //
    // Disables the block if its edges intersect with a given plane.
    // @param const toolkit::Plane & plane The reference to a plane.
    //
    void DisableIfIntersects(const toolkit::Plane & plane);

    //
    // Enabled the block. An enabled block participates in the Static Equilibrium 
    // Analysis.
    //
    void Enable();

    //
    // Returns the index of the face associated to the block.
    // @return size_t The index of the face associated to the blocl.
    //
    size_t FaceIndex() const;

    //
    // Returns the index of the first face found to be coplanar with a given plane.
    // @param const toolkit::Plane & plane The reference to a plane.
    // @param double threshold The threshold for values close to zero.
    // @return size_t The index of the face coplanar with the plane.
    //
    size_t GetCoplanarFace(const toolkit::Plane & plane, double threshold = 1e-8) const;

    // 
    // Returns the density of the block.
    // @return double
    double GetDensity() const;

    // 
    // Returns the array with the force and torque loads of the block.
    // @return const std::array<double, 6>&
    // 
    const std::array<double, 6>& GetLoads() const;

    //
    // Calculates the interface polygon between this block and another block. Keep
    // in mind an interface polygon between both blocks exists only if there 
    // exists a plane such that each block has a face coplanar to it, and there 
    // exists intersection points between such coplanar faces. Otherwise, the 
    // interface polygon doesn't exist. (MUST REVISE THIS DEFINITION)
    // @param std::shared_ptr<block> block The pointer to a block.
    // @param InterfacePolygon & interfacePolygon The reference to the interface 
    // polygon between the blocks.
    // @param double threshold The threshold for values close to zero.
    // @return VF The vertices and face of the interface polygon.
    //
    VF InterfacePolygon(const std::shared_ptr<Block> block, double threshold = 1e-8) const;

private:

    //
    // Calculates the vertices and faces of the interface polygon between this 
    // block and a given one.
    // @param const std::shared_ptr<block> block The pointer to a block.
    // @param std::list<Eigen::Vector3d> & interfaceVertices The reference to the 
    // list to store the vertices of the interface.
    // @param std::list<std::vector<size_t>> & interfaceFace The reference to the 
    // list to store the face of the interface. (Although we assume there is only 
    // one interface between two blocks, I keep this just in case we support 
    // concave blocks, which could have more than one interface).
    // @param const double threshold The threshold for values close to zero.
    //
    void InterfacePolygonElements(
        const std::shared_ptr<Block> block,
        std::list<Eigen::Vector3d> & interfaceVertices,
        std::list<std::vector<size_t>> & interfaceFace,
        double threshold = 1e-8) const;

public:

    //
    // Returns the geometry of the intersection between this block and a given one.
    // @param const std::shared_ptr<block> block The pointer to a block.
    // @param double threshold The threshold for values close to zero.
    // @return VF The vertices and faces of the intersection between the blocks.
    //
    VF Intersect(const std::shared_ptr<Block> block, double threshold = 1e-8) const;

    //
    // Returns the intersection point between the block and the ray defined by the
    // centroid and normal vector from the associated face of the block. By 
    // default it finds the intersection point at the top of the block (along the 
    // direction of the normal vector of the face); giving a negative direction 
    // value makes the function to find the intersection point at the bottom of 
    // the block (along the opposite direction of the normal vector of the face). 
    // Since the ray origins within the block then we assume there is always an 
    // intersection point.
    // @param const std::shared_ptr<Tessellation> tessellation The pointer to the 
    // tessellation.
    // @param double & t The parameter value for the intersection point along the 
    // ray.
    // @param int direction The direction of the central length. The function only
    // considers the sign of the given value.
    // @param double threshold The threshold for values close to zero.
    // @return Eigen::Vector3d The intersection point between the ray and the 
    // block.
    //
    Eigen::Vector3d Intersect(
        const std::shared_ptr<Tessellation> tessellation, 
        double & t, 
        int direction = 1, 
        double threshold = 1e-8) const;

    //
    // @param const toolkit::Ray & ray
    // @param double & t
    // @param double threshold
    //
    Eigen::Vector3d Intersect(const toolkit::Ray & ray, double & t, double threshold = 1e-8) const;

    //
    // Indicates if the block is enabled.
    // @return bool Indicates if the block is enabled.
    //
    bool IsEnabled() const;

    //
    // Returns a vector with the loads of the block. The are six values in the 
    // vector: three force loads along the X, Y and Z axes, and three momentum 
    // loads along the X, Y and Z axes. All but the Y force load are set to zero. 
    // Y force load is the volume of the block times the density value times the 
    // gravity value. Adjust the gravity value for specific units. it is usually 
    // -9.8 m/s^2 if working with meters, or -980 cm/s^2 if working with 
    // centimeters.
    // @param double density The density of the block.
    // @param double gravity The gravity value.
    // @return std::vector<double> The vector with the loads of the block.
    //
    //std::vector<double> Loads(double density, double gravity = -9.8) const;

    //
    // @param const toolkit::Plane & plane
    // @param double extrados
    // @param double intrados
    // @param double threshold
    // @return VF
    //
    VF PlaneOffsetClipping(
        const toolkit::Plane & plane, 
        double extrados, 
        double intrados, 
        double threshold = 1e-8) const;

    // 
    // 
    // 
    void SetDensity(double density);

    //
    // Sets the enabled status of the block.
    // @para bool enabled The enabled status for the block.
    //
    void SetEnabled(bool enabled);

    //
    //
    void SetForceLoad(double x, double y, double z);

    //
    //
    void SetTorqueLoad(double x, double y, double z);

    //
    // Labels the vertices of the block as located at the extrados (above) or the 
    // intrados (below) with respect of a given plane. That is, the vertices of 
    // the block located above a plane in the direction of its normal vector are 
    // labeled as extrados; otherwise, they are labeled as intrados. Vertex 
    // location labels are stored in the ATTRIB_LOCATION dynamic attribute of the
    // vertex.
    // @param const toolkit::Plane & plane The reference to a plane.
    // @param double threshold The threshold for values close to zero.
    // @return std::vector<int> The vector with the location of the vertices.
    //
    std::vector<int> VertexLocation(const toolkit::Plane & plane, double threshold = 1e-8);

    //
    // Write the information of the block.
    //
    void Write() const;

    //
    // Write the information of the block.
    // @param std::stringstream & The reference to the string stream where content
    // will be written.
    //
    void Write(std::stringstream & ss) const;

};

#endif
