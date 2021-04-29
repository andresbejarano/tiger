#ifndef _VF_H_
#define _VF_H_

#pragma once

#include <tiger/attributes.h>
#include <tiger/toolkit/plane.h>
#include <tiger/toolkit/linesegment.h>
#include <Eigen/core>

#include <list>
#include <vector>

/*
The class to represent the vertices and faces of a geometry.
*/
class VF : public AttributedClass
{

private:

    // The index to add a new face
    size_t m_addFaceIndex;

    // The index to add a new vertex
    Eigen::Index m_addVertexIndex;

    // The vector with the faces of the geometry
    std::vector<std::vector<size_t>> m_F;

    // The vector with the dynamic attributes of the faces
    std::vector<DynamicAttributes> m_facesAttributes;

    // The matrix with the vertices of the geometry. Each vertex has 
    // coordinates (x, y, z, w)
    Eigen::MatrixXd m_V;

public:

    enum FACE_CENTER {BARYCENTER, CENTROID};

    /*
    Constructor of the class.
    */
    VF();

    /*
    Constructor of the class.
    @param const size_t nV The number of vertices;
    @param const size_t nF The number of faces.
    */
    VF(const size_t nV, const size_t nF);

    /*
    Constructor of the class.
    @param const std::list<Eigen::Vector3d> & vertices The reference to the 
    vector with the vertices.
    @param const std::list<std::vector<size_t>> & faces The reference to the 
    vector with the faces.
    */
    VF(
        const std::list<Eigen::Vector3d> & vertices, 
        const std::list<std::vector<size_t>> & faces);

    /*
    Constructor of the class.
    @param const VF & vf The reference to another VF object.
    */
    VF(const VF & vf);

    /*
    Destructor of the class.
    */
    ~VF();

    /*
    Adds a face with two sides.
    @param size_t v0 The index of the first vertex.
    @param size_t v1 The index of the second vertex.
    @return size_t The index of the new face.
    */
    size_t addFace(size_t v0, size_t v1);

    /*
    Adds a face with three sides.
    @param size_t v0 The index of the first vertex.
    @param size_t v1 The index of the second vertex.
    @param size_t v2 The index of the third vertex.
    @return size_t The index of the new face.
    */
    size_t addFace(size_t v0, size_t v1, size_t v2);

    /*
    Adds a face with four sides.
    @param size_t v0 The index of the first vertex.
    @param size_t v1 The index of the second vertex.
    @param size_t v2 The index of the third vertex.
    @param size_t v3 The index of the fourth index.
    @return size_t The index of the new face.
    */
    size_t addFace(size_t v0, size_t v1, size_t v2, size_t v3);

    /*
    Adds a face.
    @param const std::vector<size_t> & indices The reference to the vector with
    the indices of the new face.
    @return size_t The index of the new face.
    */
    size_t addFace(const std::vector<size_t> & indices);

    /*
    Adds a vertex.
    @param double x The X coordinate value of the vertex.
    @param double y The Y coordinate value of the vertex.
    @param double z The Z coordinate value of the vertex.
    @return size_t The index of the new vertex.
    */
    size_t addVertex(double x, double y, double z);

    /*
    Adds a vertex. Returns the index of the vertex.
    @param const Eigen::Vector3d & P The reference to a point.
    @return size_t The index of the new vertex.
    */
    size_t addVertex(const Eigen::Vector3d & P);

    /*
    * Calculates the area of the geometry.
    * @return double The area of the geometry.
    */
    double area() const;

    /*
    * Calculates the area of the face at the given index.
    * @return double The area of the face.
    */
    double area(size_t fIdx) const;

    /*
    Calculates the Axis Aligned Bounding Box of the geometry.
    @param Eigen::Vector3d & min The coordinates of the minimum corner of the 
    box.
    @param Eigen::Vector3d & max The coordinates of the maximum corner of the 
    box.
    */
    void axisAlignedBoundingBox(Eigen::Vector3d & min, Eigen::Vector3d & max) const;

    /*
    Calculates the centroid (arithmetic mean) of the vertices of the geometry.
    @param bool fixZeros Indicates whether to fix the values close to zero.
    @param double threshold The threshold for values close to zero.
    @return Eigen::Vector3d The centroid of the geometry.
    */
    Eigen::Vector3d centroid(bool fixZeros = false, double threshold = 1e-8) const;

    /*
    Returns the centroid of the face at a given index. Centroid is calculated 
    as the arithmetic mean of the vertices.
    @param size_t index The index of the face.
    @param bool fixZeros
    @param double threshold
    @return Eigen::Vector3d The centroid of the face.
    */
    Eigen::Vector3d centroid(size_t index, bool fixZeros = false, double threshold = 1e-8) const;

    /*
    Clears the content.
    */
    void clear();

    /*
    Clones the object. Dynamic attributes are not cloned though.
    @return VF The cloned object.
    */
    VF clone() const;

    /*
    Returns the number of edges of the geometry.
    @return size_t The number of edges of the geometry.
    */
    size_t countEdges() const;

    /*
    Returns the number of faces of the geometry.
    @return size_t The number of faces.
    */
    size_t countFaces() const;

    /*
    Returns the number of sides of a face in the geometry.
    @param size_t face The index of the face.
    @return size_t The number of sides of the face.
    */
    size_t countSides(size_t face) const;

    /*
    Returns the number of triangles to describe the geometry.
    @return size_t The number of triangles to describe the geometry.
    */
    size_t countTriangles() const;

    /*
    Returns the number of triangles required to represent a face in the 
    geometry.
    @param size_t face The index of the face.
    @return size_t The number of triangles.
    */
    size_t countTriangles(size_t face) const;

    /*
    Returns the number of vertices of the geometry.
    @return size_t The number of vertices of the geometry.
    */
    size_t countVertices() const;

    /*
    Returns the number of vertices of a face in the geometry.
    @param size_t face The index of the face.
    @return size_t The number of vertices of the face.
    */
    size_t countVertices(size_t face) const;

    /*
    Returns the direction vector of an edge in a face.
    @param size_t face The index of the face.
    @param size_t edge The index of the edge.
    @param bool normalize Indicates whether to normalize or not the direction 
    vector.
    @param bool fixZeros Indicates whether to fix or not the zero values of the
    direction vector.
    @param double threshold The threshold for values close to zero.
    @return Eigen::Vector3d The direction of the edge in the face.
    */
    Eigen::Vector3d direction(
        size_t face, 
        size_t edge, 
        bool normalize = false, 
        bool fixZeros = false, 
        double threshold = 1e-8) const;

    /*
    Returns the reference to the vector with the vertex indices of a face.
    @param size_t index The index of the face.
    @return const std::vector<size_t> & The reference to the vector with the 
    vertex indices of the face.
    */
    const std::vector<size_t> & face(size_t index) const;

    /*
    @param size_t index
    @return DynamicAttributes & The reference to the dynamic attributes.
    */
    DynamicAttributes & faceAttributes(size_t index);

    /*
    Returns the coordinates of the vertices of a face in a matrix.
    @param size_t face The index of the face.
    @return Eigen::MatrixXd The matrix with the coordinates of the face.
    */
    Eigen::MatrixXd faceVertices(size_t face) const;

    /*
    Checks if all faces have an even number of sides.
    @return bool Indicates whether all faces have an even number of sides.
    */
    bool facesHaveEvenNumberOfSides() const;

    /*
    Inverts the orientation of the faces.
    */
    void flip();

    /*
    Checks if a face has an even number of sides.
    @param size_t face The index of the face.
    @return bool Indicates whether the face has an even number of sides or not.
    */
    bool HasEvenNumberOfSides(size_t face) const;

    /*
    Indicates whether the information of the geometry is complete. That is, it 
    has the number of vertices and faces indicated when constructing the object.
    */
    bool IsComplete() const;

    /*
    Checks if a face in the geometry is coplanar with a plane.
    @param size_t face The index of the face.
    @param const toolkit::Plane & plane The reference to the plane.
    @param double threshold The threshold for values close to zero.
    @return bool Indicates whether the face is coplanar with the given plane or
    not.
    */
    bool IsCoplanar(size_t face, const toolkit::Plane & plane, double threshold = 1e-8) const;

    /*
    Checks if a face is planar. That is, all vertices of the face are coplanar 
    with the plane of the face.
    @param size_t face The index of the face.
    @param double threshold The threshold value for zero values.
    @return bool Indicates whether the face is planar.
    */
    bool IsPlanar(size_t face, double threshold = 1e-8) const;

    /*
    Checks whether a point P lies either within a given face or in one of its 
    edges. We assume the face is convex and its vertices are coplanar.
    @param size_t face The index of the face.
    @param const Eigen::Vector3d & P The reference to the point.
    @param double threshold The threshold for values close to zero.
    @return bool Indicates whether the point is in the face or not.
    */
    bool IsPointIn(size_t face, const Eigen::Vector3d & P, double threshold = 1e-8) const;

    /*
    @param const Eigen::Vector3d & P The reference to the point.
    @param double threshold The threshold for values close to zero.
    @return bool Indicates whether the point is in the face or not.
    */
    bool IsPointInSurface(const Eigen::Vector3d & P, double threshold = 1e-8) const;

    /*
    Returns the normal vector of a face.
    @param size_t face The index of the face.
    @param bool normalize Indicates whether to normalize or not the normal 
    vector.
    @param bool fixZeros Indicates whether to fix or not the zeros of the 
    normal vector.
    @param double threshold The threshold for values close to zero.
    @return Eigen::Vector3d The normal vector of the face.
    */
    Eigen::Vector3d Normal(
        size_t face, 
        bool normalize = false, 
        bool fixZeros = false, 
        double threshold = 1e-8) const;

    /*
    Normalizes the vertices of the geometry.
    */
    void NormalizeVertices();

    /*
    Returns the plane that contains the information of the faces at a given 
    index.
    @param size_t index The index of the face.
    @param bool normalize Indicates whether to normalize the normal vector of 
    the plane.
    @param bool fixZeros Indicates whether to fix the zero values.
    @param double threshold The threshold for values close to zero.
    @return toolkit::Plane The plane.
    */
    toolkit::Plane Plane(
        size_t face, 
        bool normalize = false, 
        bool fixZeros = false, 
        double threshold = 1e-8) const;

    /*
    @param size_t face The index of the face.
    @param const Eigen::Vector3d & point
    @param double threshold
    @return int
    */
    int PointLocation(size_t face, const Eigen::Vector3d & point, double threshold = 1e-8) const;

    /*
    Resets the content of the vertices and faces.
    */
    void Reset();

    /*
    Rotates the vertices of the geometry.
    @param const Eigen::Vector3d & K The reference to the axis vector.
    @param double angle The rotation angle in radians.
    */
    void Rotate(const Eigen::Vector3d & K, double angle);

    /*
    Scales the geometry.
    @param double s The scale factor.
    @param const Eigen::Vector3d & C The reference to the scale pivot point.
    */
    void Scale(double s, const Eigen::Vector3d & C = Eigen::Vector3d::Zero());

    /*
    @param const VF & vf
    */
    void Set(const VF & vf);

    /*
    Returns the centroid of the vertices by calculating their arithmetic mean.
    @return Eigen::Vector3d The coordinates of the centroid of the vertices.
    */
    //Eigen::Vector3d SurfaceCentroid() const;

    /*
    @param const Eigen::Vector3d & D The reference to the translation vector.
    */
    void Translate(const Eigen::Vector3d & D);

    /*
    */
    VF TriangulateFacesByVertices() const;

    /*
    Calculates the volume of the geometry. The volumes is calculated using the 
    Divergence
        Theorem. Sources:
        - https://en.wikipedia.org/wiki/Polyhedron#Volume
        - https://www.quora.com/How-can-I-obtain-the-volume-of-an-irregular-body-if-I-just-know-the-vertex-coordinates-of-the-irregular-body
        @return double The volume of the geometry.
    */
    //double Volume() const;

    /*
    Returns a vertex.
    @param size_t index The index of the vertex.
    @return Eigen::Vector3d & The vertex.
    */
    Eigen::Vector3d Vertex(size_t index) const;

    /*
    @param size_t index
    @param const Eigen::Vector3d& P
    */
    void Vertex(size_t index, const Eigen::Vector3d& P);

    /*
    @param size_t index
    @param double x
    @param double y
    @param double z
    */
    void Vertex(size_t index, double x, double y, double z);

    /*
    Calculates the volume of the geometry. This method works on closed 
    geometries only. Volume calculation from 
    https://en.wikipedia.org/wiki/Polyhedron#Volume (formal approach in 
    https://wwwf.imperial.ac.uk/~rn/centroid.pdf
    @return double The volume of the geometry.
    */
    double Volume() const;

    /*
    Writes the vertices and faces in the console.
    */
    void Write() const;

    /*
    */
    void Write(std::stringstream & ss) const;

    /*
        @param std::stringstream & ss
        @param const std::string & prefix
        @param int r
        @param int g
        @param int b
        @param int a
        @param double threshold The threshold for values close to zero.
        */
    void WriteGeogebraJs(
        std::stringstream & ss,
        const std::string & prefix,
        int r = 128,
        int g = 128,
        int b = 128,
        int a = 1,
        double threshold = 1e-8) const;

    /*
    Saves the geometry in a OBJ file.
    @param const std::string & filepath The path to the file.
    @param double threshold The threshold for values close to zero.
    */
    void WriteObjFile(const std::string & filepath, double threshold = 1e-8) const;
};

#endif
