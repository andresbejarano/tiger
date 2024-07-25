#ifndef _DCEL_H_
#define _DCEL_H_

#pragma once

#include <tiger/ds/vf.h>
#include <tiger/toolkit/ray.h>

#include <memory>

namespace dcel
{
    // Predeclare the classes that define a Double Connected Edge List
    class Vertex;
    class Face;
    class Halfedge;

    //
    // The class to represent a vertex in a Doubly Connected Edge List.
    //
    class Vertex : public AttributedClass, public std::enable_shared_from_this<Vertex>
    {

    private:

        // The coordinates of the vertex
        Eigen::Vector3d m_coords;

    public:

        // The pointer to the incident half edge of the vertex
        std::shared_ptr<Halfedge> halfedge;

        //
        // Constructor of the class.
        //
        Vertex();

        //
        // Constructor of the class.
        // @param const Eigen::Vector3d & C The reference to a 3D vector with the coordinates.
        //
        Vertex(const Eigen::Vector3d & C);

        //
        // Constructor of the class.
        // @param double x The X coordinate of the vertex.
        // @param double y The Y coordinate of the vertex.
        // @param double z The Z coordinate of the vertex.
        //
        Vertex(double x, double y, double z);

        //
        // Destructor of the class.
        //
        ~Vertex();

        //
        // Checks the consistency of the geometric information stored in the vertex. If a test fails
        // the application throws an error.
        //
        void CheckConsistency() const;

        //
        // Clear the content of the vertex.
        //
        void Clear();

        //
        // Returns the reference to the coordinates object of the vertex.
        // @return Eigen::Vector3d & The reference to the coordinates object of the vertex.
        //
        const Eigen::Vector3d & Coords() const;

        //
        // Sets the coordinate values of the vertex.
        // @param Eigen::Vector3d & V The reference to the vector with the coordinate values of the
        // vertex.
        // @param bool fixZeros Indicates whether to fix values close to zero.
        // @param double threshold The threshold for values close to zero.
        //
        void Coords(const Eigen::Vector3d & V, bool fixZeros = false, double threshold = 1e-8);

        //
        // Sets the coordinate values of the vertex.
        // @param double x The X coordinate value of the vertex.
        // @param double y The Y coordinate value of the vertex.
        // @param double z The Z coordinate value of the vertex.
        // @param bool fixZeros Indicates whether to fix the values close to zero.
        // @param double threshold The threshold for values close to zero.
        //
        void Coords(
            double x,
            double y,
            double z,
            bool fixZeros = false,
            double threshold = 1e-8);

        //
        // Returns a String with the cooordinates of the vertex.
        // @param const std::string & sep The separator between coordinate values.
        // @param const std::string & opening The opening character.
        // @param const std::string & closing The closing character.
        // @param double threshold The threshold for values close to zero.
        // @return std::string A string with the coordinate values of the vertex.
        //
        std::string CoordsToString(
            const std::string & sep = ", ", 
            const std::string & opening = "(", 
            const std::string & closing = ")", 
            double threshold = 1e-8) const;

        //
        // Fixes the zero values of the coordinate of the vertex.
        // @param double threshold The threshold for values close to zero.
        //
        void FixZeros(double threshold = 1e-8);

		//
        // Indicates whether the vertex is fully surrounded by incident faces. That is, it does not 
        // lie at a boundary.
        // @return bool Indicates if the vertex is fully surrounded.
		//
		bool IsSurrounded() const;

        //
        // Counts the number of faces incident to the vertex.
        // @return size_t The number of faces incident to the vertex.
        //
        size_t NumberOfIncidentFaces() const;

        //
        // Writes the JavaScript commands to generate the current vertex in GeoGebra. It is assumed 
        // the vertex index is stored in the ATTRIB_INDEX dynamic attribute.
        // @param std::stringstream & ss The reference to the string stream to write the coordinates 
        // of the vertex.
        // @param const std::string prefix The prefix for the vertex.
        // @param int r The red value for the color of the vertex.
        // @param int g The green value for the color of the vertex.
        // @param int b The blue value for the color of the vertex.
        // @param double threshold The threshold for values close to zero.
        //
        void WriteGeogebraJs(
            std::stringstream & ss, 
            const std::string prefix, 
            int r = 128, 
            int g = 128, 
            int b = 128, 
            double threshold = 1e-8) const;
    };


    //
    // The class to represent a face in a Doubly Connected Edge List.
    //
    class Face : public AttributedClass, public std::enable_shared_from_this<Face>
    {

    public:

        // The pointer to the incident half edge of the face
        std::shared_ptr<Halfedge> halfedge;

        //
        // Constructor of the class.
        //
        Face();

        //
        // Destructor of the class.
        //
        ~Face();

        //
        // Calculates the area of the face.
        // @return double The area of the face.
        //
        double Area() const;

        //
        // Calculates the centroid (arithmetic mean) of the vertices of the face.
        // @param bool fixZeros Indicates whether to fix values close to zero.
        // @param double threshold The threshold for values close to zero.
        // @return Eigen::Vector3d The centroid of the face.
        //
        Eigen::Vector3d Centroid(bool fixZeros = false, double threshold = 1e-8) const;

        //
        // Checks the consistency of the geometric information stored in the face. If a test fails the
        // application throws an error.
        //
        void CheckConsistency() const;

        //
        // Clears the content of the face.
        //
        void Clear();

        //
        // Returns the number of existing face neighbors of the face.
        // @return size_t The number of existing face neighbors of the face.
        //
        size_t CountNeighbors() const;

        //
        // Returns the number of sides of the face.
        // @return size_t The number of sides of the face.
        //
        size_t CountSides() const;

        //
        // Returns the number of triangles required to represent the face. The number of triangles is 
        // the number of vertices of the face minus two.
        // @return size_t The number of triangles required to represent the face.
        //
        size_t CountTriangles() const;

        //
        // Returns the number of vertices of the face. This is an alias for the number of sides of the
        // face.
        // @return size_t The number of vertices of the face.
        //
        size_t CountVertices() const;

        //
        // @return VF
        //
        VF Flip() const;

        //
        // Checks whether the face has all of its neighbors. That is, the faces incident to the twin
        // half edges of the face exist.
        // @return bool Indicates whether all neighbors of the face exist.
        //
        bool HasAllNeighbors() const;

        //
        // Checks whether the face has an even number of sides.
        // @return bool Indicates whether the face has an even number of sides.
        //
        bool HasEvenNumberOfSides() const;

        //
        // Checks whether the face is at a boundary. That is, it misses one of its neighbors.
        // @return bool Indicates whether the face is at a boundary.
        //
        bool IsAtBoundary() const;

        //
        // Checks if the face is coplanar with a given plane. That is, the start vertices of the half
        // edges incident to the face are coplanar with the given plane.
        // @param const toolkit::Plane & plane The reference to a plane.
        // @param double threshold The threshold value for zero values.
        // @return bool Indicates whether the face is coplanar with the plane.
        //
        bool IsCoplanar(const toolkit::Plane & plane, double threshold = 1e-8) const;

        //
        // Checks if the face is planar. That is, all vertices of the face are coplanar with the plane
        // of the face.
        // @param double threshold The threshold value for zero values.
        // @return bool Indicates whether the face is planar.
        //
        bool IsPlanar(double threshold = 1e-8) const;

        //
        // Checks whether a point P lies either within the face or in one of its edges. We assume the
        // face is convex and its vertices are coplanar.
        // @param const Eigen::Vector3d & P The reference to a point.
        // @param double threshold The threshold value for zero values.
        // @return bool Indicates whether P is in the face or not.
        //
        bool IsPointIn(const Eigen::Vector3d & P, double threshold = 1e-8) const;

        //
        // Calculates the normal vector of the face.
        // @param bool normalize Indicates whether to normalize the normal vector of the face.
        // @param bool fixZeros Indicates whether to fix values close to zero.
        // @param double threshold The threshold for values close to zero.
        // @return Eigen::Vector3d The normal vector of the face.
        //
        Eigen::Vector3d Normal(
            bool normalize = false,
            bool fixZeros = false,
            double threshold = 1e-8) const;

        //
        // Calculates the plane where the face lies. The plane is defined by both the normal vector
        // and barycenter of the face.
        // @param bool normalize Indicates whether to normalize the normal vector of the plane.
        // @param bool fixZeros Indicates whether to fix values close to zero.
        // @param double threshold The threshold for values close to zero.
        // @return toolkit::Plane The plane where the face lies.
        //
        toolkit::Plane Plane(
            bool normalize = true,
            bool fixZeros = false,
            double threshold = 1e-8) const;

        //
        // Checks the location of a point with respect to the geometry and the normal vector of the
        // face. A point could be above, on or below the face.
        // @param const Eigen::Vector3d & P The reference to the point.
        // @param double threshold The threshold for values close to zero.
        // @return int The value indicating the location of the point with respecto to the face.
        //
        int PointLocation(const Eigen::Vector3d & P, double threshold = 1e-8) const;

        //
        // Sets the given attribute name and value to the vertices of the face.
        // @param std::string name The name of the attribute.
        // @param T value The value of the attribute.
        //
        template <typename T>
        void SetVerticesAttribute(const std::string & name, T value);

        //
        // Writes the vertices of the piece.
        //
        void Write() const;

        //
        // Writes the JavaScript code to generate the geometry of the face in Geogebra. It is assumed 
        // the face is planar. It is assumed the vertices of the face, and the face itself, have the 
        // ATTRIB_INDEX dynamic attribute.
        // @param std::stringstream & ss The reference to the String Stream to write the content.
        // @param const std::string & prefix The prefix of the elements.
        //
        void WritePlanarGeogebraJs(std::stringstream & ss, const std::string & prefix) const;

        //
        // Writes the JavaScript code to generate the geometry of the face in Geogebra. Face is
        // triangulated with respect of the start vertex from the incident half edge of the face. It
        // is assumed the vertices of the face, and the face itself, have the ATTRIB_INDEX dynamic 
        // attribute.
        // @param std::stringstream & ss The reference to the String Stream to write the content.
        // @param const std::string & prefix The prefix of the elements.
        //
        void WriteTriangularGeogebraJs(
            std::stringstream & ss,
            const std::string & prefix) const;

        //
        // Returns the vertex coordinates and vertex indices of the face.
        // @return VF The vertex coordinates and vertex indices of the face.
        //
        VF vf() const;

    };


    //
    // The class to represent a half edge in a Doubly Connected Edge List.
    //
    class Halfedge : public AttributedClass, public std::enable_shared_from_this<Halfedge>
    {

    public:

        // The pointer to the start vertex
        std::shared_ptr<Vertex> start;

        // The pointer to the previous half edge
        std::shared_ptr<Halfedge> previous;

        // The pointer to the next half edge
        std::shared_ptr<Halfedge> next;

        // The pointer to the twin half edge
        std::shared_ptr<Halfedge> twin;

        // The pointer to the incident face
        std::shared_ptr<Face> face;

        //
        // Constructor of the class.
        //
        Halfedge();

        //
        // Destructor of the class.
        //
        ~Halfedge();

        //
        // Calculates the direction vector of the half edge.
        // @param bool normalize Indicates whether to normalize the direction vector.
        // @param bool fixZeros Indicates whether to fix values of the direction vector close to zero.
        // @param double threshold The threshold for values close to zero.
        // @return Eigen::Vector3d The direction vector of the half edge.
        //
        Eigen::Vector3d Direction(
            bool normalize = false,
            bool fixZeros = false,
            double threshold = 1e-8) const;

        //
        // Checks the consistency of the geometric information stored in the half edge. If a test
        // fails the application throws an error.
        //
        void CheckConsistency() const;

        //
        // Clears the content of the half edge.
        //
        void Clear();

        //
        // Returns the ray defines by the coordinates of the start vertex and the direction of the 
        // half edge.
        // @param bool normalize Indicates whether to normalize the direction vector of the ray.
        // @return toolkit::Ray The ray defined by the half edge.
        //
        toolkit::Ray Ray(bool normalize = false) const;

        //
        // Calculates the intersection point of the half edge and another half edge H. If both half
        // edges are skewed then it calculates the midpoint between the closest points between both
        // half edges. The algorithm is adapted from the solution by Paul Bourke in "The shortest line
        // between two lines in 3D" found in http://paulbourke.net/geometry/pointlineplane/
        // @param std::shared_ptr<Halfedge> H The pointer to a half edge.
        // @param Eigen::Vector3d & P
        // @param double & s
        // @param double & t
        // @param double threshold The threshold for values close to zero.
        // @return bool
        //
        bool Intersect(
            std::shared_ptr<Halfedge> H,
            Eigen::Vector3d & P,
            double & s,
            double & t,
            double threshold = 1e-8) const;

        //
        // Checks if a point lies within the line segment defined by the end points of the half edge.
        // @param const Eigen::Vector3d & P The reference to a point.
        // @param double threshold The threshold for values close to zero.
        // @return bool Indicates if P lies in the half edge.
        //
        bool IsPointIn(const Eigen::Vector3d & P, double threshold = 1e-8) const;

        //
        // Calculates the length of the half edge.
        // @return double The length of the half edge.
        //
        double Length() const;

        //
        // Calculates a point along the half edge using linear interpolation with parameter t in
        // [0, 1].
        // @param double t The interpolation parameter.
        // @return Eigen::Vector3d The point along the half edge at parameter t.
        //
        Eigen::Vector3d Lerp(double t) const;

        //
        // Returns the line segment representing the half edge. It goes from the start vertex to the 
        // twin's start vertex.
        // @return toolkit::LineSegment The line segment representing the half edge.
        //
        toolkit::LineSegment LineSegment() const;

        //
        // Calculates the midpoint of the half edge.
        // @param bool fixZeros Indicates whether to fix values of the midpoint close to zero.
        // @param double threshold The threshold for values close to zero.
        // @return Eigen::Vector3d The midpoint of the half edge.
        //
        Eigen::Vector3d Midpoint(bool fixZeros = false, double threshold = 1e-8) const;

        //
        // Calculates the normal vector of the half edge.
        // @param bool normalize Indicates whether to normalize the normal vector.
        // @param bool fixZeros Indicates whether to fix the zero values of the normal vector.
        // @param double threshold The threshold for values close to zero.
        // @return Eigen::Vector3d The normal vector of the half edge.
        //
        Eigen::Vector3d Normal(
            bool normalize = false,
            bool fixZeros = false,
            double threshold = 1e-8) const;

    };


    //
    // The class representing a Doubly Connected Edge List.
    //
    class DCEL : public AttributedClass
    {

    private:

        // The vector to store the pointers to the vertices of the geometric domain
        std::vector<std::shared_ptr<Vertex>> m_vertices;

        // The vector to store the pointers to the faces of the geometric domain
        std::vector<std::shared_ptr<Face>> m_faces;

        // The vector to store the pointers to the half edges of the geometric domain
        std::vector<std::shared_ptr<Halfedge>> m_halfedges;
    
    public:

        //
        // Constructor of the class.
        //
        DCEL();

        //
        // Constructor of the class.
        // @param const VF & vf The reference to the vertex coordinates and vertex indices of the
        // geometry.
        //
        DCEL(const VF & vf);

        //
        // Constructor of the class.
        // @param const DCEL & dcel The reference to a DCEL.
        //
        DCEL(const DCEL & dcel);

        //
        // Constructor of the class.
        // @param std::string & filename The reference to the file name.
        //
        DCEL(const std::string & filename);

        //
        // Destructor of the class.
        //
        ~DCEL();

        //
        // Checks if all faces in the geometry have an even number of sides.
        // @return bool Indicates if the faces have an even number of sides.
        //
        bool AreFacesEvenSided() const;

        //
        // Calculates the Axis Aligned Bounding Box of the geometry.
        // @param Eigen::Vector3d & min The coordinates of the minimum corner of the box.
        // @param Eigen::Vector3d & max The coordinates of the maximum corner of the box.
        //
        void AxisAlignedBoundingBox(Eigen::Vector3d & min, Eigen::Vector3d & max) const;

        //
        // Calculates the centroid of the geometry. This approach is based on the divergence theorem.
        // From:  https://stackoverflow.com/questions/9325303/centroid-of-convex-polyhedron and
        // https://wwwf.imperial.ac.uk/~rn/centroid.pdf
        // @param bool fixZeros Indicates whether to fix values of the centroid close to zero.
        // @param double threshold The threshold for values close to zero.
        // @return Eigen::Vector3d The centroid of the geometry.
        //
        Eigen::Vector3d Centroid(bool fixZeros = false, double threshold = 1e-8) const;

        //
        // Checks the consistency of the geometric information stored in the elements. If a test
        // fails the application throws an error.
        //
        void CheckConsistency() const;

        //
        // Clears the content of the DCEL.
        //
        void Clear();

        //
        // Clips the geometry by a plane.
        // @param const toolkit::Plane & plane The reference to the clipping plane.
        // @param const double threshold The threshold for values close to zero.
        // @return VF The vertices and faces of the clipped geometry.
        //
        VF Clip(const toolkit::Plane & plane, const double threshold = 1e-8) const;

    private:

        //
        // @param std::list<std::vector<size_t>> & clippedFaces The vector to store the clipped faces.
        // @param bool & faceAtClippingPlane
        // @param double threshold The threshold for values close to zero.
        //
        void ClippedFaces(
            std::list<std::vector<size_t>> & clippedFaces, 
            bool & faceAtClippingPlane, 
            double threshold = 1e-8) const;

        //
        // Calculates the vertices and faces (vertex indices) of the geometry clipped by a clipping 
        // plane. The clipped geometry is the section of the original geometry that lies at the 
        // positive half space of the clipping plane.
        // @param const toolkit::Plane & plane The reference to the clipping plane.
        // @param std::list<Eigen::Vector3d> & vertices The reference to the list to store the 
        // vertices of the clipped geometry.
        // @param std::list<std::vector<size_t>> & faces The reference to the list to store the faces
        // (vertex indices) of the clipped geometry.
        // @param double threshold The threshold for values close to zero.
        //
        void ClippedElements(
            const toolkit::Plane & plane, 
            std::list<Eigen::Vector3d> & clippedVertices, 
            std::list<std::vector<size_t>> & clippedFaces, 
            double threshold = 1e-8) const;

        //
        // NOTE: The dynamic ATTRIB_INDEX attribute is set to the vertices of the geometry 
        // @param const toolkit::Plane & plane The reference to the clipping plane.
        // @param std::list<Eigen::Vector3d> & clippedVertices The reference to the list to store the
        // vertices of the clipped geometry.
        // @param size_t & pointsInPlane Returns the number of clipped vertices that lie at the 
        // clipping plane.
        // @param double threshold The threshold for values close to zero.
        //
        void ClippedVertices(
            const toolkit::Plane & plane,
            std::list<Eigen::Vector3d> & clippedVertices,
            size_t & pointsInPlane,
            double threshold = 1e-8) const;

    public:

        //
        // @return DCEL
        //
        DCEL Clone() const;

        //
        //
        void CloseLoops();

        //
        // Returns a copy of the geometry. Dynamic attributes are not copied.
        // @return DCEL A copy of the geometry.
        //
        DCEL Copy() const;

		//
        // Returns the dual of the geometry.
        // @return VF The vertices and faces of the dual geometry.
		//
		VF Dual() const;

    private:

        //
        // Calculate ther vertices and faces of the dual geometry.
        // @param std::list<Eigen::Vector3d> & dualVertices
        // @param std::list<std::vector<size_t>> & dualFaces
        //
        void DualElements(
            std::list<Eigen::Vector3d> & dualVertices, 
            std::list<std::vector<size_t>> & dualFaces) const;

    public:

        //
        // Returns the reference to the vector that stores the pointers to the faces of the geometry.
        // @return const std::vector<std::shared_ptr<Face>> & The reference to the vector that stores
        // the pointers to the faces of the geometry.
        //
        const std::vector<std::shared_ptr<Face>> & Faces() const;

        //
        // Returns the reference to the vector that stores the pointers to the half edges of the 
        // geometry.
        // @return const std::vector<std::shared_ptr<Halfedge>> & The reference to the vector that 
        // stores the half edges of the geometry.
        //
        const std::vector<std::shared_ptr<Halfedge>> & Halfedges() const;

        //
        // Checks if this geometry interstcs with a given geometry. This is a brute force algorithm, 
        // nothing fancy is hapening here.
        // @param const DCEL & dcel The reference to a geometry.
        // @param double threshold The threshold for values close to zero.
        // @return bool Indicates if the geometries intersect.
        //
        bool Intersects(const DCEL & dcel, double threshold = 1e-8) const;

        //
        // Checks if a plane intersects the geometry. A plane intersects if: 1) At least one vertex 
        // lies in the plane, or 2) At least one half edge has its end points in different half spaces
        // defined by the plane.
        // @param const toolkit::Plane & plane The reference to a plane.
        // @param double threshold The threshold for values close to zero.
        // @return bool Indicates if the plane intersects the geometry.
        //
        bool Intersects(const toolkit::Plane & plane, double threshold = 1e-8) const;

        //
        // Checks if a point lies inside a geometry. That is, the point is at the back side of all 
        // faces of the geometry, or it lies in a face. This method works for closed geometries only.
        // @param const Eigen::Vector3d & point The reference to a point.
        // @param double threshold The threshold for values close to zero.
        // @return bool Indicates if the point is inside or in the geometry.
        //
        bool IsPointIn(const Eigen::Vector3d & point, double threshold = 1e-8) const;

        //
        // Generates the vertex coordinates and vertex indices of the line segments representing the
        // edges of the geometry.
        // @return VF The vertex coordinates and vertex indices of the line segments.
        //
        VF LineSegments() const;

        //
        // Loads the content of a dcel file.
        // @param const std::string & filename The name of the file.
        //
        void LoadDcelFile(const std::string & filename);

        //
        // Normalizes the vertices of the DCEL such that they lie in a sphere of given radius. Some
        // solids such as platonic solids or archimedean solids become a sphere. Unexpected results
        // happen on free-form meshes.
        // @param double radius The radius of the sphere.
        // @param const Eigen::Vector3d & C The reference to the scale point.
        //
        void NormalizeVertices(
            double radius, 
            const Eigen::Vector3d & C = Eigen::Vector3d(0, 0, 0));

        //
        // Returns the number of internal edges of the geometry. That is, the number of edges incident
        // to two faces.
        // @return size_t The number of internal edges of the geometry.
        //
        size_t NumberOfInternalEdges() const;

        //
        // Returns the number of triangles required to represent the geometry.
        // @return size_t The number of triangles required to represent the geometry.
        //
        size_t NumberOfTriangles() const;

        //
        // Returns the number of vertices of the geometry.
        // @return size_t The number of vertices of the geometry.
        //
        size_t NumberOfVertices() const;

        //
        // Runs a quadrilateral subdivision over all faces in the DCEL. It is assumed all faces have
        // the ATTRIB_CENTER attribute with the coordinates of the selected point that subdivides each
        // face.
        //
        void QuadrangulateFaces();

        //
        // Removes an attribute from all faces in the DCEL.
        // @param std::string name The name of the attribute.
        //
        void RemoveFacesAttribute(const std::string & name) const;

        //
        // Removes an attribute from all half edges in the DCEL.
        // @param std::string & name The name of the attribute.
        //
        void RemoveHalfedgesAttribute(const std::string & name) const;

        //
        // Removes an attribute from all vertices in the DCEL.
        // @param std::string & name The name of the attribute.
        //
        void RemoveVerticesAttribute(const std::string & name) const;

        //
        //
        void Rotate(const Eigen::Vector3d & K, const double angle);

        //
        // Scales the geometry with respect of the given point and a scale factor.
        // @param double factor The scale factor.
        // @param const Eigen::Vector3d & C The reference to the scale point.
        //
        void Scale(double factor, const Eigen::Vector3d & C = Eigen::Vector3d(0.0, 0.0, 0.0));

        //
        // Sets the content of the DCEL.
        // @param const VF & vf The reference to the vertex coordinates and vertex indices.
        //
        void Set(const VF & vf);

        //
        // Sets the content of the DCEL. It copies the geometry of the given DCEL. Dynamic attributes
        // are not copied.
        // @param const DCEL & dcel The reference to a DCEL.
        //
        void Set(const DCEL & dcel);

        //
        // Sets the index of each face in the faces vector to the ATTRIB_INDEX dynamic attribute.
        //
        void SetFacesIndex() const;

        //
        // Sets an attribute on all faces in the DCEL.
        // @param std::string & name The name of the attribute.
        // @param T value The value of the attribute.
        //
        template <typename T>
        void SetFacesAttribute(const std::string & name, T value) const;

        //
        // Sets the index of each half edge in the half edges vector to the ATTRIB_INDEX dynamic
        // attribute.
        //
        void SetHalfedgesIndex() const;

        //
        // Sets an attribute on all half edges in the DCEL.
        // @param std::string & name The name of the attribute.
        // @param T value The value of the attribute.
        //
        template <typename T>
        void SetHalfedgesAttribute(const std::string & name, T value) const;

        //
        // Sets an attribute on all vertices in the DCEL.
        // @param std::string & name The name of the attribute.
        // @param T value The value of the attribute.
        //
        template <typename T>
        void SetVerticesAttribute(const std::string & name, T value) const;

        //
        // Sets the index of the vertex in the vertices vector to the ATTRIB_INDEX dynamic attribute.
        //
        void SetVerticesIndex() const;

        //
        // Sets the location of the vertices with respect of a plane. The location value of a vertex 
        // could be: 
        // 1: The vertex is in the positive half space defined by the plane.
        // 0: The vertex is in the plane.
        // -1: The vertex is in the negative half space defined by the plane.
        // Location values are stored in the ATTRIB_LOCATION dynamic attribute of each vertex.
        // @param const toolkit::Plane & plane The reference to a plane.
        // @param double threshold The threshold for values close to zero.
        //
        void SetVerticesLocation(const toolkit::Plane & plane, double threshold = 1e-8) const;

        //
        //
        Eigen::Vector3d SimpleCentroid(bool fixZeros = false, double threshold = 1e-8) const;

    private:

        //
        // Splits a half edge by its midpoint. This procedure adds a new vertex to the geometry and 
        // generates the respective half edges.
        // @param std::shared_ptr<Halfedge> The pointer to the half edge to be split.
        //
        void SplitHalfedgeByMidpoint(std::shared_ptr<Halfedge> halfedge);

    public:

        //
        // Splits each half edge by its midpoint into two half edges. All vertices are set the 
        // ATTRIB_ORIGINAL dynamic attribute indicating whether they are part of the original geometry
        // or not.
        //
        void SplitHalfedgesByMidpoint();

		//
		//
		void SubdivideFaceByMidpoints(std::shared_ptr<Face> face);

    private:

        //
        // Subdivides a face into quadrilaterals.
        // @param std::shared_ptr<Face> face The pointer to the face to be subdivided.
        //
        void SubdivideFaceIntoQuadrilaterals(std::shared_ptr<Face> face);

    public:

		//
		//
		void SubdivideFacesByMidpoints();

        //
        //
        void Translate(const Eigen::Vector3d & V);

		//
        // Triangulates a face with respect of a center point and the vertices of the face.
        // @param std::shared_ptr<Face> face The pointer to the face.
        // @param const Eigen::Vector3d & C The reference to the center point.
		//
		void TriangulateFaceByVertices(std::shared_ptr<Face> face, const Eigen::Vector3d & C);

		//
        // Subdivides the faces of the geometric domain into triangles with respect of the center of 
        // each face and its vertices. We assume each face in the geometric domain has the 
        // ATTRIB_CENTER dynamic attribute.
		//
		void TriangulateFacesByVertices();

        //
        // Returns the reference to the vector that stores the pointers to the vertices of the 
        // geometry.
        // @return const std::vector<std::shared_ptr<Vertex>> & The reference to the vector that 
        // stores the pointers to the vertices of the geometry.
        //
        const std::vector<std::shared_ptr<Vertex>> & Vertices() const;

        //
        // Returns the vertex coordinates and vertex indices of the geometry.
        // @return VF The vertex coordinates and vertex indices of the geometry.
        //
        VF vf() const;

        //
        // Calculates the volume of the geometry. The volumes is calculated using the Divergence
        // Theorem. Sources:
        // - https://en.wikipedia.org/wiki/Polyhedron#Volume
        // - https://www.quora.com/How-can-I-obtain-the-volume-of-an-irregular-body-if-I-just-know-the-vertex-coordinates-of-the-irregular-body
        // @return double The volume of the geometry.
        //
        double Volume() const;

        //
        // Writes the information of the geometry.
        //
        void Write() const;

        //
        // Writes the information of the geometry.
        // @param std::ofstream & The reference to the output file stream. It should be open.
        //
        void Write(std::ofstream & file) const;

        //
        // Writes the information of the geometry.
        // @param std::stringstream & The reference to the string stream to write the content.
        //
        void Write(std::stringstream & ss) const;

        //
        // Writes the information of the geometry in a file.
        // @param std::string & filename The reference to the file name.
        //
        void WriteDcelFile(const std::string & filename) const;

        //
        // @param std::stringstream & ss
        // @param const std::string & prefix
        // @param int r
        // @param int g
        // @param int b
        // @param int a
        // @param double threshold The threshold for values close to zero.
        //
        void WriteGeogebraJs(
            std::stringstream & ss, 
            const std::string & prefix, 
            int r = 128, 
            int g = 128, 
            int b = 128, 
            int a = 1, 
            double threshold = 1e-8) const;

    };

}

#endif
