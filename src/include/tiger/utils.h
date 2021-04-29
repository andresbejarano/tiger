#ifndef _UTILS_H_
#define _UTILS_H_

#pragma once

#include <Eigen/core>
#include <list>
#include <vector>

/*
*/
namespace utils
{
    // The constant PI
    const double PI = acos(-1.0);

    /*
    Checks if all face vertices have the ATTRIB_COPLANAR attribute set with value true.
    @param const std::shared_ptr<dcel::Face> face The pointer to a face.
    @return bool Indicates if all face vertices have the ATTRIB_COPLANAR attribute set to true.
    */
    //bool areFaceVerticesCoplanar(const std::shared_ptr<dcel::Face> face);

    /*
    Rotates vector V around an axis vector K and angle of rotation a. This method is also named as
    Rodrigues' rotation formula: (V * cos(a)) + ((K x N) * sin(a)) + (K * (K . V) * (1 - cos(a)))
    Its description is in  https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    @param const Eigen::Vector3d & V The reference to the vector to be rotated.
    @param const Eigen::Vector3d & K The reference to the rotation axis vector.
    @param double angle The rotation angle (in radians).
    @return Eigen::Vector3d The rotated vector.
    */
    Eigen::Vector3d axisAngleRotation(const Eigen::Vector3d & V, const Eigen::Vector3d & K, double angle);

    /*
    Calculates the centroid (arithmetic mean) of a list of points.
    @param const std::list<Eigen::Vector3d> & points The reference to a list with the points.
    @param bool fixZeros Indicates whether to fix the values close to zero.
    @param double threshold The threshold for values close to zero.
    @return Eigen::Vector3d The centroid of the points.
    */
    Eigen::Vector3d centroid(const std::list<Eigen::Vector3d> & points, bool fixZeros = false, double threshold = 1e-8);

    /*
    Calculates the centroid (arithmetic mean) of a vector of points.
    @param const std::vector<Eigen::Vector3d> & points The reference to a vector with the points.
    @param bool fixZeros Indicates whether to fix the values close to zero.
    @param double threshold The threshold for values close to zero.
    @return Eigen::Vector3d The centroid of the points.
    */
    Eigen::Vector3d centroid(const std::vector<Eigen::Vector3d> & points, bool fixZeros = false, double threshold = 1e-8);

    /*
    Returns the determinant of the 3x3 matrix specified as:
    |a, b, c|
    |d, e, f|
    |g, h, i|
    @param double a
    @param double b
    @param double c
    @param double d
    @param double e
    @param double f
    @param double g
    @param double h
    @param double i
    @return double The determinant value.
    */
    double det(double a, double b, double c, double d, double e, double f, double g, double h, double i);

    /*
    Calculates the parameters for the point where a face and a ray intersect.
    @param const std::shared_ptr<dcel::Face> face The pointer to a face in a DCEL.
    @param const toolkit::Ray & ray The reference to a ray.
    @param double & t The parameter value for the intersection point along the ray.
    @param double threshold The threshold for values close to zero.
    @return bool Indicates whether there is an intersection point between the face and the ray.
    */
    //bool faceRayIntersection(const std::shared_ptr<dcel::Face> face, const toolkit::Ray & ray, double & t, double threshold = 1e-8);

    /*
    Finds a path between the points in a graph represented by a directed adjacency matrix.
    @param Eigen::MatrixXd & adj The reference to the adjacency matrix.
    @param std::vector<bool> & visited The reference to the vector indicating which points have
    been visited.
    @param size_t i The start point index of the path.
    @param std::vector<size_t> & path The reference to the vector with the sequence of indices that
    form the path.
    @return bool Indicates if there is a path.
    */
    bool findPath(Eigen::MatrixXd & adj, std::vector<bool> & visited, size_t i, std::vector<size_t> & path);

    /*
    Fix the zero of a given value.
    @param double value The reference to the variable.
    @param double threshold The threshold for zero values.
    */
    void fixZero(double & val, double threshold = 1e-8);

    /*
    Fix the zero values of a given 3D vector.
    @param Eigen::Vector3d & V The reference to the vector.
    @param double threshold The threshold for zero values.
    */
    void fixZeros(Eigen::Vector3d & V, double threshold = 1e-8);

    /*
    Fix the zero values of a given vector with 3D vectors.
    @param std::vector<Eigen::Vector3d> & V The reference to the vector with the 3D vectors.
    @param double threshold The threshold value.
    */
    void fixZeros(std::vector<Eigen::Vector3d> & V, double threshold = 1e-8);

    /*
    Calculates the intersection points between Planes.
    @param const std::vector<toolkit::Plane> & Planes The reference to the vector with the Planes.
    @param std::vector<Eigen::Vector3d> & points The reference to the vector to store the
    intersection points between the Planes.
    @return bool Indicates if there are intersection points between the Planes. If at least one
    point is missed the result is false.
    */
    //bool getIntersectionPointsFromPlanes(const std::vector<toolkit::Plane> & Planes, std::vector<Eigen::Vector3d> & points);

    /*
    Populates a vector with unique points.
    @param std::vector<Eigen::Vector3d> & points The reference to the vector with the points.
    @param std::vector<Eigen::Vector3d> & unique The reference to the vector where unique points
    will be stored.
    @param double threshold The threshold for values close to zero.
    @return bool Indicates if there is at least one unique point.
    */
    bool getUniquePoints(std::vector<Eigen::Vector3d> & points, std::vector<Eigen::Vector3d> & unique, double threshold = 1e-8);

    /*
    @param const VF & vf1
    @param size_t face1
    @param const VF & vf2
    @param size_t face2
    @param std::vector<Eigen::Vector3d> & points
    @param double threshold
    @return bool
    */
    //bool intersectCoplanarFaces(const VF & vf1, size_t face1, const VF & vf2, size_t face2, std::vector<Eigen::Vector3d> & points, double threshold = 1e-8);

    /*
    Finds the intersection points of a DCEL face F1 with another DCEL face F2. Note that this
    approach does not calculate the intersection points of face F2 with face F1 since they are not
    necessarily the same.
    @param std::shared_ptr<dcel::Face> F1 The pointer to a face.
    @param std::shared_ptr<dcel::Face> F2 The pointer to a face.
    @param std::vector<Eigen::Vector3d> & points The reference to the vector where the intersection
    points will be stored.
    @param double threshold The threshold for values close to zero.
    @return double Indicates if there is at least one intersection point between the faces.
    */
    //bool intersectCoplanarFaces(std::shared_ptr<dcel::Face> F1, std::shared_ptr<dcel::Face> F2, std::vector<Eigen::Vector3d> & points, double threshold = 1e-8);

    /*
    Indicates if
    Calculate N * ((P2 - P1) x (Q - P1)), where * is dot product and x is cross product. This value
    indicates if a vector from P1 to P2 goes in CCW direction with respect to the normal vector N 
    and point Q (the resultant vector and N must point at the same half space)
    @param const Eigen::Vector3d & P1
    @param const Eigen::Vector3d & P2
    @param const toolkit::Plane & plane
    @param const double threshold
    @return int
    */
    //int isLeftTurn(const Eigen::Vector3d & P1, const Eigen::Vector3d & P2, const toolkit::Plane & plane, const double threshold = 1e-8);

    /*
    Checks if a point lies within the triangle defined by its vertices.
    From: https://stackoverflow.com/questions/995445/determine-if-a-3d-point-is-within-a-triangle
    @param const Eigen::Vector3d & P The reference to a point.
    @param const Eigen::Vector3d & v0 The reference to a vertex of the triangle.
    @param const Eigen::Vector3d & v1 The reference to a vertex of the triangle.
    @param const Eigen::Vector3d & v2 The reference to a vertex of the triangle.
    @param double threshold The threshold for values close to zero.
    @return bool Indicates whether the point lies within the triangle.
    */
    bool isPointInTriangle(const Eigen::Vector3d & P, const Eigen::Vector3d & v0, const Eigen::Vector3d & v1, const Eigen::Vector3d & v2, double threshold = 1e-8);

    /*
    Calculates the parameters for the point where a plane and a line segment intersect.
    @param const toolkit::Plane & plane The reference to a plane
    @param const toolkit::LineSegment & linesegment The reference to a line segment.
    @param double & t The parameter value for the intersection point along the line segment.
    @param double threshold The threshold for values close to zero.
    @return bool Indicates whether there is an intersection point between the face and the line
    segment.
    */
    //bool planeLineSegmentIntersection(const toolkit::Plane & plane, const toolkit::LineSegment & linesegment, double & t, double threshold = 1e-8);

    /*
    Calculates the parameter for the point where a toolkit::Plane and a ray intersect.
    @param const toolkit::Plane & toolkit::Plane The reference to a toolkit::Plane.
    @param const Ray & ray The reference to a ray.
    @param double & t The parameter value for the intersection point along the ray.
    @param double threshold The threshold for values close to zero.
    @return bool Indicates whether there is an intersection point between the toolkit::Plane and the ray.
    */
    //bool planeRayIntersection(const toolkit::Plane & Plane, const toolkit::Ray & ray, double & t, double threshold = 1e-8);
    
    /// <summary>
    /// 
    /// </summary>
    /// <param name="s"></param>
    void print(const char * s);
    
    /// <summary>
    /// 
    /// </summary>
    /// <param name="s"></param>
    /// <param name="v"></param>
    void print(const char * s, double v);
    
    /// <summary>
    /// 
    /// </summary>
    /// <param name="s"></param>
    /// <param name="V"></param>
    void print(const char * s, const Eigen::Vector3d & V);
    
    /// <summary>
    /// 
    /// </summary>
    /// <param name="s"></param>
    /// <param name="i"></param>
    /// <param name="V"></param>
    void print(const char * s, size_t i, const Eigen::Vector3d & V);
    
    /// <summary>
    /// 
    /// </summary>
    /// <param name="s1"></param>
    /// <param name="i1"></param>
    /// <param name="s2"></param>
    /// <param name="i2"></param>
    void print(const char * s1, size_t i1, const char * s2, size_t i2);
    
    /// <summary>
    /// 
    /// </summary>
    /// <param name="s1"></param>
    /// <param name="i1"></param>
    /// <param name="s2"></param>
    /// <param name="i2"></param>
    /// <param name="V"></param>
    void print(const char * s1, size_t i1, const char * s2, size_t i2, const Eigen::Vector3d & V);
    
    /// <summary>
    /// 
    /// </summary>
    /// <param name="V"></param>
    void print(const std::vector<Eigen::Vector3d> & V);

    /*
    @param const Eigen::Vector3d & P
    @param const Eigen::Vector3d & Q
    @param const Eigen::Vector3d & A
    @param const Eigen::Vector3d & B
    @param double threshold The threshold for values close to zero.
    @return bool
    */
    bool sameSide(const Eigen::Vector3d & P, const Eigen::Vector3d & Q, const Eigen::Vector3d & A, const Eigen::Vector3d & B, double threshold = 1e-8);

    /*
    Sorts a vector of points in CCW order with respect of a point P and a normal vector N.
    @param std::vector<Eigen::Vector3d> & points The reference to the vector with the points to be
    sorted.
    @param Eigen::Vector3d & N The reference to the normal vector.
    @param Eigen::Vector3d & Q The reference to the point.
    @param std::vector<Eigen::Vector3d> & sorted The reference to the vector where the points will
    be stored in CCW.
    @param double threshold The threshold for values close to zero.
    @return bool
    */
    //bool sortPointsCCW(std::vector<Eigen::Vector3d> & points, const Eigen::Vector3d & N, const Eigen::Vector3d & Q, std::vector<Eigen::Vector3d> & sorted, double threshold = 1e-8);

    /*
    Sorts a vector of points in CCW order with respect of a Plane.
    @param std::vector<Eigen::Vector3d> & points The reference to the vector with the points to be
    sorted.
    @param toolkit::Plane & plane The reference to a Plane.
    @param std::vector<Eigen::Vector3d> & sorted The reference to the vector where the points will
    be stored in CCW.
    @param double threshold
    @return bool
    */
    //bool sortPointsCCW(std::vector<Eigen::Vector3d> & points, toolkit::Plane & plane, std::vector<Eigen::Vector3d> & sorted, double threshold = 1e-8);

    /*
    Returns the intersection point between three planes. Here we assume the planes do intersect at
    one point, we don't care about the other possible cases. Intersection point is found by solving
    the system of equations:
    A.a * x + A.b * y + A.c * z - A.d = 0
    B.a * x + B.b * y + B.c * z - B.d = 0
    C.a * x + C.b * y + C.c * z - C.d = 0
    @param toolkit::Plane & A The reference to a plane.
    @param toolkit::Plane & B The reference to a plane.
    @param toolkit::Plane & C The reference to a plane.
    @param Eigen::Vector3d & P The reference to the intersection point.
    @return bool Indicates if there is an intersection point between the three Planes.
    */
    //bool threePlanesIntersection(const toolkit::Plane & A, const toolkit::Plane & B, const toolkit::Plane & C, Eigen::Vector3d & P);

    /*
    Calculates a function using the values in a vector.
    @param std::vector<double> & values The vector with the values for the function.
    @param const std::string & function The name of the function.
    @return double The value of the function over the values in the vector.
    */
    double vectorFunction(std::vector<double> & values, const std::string & function);

    /*
    Writes the content of a 3D vector.
    @param const Eigen::Vector3d & V The reference to the vector.
    @param bool newline Indicates whether to write a new line character after the vector or not.
    */
    void write(const Eigen::Vector3d & V, bool newline = false);
    
    /// <summary>
    /// 
    /// </summary>
    /// <param name="points"></param>
    void write(const std::list<Eigen::Vector3d> & points);
    
    /// <summary>
    /// 
    /// </summary>
    /// <param name="L"></param>
    void write(const std::list<std::vector<size_t>> & L);
    
    /// <summary>
    /// 
    /// </summary>
    /// <param name="V"></param>
    void write(const std::vector<Eigen::Vector3d> & V);
    
    /// <summary>
    /// 
    /// </summary>
    /// <param name="V"></param>
    //void write(const std::vector<dcel::DCEL> & V);
    
    /// <summary>
    /// 
    /// </summary>
    /// <param name="V"></param>
    void write(const std::vector<double> & V);

    /*
    @param const std::vector<dcel::DCEL> & blocks The reference to the vector with the geometry of
    the blocks.
    @param const std::vector<dcel::DCEL> & interfaces The reference to the vector with the geometry
    of the interface polygons.
    */
    //void writeGeogebraJS(const std::vector<dcel::DCEL> & blocks, const std::vector<dcel::DCEL> & interfaces, const std::string filename = "geogebra.js", bool writeReference = false, const Eigen::Vector3d & O = Eigen::Vector3d(0, 0, 0), const Eigen::Vector3d & X = Eigen::Vector3d(1, 0, 0), const Eigen::Vector3d & Y = Eigen::Vector3d(0, 1, 0), const Eigen::Vector3d & Z = Eigen::Vector3d(0, 0, 1));

}

#endif
