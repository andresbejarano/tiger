#ifndef _UTILS_H_
#define _UTILS_H_

#pragma once

#include <Eigen/core>
#include <list>
#include <vector>

namespace utils
{
    // The constant PI
    const double PI = acos(-1.0);
    const double TWO_PI = PI * 2.0;

    // 
    // Rotates vector V around an axis vector K and angle of rotation a. This method is also named as
    // Rodrigues' rotation formula: (V * cos(a)) + ((K x N) * sin(a)) + (K * (K . V) * (1 - cos(a)))
    // Its description is in  https:// en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    // @param const Eigen::Vector3d & V The reference to the vector to be rotated.
    // @param const Eigen::Vector3d & K The reference to the rotation axis vector.
    // @param double angle The rotation angle (in radians).
    // @return Eigen::Vector3d The rotated vector.
    // 
    Eigen::Vector3d axisAngleRotation(const Eigen::Vector3d & V, const Eigen::Vector3d & K, double angle);

    // 
    // Calculates the centroid (arithmetic mean) of a list of points.
    // @param const std::list<Eigen::Vector3d> & points The reference to a list with the points.
    // @param bool fixZeros Indicates whether to fix the values close to zero.
    // @param double threshold The threshold for values close to zero.
    // @return Eigen::Vector3d The centroid of the points.
    // 
    Eigen::Vector3d centroid(const std::list<Eigen::Vector3d> & points, bool fixZeros = false, double threshold = 1e-8);

    // 
    // Calculates the centroid (arithmetic mean) of a vector of points.
    // @param const std::vector<Eigen::Vector3d> & points The reference to a vector with the points.
    // @param bool fixZeros Indicates whether to fix the values close to zero.
    // @param double threshold The threshold for values close to zero.
    // @return Eigen::Vector3d The centroid of the points.
    // 
    Eigen::Vector3d centroid(const std::vector<Eigen::Vector3d> & points, bool fixZeros = false, double threshold = 1e-8);

    // 
    // Returns the determinant of the 3x3 matrix specified as:
    // |a, b, c|
    // |d, e, f|
    // |g, h, i|
    // @param double a
    // @param double b
    // @param double c
    // @param double d
    // @param double e
    // @param double f
    // @param double g
    // @param double h
    // @param double i
    // @return double The determinant value.
    // 
    double det(double a, double b, double c, double d, double e, double f, double g, double h, double i);

    // 
    // Finds a path between the points in a graph represented by a directed adjacency matrix.
    // @param Eigen::MatrixXd & adj The reference to the adjacency matrix.
    // @param std::vector<bool> & visited The reference to the vector indicating which points have
    // been visited.
    // @param size_t i The start point index of the path.
    // @param std::vector<size_t> & path The reference to the vector with the sequence of indices that
    // form the path.
    // @return bool Indicates if there is a path.
    // 
    bool findPath(Eigen::MatrixXd & adj, std::vector<bool> & visited, size_t i, std::vector<size_t> & path);

    // 
    // Fix the zero of a given value.
    // @param double value The reference to the variable.
    // @param double threshold The threshold for zero values.
    // 
    void fixZero(double & val, double threshold = 1e-8);

    // 
    // Fix the zero values of a given 3D vector.
    // @param Eigen::Vector3d & V The reference to the vector.
    // @param double threshold The threshold for zero values.
    // 
    void fixZeros(Eigen::Vector3d & V, double threshold = 1e-8);

    // 
    // Fix the zero values of a given vector with 3D vectors.
    // @param std::vector<Eigen::Vector3d> & V The reference to the vector with the 3D vectors.
    // @param double threshold The threshold value.
    // 
    void fixZeros(std::vector<Eigen::Vector3d> & V, double threshold = 1e-8);

    // 
    // Populates a vector with unique points.
    // @param std::vector<Eigen::Vector3d> & points The reference to the vector with the points.
    // @param std::vector<Eigen::Vector3d> & unique The reference to the vector where unique points
    // will be stored.
    // @param double threshold The threshold for values close to zero.
    // @return bool Indicates if there is at least one unique point.
    // 
    bool getUniquePoints(std::vector<Eigen::Vector3d> & points, std::vector<Eigen::Vector3d> & unique, double threshold = 1e-8);

    // 
    // Checks if a point lies within the triangle defined by its vertices.
    // From: https:// stackoverflow.com/questions/995445/determine-if-a-3d-point-is-within-a-triangle
    // @param const Eigen::Vector3d & P The reference to a point.
    // @param const Eigen::Vector3d & v0 The reference to a vertex of the triangle.
    // @param const Eigen::Vector3d & v1 The reference to a vertex of the triangle.
    // @param const Eigen::Vector3d & v2 The reference to a vertex of the triangle.
    // @param double threshold The threshold for values close to zero.
    // @return bool Indicates whether the point lies within the triangle.
    // 
    bool isPointInTriangle(const Eigen::Vector3d & P, const Eigen::Vector3d & v0, const Eigen::Vector3d & v1, const Eigen::Vector3d & v2, double threshold = 1e-8);
    
    //
    void print(const char * s);
    
    //
    void print(const char * s, double v);
    
    //
    void print(const char * s, const Eigen::Vector3d & V);
    
    //
    void print(const char * s, size_t i, const Eigen::Vector3d & V);
    
    //
    void print(const char * s1, size_t i1, const char * s2, size_t i2);
    
    //
    void print(const char * s1, size_t i1, const char * s2, size_t i2, const Eigen::Vector3d & V);
    
    //
    void print(const std::vector<Eigen::Vector3d> & V);

    // 
    // @param const Eigen::Vector3d & P
    // @param const Eigen::Vector3d & Q
    // @param const Eigen::Vector3d & A
    // @param const Eigen::Vector3d & B
    // @param double threshold The threshold for values close to zero.
    // @return bool
    // 
    bool sameSide(const Eigen::Vector3d & P, const Eigen::Vector3d & Q, const Eigen::Vector3d & A, const Eigen::Vector3d & B, double threshold = 1e-8);

    // 
    // Calculates a function using the values in a vector.
    // @param std::vector<double> & values The vector with the values for the function.
    // @param const std::string & function The name of the function.
    // @return double The value of the function over the values in the vector.
    // 
    double vectorFunction(std::vector<double> & values, const std::string & function);

    // 
    // Writes the content of a 3D vector.
    // @param const Eigen::Vector3d & V The reference to the vector.
    // @param bool newline Indicates whether to write a new line character after the vector or not.
    // 
    void write(const Eigen::Vector3d & V, bool newline = false);
    
    //
    void write(const std::list<Eigen::Vector3d> & points);
    
    //
    void write(const std::list<std::vector<size_t>> & L);
    
    //
    void write(const std::vector<Eigen::Vector3d> & V);
    
    //
    void write(const std::vector<double> & V);

}

#endif
