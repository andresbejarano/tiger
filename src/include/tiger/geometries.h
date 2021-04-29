#ifndef _GEOMETRIES_H_
#define _GEOMETRIES_H_

#pragma once

#include <tiger/ds/vf.h>
#include <math.h>
#include <functional>

namespace geometries
{
    /*
    @param double radius
    @param double length
    @param size_t rs
    @param size_t ls
    @param int direction
    @param const Eigen::Vector3d & U
    @param const Eigen::Vector3d & V
    @return VF The vertices and faces of the geometry.
    */
    /*VF BendedSquare(
        double radius, 
        double length, 
        size_t rs, 
        size_t ls, 
        int direction, 
        const Eigen::Vector3d & U = Eigen::Vector3d(1, 0, 0), 
        const Eigen::Vector3d & V = Eigen::Vector3d(0, 1, 0));*/

    /*
    Get the vertices and faces of a bended square (centered at the origin) along
    the Z axis.
    @param double radius The bending radius. Must be non-negative.
    @param double length The length of the square. Must be non-negative.
    @param size_t rs The number of radial segments. Must be non-negative.
    @param size_t ls The number of length segments. Must be greater or equal to 3.
    @param int bending Indicates whether the bending goes upward (1) or downward (-1).
    @return VF The vertices and faces of the geometry.
    */
    VF BendedSquare(double radius, double length, size_t rs, size_t ls, int bending);

    /*
    Get the vertex coordinates and the vertex indices of a cone
    @param const Eigen::Vector3d & A The reference to the apex (Singularity) point.
    @param const Eigen::Vector3d & B The reference to the point at the end of the cone.
    @param double radius The radius of the cone.
    @param size_t ls The number of length sections along vector AB.
    @param size_t rs The number of radial sections.
    @param bool open Indicates whether to open or close the cap of the cone.
    @return VF The vertices and faces of the geometry.
    */
    VF Cone(
        const Eigen::Vector3d & A, 
        const Eigen::Vector3d & B, 
        double radius, 
        size_t ls, 
        size_t rs, 
        bool open = true);

    /// <summary>
    /// 
    /// </summary>
    /// <param name="a"></param>
    /// <param name="b"></param>
    /// <param name="c"></param>
    /// <param name="d"></param>
    /// <returns></returns>
    VF Cyclide(double a, double b, double c, double d, size_t Rs, size_t rs);

    /*
    Get the vertices and faces of a cylinder (centered at a given point) along 
    the V axis vector. The parametric equation of the cylinder is defined as:
    Cylinder(C, K, r, l) = [(A + t * (B - A)) + (r * R * cos(a)) + (r * S * sin(a))]
    where C is the center point of the cylinder, K is the axis vector of the cylinder, r is the
    radius of the cylinder, l is the length of the cylinder, A = C - (l / 2) * ||K||,
    B = C + (l / 2) * ||K||, P is a random point but along line segment AB, R = ||K X (P - C)||,
    S = ||K x R||, t in [0, 1] and a in [0, 2 pi)
    Reference: http://paulbourke.net/geometry/circlesphere/ in Creating a plane/disk perpendicular
    to a line segment.
    @param const Eigen::Vector3d & C The reference to the center point of the cylinder.
    @param const Eigen::Vector3d & K The reference to the axis vector of the cylinder.
    @param double radius The radius of the cylinder. Must be non-negative.
    @param double length The length of the cylinder. Must be non-negative.
    @param size_t rs The number of radial segments. Must be non-negative.
    @param size_t ls The number of length segments. Must be greater or equal to 3.
    @param bool open
    @return VF The vertices and faces of the geometry.
    */
    VF Cylinder(
        const Eigen::Vector3d & C, 
        const Eigen::Vector3d & K, 
        double radius, 
        double length, 
        size_t rs, 
        size_t ls, 
        bool open = true);

    /*
    Get the vertices and faces of a cylinder between points A and B. 
    @param const Eigen::Vector3d & A The reference to one of the end points of the cylinder.
    @param const Eigen::Vector3d & B The reference to one of the end points of the cylinder.
    @param double radius The radius of the cylinder. Must be non-negative.
    @param size_t rs The number of radial segments. Must be non-negative.
    @param size_t ls The number of length segments. Must be greater or equal to 3.
    @param bool open
    @return VF The vertices and faces of the geometry.
    */
    VF Cylinder(
        const Eigen::Vector3d & A,
        const Eigen::Vector3d & B,
        double radius,
        size_t rs,
        size_t ls, 
        bool open = true);

    /*
    Get the vertices and faces of a regular dodecahedron or cube (centered at 
    the origin) as described in https://en.wikipedia.org/wiki/Dodecahedron
    @param double radius The radius of the circumscribed sphere.
    @return VF The vertices and faces of the geometry.
    */
    VF Dodecahedron(double radius = 1.0);

    /*
    @param double a
    @param double b
    @param double width
    @param double height
    @param size_t ws
    @param size_t hs
    @param const Eigen::Vector3d& X = Eigen::Vector3d(1, 0, 0)
    @param const Eigen::Vector3d& Y = Eigen::Vector3d(0, 1, 0)
    @return VF The vertices and faces of the geometry.
    */
    VF EllipticParaboloid(
        double a, 
        double b, 
        double width, 
        double height, 
        size_t ws, 
        size_t hs, 
        const Eigen::Vector3d& X = Eigen::Vector3d(1, 0, 0), 
        const Eigen::Vector3d& Y = Eigen::Vector3d(0, 1, 0));

    /*
    Generates an equilateral triangle with a given length sides. Triangle centroid is at the 
    origin.
    @param double length The length of the triangle sides.
    @param const Eigen::Vector3d & U The reference to the vector with the U direction.
    @param const Eigen::Vector3d & V The reference to the vector with the V direction.
    @return VF The vertices and faces of the geometry.
    */
    VF EquilateralTriangle(
        double length,
        const Eigen::Vector3d & U = Eigen::Vector3d(1, 0, 0),
        const Eigen::Vector3d & V = Eigen::Vector3d(0, 1, 0));

    /*
    Get the vertices and faces of a regular hexahedron or cube (centered at the 
    origin) as described in https://en.wikipedia.org/wiki/Cube
    @param double radius The radius of the circumscribed sphere.
    @return VF The vertices and faces of the geometry.
    */
    VF Hexahedron(double radius);

    /*
    @param double width
    @param double height
    @param size_t ws
    @param size_t hs
    @param const std::function<double(double& x, double& y)>& Z
    @param const Eigen::Vector3d& X = Eigen::Vector3d(1, 0, 0)
    @param const Eigen::Vector3d& Y = Eigen::Vector3d(0, 1, 0)
    @return VF The vertices and faces of the geometry.
    */
    VF Grid(
        double width,
        double height,
        size_t ws,
        size_t hs,
        const std::function<double(double& x, double& y)>& Z,
        const Eigen::Vector3d& X = Eigen::Vector3d(1, 0, 0),
        const Eigen::Vector3d& Y = Eigen::Vector3d(0, 1, 0));

    /*
    @param double a
    @param double b
    @param double width
    @param double height
    @param size_t ws
    @param size_t hs
    @param const Eigen::Vector3d& X = Eigen::Vector3d(1, 0, 0)
    @param const Eigen::Vector3d& Y = Eigen::Vector3d(0, 1, 0)
    @return VF The vertices and faces of the geometry.
    */
    VF HyperbolicParaboloid(
        double a,
        double b,
        double width,
        double height,
        size_t ws,
        size_t hs,
        const Eigen::Vector3d& X = Eigen::Vector3d(1, 0, 0),
        const Eigen::Vector3d& Y = Eigen::Vector3d(0, 1, 0));

    /*
    Get the vertices and faces of a regular icosahedron (centered at the origin)
    as described in https://en.wikipedia.org/wiki/Icosahedron
    @param double radius The radius of the circumscribed sphere.
    @return VF The vertices and faces of the geometry.
    */
    VF Icosahedron(double radius);

    /*
    Get the vertices and faces of a regular octahedron (centered at the origin) 
    as described in https://en.wikipedia.org/wiki/Octahedron
    @param double radius The radius of the circumscribed sphere.
    @return VF The vertices and faces of the geometry.
    */
    VF Octahedron(double radius);

    /*
    Get the vertices and faces of a regular polygon in the plane spanned by two UV vectors.
    @param size_t sides The number of sides of the polygon.
    @param double length The length of the sides of the polygon.
    @param const Eigen::Vector3d & U The reference to the vector with the U direction.
    @param const Eigen::Vector3d & V The reference to the vector with the V direction.
    @return VF The vertices and faces of the geometry.
    */
    VF Polygon(
        size_t sides,
        double length,
        const Eigen::Vector3d & U = Eigen::Vector3d(1, 0, 0),
        const Eigen::Vector3d & V = Eigen::Vector3d(0, 1, 0));

    /*
    Returns a random number between 0.0 and 1.0.
    @return double A random number between 0.0 and 1.0.
    */
    double Random01();

    /*
    Generates the vertices and face of a rectangle in the plane spanned by two UV vectors.
    @param double width The width of the rectangle.
    @param double height The height of the rectangle.
    @param const Eigen::Vector3d & W The reference to the vector with the width direction.
    @param const Eigen::Vector3d & H The reference to the vector with the height direction.
    @return VF The vertices and faces of the geometry.
    */
    VF Rectangle(
        double width, 
        double height, 
        const Eigen::Vector3d & W = Eigen::Vector3d(1, 0, 0),
        const Eigen::Vector3d & H = Eigen::Vector3d(0, 1, 0));

    /*
    @param double width
    @param double height
    @param size_t ws
    @param size_t hs
    @param const Eigen::Vector3d & X
    @param const Eigen::Vector3d & Y
    @return VF The vertices and faces of the geometry.
    */
    VF Saddle(
        double width, 
        double height, 
        size_t ws, 
        size_t hs, 
        const Eigen::Vector3d& X = Eigen::Vector3d(1, 0, 0), 
        const Eigen::Vector3d& Y = Eigen::Vector3d(0, 1, 0));

    /*
    Generates the vertices and face of a square.
    @param double length The length of the square sides.
    @param const Eigen::Vector3d & U The reference to the vector with the U direction.
    @param const Eigen::Vector3d & V The reference to the vector with the V direction.
    @return VF The vertices and faces of the geometry.
    */
    VF Square(
        double length, 
        const Eigen::Vector3d & U = Eigen::Vector3d(1, 0, 0), 
        const Eigen::Vector3d & V = Eigen::Vector3d(0, 1, 0));

    /*
    Generates a grid of qudrilaterals in the plane spanned by two vectors.
    @param double width The width of the grid.
    @param double height The height of the grid.
    @param size_t ws The number of width segments.
    @param size_t hs The number of height segments.
    @param const Eigen::Vector3d& X The reference to the vector with the width direction.
    @param const Eigen::Vector3d& Y The reference to the vector with the height direction.
    @return VF The vertices and faces of the geometry.
    */
    VF SquareGrid(
        double width,
        double height,
        size_t ws,
        size_t hs,
        const Eigen::Vector3d& X = Eigen::Vector3d(1.0, 0.0, 0.0),
        const Eigen::Vector3d& Y = Eigen::Vector3d(0.0, 1.0, 0.0));

    /*
    Get the vertices and faces of a regular tetrahedron (centered at the origin)
    as described in https://en.wikipedia.org/wiki/Tetrahedron
    @param double radius The radius of the circumscribed sphere.
    @return VF The vertices and faces of the geometry.
    */
    VF Tetrahedron(double radius);

    /*
    Get the vertices and faces of a torus centered at a point C and along a V 
    axis vector. Torus(C, K, r1, r2) =
        C + (r1 * (R * cos(u) + S * sin(u))) + (r2 * ((R * cos(u) + S * sin(u)) * cos(v) + D * sin(v)))
    where C is the center point, K is the axis vector, r1 is the major radius, r2 is the minor
    radius, P a random point but collinear with C + t * K, D = ||K||, R = ||D X (P - C)||, 
    S = ||D x R||, u in [0, 2 pi) and v in [0, 2 pi]
    @param const Eigen::Vector3d C The center point of the torus.
    @param const Eigen::Vector3d K The axis vector of the torus.
    @param double r1 The major radius of the torus.
    @param double r2 The minor radius of the torus.
    @param int r1s The number of segments along the major section.
    @param int r2s The number of segments along the minor section.
    @return VF The vertices and faces of the geometry.
    */
    VF Torus(
        double R, 
        double r, 
        size_t Rs, 
        size_t rs, 
        const Eigen::Vector3d& C = Eigen::Vector3d::Zero(), 
        const Eigen::Vector3d& X = Eigen::Vector3d(1.0, 0.0, 0.0),
        const Eigen::Vector3d& Y = Eigen::Vector3d(0.0, 1.0, 0.0));

    /*
    Get the vertices and faces of a truncated cone centered at C and along the 
    axis vector K. The parametric equation of the truncated cone is defined as:
    TruncatedCone(C, K, br, tr, l) =
        [(A + t * (B - A)) + ((br + t * (tr - br)) * R * cos(a)) + ((br + t * (tr - br)) * S * sin(a))]
    where C is the center point, K is the axis vector, br is the bottom radius, tr is the top
    radius, l is the length, A = C - (l / 2) * ||K||, B = C + (l / 2) * ||K||, P is a random point
    but along line segment AB, R = ||K X (P - C)||, S = ||K x R||, t in [0, 1] and a in [0, 2 pi)
    @param const Eigen::Vector3d C The center point of the truncated cone.
    @param const Eigen::Vector3d K The axis vector of the truncated cone.
    @param double br The bottom radius of the truncated cone. Must be non-negative.
    @param double tr The top radius of the truncated cone. Must be non-negative.
    @param double length The length of the truncated cone. Must be non-negative.
    @param int rs The number of radial segments. Must be non-negative.
    @param int ls The number of length segments. Must be greater or equal to 3.
    @return VF The vertices and faces of the geometry.
    */
    VF TruncatedCone(
        const Eigen::Vector3d & C, 
        const Eigen::Vector3d & K, 
        double br, 
        double tr, 
        double length, 
        size_t rs, 
        size_t ls);
}

#endif
