#ifndef _TOOLKIT_ALGORITHMS_H_
#define _TOOLKIT_ALGORITHMS_H_

#pragma once

#include <tiger/toolkit/ray.h>
#include <tiger/toolkit/plane.h>
#include <tiger/toolkit/linesegment.h>

#include <vector>

namespace algorithms 
{
    //
    // Calculates the intersection points between Planes.
    // @param const std::vector<toolkit::Plane> & Planes The reference to the 
    // vector with the Planes.
    // @param std::vector<Eigen::Vector3d> & points The reference to the vector
    // to store the intersection points between the Planes.
    // @return bool Indicates if there are intersection points between the 
    // Planes. If at least one point is missed the result is false.
    //
    bool getIntersectionPointsFromPlanes(
        const std::vector<toolkit::Plane>& Planes,
        std::vector<Eigen::Vector3d>& points);

    //
    // Calculates the parameters for the point where a plane and a line segment
    // intersect.
    // @param const toolkit::Plane & plane The reference to a plane
    // @param const toolkit::LineSegment & linesegment The reference to a line 
    // segment.
    // @param double & t The parameter value for the intersection point along 
    // the line segment.
    // @param double threshold The threshold for values close to zero.
    // @return bool Indicates whether there is an intersection point between 
    // the face and the line segment.
    //
    bool planeLineSegmentIntersection(
        const toolkit::Plane & plane,
        const toolkit::LineSegment & linesegment,
        double& t,
        double threshold = 1e-8);

    //
    // Calculate N * ((P2 - P1) x (Q - P1)), where * is dot product and x is 
    // cross product. This value indicates if a vector from P1 to P2 goes in 
    // CCW direction with respect to the normal vector N and point Q (the 
    // resultant vector and N must point at the same half space)
    // @param const Eigen::Vector3d & P1
    // @param const Eigen::Vector3d & P2
    // @param const toolkit::Plane & plane
    // @param const double threshold
    // @return int
    //
    int isLeftTurn(
        const Eigen::Vector3d& P1,
        const Eigen::Vector3d& P2,
        const toolkit::Plane& plane,
        const double threshold = 1e-8);

    //
    // Calculates the parameter for the point where a toolkit::Plane and a ray 
    // intersect.
    // @param const toolkit::Plane & toolkit::Plane The reference to a 
    // toolkit::Plane.
    // @param const Ray & ray The reference to a ray.
    // @param double & t The parameter value for the intersection point along 
    // the ray.
    // @param double threshold The threshold for values close to zero.
    // @return bool Indicates whether there is an intersection point between 
    // the toolkit::Plane and the ray.
    //
    bool planeRayIntersection(
        const toolkit::Plane & Plane,
        const toolkit::Ray & ray,
        double& t,
        double threshold = 1e-8);

    //
    // Sorts a vector of points in CCW order with respect of a Plane.
    // @param std::vector<Eigen::Vector3d> & points The reference to the vector
    // with the points to be sorted.
    // @param toolkit::Plane & plane The reference to a Plane.
    // @param std::vector<Eigen::Vector3d> & sorted The reference to the vector
    // where the points will be stored in CCW.
    // @param double threshold
    // @return bool
    //
    bool sortPointsCCW(
        std::vector<Eigen::Vector3d> & points,
        toolkit::Plane & plane,
        std::vector<Eigen::Vector3d> & sorted,
        double threshold = 1e-8);

    //
    bool sortPointsCCW(
        std::vector<Eigen::Vector3d>& points,
        const Eigen::Vector3d& N,
        const Eigen::Vector3d& Q,
        std::vector<Eigen::Vector3d>& sorted,
        double threshold = 1e-8);

    //
    // Returns the intersection point between three planes. Here we assume the 
    // planes do intersect at one point, we don't care about the other possible
    // cases. Intersection point is found by solving the system of equations:
    // A.a * x + A.b * y + A.c * z - A.d = 0
    // B.a * x + B.b * y + B.c * z - B.d = 0
    // C.a * x + C.b * y + C.c * z - C.d = 0
    // @param toolkit::Plane & A The reference to a plane.
    // @param toolkit::Plane & B The reference to a plane.
    // @param toolkit::Plane & C The reference to a plane.
    // @param Eigen::Vector3d & P The reference to the intersection point.
    // @return bool Indicates if there is an intersection point between the 
    // three Planes.
    //
    bool threePlanesIntersection(
        const toolkit::Plane& A,
        const toolkit::Plane& B,
        const toolkit::Plane& C,
        Eigen::Vector3d& P);
}

#endif
