#ifndef _DS_ALGORITHMS_H_
#define _DS_ALGORITHMS_H_

#pragma once

#include <tiger/ds/vf.h>
#include <tiger/ds/dcel.h>
#include <tiger/toolkit/ray.h>

namespace algorithms 
{
    /*
    Checks if all face vertices have the ATTRIB_COPLANAR attribute set with value true.
    @param const std::shared_ptr<dcel::Face> face The pointer to a face.
    @return bool Indicates if all face vertices have the ATTRIB_COPLANAR attribute set to true.
    */
    bool areFaceVerticesCoplanar(const std::shared_ptr<dcel::Face> face);

    /*
    Calculates the parameters for the point where a face and a ray intersect.
    @param const std::shared_ptr<dcel::Face> face The pointer to a face in a DCEL.
    @param const toolkit::Ray & ray The reference to a ray.
    @param double & t The parameter value for the intersection point along the ray.
    @param double threshold The threshold for values close to zero.
    @return bool Indicates whether there is an intersection point between the face and the ray.
    */
    bool faceRayIntersection(
        const std::shared_ptr<dcel::Face> face,
        const toolkit::Ray & ray,
        double& t,
        double threshold = 1e-8);

    /*
    @param const VF & vf1
    @param size_t face1
    @param const VF & vf2
    @param size_t face2
    @param std::vector<Eigen::Vector3d> & points
    @param double threshold
    @return bool
    */
    bool intersectCoplanarFaces(
        const VF & vf1,
        size_t face1,
        const VF & vf2,
        size_t face2,
        std::vector<Eigen::Vector3d> & points,
        double threshold = 1e-8);

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
    bool intersectCoplanarFaces(
        std::shared_ptr<dcel::Face> F1,
        std::shared_ptr<dcel::Face> F2,
        std::vector<Eigen::Vector3d>& points,
        double threshold = 1e-8);

    /// <summary>
    /// 
    /// </summary>
    /// <param name="V"></param>
    void write(const std::vector<dcel::DCEL> & V);

    /*
    @param const std::vector<dcel::DCEL> & blocks The reference to the vector with the geometry of
    the blocks.
    @param const std::vector<dcel::DCEL> & interfaces The reference to the vector with the geometry
    of the interface polygons.
    */
    void writeGeogebraJS(
        const std::vector<dcel::DCEL>& blocks,
        const std::vector<dcel::DCEL>& interfaces,
        const std::string filename = "geogebra.js",
        bool writeReference = false,
        const Eigen::Vector3d& O = Eigen::Vector3d(0, 0, 0),
        const Eigen::Vector3d& X = Eigen::Vector3d(1, 0, 0),
        const Eigen::Vector3d& Y = Eigen::Vector3d(0, 1, 0),
        const Eigen::Vector3d& Z = Eigen::Vector3d(0, 0, 1));
}

#endif