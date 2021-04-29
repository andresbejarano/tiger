#ifndef _RAY_H_
#define _RAY_H_

#pragma once

#include <Eigen/Core>

namespace toolkit 
{
    /// <summary>
    /// The class representing a ray in 3D. A ray is defined by a point A and a direction vector D,
    /// a point in the ray is P(t) = A + tD, where t is a scalar parameter.
    /// </summary>
    class Ray
    {
    
    private:
        
        // The start point of the ray
        Eigen::Vector3d m_A;
        
        // The direction vector of the ray
        Eigen::Vector3d m_D;
    
    public:
        
        /// <summary>
        /// Constructor of the class. The ray is initialize with point (0, 0, 0) and direction 
        /// vector along the X axis(1, 0, 0).
        /// </summary>
        Ray();
        
        /// <summary>
        /// Constructor of the class.
        /// </summary>
        /// <param name="P">The reference to the coordinates of the start point.</param>
        /// <param name="D">The reference to the coordinates of the direction vector.</param>
        /// <param name="normalize">Indicates whether to normalize the direction vector.</param>
        Ray(const Eigen::Vector3d & P, const Eigen::Vector3d & D, bool normalize = false);
        
        /// <summary>
        /// Returns the point at parameter t along the ray.
        /// </summary>
        /// <param name="t">The parameter value.</param>
        /// <returns>The point at parameter t along the ray.</returns>
        Eigen::Vector3d at(double t) const;
        
        /// <summary>
        /// Returns a clone of the ray.
        /// </summary>
        /// <returns>Ray A clone of the ray.</returns>
        Ray clone() const;
        
        /// <summary>
        /// Returns the reference to the direction vector of the ray.
        /// </summary>
        /// <returns>The reference to the direction vector.</returns>
        const Eigen::Vector3d & direction() const;
        
        /// <summary>
        /// Indicates whether a point P lies along the ray.
        /// </summary>
        /// <param name="P">The reference to the point.</param>
        /// <param name="t">The parameter for the point if it lies along the ray.</param>
        /// <param name="threshold">The threshold for values close to zero.</param>
        /// <returns>Indicates whether the point lies along the ray.</returns>
        bool isPointIn(const Eigen::Vector3d & P, double & t, double threshold = 1e-8) const;
        
        /// <summary>
        /// Normalizes the direction vector of the ray.
        /// </summary>
        void normalize();
        
        /// <summary>
        /// Sets the start point and the direction vector of the ray.
        /// /// </summary>
        /// <param name="P">The reference to the coordinates of the start point.</param>
        /// <param name="D">The reference to the coordinates of the direction vector.</param>
        /// <param name="normalize">Indicates whether to normalize the direction vector.</param>
        void set(const Eigen::Vector3d & P, const Eigen::Vector3d & D, bool normalize = false);
        
        /// <summary>
        /// Sets the direction vector of the ray.
        /// </summary>
        /// <param name="D">The reference to the coordinates of the direction vector.</param>
        /// <param name="normalize">Indicates whether to normalize the direction vector.</param>
        void setDirection(const Eigen::Vector3d & D, bool normalize = false);
        
        /// <summary>
        /// Sets the start point of the ray.
        /// </summary>
        /// <param name="P">The reference to the coordinates of the start point.</param>
        void setStart(const Eigen::Vector3d & P);
        
        /// <summary>
        /// Returns the reference to the start point of the ray.
        /// </summary>
        /// <returns>The reference to the start point.</returns>
        const Eigen::Vector3d & start() const;
        
        /// <summary>
        /// Writes the information of the ray.
        /// </summary>
        void write() const;

    };
}

#endif
