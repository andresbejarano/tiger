#ifndef _PLANE_H_
#define _PLANE_H_

#pragma once

#include <Eigen/Core>

namespace toolkit 
{
    //
    // The class to represent a plane in 3D. A plane is defined by a normal vector N and a point P, 
    // the implicit definition of the plane is N * (Q - P) = 0, where Q is any point in the plane and 
    // * is the dot product.
    //
    class Plane
    {

    private:

        // The normal vector of the plane
        Eigen::Vector3d m_N;

        // A point in the plane
        Eigen::Vector3d m_P;

    public:

        //
        // Constructor of the class. 
        // The plane initializes with normal vector (0, 1, 0) and point (0, 0, 0).
        //
        Plane();

        //
        // Constructor of the class. Represents a plane defined by the normal vector N and a point P 
        // in the plane.
        // @param const Eigen::Vector3d & N The reference to the normal vector of the plane.
        // @param const Eigen::Vector3d & P The reference to a point in the plane.
        // @param bool normalize Indicates whether to normalize the normal vector.
        // @param bool fixZeros Indicates whether to fix the zero values.
        // @param double threshold The threshold for values close to zero.
        //
        Plane(
            const Eigen::Vector3d & N, 
            const Eigen::Vector3d & P, 
            bool normalize = false, 
            bool fixZeros = false, 
            double threshold = 1e-8);

        //
        // Constructor of the class.
        // @param double nx The X coordinate value of the normal vector.
        // @param double ny The Y coordinate value of the normal vector
        // @param double nz The Z coordinate value of the normal vector.
        // @param double px The X coordinate value of the point in the plane.
        // @param double py The Y coordinate value of the point in the plane.
        // @param double pz The Z coordinate value of the point in the plane.
        // @param bool normalize Indicates whether to normalize the normal vector.
        // @param bool fixZeros Indicates whether to fix the zero values.
        // @param double threshold The threshold for values close to zero.
        //
        Plane(
            double nx,
            double ny,
            double nz,
            double px,
            double py,
            double pz,
            bool normalize = false, 
            bool fixZeros = false,
            double threshold = 1e-8);

        //
        // Destructor of the class.
        //
        ~Plane();

        //
        // Returns a clone of the plane.
        // @return Plane A clone of the plane.
        //
        Plane Clone() const;

        //
        // Checks if a point Q is in the plane. Point Q must satisfy Dot(N, Q - P) = 0.
        // @param const Eigen::Vector3d & Q The reference to a point.
        // @param double threshold The threshold for values close to zero.
        // @return bool Indicates if the point is in the plane.
        //
        bool IsPointInPlane(const Eigen::Vector3d & Q, double threshold = 1e-8) const;

        //
        // Returns the reference to the normal vector of the plane.
        // @return const Eigen::Vector3d The normal vector of the plane.
        //
        const Eigen::Vector3d & Normal() const;

        //
        // Sets the normal vector of the plane.
        // @param const Eigen::Vector3d & N The reference to the normal vector of the plane.
        // @param bool normalize Indicates whether to normalize or not the normal vector.
        //
        void Normal(const Eigen::Vector3d & N, bool normalize = false);

        //
        // Sets the normal vector of the plane.
        // @param double x The X coordinate value of the normal vector.
        // @param double y The Y coordinate value of the normal vector.
        // @param double z The Z coordinate value of the normal vector.
        // @param bool normalize Indicates whether to normalize or not the normal vector.
        //
        void Normal(double x, double y, double z, bool normalize = false);

        //
        // Returns the reference to the point of the plane.
        // @return const Eigen::Vector3d The reference to the point of the plane.
        //
        const Eigen::Vector3d & Point() const;

        //
        // Sets the point of the plane.
        // @param const Eigen::VBector3d & P The reference to the point in the plane.
        //
        void Point(const Eigen::Vector3d & P);

        //
        // Sets the point of the plane.
        // @param double x The X coordinate value of the point in the plane.
        // @param double y The Y coordinate value of the point in the plane.
        // @param double z The Z coordinate value of the point in the plane.
        //
        void Point(double x, double y, double z);

        //
        // Checks the location of a point Q with respect of the plane. A 1 value indicates the point 
        // is in the positive half space of the plane (in the direction of the normal vector), a 0 
        // value indicates the point is at the plane, a -1 value indicates the point is at the 
        // negative half space of the plane (in the opposite direction of the normal vector). The 
        // indicator comes from Dot(N, Q - P).
        // @param const Eigen::Vector3d & Q The reference to a point.
        // @param double threshold The threshold for values close to zero.
        // @return int The indicator of the location of the point.
        //
        int PointLocation(const Eigen::Vector3d & Q, double threshold = 1e-8) const;

        //
        // Sets the normal vector and point of the plane.
        // @param const Eigen::Vector3d & N The reference to the normal vector of the plane.
        // @param const Eigen::Vector3d & P The reference to a point in the plane.
        // @param bool normalize Indicates whether to normalize or not the normal vector.
        //
        void Set(const Eigen::Vector3d & N, const Eigen::Vector3d & P, bool normalize = false);

        //
        // Sets the normal vector and point of the plane.
        // @param double nx The X coordinate value of the normal vector.
        // @param double ny The Y coordinate value of the normal vector
        // @param double nz The Z coordinate value of the normal vector.
        // @param double px The X coordinate value of the point in the plane.
        // @param double py The Y coordinate value of the point in the plane.
        // @param double pz The Z coordinate value of the point in the plane.
        // @param bool normalize Indicates whether to normalize the normal vector.
        //
        void Set(
            double nx,
            double ny,
            double nz,
            double px,
            double py,
            double pz,
            bool normalize = false);

        //
        // Writes the information of the plane.
        //
        void Write() const;

        //
        // Returns the coefficient of the X variable of the plane equation ax + by + cz = d.
        // @return double The coefficient of the X variable.
        //
        double a() const;

        //
        // Returns the coefficient of the Y variable of the plane equation ax + by + cz = d.
        // @return double The coefficient of the Y variable.
        //
        double b() const;

        //
        // Returns the coefficient of the Z variable of the plane equation ax + by + cz = d.
        // @return double The coefficient of the Z variable.
        //
        double c() const;

        //
        // Returns the d value of the plane equation ax + by + cz = d.
        // @return double The d value.
        //
        double d() const;

    };
}

#endif
