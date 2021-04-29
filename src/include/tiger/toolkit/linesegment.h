#ifndef _LINESEGMENT_H_
#define _LINESEGMENT_H_

#pragma once

#include <Eigen/Core>

namespace toolkit 
{
    /*
    The class to represent a line segment in 3D. A line segment is defined by its end points A and 
    B. A point in the line segment is P(t) = (1 - t)A + tB, where t is a scalar in range [0, 1].
    */
    class LineSegment
    {

    private:

        // The first end point of the line segment
        Eigen::Vector3d m_A;

        // The second end point of the line segment
        Eigen::Vector3d m_B;

    public:

        /*
        Constructor of the class.
        The line segment initializes with end points (0, 0, 0) and (1, 0, 0).
        */
        LineSegment();

        /*
        Constructor of the class.
        @param const Eigen::Vector3d & A The reference to the coordinates of the first end point.
        @param const Eigen::Vector3d & B The reference to the coordinates of the second end point.
        */
        LineSegment(const Eigen::Vector3d & A, const Eigen::Vector3d & B);

        /*
        Constructor of the class.
        @param const LineSegment & L The reference to a line segment.
        */
        LineSegment(const LineSegment & L);

        /*
        Constructor of the class.
        @param double Ax The X coordinate of the first end point.
        @param double Ay The Y coordinate of the first end point.
        @param double Az The Z coordinate of the first end point.
        @param double Bx The X coordinate of the second end point.
        @param double By The Y coordinate of the second end point.
        @param double Bz The Z coordinate of the second end point.
        */
        LineSegment(double Ax, double Ay, double Az, double Bx, double By, double Bz);

        /*
        Destructor of the class.
        */
        ~LineSegment();

        /*
        Returns the reference to the first end point of the line segment.
        @return const Eigen::Vector3d % The reference to the first end point of the line segment.
        */
        const Eigen::Vector3d & A() const;

        /*
        Sets the coordinates of the first end point of the line segment.
        @param const Eigen::Vector3d & A The coordinates for the first end point.
        */
        void A(const Eigen::Vector3d & A);

        /*
        Returns the point at parameter t along the line segment.
        @param double t The parameter value.
        @return Eigen::Vector3d The point at parameter t along the line segment.
        */
        Eigen::Vector3d At(double t) const;

        /*
        Returns the reference to the second end point of the line segment.
        @return const Eigen::Vector3d % The reference to the second end point of the line segment.
        */
        const Eigen::Vector3d & B() const;

        /*
        Sets the coordinates of the second end point of the line segment.
        @param const Eigen::Vector3d & B The coordinates for the second end point.
        */
        void B(const Eigen::Vector3d & B);

        /*
        Returns a clone of the line segment.
        @return Ray A clone of the line segment.
        */
        LineSegment Clone() const;

        /*
        Returns the vector representing the direction of the line segment.
        @param double normalize Indicates whether to normalize the direction vector or not.
        @return const Eigen::Vector3d The direction vector of the line segment.
        */
        Eigen::Vector3d Direction(bool normalize = false) const;

        /*
        Fix the zeros of the end points of the line segment.
        @param double threshold The threshold for values close to zero.
        */
        void FixZeros(double threshold = 1e-8);

        /*
        Calculates the intersection between the line segment and another line segment. If both line
        segments are skewed then it calculates the midpoint between the closest points from both
        line segments. The algorithm is adapted from the solution by Paul Bourke in "The shortest 
        line between two lines in 3D" found in http://paulbourke.net/geometry/pointlineplane/
        @param const toolkit::LineSegment & L The reference to a line segment.
        @param Eigen::Vector3d & P The reference to store the intersection point.
        @param double & s The parameter of the intersection point along the line segment.
        @param double & t The parameter of the intersection point along the given line segment.
        @param double threshold The threshold for values close to zero.
        @return bool Indicates whether there is an intersection point or not.
        */
        bool Intersect(
            const toolkit::LineSegment & L, 
            Eigen::Vector3d & P, 
            double & s, 
            double & t, 
            double threshold = 1e-8) const;

        /*
        Indicates whether a point P lies along the line segment.
        @param const Eigen::Vector3d & P The reference to the point.
        @param double & t The parameter for the point if it lies along the line segment.
        @param double threshold The threshold for values close to zero.
        @return Indicates whether the point lies along the line segment.
        */
        bool IsPointIn(const Eigen::Vector3d & P, double & t, double threshold = 1e-8) const;

        /*
        Returns the midpoint of the line segment.
        @return Eigen::Vector3d The midpoint of the line segment.
        */
        Eigen::Vector3d Midpoint() const;

        /*
        Sets the end points of the line segment.
        @param const Eigen::Vector3d & A The reference to the coordinates of an end point.
        @param const Eigen::Vector3d & B The reference to the coordinates of an end point.
        */
        void Set(const Eigen::Vector3d & A, const Eigen::Vector3d & B);

        /*
        Sets the end points of a line segment.
        @param const LineSegment & L The reference to a line segment.
        */
        void Set(const LineSegment & L);

        /*
        Sets the end points of the line segment.
        @param double Ax
        @param double Ay
        @param double Az
        @param double Bx
        @param double By
        @param double Bz
        */
        void Set(double Ax, double Ay, double Az, double Bx, double By, double Bz);

        /*
        Writes the information of the line segment.
        */
        void Write() const;

    };
}

#endif
