#include "linesegment.h"
#include "utils.h"
#include <iostream>

toolkit::LineSegment::LineSegment() : m_A(0.0, 0.0, 0.0), m_B(1.0, 0.0, 0.0)
{
}

toolkit::LineSegment::LineSegment(const Eigen::Vector3d & A, const Eigen::Vector3d & B) : 
    m_A(A), 
    m_B(B)
{
}

toolkit::LineSegment::LineSegment(const LineSegment & L) : m_A(L.A()), m_B(L.B())
{
}

toolkit::LineSegment::LineSegment(
    double Ax, 
    double Ay, 
    double Az, 
    double Bx, 
    double By, 
    double Bz) : m_A(Ax, Ay, Az), m_B(Bx, By, Bz)
{
}

toolkit::LineSegment::~LineSegment()
{
}

const Eigen::Vector3d & toolkit::LineSegment::A() const
{
    return m_A;
}

void toolkit::LineSegment::A(const Eigen::Vector3d & A)
{
    m_A << A;
}

Eigen::Vector3d toolkit::LineSegment::At(double t) const
{
    return ((1.0 - t) * m_A) + (t * m_B);
}

const Eigen::Vector3d & toolkit::LineSegment::B() const
{
    return m_B;
}

void toolkit::LineSegment::B(const Eigen::Vector3d & B)
{
    m_B << B;
}

toolkit::LineSegment toolkit::LineSegment::Clone() const
{
    return LineSegment(m_A, m_B);
}

Eigen::Vector3d toolkit::LineSegment::Direction(bool normalize) const
{
    return (normalize) ? (m_B - m_A).normalized() : m_B - m_A;
}

void toolkit::LineSegment::FixZeros(double threshold)
{
    utils::fixZeros(m_A, threshold);
    utils::fixZeros(m_B, threshold);
}

bool toolkit::LineSegment::Intersect(
    const toolkit::LineSegment & L, 
    Eigen::Vector3d & P, 
    double & s, 
    double & t, 
    double threshold) const
{
    // 
    Eigen::Vector3d t1 = m_A - L.A(), t2 = L.Direction();

    // 
    if (abs(t2.x()) < threshold && abs(t2.y()) < threshold && abs(t2.z()) < threshold)
    {
        return false;
    }

    // 
    Eigen::Vector3d t3 = Direction();

    //
    if (abs(t3.x()) < threshold && abs(t3.y()) < threshold && abs(t3.z()) < threshold)
    {
        return false;
    }

    //
    double d1343 = (t1.x() * t2.x()) + (t1.y() * t2.y()) + (t1.z() * t2.z());
    double d4321 = (t2.x() * t3.x()) + (t2.y() * t3.y()) + (t2.z() * t3.z());
    double d1321 = (t1.x() * t3.x()) + (t1.y() * t3.y()) + (t1.z() * t3.z());
    double d4343 = (t2.x() * t2.x()) + (t2.y() * t2.y()) + (t2.z() * t2.z());
    double d2121 = (t3.x() * t3.x()) + (t3.y() * t3.y()) + (t3.z() * t3.z());

    //
    double denom = (d2121 * d4343) - (d4321 * d4321);

    // 
    if (abs(denom) < threshold)
    {
        return false;
    }

    //
    double numer = (d1343 * d4321) - (d1321 * d4343);

    // Calculate the parameters for points Pa and Pb
    s = numer / denom;
    t = (d1343 + (d4321 * s)) / d4343;

    // Calculate points Pa and Pb
    Eigen::Vector3d pa = (t3 * s) + m_A, pb = (t2 * t) + L.A();

    // Calculate the midpoint between points Pa and Pb
    P = (pa + pb) / 2.0;

    // Indicate there is a result and set M as the intersection point between the line segments
    return true;
}

bool toolkit::LineSegment::IsPointIn(
    const Eigen::Vector3d & P, 
    double & t, 
    double threshold) const
{
    // Get the vector from the start vertex of the line segment and point P
    Eigen::Vector3d AP = P - m_A;

    // Get the direction vector of the line segment
    Eigen::Vector3d AB = Direction();

    // Get the cross product between AP and the direction vector of the line segment
    Eigen::Vector3d test = AP.cross(AB);

    // Fix the zeros of the test vector
    utils::fixZeros(test, threshold);

    // If the test result is not the zero vector then P does not lie along the line defined by the 
    // end points of the line segment
    if (!(test(0) == 0.0 && test(1) == 0.0 && test(2) == 0.0))
    {
        return false;
    }

    // Get the dot product between AP and AB, and AB with itself
    double dotAP = AP.dot(AB);
    double dotAB = AB.dot(AB);

    // Fix zeros
    utils::fixZero(dotAP, threshold);
    utils::fixZero(dotAB, threshold);

    // If dotAP is equal to 0 then P matches with the start vertex of the line segment. If it is equal 
    // to dotAB then P matches with the end vertex of the line segment. If it is greater than zero but
    // less than dotAB then P lies in the line segment between the end points of the line segment. In 
    // such cases return true. Otherwise, return false
    return (dotAP == 0.0 || dotAP == dotAB || (dotAP > 0.0 && dotAP < dotAB));
}

Eigen::Vector3d toolkit::LineSegment::Midpoint() const
{
    return (m_A + m_B) / 2.0;
}

void toolkit::LineSegment::Set(const Eigen::Vector3d & A, const Eigen::Vector3d & B)
{
    m_A << A;
    m_B << B;
}

void toolkit::LineSegment::Set(const LineSegment & L)
{
    m_A << L.A();
    m_B << L.B();
}

void toolkit::LineSegment::Set(double Ax, double Ay, double Az, double Bx, double By, double Bz)
{
    m_A << Ax, Ay, Az;
    m_B << Bx, By, Bz;
}

void toolkit::LineSegment::Write() const
{
    std::cout << "A = (" << m_A.x() << ", " << m_A.y() << ", " << m_A.z() << ")" << std::endl;
    std::cout << "B = (" << m_B.x() << ", " << m_B.y() << ", " << m_B.z() << ")" << std::endl;
}

