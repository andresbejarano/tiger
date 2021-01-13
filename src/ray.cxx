#include "ray.h"
#include "utils.h"
#include <iostream>

toolkit::Ray::Ray() :
    m_A(0.0, 0.0, 0.0),
    m_D(1.0, 0.0, 0.0) 
{
}

toolkit::Ray::Ray(const Eigen::Vector3d & P, const Eigen::Vector3d & D, bool normalize) : 
    m_A(P), 
    m_D(D)
{
    if (normalize) 
    {
        m_D.normalize();
    }
}

Eigen::Vector3d toolkit::Ray::at(double t) const
{
    return m_A + (m_D * t);
}

toolkit::Ray toolkit::Ray::clone() const
{
    return Ray(m_A, m_D);
}

const Eigen::Vector3d & toolkit::Ray::direction() const
{
    return m_D;
}

bool toolkit::Ray::isPointIn(const Eigen::Vector3d & P, double & t, double threshold) const
{
    // Get the normalized vector between the start point of the ray and P. Fix its zeros
    Eigen::Vector3d NP = (P - m_A).normalized();
    utils::fixZeros(NP, threshold);

    // Get the normalized direction vector of the ray. Fix its zeros
    Eigen::Vector3d ND = m_D.normalized();
    utils::fixZeros(ND, threshold);

    // Calculate the dot product between the normalized vectors
    double dot = NP.dot(ND);

    // If the absolute value of the dot product is not equal to 1 then return false since P does 
    // not lie along the ray (if it is -1 the it lies in the opposite direction with respecto to 
    // the direction vector)
    if (abs(dot) != 1.0) 
    {
        return false;
    }

    // Calculate the parameter for the point along the ray. Try different coordinates if the 
    // denominator is zero for any of them
    if (abs(m_D.x()) > threshold) 
    {
        t = (P.x() - m_A.x()) / m_D.x();
    }
    else if (abs(m_D.y()) > threshold) 
    {
        t = (P.y() - m_A.y()) / m_D.y();
    }
    else if (abs(m_D.z()) > threshold) 
    {
        t = (P.z() - m_A.z()) / m_D.z();
    }
    else 
    {
        assert(false);
    }
    
    // Return true indicating that P lies along the ray
    return true;
}

void toolkit::Ray::normalize()
{
    m_D.normalize();
}

void toolkit::Ray::set(const Eigen::Vector3d & P, const Eigen::Vector3d & D, bool normalize)
{
    m_A << P;

    setDirection(D, normalize);
}

void toolkit::Ray::setDirection(const Eigen::Vector3d & D, bool normalize)
{
    m_D << D;

    // Normalize the direction vector if indicated
    if (normalize)
    {
        m_D.normalize();
    }
}

void toolkit::Ray::setStart(const Eigen::Vector3d & P)
{
    m_A << P;
}

const Eigen::Vector3d & toolkit::Ray::start() const
{
    return m_A;
}

void toolkit::Ray::write() const
{
    std::cout << "A = (" << m_A.x() << ", " << m_A.y() << ", " << m_A.z() << ")" << std::endl;
    std::cout << "D = (" << m_D.x() << ", " << m_D.y() << ", " << m_D.z() << ")" << std::endl;
}
