#include <tiger/toolkit/plane.h>
#include <tiger/utils.h>
#include <iostream>

toolkit::Plane::Plane() : m_N(0.0, 1.0, 0.0), m_P(0.0, 0.0, 0.0) 
{
}

toolkit::Plane::Plane(
    const Eigen::Vector3d & N, 
    const Eigen::Vector3d & P, 
    bool normalize, 
    bool fixZeros, 
    double threshold) : m_N(N), m_P(P)
{
    if (normalize)
    {
        m_N.normalize();
    }

    if (fixZeros) 
    {
        utils::fixZeros(m_N, threshold);
        utils::fixZeros(m_P, threshold);
    }
}

toolkit::Plane::Plane(
    double nx, 
    double ny, 
    double nz, 
    double px, 
    double py, 
    double pz, 
    bool normalize, 
    bool fixZeros, 
    double threshold) : m_N(nx, ny, nz), m_P(px, py, pz)
{
    if (normalize)
    {
        m_N.normalize();
    }

    if (fixZeros) 
    {
        utils::fixZeros(m_N, threshold);
        utils::fixZeros(m_P, threshold);
    }
}

toolkit::Plane::~Plane() 
{
}

toolkit::Plane toolkit::Plane::Clone() const
{
    return Plane(m_N, m_P);
}

bool toolkit::Plane::IsPointInPlane(const Eigen::Vector3d & Q, double threshold) const 
{
    return PointLocation(Q, threshold) == 0;
}

const Eigen::Vector3d & toolkit::Plane::Normal() const 
{
    return m_N;
}

void toolkit::Plane::Normal(const Eigen::Vector3d & N, bool normalize)
{
    m_N << N;

    if (normalize)
    {
        m_N.normalize();
    }
}

void toolkit::Plane::Normal(double x, double y, double z, bool normalize)
{
    m_N << x, y, z;

    if (normalize)
    {
        m_N.normalize();
    }
}

const Eigen::Vector3d & toolkit::Plane::Point() const 
{
    return m_P;
}

void toolkit::Plane::Point(const Eigen::Vector3d & P)
{
    m_P << P;
}

void toolkit::Plane::Point(double x, double y, double z)
{
    m_P << x, y, z;
}

int toolkit::Plane::PointLocation(const Eigen::Vector3d & Q, double threshold) const
{
    // Calculate the dot product between the normal vector and the PQ vector. Then, return its 
    // value (check if it is close to zero)
    double dot = m_N.dot(Q - m_P);
    return (abs(dot) <= threshold) ? 0 : ((dot > 0.0) ? 1 : -1);
}

void toolkit::Plane::Set(const Eigen::Vector3d & N, const Eigen::Vector3d & P, bool normalize)
{
    m_N << N;
    m_P << P;

    if (normalize)
    {
        m_N.normalize();
    }
}

void toolkit::Plane::Set(
    double nx, 
    double ny, 
    double nz, 
    double px, 
    double py, 
    double pz, 
    bool normalize) 
{
    m_N << nx, ny, nz;
    m_P << px, py, pz;

    if (normalize)
    {
        m_N.normalize();
    }
}

void toolkit::Plane::Write() const 
{
    std::cout << "N = (" << m_N.x() << ", " << m_N.y() << ", " << m_N.z() << ")" << std::endl;
    std::cout << "P = (" << m_P.x() << ", " << m_P.y() << ", " << m_P.z() << ")" << std::endl;
}

double toolkit::Plane::a() const
{
    return m_N.x();
}

double toolkit::Plane::b() const
{
    return m_N.y();
}

double toolkit::Plane::c() const
{
    return m_N.z();
}

double toolkit::Plane::d() const
{
    return m_N.dot(m_P);
}