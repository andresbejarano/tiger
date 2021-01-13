#include "geometricclass.h"

GeometricClass::GeometricClass() : m_dcel(nullptr) 
{
}

GeometricClass::GeometricClass(const VF & vf) : m_dcel(nullptr)
{
    Geometry(vf);
}

GeometricClass::~GeometricClass()
{
}

void GeometricClass::Clear()
{
    ClearDCEL();

    m_vf.clear();
}

void GeometricClass::ClearDCEL()
{
    if (m_dcel) 
    {
        m_dcel->Clear();
        m_dcel = nullptr;
    }
}

std::shared_ptr<dcel::DCEL> GeometricClass::DCEL(bool load) 
{
    if (load) 
    {
        LoadDCEL();
    }

    return m_dcel;
}

VF & GeometricClass::Geometry() 
{
    return m_vf;
}

void GeometricClass::Geometry(const VF & vf)
{
    Clear();

    m_vf.Set(vf);
}

bool GeometricClass::IsDCELLoaded()
{
    return m_dcel != nullptr;
}

void GeometricClass::LoadDCEL() 
{
    ClearDCEL();

    m_dcel = std::make_shared<dcel::DCEL>(m_vf);
}

void GeometricClass::UpdateFromDCEL()
{
    assert(m_dcel);
    m_vf = m_dcel->vf();
}
