#include <tiger/tic/interfaces.h>
#include <iostream>

InterfacePolygons::InterfacePolygons()
{
}

InterfacePolygons::~InterfacePolygons()
{
    Clear();
}

void InterfacePolygons::Clear()
{
    // Traverse through the interface polygons and delete them properly. Keep in mind there could 
    // be null pointers in the interfaces vector (they occur when two incident pieces do not have 
    // coplanar faces between them)
    for (auto it = m_interfaces.begin(); it != m_interfaces.end(); ++it) 
    {
        if (*it) 
        {
            *it = nullptr;
        }
    }

    m_interfaces.clear();
}

size_t InterfacePolygons::CountInterfaces() const
{
    // Initialize the count of existing interface polygons
    size_t nInterfaces = 0;

    // Traverse through the vector with the pointers to the interface polygons, increment the 
    // counter when there is a non-null pointer
    for (auto it = m_interfaces.begin(); it != m_interfaces.end(); ++it)
    {
        if (*it)
        {
            nInterfaces += 1;
        }
    }

    // Return the number of existing interface polygons
    return nInterfaces;
}

size_t InterfacePolygons::CountTriangles() const
{
    // Initialize the count of triangles required for representing the interface polygons
    size_t nTriangles = 0;

    // Traverse through the interfaces and count the number of triangles required for representing 
    // each one of them
    for (auto it = m_interfaces.begin(); it != m_interfaces.end(); ++it)
    {
        if (*it)
        {
            nTriangles += (*it)->countTriangles();
        }
    }

    // Return the number of triangles
    return nTriangles;
}

size_t InterfacePolygons::CountVertices() const
{
    // Initialize the number of vertices of the interface polygons
    size_t nVertices = 0;

    // Traverse through the interface polygons and count the number of vertices of each interface
    for (auto it = m_interfaces.begin(); it != m_interfaces.end(); ++it)
    {
        if (*it)
        {
            nVertices += (*it)->countVertices();
        }
    }

    // Return the number of vertices
    return nVertices;
}

const std::vector<std::shared_ptr<VF>> & InterfacePolygons::GetInterfaces() const
{
    return m_interfaces;
}

void InterfacePolygons::Set(const std::vector<std::shared_ptr<VF>> & interfaces)
{
    // Clear the interfaces and resize the vector
    Clear();

    // Get the number of given interface polygons
    size_t nInterfaces = interfaces.size();

    // Resize the vector to store the pointers of the interface polygons
    m_interfaces.resize(nInterfaces, nullptr);

    // Traverse through the vector with the given interface polygons and copy the pointers
    for (size_t i = 0; i < nInterfaces; i += 1) 
    {
        m_interfaces[i] = interfaces[i];
    }
}

void InterfacePolygons::Write() const
{
    // Get the number of interface polygons
    size_t nInterfaces = m_interfaces.size();

    // Write the label indicating the start of the interfaces content
    std::cout << "interfaces " << nInterfaces << std::endl;

    // Traverse through the interfaces and write their content
    for (size_t i = 0; i < nInterfaces; i += 1)
    {
        if (m_interfaces[i]) 
        {
            std::cout << i << std::endl;
            m_interfaces[i]->Write();
        }
    }

    // Write the label indicating the end of the interfaces content
    std::cout << "/interfaces" << std::endl;
}

void InterfacePolygons::Write(std::stringstream & ss) const
{
    // Get the number of interface polygons
    size_t nInterfaces = m_interfaces.size();

    // Write the label indicating the start of the interfaces content
    ss << "interfaces " << nInterfaces << std::endl;

    // Traverse through the interfaces and write their content
    for (size_t i = 0; i < nInterfaces; i += 1) 
    {
        if (m_interfaces[i]) 
        {
            ss << i << std::endl;
            m_interfaces[i]->Write(ss);
        }
    }

    // Write the label indicating the end of the interfaces content
    ss << "/interfaces" << std::endl;
}