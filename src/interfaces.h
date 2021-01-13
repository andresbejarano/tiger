#ifndef _INTERFACES_H_
#define _INTERFACES_H_

#pragma once

#include "vf.h"

/*
The class representing a collection of interface polygons from an assembly of blocks.
*/
class InterfacePolygons
{

private:

    // The vector to store the pointers to the geometries of the interface polygons. It is possible
    // an interface between two blocks doesn't exist. In such cases we represent that by having a 
    // null pointer
    std::vector<std::shared_ptr<VF>> m_interfaces;

public:
    
    /*
    Constructor of the class.
    */
    InterfacePolygons();

    /*
    Destructor of the class.
    */
    ~InterfacePolygons();

    /*
    Clears any content stored in the object.
    */
    void Clear();

    /*
    Returns the number of existing interface polygons. It might not be the same as
    GetInterfaces().size() which is the length of the vector storing the pointers to the interface
    polygons (some of them could be null pointers).
    @return size_t The number of interface polygons.
    */
    size_t CountInterfaces() const;

    /*
    Returns the number of triangles required to represent the interface polygons.
    @return size_t The number of triangles required to represent the interface polygons.
    */
    size_t CountTriangles() const;

    /*
    Returns the number of vertices of the interface polygons.
    @return size_t The number of vertices of the interface polygons.
    */
    size_t CountVertices() const;

    /*
    Returns the reference to the vector with the interface polygons.
    @return const std::vector<VF> & The reference to the vector with the pointers to the geometries
    of the interface polygons.
    */
    const std::vector<std::shared_ptr<VF>> & GetInterfaces() const;

    /*
    Sets the interface polygons. This method clears any previous content and copies the pointers 
    given in the vector. No new objects are generated.
    @param const std::vector<std::shared_ptr<VF>> & interfaces The reference to the vector with the
    pointers to the geometries of the interface polygons.
    */
    void Set(const std::vector<std::shared_ptr<VF>> & interfaces);

    /*
    Writes the content of the interface polygons.
    */
    void Write() const;

    /*
    Writes the content of the interface polygons.
    @param std::stringstream & ss The reference to the string stream for writing the interfaces 
    content.
    */
    void Write(std::stringstream & ss) const;

};

#endif
