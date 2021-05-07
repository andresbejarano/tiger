#ifndef _GEOMETRIC_CLASS_H_
#define _GEOMETRIC_CLASS_H_

#pragma once

#include <tiger/ds/vf.h>
#include <tiger/ds/dcel.h>

class GeometricClass
{

protected:

    // The pointer to the DCEL
    std::shared_ptr <dcel::DCEL> m_dcel;

    // The vertices and faces of the geometry
    VF m_vf;

public:

    //
    // Constructor of the class.
    //
    GeometricClass();

    //
    // Constructor of the class.
    //
    GeometricClass(const VF & vf);

    //
    // Destructor of the class.
    //
    virtual ~GeometricClass();

    //
    // Clears the content of the geometric object.
    //
    virtual void Clear();

    //
    // Clears the pointer to the DCEL object with the information of the geometry.
    //
    void ClearDCEL();

    //
    // Returns the pointer to the DCEL that contains the information of the 
    // geometry.
    // @param bool load Indicates whether to load the current information of the 
    // VF object in the DCEL object before returning its pointer.
    // @return std::shared_ptr<dcel::DCEL> The pointer to the DCEL object with the
    // information of the geometry.
    //
    std::shared_ptr<dcel::DCEL> DCEL(bool load = false);

    //
    // Returns the reference to the VF object with the information of the geometry.
    // @return VF & The reference to the object with the information of the 
    // geometry.
    //
    VF & Geometry();

    //
    // Sets the geometry of the object. This function clear the content of the 
    // object before setting the given geometry.
    // @param const VF & vf The reference to the vertex coordinates and vertex 
    // indices of the geometry.
    //
    void Geometry(const VF & vf);

    //
    //
    bool IsDCELLoaded();

    //
    // Loads the content of the VF object into the DCEL object.
    //
    void LoadDCEL();

    //
    //
    void UpdateFromDCEL();

};

#endif