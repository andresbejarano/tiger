#ifndef _TESSELLATION_H_
#define _TESSELLATION_H_

#pragma once

#include <tiger/tic/geometricclass.h>

//
// The class representing a tessellation.
//
class Tessellation : public GeometricClass
{

public:

    // 
    enum TILE_SUBDIVISION_TYPE
    {
        QUADRILATERALS, 
		MIDPOINTS, 
		TRIANGLES,
    };

    //
    // Constructor of the class.
    //
    Tessellation();

    //
    // Constructor of the class. Copies the given vertex coordinates and vertex indices.
    // @param const VF & vf The reference to the vertex coordinates and vertex indices of the 
    // geometry.
    //
    Tessellation(const VF & vf);

    //
    // Destructor of the class.
    //
    ~Tessellation();

    //
    // Checks if all tiles in the tessellation have a center point stored in the ATTRIB_CENTER dynamic
    // attribute. NOTE: call LoadDCEL before if the tessellation has changed.
    // @return bool Indicates if all tiles have a center point.
    //
    bool CheckTileCenters() const;

    //
    // Checks if all tiles in the tessellation have an even number of sides.
    // @return bool Indicates if all tiles have an even number of sides.
    //
    bool CheckTilesAreEvenSided() const;

    //
    // Checks if all half edges in the tessellation have both direction value and direction vector
    // stored in the respective ATTRIB_DIRECTION_VALUE and ATTRIB_DIRECTION_VECTOR dynamic attributes.
    // NOTE: call LoadDCEL before if the tessellation has changed.
    // @return bool Indicates if all half edges have a direction value and a direction vector.
    //
    bool CheckHalfedgeDirections() const;

    //
    // Clears the content of the tessellation.
    //
    void Clear();

    //
    // Returns the tile subdivision value based on the given subdivision name. If the given 
    // subdivision name is not valid then it returns the quadrilaterals type.
    // @param const std::string subdivision The subdivision name.
    // @return tile_SUBDIVISION_TYPE The respective tile subdivision type.
    //
    static TILE_SUBDIVISION_TYPE GetTileSubdivisionType(const std::string subdivision);

    //
    // Indicates if the half edges of the tessellation have the direction values. NOTE: call 
    // LoadDCEL before if the tessellation has changed.
    // @return bool Indicates if all half edges in the tessellation have the direction values.
    //
    bool HasEdgeDirections() const;

    //
    //
    void QuadrangulateTiles();

private:

    //
    // Sets the direction values on the edges of the tessellation starting on the given half edge.
    // Direction values are stored in the ATTRIB_DIRECTION_VALUE dynamic attribute of each half edge.
    // A positive direction value indicates it goes inward the tile, a negative direction value
    // indicates it goes outwards the tile.
    // @param std::shared_ptr<dcel::Halfedge> halfedge The pointer to the initial half edge in the
    // tessellation.
    // @param int dir The initial direction value.
    //
    void SetEdgeDirectionValueFromHalfedge(std::shared_ptr<dcel::Halfedge> halfedge, int dir);

public:

    //
    // Sets the direction values on the edges of the tessellation. Direction values are stored in
    // the ATTRIB_DIRECTION_VALUE dynamic attribute of each half edge. A positive direction value
    // indicates it goes inward the tile, a negative direction value indicates it goes outwards the
    // tile. NOTE: call LoadDCEL before if the tessellation has changed.
    // @param int dir The initial direction value.
    //
    void SetEdgeDirectionValues(int dir);

private:

    //
    // Sets the direction vectors on the edges of the tessellation. Direction vectors are stored
    // in the ATTRIB_DIRECTION_VECTOR dynamic attribute on each half edge. We assume the half edges
    // have the ATTRIB_DIRECTION_VALUE dynamic attribute and tiles have the ATTRIB_CENTER dynamic
    // attribute. NOTE: call LoadDCEL before if the tessellation has changed.
    //
    void SetEdgesDirectionVector();

public:

    //
    // Sets the directions on the edges of the tessellation. Direction values are stored in the
    // ATTRIB_DIRECTION_VALUE dynamic attribute on each half edge. Direction vectors are stored in the
    // ATTRIB_DIRECTION_VECTOR dynamic attribute on each half edge. We assume there is a geometric
    // domain loaded in workspace, and all tiles have the ATTRIB_CENTER dynamic attribute set.
    //
    void SetEdgesDirection(int dir);

    //
    // Sets the midpoints on the edges of the tessellation. Midpoints are stored in the
    // ATTRIB_MIDPOINT dynamic attribute of the half edges. NOTE: call LoadDCEL before if the 
    // tessellation has changed.
    //
    void SetEdgesMidpoint();

    //
    // Sets the normal vector on the edges of the tessellation. Normal vectors are stored in the
    // ATTRIB_NORMAL dynamic attribute of the half edges. NOTE: call LoadDCEL before if the geometric 
    // domain has changed.
    //
    void SetEdgesNormal();

    //
    // Sets the rotated vector on the edges of the tessellation using the given angle value. We
    // assume the half edges have the ATTRIB_DIRECTION_VALUE dynamic attribute. Rotated vectors are
    // stored in the ATTRIB_ROTATED_VECTOR dynamic attributes of the half edges. NOTE: call LoadDCEL 
    // before if the tessellation has changed.
    // @param double angle The rotation angle for all rotated vectors.
    //
    void SetEdgeRotatedVectorsUsingAngle(double angle);

    //
    // Sets the rotated vectors on the edges of the tessellation using the top and bottom section
    // points of the tiles. We assume each tile in the tessellation has the ATTRIB_TOP_POINT and
    // ATTRIB_BOTTOM_POINT dynamic attributes, and each half edge has the ATTRIB_DIRECTION_VALUE and
    // ATTRIB_MIDPOINT dynamic attributes. Rotated vectors are stored in the ATTRIB_ROTATED_VECTOR
    // dynamic attribute of the half edges. NOTE: call LoadDCEL before if the tessellation has 
    // changed.
    // @param bool boundary Indicates whether to calculate the pieces at the boundary of the geometric
    // domain.
    //
    void SetEdgeRotatedVectorsUsingTopBottomSectionPoints(bool boundary);

public:

    //
    // Checks if each tile in the tessellation is at the boundary. Boundary indicator is stored in
    // the ATTRIB_BOUNDARY dynamic attribute of each tile. NOTE: call LoadDCEL before if the geometric
    // domain has changed.
    //
    void SetTileBoundaries();

    //
    // Sets the centroids as a dynamic attribute of each tile in the tessellation. Centers 
    // are stored in the ATTRIB_tile dynamic attribute. NOTE: call LoadDCEL before if the geometric 
    // domain has changed.
    //
    void SetTilesCenter();

    //
    // Calculates the top and bottom section points of each tile in the tessellation. We assume
    // all tiles have the ATTRIB_CENTER dynamic attribute. Top and bottom points are stored in the
    // ATTRIB_TOP_POINT and ATTRIB_BOTTOM_POINT dynamic attributes of the tile respectively. 
    // NOTE: call LoadDCEL before if the tessellation has changed.
    // @param double height The height parameter from the center of the tile to the respective section
    // points.
    //
    void SetTileTopBottomSectionPoints(double height);

    //
    //
    void SubdivideTilesByMidpoints();

    //
    //
    void TriangulateTilesByVertices();

    //
    // Writes the content of the tessellation.
    //
    void Write() const;

    //
    // Writes the content of the tessellation.
    // @param std::stringstream & ss The reference to the string stream for writing the content.
    //
    void Write(std::stringstream & ss) const;

};

#endif
