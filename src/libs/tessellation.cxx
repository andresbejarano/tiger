#include <tiger/tic/tessellation.h>
#include <Eigen/Geometry>
#include <queue>
#include <tiger/utils.h>
#include <iostream>

Tessellation::Tessellation() 
{
}

Tessellation::Tessellation(const VF & vf) : GeometricClass(vf)
{
}

Tessellation::~Tessellation() 
{
    Clear();
}

bool Tessellation::CheckTileCenters() const 
{
    // If there is no DCEL then return false since the centers are stored at each tile
    if (!m_dcel) 
    {
        return false;
    }

    // Traverse through the tiles of the geometry and check if they have the ATTRIB_CENTER 
    // attribute
    for (auto it = m_dcel->Faces().begin(); it != m_dcel->Faces().end(); ++it)
    {
        // Get the pointer to the current tile
        std::shared_ptr<dcel::Face> tile = *it;

        // If the tile doesn't have the center attribute then return false
        if (!tile->Attributes().Has(ATTRIB_CENTER))
        {
            return false;
        }
    }

    // Return true since all tiles in the geometric domain have the center attribute
    return true;
}

bool Tessellation::CheckTilesAreEvenSided() const 
{
    return m_vf.facesHaveEvenNumberOfSides();
}

bool Tessellation::CheckHalfedgeDirections() const 
{
    assert(m_dcel);

    // Traverse through the half edges of the DCEL and check they have the direction attribute
    for (auto itH = m_dcel->Halfedges().begin(); itH != m_dcel->Halfedges().end(); ++itH)
    {
        // Get the pointer to the current half edge
        std::shared_ptr<dcel::Halfedge> halfedge = *itH;

        // If the half edge doesn't have both direction value and direction vector attributes then 
        // return false
        if (!(halfedge->Attributes().Has(ATTRIB_DIRECTION_VALUE) &&
            halfedge->Attributes().Has(ATTRIB_DIRECTION_VECTOR)))
        {
            return false;
        }
    }

    // Return true since all half edges have both direction value and direction vector attributes
    return true;
}

void Tessellation::Clear()
{
    if (m_dcel) 
    {
        m_dcel->Clear();
        m_dcel = nullptr;
    }

    m_vf.clear();
}

Tessellation::TILE_SUBDIVISION_TYPE Tessellation::GetTileSubdivisionType(const std::string subdivision)
{
	return (subdivision.compare("Quadrilaterals") == 0) ? TILE_SUBDIVISION_TYPE::QUADRILATERALS :
		((subdivision.compare("Midpoints") == 0) ? TILE_SUBDIVISION_TYPE::MIDPOINTS : 
			TILE_SUBDIVISION_TYPE::TRIANGLES);
}

bool Tessellation::HasEdgeDirections() const 
{
    assert(m_dcel);

    // Traverse through the half edges of the geometric domain
    for (auto it = m_dcel->Halfedges().begin(); it != m_dcel->Halfedges().end(); ++it)
    {
        // Get the pointer to the current half edge
        std::shared_ptr<dcel::Halfedge> halfedge = *it;

        // If the current half edge does not have the direction value then return false
        if (!halfedge->Attributes().Has(ATTRIB_DIRECTION_VALUE))
        {
            return false;
        }
    }

    // Return true since all half edges have the direction value dynamic attribute
    return true;
}

void Tessellation::QuadrangulateTiles()
{
    LoadDCEL();

    SetTilesCenter();

    m_dcel->QuadrangulateFaces();

    m_vf.Set(m_dcel->vf());

    ClearDCEL();
}

void Tessellation::SetEdgeDirectionValueFromHalfedge(
    std::shared_ptr<dcel::Halfedge> halfedge,
    int direction)
{
    // Initialize a queue for the half edges and their respective direction value
    std::queue<std::tuple<std::shared_ptr<dcel::Halfedge>, int>> queue;

    // Enqueue the given half edge along with the given direction value
    queue.push(std::make_tuple(halfedge, direction));

    // Repeat while the queue is not empty
    while (!queue.empty())
    {
        // Get the next tuple in the queue
        std::tuple<std::shared_ptr<dcel::Halfedge>, int> tuple = queue.front();
        queue.pop();

        // Get the pointer to the current half edge and its respective direction value
        std::shared_ptr<dcel::Halfedge> currentHalfedge = std::get<0>(tuple);
        int dir = std::get<1>(tuple);

        // If the current half edge already has a direction value then continue with the next half 
        // edge in the queue
        if (currentHalfedge->Attributes().Has(ATTRIB_DIRECTION_VALUE))
        {
            continue;
        }

        // Set the direction value to the current half edge
        currentHalfedge->Attributes().Set<int>(ATTRIB_DIRECTION_VALUE, dir);

        // If the next half edge to the current half edge doesn't have a direction value then 
        // enqueue it
        if (!currentHalfedge->next->Attributes().Has(ATTRIB_DIRECTION_VALUE))
        {
            queue.push(std::make_tuple(currentHalfedge->next, dir * -1));
        }

        // If the previous half edge to the current half edge doesn't have a direction value then 
        // enqueue it
        if (!currentHalfedge->previous->Attributes().Has(ATTRIB_DIRECTION_VALUE))
        {
            queue.push(std::make_tuple(currentHalfedge->previous, dir * -1));
        }

        // If the twin half edge to the current half edge doesn't have a direction value then 
        // enqueue it
        if (!currentHalfedge->twin->Attributes().Has(ATTRIB_DIRECTION_VALUE))
        {
            queue.push(std::make_tuple(currentHalfedge->twin, dir * -1));
        }
    }
}

void Tessellation::SetEdgeDirectionValues(int dir)
{
    assert(m_dcel);

    // Remove the direction value attribute from all half edges
    m_dcel->RemoveHalfedgesAttribute(ATTRIB_DIRECTION_VALUE);

    // Traverse through the half edges of the geometry
    for (auto it = m_dcel->Halfedges().begin(); it != m_dcel->Halfedges().end(); ++it)
    {
        // Get the pointer to the current half edge
        std::shared_ptr<dcel::Halfedge> halfedge = *it;

        // If the current half edge doesn't have the direction attribute then start the direction 
        // setup on such half edge
        if (!halfedge->Attributes().Has(ATTRIB_DIRECTION_VALUE))
        {
            SetEdgeDirectionValueFromHalfedge(halfedge, dir);
        }
    }
}

void Tessellation::SetEdgesDirectionVector()
{
    assert(m_dcel);

    // Remove the direction vector attribute from all half edges in the geometric domain
    m_dcel->RemoveHalfedgesAttribute(ATTRIB_DIRECTION_VECTOR);

    // Traverse through the half edges and generate the direction vectors
    for (auto it = m_dcel->Halfedges().begin(); it != m_dcel->Halfedges().end(); ++it)
    {
        // Get the pointer to the current half edge
        std::shared_ptr<dcel::Halfedge> halfedge = *it;

        // If the current half edge already has the direction vector attribute then continue to the
        // next half edge
        if (halfedge->Attributes().Has(ATTRIB_DIRECTION_VECTOR))
        {
            continue;
        }

        // Declare the variables for points A, B and C
        Eigen::Vector3d A, B, C;

        // Declare the variable for the correct direction value of the vector
        int dir;

        // If the current half edge has an incident tile use C from the tile. Otherwise, use the 
        // twin's tile
        if (halfedge->face)
        {
            // Use the start and end vertices of the half edge as A and B respectively
            A << halfedge->start->Coords();
            B << halfedge->twin->start->Coords();

            // Define the direction value based on the direction attribute of the half edge. 
            // Remember that positive values point inwards the tile while negative values point 
            // outwards the tile
            assert(halfedge->Attributes().Get<int>(ATTRIB_DIRECTION_VALUE, dir));
            dir = dir > 0 ? 1 : -1;

            // Use the center of the incident tile of the half edge
            assert(halfedge->face->Attributes().Get<Eigen::Vector3d>(ATTRIB_CENTER, C));
        }
        else
        {
            // Use the start and end vertices of the twin half edge as A and B respectively
            A << halfedge->twin->start->Coords();
            B << halfedge->twin->twin->start->Coords();

            // Define the direction value based on the direction attribute of the twin half edge. 
            // Remember that positive values point inwards the tile while negative values point 
            // outwards the tile
            assert(halfedge->twin->Attributes().Get<int>(ATTRIB_DIRECTION_VALUE, dir));
            dir = dir > 0 ? 1 : -1;

            // Use the center of the incident tile of the twin half edge
            assert(halfedge->twin->Attributes().Get<Eigen::Vector3d>(ATTRIB_CENTER, C));
        }

        // Calculate the vector from A to B
        Eigen::Vector3d AB = B - A;

        // Calculate the vector from A to C
        Eigen::Vector3d AC = C - A;

        // Calculate the parameter for the point P along AB where AB and AC are orthogonal
        double t = AB.dot(AC) / AB.dot(AB);

        // Calculate P(t) = A + t(B - A)
        Eigen::Vector3d P = A + (AB * t);

        // Calculate the vector from P to C and normalize it. Multiply it by the direction value
        Eigen::Vector3d D = (C - P).normalized() * double(dir);

        // Set the ATTRIB_DIRECTION_VECTOR attribute to both half edge and its twin
        halfedge->Attributes().Set<Eigen::Vector3d>(ATTRIB_DIRECTION_VECTOR, D);
        halfedge->twin->Attributes().Set<Eigen::Vector3d>(ATTRIB_DIRECTION_VECTOR, D);
    }
}

void Tessellation::SetEdgesDirection(int dir)
{
    // Set the direction values of the edges
    SetEdgeDirectionValues(dir);

    // Set the direction vectors of the edges
    SetEdgesDirectionVector();

    // Set the midpoints of the edges
    SetEdgesMidpoint();
}

void Tessellation::SetEdgesMidpoint()
{
    assert(m_dcel);

    // Remove the midpoint attribute from the half edges of the geometry
    m_dcel->RemoveHalfedgesAttribute(ATTRIB_MIDPOINT);

    // Traverse through the half edges and set their midpoints
    for (auto it = m_dcel->Halfedges().begin(); it != m_dcel->Halfedges().end(); ++it)
    {
        // Get the pointer to the current half edge
        std::shared_ptr<dcel::Halfedge> halfedge = *it;

        // If the current half edge has the midpoint then move to the next one
        if (halfedge->Attributes().Has(ATTRIB_MIDPOINT))
        {
            continue;
        }

        // Calculate the midpoint of the half edge
        Eigen::Vector3d midpoint = halfedge->Midpoint();

        // Set the midpoint as an attribute to the half edge and its twin
        halfedge->Attributes().Set<Eigen::Vector3d>(ATTRIB_MIDPOINT, midpoint);
        halfedge->twin->Attributes().Set<Eigen::Vector3d>(ATTRIB_MIDPOINT, midpoint);
    }
}

void Tessellation::SetEdgesNormal()
{
    assert(m_dcel);

    // Remove the normal attribute from all half edges of the DCEL
    m_dcel->RemoveHalfedgesAttribute(ATTRIB_NORMAL);

    // Traverse through the half edges of the DCEL and set their normal vector as an attribute
    for (auto it = m_dcel->Halfedges().begin(); it != m_dcel->Halfedges().end(); ++it)
    {
        // Get the pointer to the current half edge
        std::shared_ptr<dcel::Halfedge> halfedge = *it;

        // If the current half edge already has the normal attribute then continue to the next half
        // edge
        if (halfedge->Attributes().Has(ATTRIB_NORMAL))
        {
            continue;
        }

        // Calculate the normalized normal vector of the half edge
        Eigen::Vector3d N = halfedge->Normal(true);

        // Set the normal vector to both half edge and its twin
        halfedge->Attributes().Set<Eigen::Vector3d>(ATTRIB_NORMAL, N);
        halfedge->twin->Attributes().Set<Eigen::Vector3d>(ATTRIB_NORMAL, N);
    }
}

void Tessellation::SetEdgeRotatedVectorsUsingAngle(double angle)
{
    assert(m_dcel);

    // Remove the rotated vector attributes from all half edges of the DCEL
    m_dcel->RemoveHalfedgesAttribute(ATTRIB_ROTATED_VECTOR);

    // Traverse through the half edges of the DCEL and calculated the rotated vectors
    for (auto it = m_dcel->Halfedges().begin(); it != m_dcel->Halfedges().end(); ++it)
    {
        // Get the pointer to the current half edge
        std::shared_ptr<dcel::Halfedge> halfedge = *it;

        // If the current half edge already has the rotated vector attribute then continue to the 
        // next half edge
        if (halfedge->Attributes().Has(ATTRIB_ROTATED_VECTOR))
        {
            continue;
        }

        // Get the direction vector of the current half edge and normalize it. It will be used as 
        // the rotation axis
        Eigen::Vector3d K = halfedge->Direction(true);

        // Get the normal vector of the half edge and normalize it
        Eigen::Vector3d V = halfedge->Normal(true);

        // Calculate the sign of the rotation. It is determined by the direction value of the half 
        // edge.
        // IMPORTANT: The Axis-Angle rotation follows the right-hand rule, so we need to invert the
        // direction of the half edge here for the correct angle value.
        int sign;
        assert(halfedge->Attributes().Get<int>(ATTRIB_DIRECTION_VALUE, sign));
        sign = sign > 0 ? -1 : 1;

        // Calculate the rotation angle in radians and with the respective sign
        double a = double(sign) * angle * utils::PI / 180.0;

        // Calculate the rotated vector
        Eigen::Vector3d R = utils::axisAngleRotation(V, K, a);

        // Set the rotated vector to the half edge and its twin
        halfedge->Attributes().Set<Eigen::Vector3d>(ATTRIB_ROTATED_VECTOR, R);
        halfedge->twin->Attributes().Set<Eigen::Vector3d>(ATTRIB_ROTATED_VECTOR, R);
    }
}

void Tessellation::SetEdgeRotatedVectorsUsingTopBottomSectionPoints(bool boundary)
{
    assert(m_dcel);

    // Remove the rotated vector attributes from all half edges of the DCEL
    m_dcel->RemoveHalfedgesAttribute(ATTRIB_ROTATED_VECTOR);

    // Traverse through the half edges of the geometric domain and calculated their rotated vectors
    for (auto it = m_dcel->Halfedges().begin(); it != m_dcel->Halfedges().end(); ++it)
    {
        // Get the pointer to the current half edge
        std::shared_ptr<dcel::Halfedge> halfedge = *it;

        // If the current half edge already has the rotated vector attribute then continue to the 
        // next half edge
        if (halfedge->Attributes().Has(ATTRIB_ROTATED_VECTOR))
        {
            continue;
        }

        // If either the current half edge or its twin doesn't have an incident tile and the pieces
        // at the boundary of the geometric domain should not be generated then move to the next 
        // half edge
        if (!(halfedge->face && halfedge->twin->face) && !boundary)
        {
            continue;
        }

        // Get the direction value of the current half edge
        int direction;
        assert(halfedge->Attributes().Get<int>(ATTRIB_DIRECTION_VALUE, direction));

        // Get the references to coordinates of the the start and end vertices of the current half 
        // edge
        const Eigen::Vector3d & A = halfedge->start->Coords();
        const Eigen::Vector3d & B = halfedge->twin->start->Coords();

        // Declare the variables for the respective points of interest from the twin tile (or the 
        // current half edge if boundary is enabled)
        Eigen::Vector3d D, E;

        // If the current half edge has an incident tile then use its top and bottom points 
        // according to the current half edge direction value. Otherwise, use the top and bottom 
        // points of the incident tile of the twin half edge
        if (halfedge->face)
        {
            // Use the bottom point of the tile if the direction of the half edge is positive; 
            // otherwise, use its top point
            if (direction > 0)
            {
                assert(halfedge->face->Attributes().Get<Eigen::Vector3d>(ATTRIB_BOTTOM_POINT, D));
            }
            else
            {
                assert(halfedge->face->Attributes().Get<Eigen::Vector3d>(ATTRIB_TOP_POINT, D));
            }
        }
        else
        {
            // Use the top point of the twin's tile if the direction of the half edge is 
            // positive; otherwise, use its bottom point
            if (direction > 0)
            {
                assert(halfedge->twin->face->Attributes().Get<Eigen::Vector3d>(ATTRIB_TOP_POINT, D));
            }
            else
            {
                assert(halfedge->twin->face->Attributes().Get<Eigen::Vector3d>(ATTRIB_BOTTOM_POINT, D));
            }
        }

        // If the twin half edge of the current half edge has an incident tile then use its top and 
        // bottom points according to the current half edge direction value. Otherwise, use the top 
        // and bottom points of the incident tile of the current half edge
        if (halfedge->twin->face)
        {
            if (direction > 0)
            {
                assert(halfedge->twin->face->Attributes().Get<Eigen::Vector3d>(ATTRIB_TOP_POINT, E));
            }
            else
            {
                assert(halfedge->twin->face->Attributes().Get<Eigen::Vector3d>(ATTRIB_BOTTOM_POINT, E));
            }
        }
        else
        {
            if (direction > 0)
            {
                assert(halfedge->face->Attributes().Get<Eigen::Vector3d>(ATTRIB_BOTTOM_POINT, E));
            }
            else
            {
                assert(halfedge->face->Attributes().Get<Eigen::Vector3d>(ATTRIB_TOP_POINT, E));
            }
        }

        // Initialize the variable for the rotated vector of the half edge. This vector is the 
        // bisector vector between the planes of interest defined by points A, B, D and E (read the
        // paper)
        Eigen::Vector3d R;

        // If the points of interest are not the same (which means there exist two tiles incident 
        // to the edge represented by the current half edge) then calculate the average direction 
        // vectors. Otherwise, ...
        if (D != E)
        {
            // If the direction is positive then D is the bottom point of the incident tile of the 
            // current half edge and E is the top point of the twin's incident tile. Otherwise, D 
            // is the top point of the incident tile of the current half edge and E is the bottom 
            // point of the twin's incident tile. Calculate the two normals corresponding to the 
            // planes of interest (read paper for further details) and calculate their average.
            // NOTE: Here I have some redundant code but I don't want to factorize it since the 
            // order the vectors are defined on each case helps the reader to grasp the geometric 
            // intuition concerning the directions of the normal vectors from the planes of 
            // interest
            if (direction > 0)
            {
                // 
                Eigen::Vector3d DA = (A - D).normalized();
                Eigen::Vector3d DB = (B - D).normalized();

                // Calculate the normal vector of the plane for points A, D and B
                Eigen::Vector3d N1 = DA.cross(DB).normalized();

                // 
                Eigen::Vector3d EB = (B - E).normalized();
                Eigen::Vector3d EA = (A - E).normalized();

                // Calculate the normal vector of the plane for points A, B and E
                Eigen::Vector3d N2 = EB.cross(EA).normalized();

                // Calculate the average between both normal vectors and normalize it
                R = ((N1 + N2) / 2.0).normalized();
            }
            else
            {
                // 
                Eigen::Vector3d EB = (B - E).normalized();
                Eigen::Vector3d EA = (A - E).normalized();

                // Calculate the normal vector of the plane for points A, E and B
                Eigen::Vector3d N1 = EB.cross(EA).normalized();

                // 
                Eigen::Vector3d DA = (A - D).normalized();
                Eigen::Vector3d DB = (B - D).normalized();

                // Calculate the normal vector of the plane for points A, B and D
                Eigen::Vector3d N2 = DA.cross(DB).normalized();

                // Calculate the average between both normal vectors and normalize it
                R = ((N1 + N2) / 2.0).normalized();
            }
        }
        else
        {
            // If the direction is positive then both D and E are the bottom point of the incident 
            // tile of the current half edge. Otherwise, both D and E are the top point of the 
            // incident tile of the current half edge. Calculate the normal vector of the only 
            // plane of interest
            Eigen::Vector3d DA = (A - D).normalized();
            Eigen::Vector3d DB = (B - D).normalized();

            // Calculate the normal vector of the plane for points A, D and B and normalize it
            R = DA.cross(DB).normalized();
        }

        // Set the rotated vector to the half edge and its twin
        halfedge->Attributes().Set<Eigen::Vector3d>(ATTRIB_ROTATED_VECTOR, R);
        halfedge->twin->Attributes().Set<Eigen::Vector3d>(ATTRIB_ROTATED_VECTOR, R);
    }
}

void Tessellation::SetTileBoundaries()
{
    assert(m_dcel);

    std::shared_ptr<dcel::Face> tile = nullptr;

    // Traverse through the tiles of the geometric domain
    for (auto it = m_dcel->Faces().begin(); it != m_dcel->Faces().end(); ++it)
    {
        // Get the pointer to the current tile
        tile = *it;

        // Check if the tile has all of its neighbors. If it has all of its neighbors then it is 
        // not at the boundary of the geometric domain; otherwise, it does. Store the indicator in 
        // the ATTRIB_BOUNDARY dynamic attribute of the tile
        tile->Attributes().Set<bool>(ATTRIB_BOUNDARY, !tile->HasAllNeighbors());
    }
}

void Tessellation::SetTilesCenter() 
{
    assert(m_dcel);

    std::shared_ptr<dcel::Face> tile = nullptr;

    // Traverse through the tiles of the geometry
    for (auto it = m_dcel->Faces().begin(); it != m_dcel->Faces().end(); ++it)
    {
        // Get the pointer to the current tile
        tile = *it;

        // Set the center of the tile as a dynamic attribute
        tile->Attributes().Set<Eigen::Vector3d>(ATTRIB_CENTER, tile->Centroid());
    }
}

void Tessellation::SetTileTopBottomSectionPoints(double height)
{
    assert(m_dcel);

    std::shared_ptr<dcel::Face> tile = nullptr;

    // Traverse through the tiles of the geometry
    for (auto it = m_dcel->Faces().begin(); it != m_dcel->Faces().end(); ++it)
    {
        // Get the pointer to the current tile
        tile = *it;

        // Get the normalized normal vector of the tile. Then, multiply it by the height value
        Eigen::Vector3d N = tile->Normal(true) * height;

        // Get the center point associated to the tile
        Eigen::Vector3d C;
        assert(tile->Attributes().Get<Eigen::Vector3d>(ATTRIB_CENTER, C));

        // Calculate the top and bottom points of the tile
        Eigen::Vector3d T = C + N;
        Eigen::Vector3d B = C - N;

        // Store the top and bottom points in the ATTRIB_TOP_POINT and ATTRIB_BOTTOM_POINT dynamic 
        // attributes of the tile
        tile->Attributes().Set<Eigen::Vector3d>(ATTRIB_TOP_POINT, T);
        tile->Attributes().Set<Eigen::Vector3d>(ATTRIB_BOTTOM_POINT, B);
    }
}

void Tessellation::SubdivideTilesByMidpoints()
{
    LoadDCEL();

    m_dcel->SubdivideFacesByMidpoints();

    m_vf.Set(m_dcel->vf());

    ClearDCEL();
}

void Tessellation::TriangulateTilesByVertices()
{
    m_vf.Set(m_vf.TriangulateFacesByVertices());

    ClearDCEL();
}

void Tessellation::Write() const
{
    std::cout << "Tessellation" << std::endl;
    m_vf.Write();
    std::cout << "/Tessellation" << std::endl;
}

void Tessellation::Write(std::stringstream & ss) const 
{
    ss << "Tessellation" << std::endl;
    m_vf.Write(ss);
    ss << "/Tessellation" << std::endl;
}