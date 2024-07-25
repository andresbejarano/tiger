#include <tiger/toolkit/algorithms.h>
#include <tiger/ds/algorithms.h>
#include <tiger/utils.h>

#include <Eigen/Geometry>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <numeric>

bool algorithms::areFaceVerticesCoplanar(const std::shared_ptr<dcel::Face> face)
{
    // There must be a face
    assert(face);

    // Get the pointer to the incident half edge of the face
    std::shared_ptr<dcel::Halfedge> halfedge = face->halfedge;

    bool coplanar = false;

    // Traverse through the half edges of the face
    do
    {
        // Check if the start vertex of the current half edge has the 
        // ATTRIB_COPLANAR attribute
        assert(halfedge->start->Attributes().Get<bool>(ATTRIB_COPLANAR, coplanar));

        // If the coplanar attribute is not true then return false
        if (!coplanar)
        {
            return false;
        }

        // Move to the next half edge
        halfedge = halfedge->next;

    } while (halfedge != face->halfedge);

    // Return true since all vertices in the face have the coplanar attribute
    return true;
}

bool algorithms::faceRayIntersection(
    const std::shared_ptr<dcel::Face> face,
    const toolkit::Ray & ray,
    double & t,
    double threshold)
{
    // Exit the function if no face
    if (!face)
    {
        return false;
    }

    Eigen::Vector3d P = Eigen::Vector3d::Zero();

    // Get the plane of the face
    toolkit::Plane plane = face->Plane();

    // If the plane and the ray intersect then determine if the intersection 
    // point lies within the face
    if (planeRayIntersection(plane, ray, t, threshold))
    {
        // Get the intersection point between the plane and the ray. Then, fix 
        // its zeros
        P << ray.at(t);
        utils::fixZeros(P, threshold);

        // Indicate if the point is within the face (including the edges)
        return face->IsPointIn(P, threshold);
    }

    // Return false since the plane and the ray do not intersect
    return false;
}

bool algorithms::getIntersectionPointsFromPlanes(
    const std::vector<toolkit::Plane> & planes,
    std::vector<Eigen::Vector3d> & points)
{
    // Get the number of planes
    size_t nPlanes = planes.size(), _nPlanes = nPlanes - 1, i = 0, prevIndex = 0, currIndex = 0, 
        nextIndex = 0;

    Eigen::Vector3d P = Eigen::Vector3d::Zero();

    bool isP = false;

    // Clear the vector for the intersection points between the planes
    points.clear();
    points.resize(nPlanes);

    // Traverse through the planes and calculate their intersection points
    for (i = 0; i < nPlanes; i += 1)
    {
        // Get the indices to the previous plane, the current plane and the 
        // next plane
        prevIndex = (i == 0) ? _nPlanes : i - 1;
        currIndex = i;
        nextIndex = (i == _nPlanes) ? 0 : i + 1;

        // Calculate the intersection point between the previous, current and 
        // next planes
        isP = algorithms::threePlanesIntersection(
            planes.at(prevIndex),
            planes[currIndex],
            planes[nextIndex], P);

        // If no intersection point then return false
        if (!isP)
        {
            points.clear();
            return false;
        }

        // Insert the intersection point in the points vector
        points[i] << P;
    }

    // Return the vector with the intersection points between the planes
    return true;
}

bool algorithms::intersectCoplanarFaces(
    const VF & vf1, 
    size_t face1, 
    const VF & vf2, 
    size_t face2, 
    std::vector<Eigen::Vector3d> & points, 
    double threshold)
{
    // Clear the vector where the intersection points will be stored
    points.clear();

    // Get the references to the vectors with the vertex indices of each face
    const std::vector<size_t> & face1Indices = vf1.face(face1), & face2Indices = vf2.face(face2);

    // Get the number of vertices from each face
    size_t nVertices1 = face1Indices.size(), nVertices2 = face2Indices.size(), e1 = 0, e2 = 0;

    // 
    Eigen::Vector3d P = Eigen::Vector3d::Zero();

    // 
    toolkit::LineSegment edge1, edge2;

    double s = 0, t = 0;

    bool intersect = false;

    // Traverse through the edges from the first face
    for (e1 = 0; e1 < nVertices1; e1 += 1) 
    {
        // Get the end points of the current edge in the first face
        edge1.Set(
            vf1.Vertex(face1Indices[e1]), 
            vf1.Vertex(face1Indices[(e1 == nVertices1 - 1) ? 0 : e1 + 1]));

        // Step 1: Check if the start vertex of the current edge from the first
        // face is within the second face. If so, make a copy of the point and 
        // insert it into the intersection points set
        if (vf2.IsPointIn(face2, edge1.A(), threshold)) 
        {
            points.emplace_back(edge1.A());
        }

        // Traverse through the edges from the second face
        for (e2 = 0; e2 < nVertices2; e2 += 1) 
        {
            // Get the end points of the current edge in the second face
            edge2.Set(
                vf2.Vertex(face2Indices[e2]), 
                vf2.Vertex(face2Indices[(e2 == nVertices2 - 1) ? 0 : e2 + 1]));

            // 
            intersect = edge1.Intersect(edge2, P, s, t, threshold);

            // Step 2: Check if the current edge from the first face intersects
            // with the current edge from the second face. If so, check if the 
            // intersection parameters s and t lie between 0 and 1, then insert
            // the intersection point into the points set
            if (intersect && s >= 0.0 && s <= 1.0 && t >= 0.0 && t <= 1.0) 
            {
                points.emplace_back(P);
            }
        }

    }

    // Indicate at least one intersection point was found
    return points.size() > 0;
}

bool algorithms::intersectCoplanarFaces(
    std::shared_ptr<dcel::Face> F1,
    std::shared_ptr<dcel::Face> F2,
    std::vector<Eigen::Vector3d> & points,
    double threshold)
{
    // 
    assert(F1 && F2);

    // Clear the vector where the intersection points will be stored
    points.clear();

    // Get the pointer to the incident half edge from face F1
    std::shared_ptr<dcel::Halfedge> currentF1halfedge = F1->halfedge, currentF2Halfedge = nullptr;

    double s = 0, t = 0;
    bool intersect = false;

    Eigen::Vector3d point = Eigen::Vector3d::Zero();

    // Traverse through the half edges of face F1
    do
    {
        // Step 1: Check if the start vertex of the current half edge from face
        // F1 is within face F2. If so, make a copy of the point and insert it 
        // into the intersection points set
        if (F2->IsPointIn(currentF1halfedge->start->Coords(), threshold))
        {
            points.emplace_back(currentF1halfedge->start->Coords());
        }

        // Get the incident half edge from face F2
        currentF2Halfedge = F2->halfedge;

        // Traverse through the half edges of face F2
        do
        {
            // Step 2: Check if the current half edge in face F1 intersects 
            // with the current half edge in face F2. If so, check if the 
            // intersection parameters s and t lie between 0 and 1, then insert
            // the intersection point into the intersection points set
            intersect = currentF1halfedge->Intersect(currentF2Halfedge, point, s, t, threshold);

            // 
            if (intersect && s >= 0.0 && s <= 1.0 && t >= 0.0 && t <= 1.0)
            {
                points.emplace_back(point);
            }

            // Move to the next half edge in face F2
            currentF2Halfedge = currentF2Halfedge->next;

        } while (currentF2Halfedge != F2->halfedge);

        // Move to the next half edge in face F1
        currentF1halfedge = currentF1halfedge->next;

    } while (currentF1halfedge != F1->halfedge);

    // Indicate at least one intersection point was found
    return points.size() > 0;
}

int algorithms::isLeftTurn(
    const Eigen::Vector3d & P1, 
    const Eigen::Vector3d & P2, 
    const toolkit::Plane & plane, 
    const double threshold)
{
    double dot = (P2 - P1).cross(plane.Point() - P1).dot(plane.Normal());
    return (abs(dot) <= threshold) ? 0 : ((dot > 0) ? 1 : -1);
}

bool algorithms::planeLineSegmentIntersection(
    const toolkit::Plane & plane,
    const toolkit::LineSegment & linesegment,
    double & t,
    double threshold)
{
    // Generate a ray using the line segment, then calculate its intersection 
    // with the plane
    const toolkit::Ray ray = toolkit::Ray(linesegment.A(), linesegment.Direction());
    return algorithms::planeRayIntersection(plane, ray, t, threshold);
}

bool algorithms::planeRayIntersection(
    const toolkit::Plane & plane,
    const toolkit::Ray & ray,
    double & t,
    double threshold)
{
    // Calculate the dot product between the direction vector of the ray and 
    // the normal vector of the plane
    double denom = ray.direction().dot(plane.Normal());

    // If the absolute value of the dot product is zero then return false since
    // there is no intersection between the plane and the ray
    if (abs(denom) <= threshold)
    {
        return false;
    }

    // 
    double num = (plane.Point() - ray.start()).dot(plane.Normal());

    // Calculate the parameter value along the ray for the intersection point 
    // between the plane and the ray. Then return true since there is an 
    // intersection point
    t = num / denom;
    return true;
}

bool algorithms::sortPointsCCW(
    std::vector<Eigen::Vector3d> & points,
    const Eigen::Vector3d & N,
    const Eigen::Vector3d & Q,
    std::vector<Eigen::Vector3d> & sorted,
    double threshold)
{
    // Clear the content of the vector where the sorted points will be stored
    sorted.clear();

    // Get the number of points to be sorted
    size_t nPoints = points.size(), i = 0, j = 0, k = 0;
    assert(nPoints > 0);

    // Initialize the adjacency matrix between the points. This matrix will 
    // indicate whether two points i and j are in CCW with respect of the 
    // normal vector and the given point. A 1 indicates CCW order.
    Eigen::MatrixXd M(nPoints, nPoints);

    double dot = 0, r = 0;

    // Populate the adjacency matrix
    for (i = 0; i < nPoints; i += 1)
    {
        // Set the diagonal value with -1 (we do not allow loops)
        M(i, i) = -1.0;

        // Traverse through the remaining points in the vector
        for (j = i + 1; j < nPoints; j += 1)
        {
            // Calculate N * ((Vj - Vi) x (Q - Vi)), where * is dot product and
            // x is cross product. This value indicates if a vector from point 
            // i to point j goes in CCW direction with respect of the normal 
            // vector N and point Q (the resultant vector and N must point at
            // the same half space)
            dot = (points[j] - points[i]).cross(Q - points[i]).dot(N);

            // Check the value of the dot product. If it is zero (or close to 
            // zero) then the CCW between both points is not possible in both 
            // directions (e.g., points i, j and Q are collinear). Otherwise, 
            // process it as regular
            if (abs(dot) < threshold)
            {
                M(i, j) = -1.0;
                M(j, i) = -1.0;
            }
            else
            {
                // Determine the direction between both points based on the dot
                // product value. If it is possible to move from point i to 
                // point j in CCW then the adjacency between both points is 
                // valid; otherwise, it is not (and viceversa). Set these 
                // values in the adjacency matrix
                r = (dot > 0.0) ? 1.0 : -1.0;
                M(i, j) = r;
                M(j, i) = -r;
            }
        }
    }

    // Let's refine the adjacency matrix so we have only valid edges. It is the
    // same logic as before but rather between points in valid directions 
    // instead of comparing with respect of point Q
    for (i = 0; i < nPoints; i += 1)
    {
        for (j = 0; j < nPoints; j += 1)
        {
            // Avoid invalid directions
            if (M(i, j) == -1.0)
            {
                continue;
            }

            // 
            for (k = j + 1; k < nPoints; k += 1)
            {
                // Avoid invalid directions
                if (M(i, k) == -1.0)
                {
                    continue;
                }

                // 
                dot = (points[j] - points[i]).cross(points[k] - points[i]).dot(N);

                // 
                if (dot < 0.0)
                {
                    M(i, j) = -1.0;
                    break;
                }
                else if (dot > 0.0)
                {
                    M(i, k) = -1.0;
                }
            }
        }
    }

    // Initialize the vector indicating which points are visited
    std::vector<bool> visited(nPoints, false);

    // Initialize the vector for storing the indices of the path
    std::vector<size_t> path;

    // If there is a path through the points using M such that all points are 
    // visited then store the path in S
    if (utils::findPath(M, visited, 0, path))
    {
        // Traverse the points and clone them in the correct ccw order
        for (i = 0; i < nPoints; i += 1)
        {
            sorted.push_back(points[path[i]]);
        }
    }

    // Indicate if there are sorted points
    return sorted.size() > 0;
}

bool algorithms::sortPointsCCW(
    std::vector<Eigen::Vector3d> & points,
    toolkit::Plane & plane,
    std::vector<Eigen::Vector3d> & sorted,
    double threshold)
{
    return algorithms::sortPointsCCW(points, plane.Normal(), plane.Point(), sorted, threshold);
}

bool algorithms::threePlanesIntersection(
    const toolkit::Plane & A,
    const toolkit::Plane & B,
    const toolkit::Plane & C, 
    Eigen::Vector3d & P)
{
    // Get the components from plane A
    double Aa = A.a();
    double Ab = A.b();
    double Ac = A.c();
    double Ad = A.d();

    // Get the components from plane B
    double Ba = B.a();
    double Bb = B.b();
    double Bc = B.c();
    double Bd = B.d();

    // Get the components from plane C
    double Ca = C.a();
    double Cb = C.b();
    double Cc = C.c();
    double Cd = C.d();

    // Get the determinant of the system
    double D = utils::det(Aa, Ab, Ac, Ba, Bb, Bc, Ca, Cb, Cc);

    // If the determinant is zero return false since there is no intersection 
    // point between the planes
    if (abs(D) <= 1e-8)
    {
        return false;
    }

    // Use the Cramer's rule for obtaining the coordinates of the intersection 
    // point
    double D1 = utils::det(Ad, Ab, Ac, Bd, Bb, Bc, Cd, Cb, Cc);
    double D2 = utils::det(Aa, Ad, Ac, Ba, Bd, Bc, Ca, Cd, Cc);
    double D3 = utils::det(Aa, Ab, Ad, Ba, Bb, Bd, Ca, Cb, Cd);

    // Set the coordinates of the intersection point and return true
    P << D1 / D, D2 / D, D3 / D;
    return true;
}

void algorithms::write(const std::vector<dcel::DCEL> & V)
{
    for (auto it = V.begin(); it != V.end(); ++it)
    {
        (*it).Write();
    }
}

void algorithms::writeGeogebraJS(
    const std::vector<dcel::DCEL> & blocks,
    const std::vector<dcel::DCEL> & interfaces,
    const std::string fileprefix,
    bool writeReference,
    const Eigen::Vector3d & O,
    const Eigen::Vector3d & X,
    const Eigen::Vector3d & Y,
    const Eigen::Vector3d & Z)
{
    std::stringstream ss;
    ss << fileprefix << ".js";

    // Open the file
    std::ofstream file(ss.str(), std::ios::out);

    // Get the number of blocks
    size_t nBlocks = blocks.size();

    // Traverse through the blocks and write their geometric information
    for (size_t bIdx = 0; bIdx < nBlocks; bIdx += 1)
    {
        // Get the number of vertices of the current block
        size_t nVertices = blocks[bIdx].Vertices().size();

        // Traverse through the vertices of the current block and write the 
        // Geogebra commands that define them
        for (size_t vIdx = 0; vIdx < nVertices; vIdx += 1)
        {
            // Get a copy of the coordinates of the current vertex of the 
            // block. Then, fix its zeros
            Eigen::Vector3d V = blocks[bIdx].Vertices()[vIdx]->Coords();
            utils::fixZeros(V);

            // Write the command that defines the current vertex
            file << "ggbApplet.evalCommand(\"b" << bIdx << "v" << vIdx << " = Point({" <<
                V.x() << ", " << V.y() << ", " << V.z() << "})\");" << std::endl;

            // Set the index of the attribute in the ATTRIB_INDEX dynamic 
            // attribute
            blocks[bIdx].Vertices()[vIdx]->Attributes().Set<size_t>(ATTRIB_INDEX, vIdx);
        }

        // Get the coordinates of the centroid of the current block. Then, fix 
        // its zeros
        Eigen::Vector3d C = blocks[bIdx].Centroid();
        utils::fixZeros(C);

        // Write the commands that define the centroid of the current block, 
        // size and color
        file << "ggbApplet.evalCommand(\"b" << bIdx << "centroid = Point({" << C.x() << ", " <<
            C.y() << ", " << C.z() << "})\");" << std::endl;
        file << "ggbApplet.setColor(\"b" << bIdx << "centroid\", 0, 255, 255);" << std::endl;
        file << "ggbApplet.setPointSize(\"b" << bIdx << "centroid\", 9);" << std::endl;

        // Write the command that defines the volume of the current block
        file << "ggbApplet.evalCommand(\"b" << bIdx << "volume = " << blocks[bIdx].Volume() <<
            "\");" << std::endl;

        // Get the number of faces of the current block
        size_t nFaces = blocks[bIdx].Faces().size();

        // Traverse through the faces of the current block and generate the 
        // Geogebra commands that define them
        for (size_t fIdx = 0; fIdx < nFaces; fIdx += 1)
        {
            // Get the pointer to the incident half edge of the face
            std::shared_ptr<dcel::Halfedge> halfedge = blocks[bIdx].Faces()[fIdx]->halfedge;

            // Start the line for the current interface face
            file << "ggbApplet.evalCommand(\"b" << bIdx << "f" << fIdx << " = Polygon(";

            // 
            bool firstVertex = true;

            // Traverse through the half edges of the face
            do
            {
                // Get the index of the start vertex of the half edge
                size_t vIdx;
                assert(halfedge->start->Attributes().Get<size_t>(ATTRIB_INDEX, vIdx));

                if (firstVertex)
                {
                    file << "b" << bIdx << "v" << vIdx;

                    firstVertex = false;
                }
                else
                {
                    file << ", b" << bIdx << "v" << vIdx;
                }

                // Move to the next half edge
                halfedge = halfedge->next;

            } while (halfedge != blocks[bIdx].Faces()[fIdx]->halfedge);

            // Close the line for the current interface face
            file << ")\");" << std::endl;
        }

        // Remove the ATTRIB_INDEX dynamic attribute from the vertices
        blocks[bIdx].RemoveVerticesAttribute(ATTRIB_INDEX);
    }

    // Get the number of interface polygons
    size_t nInterfaces = interfaces.size();

    // Traverse through the interface polygons and write their geometric 
    // information
    for (size_t iIdx = 0; iIdx < nInterfaces; iIdx += 1)
    {
        // Get the number of vertices of the current interface polygon
        size_t nVertices = interfaces[iIdx].Vertices().size();

        // Traverse through the vertices of the current interface polygon and 
        // write the Geogebra commands that define them
        for (size_t vIdx = 0; vIdx < nVertices; vIdx += 1)
        {
            // Get a copy of the coordinates of the current vertex of the 
            // interface polygon. Then, fix its zeros
            Eigen::Vector3d V = interfaces[iIdx].Vertices()[vIdx]->Coords();
            utils::fixZeros(V);

            // Write the commands that define the current vertex, size and color
            file << "ggbApplet.evalCommand(\"i" << iIdx << "v" << vIdx << " = Point({" <<
                V.x() << ", " << V.y() << ", " << V.z() << "})\");" << std::endl;
            file << "ggbApplet.setColor(\"i" << iIdx << "v" << vIdx << "\", 255, 0, 255);" <<
                std::endl;
            file << "ggbApplet.setPointSize(\"i" << iIdx << "v" << vIdx << "\", 9);" << std::endl;

            // Set the index of the attribute in the ATTRIB_INDEX dynamic 
            // attribute
            interfaces[iIdx].Vertices()[vIdx]->Attributes().Set<size_t>(ATTRIB_INDEX, vIdx);
        }

        // Get the coordinates of the centroid of the current interface 
        // polygon. Then, fix its zeros
        Eigen::Vector3d C = interfaces[iIdx].Faces()[0]->Centroid();
        utils::fixZeros(C);

        // Write the command that defines the centroid of the current interface
        // polygon
        file << "ggbApplet.evalCommand(\"i" << iIdx << "center = Point({" << C.x() << ", " <<
            C.y() << ", " << C.z() << "})\");" << std::endl;

        // Get the pointer to the incident half edge of the face
        std::shared_ptr<dcel::Halfedge> halfedge = interfaces[iIdx].Faces()[0]->halfedge;

        // Start the line for the current interface face
        file << "ggbApplet.evalCommand(\"i" << iIdx << "f = Polygon(";

        // 
        bool firstVertex = true;

        // Traverse through the half edges of the face
        do
        {
            // Get the index of the start vertex of the half edge
            size_t vIdx;
            assert(halfedge->start->Attributes().Get<size_t>(ATTRIB_INDEX, vIdx));

            if (firstVertex)
            {
                file << "i" << iIdx << "v" << vIdx;

                firstVertex = false;
            }
            else
            {
                file << ", i" << iIdx << "v" << vIdx;
            }

            // Move to the next half edge
            halfedge = halfedge->next;

        } while (halfedge != interfaces[iIdx].Faces()[0]->halfedge);

        // Close the line for the current interface face
        file << ")\");" << std::endl;

        // Color the current interface polygon
        file << "ggbApplet.setColor(\"i" << iIdx << "f\", 255, 0, 255);" << std::endl;

        // Get the normal and tangential vectors of the current interface 
        // polygon. Then, fix their zeros
        Eigen::Vector3d N = interfaces[iIdx].Faces()[0]->Normal(true);
        Eigen::Vector3d U = interfaces[iIdx].Faces()[0]->halfedge->Direction(true);
        Eigen::Vector3d V = N.cross(U).normalized();
        utils::fixZeros(N);
        utils::fixZeros(U);
        utils::fixZeros(V);

        // Write the commands that define the normal and tangential vectors of 
        // the current interface polygon
        file << "ggbApplet.evalCommand(\"i" << iIdx << "N = UnitVector(Segment(i" << iIdx <<
            "center, i" << iIdx << "center + (" << N.x() << ", " << N.y() << ", " << N.z() <<
            ")))\");" << std::endl;
        file << "ggbApplet.setColor(\"i" << iIdx << "N\", 0, 255, 255);" << std::endl;

        file << "ggbApplet.evalCommand(\"i" << iIdx << "U = UnitVector(Segment(i" << iIdx <<
            "center, i" << iIdx << "center + (" << U.x() << ", " << U.y() << ", " << U.z() <<
            ")))\");" << std::endl;
        file << "ggbApplet.setColor(\"i" << iIdx << "U\", 255, 0, 255);" << std::endl;

        file << "ggbApplet.evalCommand(\"i" << iIdx << "V = UnitVector(Segment(i" << iIdx <<
            "center, i" << iIdx << "center + (" << V.x() << ", " << V.y() << ", " << V.z() <<
            ")))\");" << std::endl;
        file << "ggbApplet.setColor(\"i" << iIdx << "V\", 255, 255, 0);" << std::endl;

        // Remove the ATTRIB_INDEX dynamic attribute from the vertices
        interfaces[iIdx].RemoveVerticesAttribute(ATTRIB_INDEX);
    }

    // 
    if (writeReference)
    {
        Eigen::Vector3d origin(O);
        Eigen::Vector3d xref(X);
        Eigen::Vector3d yref(Y);
        Eigen::Vector3d zref(Z);
        utils::fixZeros(origin);
        utils::fixZeros(xref);
        utils::fixZeros(yref);
        utils::fixZeros(zref);

        // Write the command that defines the centroid of the current interface
        // polygon
        file << "ggbApplet.evalCommand(\"O = Point({" << origin.x() << ", " << origin.y() <<
            ", " << origin.z() << "})\");" << std::endl;

        file << "ggbApplet.evalCommand(\"X = UnitVector(Segment(O, O + (" << xref.x() << ", " <<
            xref.y() << ", " << xref.z() << ")))\");" << std::endl;
        file << "ggbApplet.setColor(\"X\", 255, 0, 0);" << std::endl;

        file << "ggbApplet.evalCommand(\"Y = UnitVector(Segment(O, O + (" << yref.x() << ", " <<
            yref.y() << ", " << yref.z() << ")))\");" << std::endl;
        file << "ggbApplet.setColor(\"Y\", 0, 255, 0);" << std::endl;

        file << "ggbApplet.evalCommand(\"Z = UnitVector(Segment(O, O + (" << zref.x() << ", " <<
            zref.y() << ", " << zref.z() << ")))\");" << std::endl;
        file << "ggbApplet.setColor(\"Z\", 0, 0, 255);" << std::endl;

        // Write the command that defines the reference plane
        file << "ggbApplet.evalCommand(\"PerpendicularPlane(O, Y)\");" << std::endl;
    }

    // Close the file
    file.close();
}
