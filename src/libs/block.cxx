#include <tiger/tic/block.h>
#include <tiger/toolkit/ray.h>
#include <tiger/utils.h>
#include <tiger/ds/algorithms.h>
#include <tiger/toolkit/algorithms.h>
#include <iostream>

Block::Block() : 
    m_enabled(true), 
    m_faceIndex(0), 
    m_loads({ 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 })
{
}

Block::Block(const VF & vf, size_t faceIndex) : 
    GeometricClass(vf), 
    m_enabled(true), 
    m_faceIndex(faceIndex), 
    m_loads({ 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 })
{
}

Block::~Block() 
{
    Clear();
}

void Block::AddForceLoad(double x, double y, double z)
{
    m_loads[0] += x;
    m_loads[1] += y;
    m_loads[2] += z;
}

void Block::AddTorqueLoad(double x, double y, double z)
{
    m_loads[3] += x;
    m_loads[4] += y;
    m_loads[5] += z;
}

double Block::CentralLength(
    const std::shared_ptr<Tessellation> tessellation, 
    int direction, 
    double threshold) const
{
    assert(tessellation);

    double sign = direction >= 0 ? 1.0 : -1.0;

    // Get the centroid and normal vector for the ray
    Eigen::Vector3d C = tessellation->DCEL()->Faces()[m_faceIndex]->Centroid();
    Eigen::Vector3d N = tessellation->DCEL()->Faces()[m_faceIndex]->Normal() * sign;

    // Calculate the intersection point between the ray and the geometry of the block
    double t;
    Eigen::Vector3d P = Intersect(toolkit::Ray(C, N), t, threshold);

    // Return the distance from the intersection point and the centroid
    return (P - C).norm();
}

void Block::Clear()
{
    m_attributes.Clear();
}

std::vector<bool> Block::CoplanarVertices(
    const Eigen::Vector3d & N,
    const Eigen::Vector3d & P,
    double threshold) const
{
    return CoplanarVertices(toolkit::Plane(N, P), threshold);
}

std::vector<bool> Block::CoplanarVertices(const toolkit::Plane & plane, double threshold) const
{
    size_t nVertices = m_vf.countVertices(), i;

    std::vector<bool> coplanar(nVertices, false);

    for (i = 0; i < nVertices; i += 1)
    {
        coplanar[i] = plane.IsPointInPlane(m_vf.Vertex(i), threshold);
    }

    return coplanar;
}

void Block::Disable()
{
    m_enabled = false;
}

void Block::DisableIfIntersects(const toolkit::Plane & plane)
{
    size_t nFaces = m_vf.countFaces(), nVertices = 0, fIdx, vIdx, jIdx;

    Eigen::Vector3d V0 = Eigen::Vector3d::Zero(), V1 = Eigen::Vector3d::Zero();

    for (fIdx = 0; fIdx < nFaces; fIdx += 1) 
    {
        const std::vector<size_t> & face = m_vf.face(fIdx);

        nVertices = face.size();

        for (vIdx = 0; vIdx < nVertices; vIdx += 1) 
        {
            jIdx = (vIdx == nVertices - 1) ? face[0] : face[vIdx + 1];

            V0 << m_vf.Vertex(face[vIdx]);
            V1 << m_vf.Vertex(jIdx);

            // If the edge vertices lie at different half spaces then the edge intersects the plane
            if (plane.PointLocation(V0) * plane.PointLocation(V1) == -1) 
            {
                m_enabled = false;
                return;
            }
        }
    }
}

void Block::Enable()
{
    m_enabled = true;
}

size_t Block::FaceIndex() const
{
    return m_faceIndex;
}

size_t Block::GetCoplanarFace(const toolkit::Plane & plane, double threshold) const
{
    size_t nFaces = m_vf.countFaces(), i;

    for (i = 0; i < nFaces; i += 1)
    {
        if (m_vf.IsCoplanar(i, plane, threshold))
        {
            return i;
        }
    }

    return nFaces;
}

VF Block::InterfacePolygon(const std::shared_ptr<Block> block, double threshold) const
{
    // Initialize the lists to store the vertices and faces of the interface polygon
    std::list<Eigen::Vector3d> interfaceVertices;
    std::list<std::vector<size_t>> interfaceFace;

    // 
    InterfacePolygonElements(block, interfaceVertices, interfaceFace, threshold);

    // 
    VF vf(interfaceVertices, interfaceFace);

    // 
    interfaceVertices.clear();
    interfaceFace.clear();

    return vf;
}

void Block::InterfacePolygonElements(
    const std::shared_ptr<Block> block,
    std::list<Eigen::Vector3d> & interfaceVertices,
    std::list<std::vector<size_t>> & interfaceFace,
    double threshold) const
{
    // 
    assert(block);

    // 
    interfaceVertices.clear();
    interfaceFace.clear();

    size_t nBlock1Faces = m_vf.countFaces(), otherFaceIdx = 0;

    //std::shared_ptr<dcel::Face> thisFace = nullptr;

    toolkit::Plane plane;

    // Traverse through the faces of the geometry of this block
    for(size_t block1FaceIdx = 0; block1FaceIdx < nBlock1Faces; block1FaceIdx += 1)
    {
        // Get the plane where the current face lies. Remember, the plane is defined by the normal 
        // vector and barycenter of the face
        plane = m_vf.Plane(block1FaceIdx, true, true, threshold);

        // Get the pointer to the face in the given block that is coplanar with the plane
        otherFaceIdx = block->GetCoplanarFace(plane, threshold);

        // If no face in the given block is coplanar with the plane then continue to the next face 
        // of the geometry from this block
        if (otherFaceIdx == block->Geometry().countFaces())
        {
            continue;
        }

        // Initialize the vectors for storing the intersection points between the faces
        std::vector<Eigen::Vector3d> V1, V2;

        algorithms::intersectCoplanarFaces(m_vf, block1FaceIdx, block->Geometry(), otherFaceIdx, V1, threshold);
        algorithms::intersectCoplanarFaces(block->Geometry(), otherFaceIdx, m_vf, block1FaceIdx, V2, threshold);

        // If at least one of the points vectors is empty then continue with the next face in this
        // block
        if (V1.size() == 0 || V2.size() == 0)
        {
            continue;
        }

        // Merge the vectors with the intersection points (this may generate duplicated points)
        V1.insert(V1.end(), V2.begin(), V2.end());

        // Fix the zeros of the intersection points (some points may have values very close to zero
        // due to numerical precision)
        utils::fixZeros(V1, threshold);

        // Get a vector with unique intersection points. If there is none then return a null 
        // pointer since the interace polygon cannot be calculated
        std::vector<Eigen::Vector3d> uniquePoints;
        if (!utils::getUniquePoints(V1, uniquePoints, threshold))
        {
            return;
        }

        // Initialize a vector for storing the middle point between the unique intersection points.
        // This point will be used for calculating their CCW order
        Eigen::Vector3d uniqueCenter;
        uniqueCenter << 0.0, 0.0, 0.0;

        // traverse through the unique points and add them up
        for (auto uit = uniquePoints.begin(); uit != uniquePoints.end(); ++uit)
        {
            uniqueCenter += (*uit);
        }

        // Calculate the average of the unique points and set it as the reference point of the 
        // plane
        uniqueCenter /= (double)uniquePoints.size();
        plane.Point(uniqueCenter);

        // Initialize the vertex coordinates and vertex indices representing the geometry of the 
        // interface polygon
        //VF vf;

        // 
        std::vector<Eigen::Vector3d> ccwPoints;

        // Sort the interface polygon points in CCW with respect of the plane. If sorting is not 
        // possible then return a null pointer since the interface polygon cannot be calculated
        if (!algorithms::sortPointsCCW(uniquePoints, plane, ccwPoints, threshold))
        {
            return;
        }

        // Get the number of vertices of the interface polygon
        size_t nVertices = ccwPoints.size();

        // Initialize the vertex indices of the interface polygon. Its size is the number of 
        // vertices of the interface polygon
        std::vector<size_t> indices(nVertices);

        // Define the vertex indices of the interface polygon
        for (size_t i = 0; i < nVertices; i += 1)
        {
            indices[i] = i;
        }

        // 
        for (auto jt = ccwPoints.begin(); jt != ccwPoints.end(); ++jt)
        {
            interfaceVertices.push_back(*jt);
        }

        // 
        interfaceFace.push_back(indices);

        return;
    }
}

VF Block::Intersect(const std::shared_ptr<Block> block, double threshold) const
{
    assert(block);

    // Make a copy of the geometry from the current block and get the reference to the geometry 
    // from the given block
    VF clipped(m_vf), & otherVf = block->Geometry();

    dcel::DCEL geom;

    size_t nFaces = otherVf.countFaces(), fIdx;

    toolkit::Plane clippingPlane;

    // Clip the cloned geometry using the planes from the faces of the given block
    for (fIdx = 0; fIdx < nFaces; fIdx += 1) 
    {
        clippingPlane.Set(otherVf.Normal(fIdx, true, threshold), otherVf.centroid(fIdx));

        geom.Set(clipped);

        clipped.Set(geom.Clip(clippingPlane, threshold));
    }

    return clipped;
}

Eigen::Vector3d Block::Intersect(
    const std::shared_ptr<Tessellation> tessellation, 
    double & t, 
    int direction, 
    double threshold) const
{
    assert(tessellation);

    std::shared_ptr<dcel::Face> face = tessellation->DCEL()->Faces()[m_faceIndex];

    double sign = (direction >= 0) ? 1.0 : -1.0;

    // Define a ray using the centroid and normal vector of the associated face of the block
    toolkit::Ray ray(face->Centroid(), face->Normal() * sign, true);

    // First, check if the ray contains a vertex of the geometry

    Eigen::Vector3d V = Eigen::Vector3d::Zero();
    
    size_t vIdx, nVertices = m_vf.countVertices();

    for (vIdx = 0; vIdx < nVertices; vIdx += 1)
    {
        V << m_vf.Vertex(vIdx);

        // Check if the current vertex lies in the "positive" section of the ray
        if (ray.isPointIn(V, t, threshold) && t > 0)
        {
            return V;
        }
    }

    // Second, check if the ray intersects a triangle of the geometry

    size_t fIdx, nFaceVertices = 0, i, j, nFaces = m_vf.countFaces();

    double denom = 0, s = 0;

    Eigen::Vector3d
        Vi = Eigen::Vector3d::Zero(),
        Vj = Eigen::Vector3d::Zero(),
        N = Eigen::Vector3d::Zero(), 
        P = Eigen::Vector3d::Zero();

    toolkit::LineSegment raySegment(ray.start(), ray.start() + ray.direction()), segment;

    for (fIdx = 0; fIdx < nFaces; fIdx += 1)
    {
        const std::vector<size_t> & face = m_vf.face(fIdx);

        V << m_vf.Vertex(face[0]);

        nFaceVertices = face.size() - 1;

        for (i = 1; i < nFaceVertices; i += 1)
        {
            j = i + 1;

            Vi << m_vf.Vertex(face[i]);
            Vj << m_vf.Vertex(face[j]);

            // Let's have a plane defined by point C and normal vector N, and a ray defined by 
            // point A and direction vector D. The plane and ray intesect in a point P if 
            // D * N != 0. The parameter for P along the ray is t = ((C - A) * N) / (D * N)

            N << (Vj - Vi).cross(V - Vi);

            // Check if the plane and ray are orthogonal to each other
            denom = ray.direction().dot(N);

            // Discard values close to zero
            if (abs(denom) <= threshold)
            {
                continue;
            }

            // Calculate the parameter for the intersection point along the ray
            t = (V - ray.start()).dot(N) / denom;

            // Discard intersection points along the "negative" section of the ray
            if (t < 0)
            {
                continue;
            }

            P = ray.at(t);

            // Check if the intersection point is in the triangle defined by vertices V0, Vi, Vj
            if (utils::isPointInTriangle(P, V, Vi, Vj, threshold))
            {
                return P;
            }

            // Third, check if there is an intersection point with any of the edges from the 
            // current triangle

            segment.Set(V, Vi);

            if (raySegment.Intersect(segment, P, s, t, threshold) && ray.isPointIn(P, t, threshold) && t > 0) 
            {
                return P;
            }

            segment.Set(Vi, Vj);

            if (raySegment.Intersect(segment, P, s, t, threshold) && ray.isPointIn(P, t, threshold) && t > 0)
            {
                return P;
            }

            segment.Set(Vj, V);

            if (raySegment.Intersect(segment, P, s, t, threshold) && ray.isPointIn(P, t, threshold) && t > 0)
            {
                return P;
            }
        }
    }



    std::cout << "Face: " << std::endl;
    face->Write();
    std::cout << "Ray: " << std::endl;
    ray.write();
    std::cout << std::endl;
    std::cout << "Geometry: " << std::endl;
    m_vf.Write();

    // There must exist an intersection point, not having one here is an error
    assert(false);
    return Eigen::Vector3d::Zero();
}

Eigen::Vector3d Block::Intersect(const toolkit::Ray & ray, double & t, double threshold) const
{
    Eigen::Vector3d V = Eigen::Vector3d::Zero();

    size_t vIdx, nVertices = m_vf.countVertices();

    // Check if the ray intersects any vertex of the block. The intersection parameter must be 
    // positive to be considered a valid intersection
    for (vIdx = 0; vIdx < nVertices; vIdx += 1)
    {
        V << m_vf.Vertex(vIdx);

        if (ray.isPointIn(V, t, threshold) && t > 0)
        {
            return V;
        }
    }

    // Second, check if the ray intersects a triangle of the geometry

    size_t fIdx, nFaceVertices = 0, i, j, nFaces = m_vf.countFaces();

    double denom = 0, s = 0;

    Eigen::Vector3d
        Vi = Eigen::Vector3d::Zero(),
        Vj = Eigen::Vector3d::Zero(),
        N = Eigen::Vector3d::Zero(),
        P = Eigen::Vector3d::Zero();

    toolkit::LineSegment raySegment(ray.start(), ray.start() + ray.direction()), segment;

    for (fIdx = 0; fIdx < nFaces; fIdx += 1)
    {
        const std::vector<size_t> & face = m_vf.face(fIdx);

        // First try, start with V at vertex 0

        V << m_vf.Vertex(face[0]);

        nFaceVertices = face.size() - 1;

        for (i = 1; i < nFaceVertices; i += 1)
        {
            // 
            Vi << m_vf.Vertex(face[i]);
            Vj << m_vf.Vertex(face[i + 1]);

            // Let's have a plane defined by point C and normal vector N, and a ray defined by 
            // point A and direction vector D. The plane and ray intesect in a point P if 
            // D * N != 0. The parameter for P along the ray is t = ((C - A) * N) / (D * N)

            N << (Vj - Vi).cross(V - Vi);

            // Check if the plane and ray are orthogonal to each other
            denom = ray.direction().dot(N);

            // Discard values close to zero
            if (abs(denom) <= threshold)
            {
                continue;
            }

            // Calculate the parameter for the intersection point along the ray
            t = (V - ray.start()).dot(N) / denom;

            // Discard intersection points along the "negative" section of the ray
            if (t < 0)
            {
                continue;
            }

            P = ray.at(t);

            // Check if the intersection point is in the triangle defined by vertices V0, Vi, Vj
            if (utils::isPointInTriangle(P, V, Vi, Vj, threshold))
            {
                return P;
            }

            // Third, check if there is an intersection point with any of the edges from the 
            // current triangle

            segment.Set(V, Vi);

            if (raySegment.Intersect(segment, P, s, t, threshold) && ray.isPointIn(P, t, threshold) && t > 0)
            {
                return P;
            }

            segment.Set(Vi, Vj);

            if (raySegment.Intersect(segment, P, s, t, threshold) && ray.isPointIn(P, t, threshold) && t > 0)
            {
                return P;
            }

            segment.Set(Vj, V);

            if (raySegment.Intersect(segment, P, s, t, threshold) && ray.isPointIn(P, t, threshold) && t > 0)
            {
                return P;
            }
        }

        // Second try, start with V at 1
        // Why do we need a second try? Well, thank numerical precision. Starting with 1 rather 
        // than 0 changes the triangulation of the faces. It could happen a different triangulation
        // is better for numerical precision. It happened to me multiple times and this is the 
        // simplest yet most effective way to handle intersections between a ray and a geometry. 
        // For sure there are better approaches, but this one works as expected and it is easier to
        // both maintain and understand

        V << m_vf.Vertex(face[1]);

        nFaceVertices = face.size();

        for (i = 2; i < nFaceVertices; i += 1)
        {
            j = (i == nFaceVertices - 1) ? 0 : i + 1;

            // 
            Vi << m_vf.Vertex(face[i]);
            Vj << m_vf.Vertex(face[j]);

            // Let's have a plane defined by point C and normal vector N, and a ray defined by 
            // point A and direction vector D. The plane and ray intesect in a point P if 
            // D * N != 0. The parameter for P along the ray is t = ((C - A) * N) / (D * N)

            N << (Vj - Vi).cross(V - Vi);

            // Check if the plane and ray are orthogonal to each other
            denom = ray.direction().dot(N);

            // Discard values close to zero
            if (abs(denom) <= threshold)
            {
                continue;
            }

            // Calculate the parameter for the intersection point along the ray
            t = (V - ray.start()).dot(N) / denom;

            // Discard intersection points along the "negative" section of the ray
            if (t < 0)
            {
                continue;
            }

            P = ray.at(t);

            // Check if the intersection point is in the triangle defined by vertices V0, Vi, Vj
            if (utils::isPointInTriangle(P, V, Vi, Vj, threshold))
            {
                return P;
            }

            // Third, check if there is an intersection point with any of the edges from the 
            // current triangle

            segment.Set(V, Vi);

            if (raySegment.Intersect(segment, P, s, t, threshold) && ray.isPointIn(P, t, threshold) && t > 0)
            {
                return P;
            }

            segment.Set(Vi, Vj);

            if (raySegment.Intersect(segment, P, s, t, threshold) && ray.isPointIn(P, t, threshold) && t > 0)
            {
                return P;
            }

            segment.Set(Vj, V);

            if (raySegment.Intersect(segment, P, s, t, threshold) && ray.isPointIn(P, t, threshold) && t > 0)
            {
                return P;
            }
        }
    }



    std::cout << "Face: " << std::endl;
    //face->Write();
    std::cout << "Ray: " << std::endl;
    ray.write();
    std::cout << std::endl;
    std::cout << "Geometry: " << std::endl;
    m_vf.Write();

    // There must exist an intersection point, not having one here is an error
    assert(false);
    return Eigen::Vector3d::Zero();
}

bool Block::IsEnabled() const
{
    return m_enabled;
}

std::vector<double> Block::Loads(double density, double gravity) const
{
    return { 0, m_vf.Volume() * density * gravity, 0, 0, 0, 0 };
}

VF Block::PlaneOffsetClipping(
    const toolkit::Plane & plane, 
    double extrados, 
    double intrados, 
    double threshold) const
{
    Eigen::Vector3d N = plane.Normal().normalized();

    // Set the plane to clip the extrados of the block
    toolkit::Plane clipPlane(N, plane.Point() + (N * extrados));

    dcel::DCEL geom(m_vf);

    // Clip the extrados
    VF clipped = geom.Clip(clipPlane, threshold);
    
    geom.Set(clipped);

    // Update the plane to clip the intrados of the block
    clipPlane.Set(-N, plane.Point() - (N * intrados));

    // Clip the intrados
    clipped = geom.Clip(clipPlane, threshold);

    // Clear the temporal geometry
    geom.Clear();

    // Return the pointer to the clipped block
    return clipped;
}

void Block::SetEnabled(bool enabled)
{
    m_enabled = enabled;
}

void Block::setForceLoad(double x, double y, double z) 
{
    m_loads[0] = x;
    m_loads[1] = y;
    m_loads[2] = z;
}

void Block::SetTorqueLoad(double x, double y, double z) 
{
    m_loads[3] = x;
    m_loads[4] = y;
    m_loads[5] = z;
}

std::vector<int> Block::VertexLocation(const toolkit::Plane & plane, double threshold)
{
    size_t nVertices = m_vf.countVertices();

    std::vector<int> locations(nVertices, 0);

    for (size_t i = 0; i < nVertices; i += 1)
    {
        locations[i] = plane.PointLocation(m_vf.Vertex(i), threshold);
    }

    return locations;
}

void Block::Write() const 
{
    // Write the label indicating the start of the block content
    std::cout << "block" << std::endl;

    // Write the information of the geometry of the block
    m_vf.Write();

    // Write the label indicating the end of the block content
    std::cout << "/block" << std::endl;
}

void Block::Write(std::stringstream & ss) const 
{
    // Write the label indicating the start of the block content
    ss << "block" << std::endl;

    // Write the information of the geometry of the block
    m_vf.Write(ss);

    // Write the label indicating the end of the block content
    ss << "/block" << std::endl;
}