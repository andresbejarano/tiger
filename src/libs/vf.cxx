#include <tiger/ds/vf.h>
#include <tiger/utils.h>

#include <iostream>
#include <fstream>
#include <set>

VF::VF() : 
    m_F(1), 
    m_facesAttributes(1), 
    m_V(1, 4),
    m_addFaceIndex(0),
    m_addVertexIndex(0)
{
}

VF::VF(const size_t nV, const size_t nF) :
    m_F(nF), 
    m_facesAttributes(nF), 
    m_V(nV, 4),
    m_addFaceIndex(0), 
    m_addVertexIndex(0)
{   
}

VF::VF(
    const std::list<Eigen::Vector3d> & vertices, 
    const std::list<std::vector<size_t>> & faces) : 
    m_F(faces.size()), 
    m_facesAttributes(faces.size()), 
    m_V(vertices.size(), 4),
    m_addFaceIndex(0),
    m_addVertexIndex(0)
{
    for (auto it = vertices.begin(); it != vertices.end(); ++it) 
    {
        addVertex(*it);
    }

    for (auto it = faces.begin(); it != faces.end(); ++it) 
    {
        addFace(*it);
    }
}

VF::VF(const VF & vf) : 
    m_F(vf.countFaces()), 
    m_facesAttributes(vf.countFaces()), 
    m_V(vf.countVertices(), 4),
    m_addFaceIndex(0),
    m_addVertexIndex(0)
{
    size_t nVertices = vf.countVertices(), nFaces = vf.countFaces(), i;

    for (i = 0; i < nVertices; i += 1) 
    {
        addVertex(vf.Vertex(i));
    }

    for (i = 0; i < nFaces; i += 1) 
    {
        addFace(vf.face(i));
    }
}

VF::~VF()
{
    Reset();
}

size_t VF::addFace(size_t v0, size_t v1)
{
    assert(m_addFaceIndex < m_F.size());

    m_F[m_addFaceIndex] = { v0, v1 };
    return m_addFaceIndex++;
}

size_t VF::addFace(size_t v0, size_t v1, size_t v2)
{
    assert(m_addFaceIndex < m_F.size());

    m_F[m_addFaceIndex] = { v0, v1, v2 };
    return m_addFaceIndex++;
}

size_t VF::addFace(size_t v0, size_t v1, size_t v2, size_t v3)
{
    assert(m_addFaceIndex < m_F.size());

    m_F[m_addFaceIndex] = { v0, v1, v2, v3 };
    return m_addFaceIndex++;
}

size_t VF::addFace(const std::vector<size_t> & indices)
{
    assert(m_addFaceIndex < m_F.size());

    m_F[m_addFaceIndex] = indices;
    return m_addFaceIndex++;
}

size_t VF::addVertex(double x, double y, double z) 
{
    assert(m_addVertexIndex < m_V.rows());

    m_V.row(m_addVertexIndex) << x, y, z, 1;
    return m_addVertexIndex++;
}

size_t VF::addVertex(const Eigen::Vector3d & P)
{
    assert(m_addVertexIndex < m_V.rows());

    m_V.row(m_addVertexIndex) << P.x(), P.y(), P.z(), 1;
    return m_addVertexIndex++;
}

double VF::area() const
{
    double sumAreas = 0;

    size_t nFaces = m_F.size(), i;

    for (i = 0; i < nFaces; i += 1)
    {
        sumAreas += area(i);
    }

    return sumAreas;
}

double VF::area(size_t fIdx) const
{
    assert(fIdx >= 0 && fIdx < m_addFaceIndex);

    double sumAreas = 0;

    Eigen::Vector3d
        V0 = Eigen::Vector3d::Zero(),
        Vi = Eigen::Vector3d::Zero(),
        Vj = Eigen::Vector3d::Zero();

    V0 << m_V(m_F[fIdx][0], 0), m_V(m_F[fIdx][0], 1), m_V(m_F[fIdx][0], 2);

    size_t nVertices = m_F[fIdx].size() - 1, i, j;

    for (i = 1; i < nVertices; i += 1) 
    {
        j = i + 1;

        Vi << m_V(m_F[fIdx][i], 0), m_V(m_F[fIdx][i], 1), m_V(m_F[fIdx][i], 2);
        Vj << m_V(m_F[fIdx][j], 0), m_V(m_F[fIdx][j], 1), m_V(m_F[fIdx][j], 2);

        sumAreas += abs((Vj - Vi).cross(V0 - Vi).norm());
    }

    return sumAreas / 2.0;
}

void VF::axisAlignedBoundingBox(Eigen::Vector3d & min, Eigen::Vector3d & max) const
{
    // Initialize both min and max points using the coordinates of the first vertex of the geometry
    min << m_V(0, 0), m_V(0, 1), m_V(0, 2);
    max << min;

     // Traverse through the vertices of the geometry and update the corner points respectively
    for (Eigen::Index i = 1; i < m_addVertexIndex; i += 1)
    {
        if (m_V(i, 0) < min.x())
        {
            min(0) = m_V(i, 0);
        }

        if (m_V(i, 1) < min.y())
        {
            min(1) = m_V(i, 1);
        }

        if (m_V(i, 2) < min.z())
        {
            min(2) = m_V(i, 2);
        }

        if (m_V(i, 0) > max.x())
        {
            max(0) = m_V(i, 0);
        }

        if (m_V(i, 1) > max.y())
        {
            max(1) = m_V(i, 1);
        }

        if (m_V(i, 2) > max.z())
        {
            max(2) = m_V(i, 2);
        }
    }
}

Eigen::Vector3d VF::centroid(bool fixZeros, double threshold) const
{
    Eigen::Vector3d C = Eigen::Vector3d::Zero(), P = Eigen::Vector3d::Zero();

    Eigen::Index i;

    for (i = 0; i < m_addVertexIndex; i += 1)
    {
        P << m_V(i, 0), m_V(i, 1), m_V(i, 2);

        C += P;
    }

    C /= (double)m_addVertexIndex;

    if (fixZeros) 
    {
        utils::fixZeros(C, threshold);
    }

    return C;
}

Eigen::Vector3d VF::centroid(size_t index, bool fixZeros, double threshold) const
{
    assert(index >= 0 && index < m_F.size());
    
    Eigen::Vector3d C = Eigen::Vector3d::Zero(), V = Eigen::Vector3d::Zero();

    size_t i, nVertices = m_F[index].size();

    for (i = 0; i < nVertices; i += 1) 
    {
        V << m_V(m_F[index][i], 0), m_V(m_F[index][i], 1), m_V(m_F[index][i], 2);

        C += V;
    }

    C /= (double)nVertices;

    if (fixZeros) 
    {
        utils::fixZeros(C, threshold);
    }

    return C;
}

void VF::clear()
{
    m_F.clear();
    m_attributes.Clear();
}

VF VF::clone() const
{
    VF clone(m_addVertexIndex, m_addFaceIndex);

    Eigen::Index i;

    for (i = 0; i < m_addVertexIndex; i += 1) 
    {
        clone.addVertex(m_V(i, 0), m_V(i, 1), m_V(i, 2));
    }

    for (auto it = m_F.begin(); it != m_F.end(); ++it) 
    {
        clone.addFace(*it);
    }

    return clone;
}

size_t VF::countEdges() const
{
    std::set<std::tuple<size_t, size_t>> edges;

    size_t nFaces = m_F.size(), nVertices, fIdx, vIdx, i, j;

    std::tuple<size_t, size_t> edge;

    // 
    for (fIdx = 0; fIdx < nFaces; fIdx += 1)
    {
        nVertices = m_F[fIdx].size() - 1;

        // 
        for (vIdx = 0; vIdx < nVertices; vIdx += 1) 
        {
            i = m_F[fIdx][vIdx];
            j = m_F[fIdx][vIdx + 1];
            
            edge = std::make_tuple(std::min(i, j), std::max(i, j));
            
            if (edges.find(edge) == edges.end()) 
            {
                edges.insert(edge);
            }
        }

        // Check for the last edge between the first and last vertices of the face
        i = m_F[fIdx][0], j = m_F[fIdx][nVertices];

        edge = std::make_tuple(std::min(i, j), std::max(i, j));

        // 
        if (edges.find(edge) == edges.end())
        {
            edges.insert(edge);
        }
    }

    return edges.size();
}

size_t VF::countFaces() const
{
    return m_addFaceIndex;
}

size_t VF::countSides(size_t face) const
{
    assert(face >= 0 && face < m_addFaceIndex);

    return m_F[face].size();
}

size_t VF::countTriangles() const
{
    size_t nTriangles = 0;

    for (auto it = m_F.begin(); it != m_F.end(); ++it)
    {
        nTriangles += (*it).size() - 2;
    }

    return nTriangles;
}

size_t VF::countTriangles(size_t face) const
{
    assert(face >= 0 && face < m_addFaceIndex);

    return m_F[face].size() - 2;
}

size_t VF::countVertices() const
{
    return m_addVertexIndex;
}

size_t VF::countVertices(size_t face) const
{
    assert(face >= 0 && face < m_addFaceIndex);

    return m_F[face].size();
}

Eigen::Vector3d VF::direction(
    size_t face, 
    size_t edge, 
    bool normalize, 
    bool fixZeros, 
    double threshold) const
{
    assert(face >= 0 && face < m_addFaceIndex);
    assert(edge >= 0 && edge < m_F[face].size());

    // 
    size_t next = (edge < m_F[face].size() - 1) ? edge + 1 : 0;

    Eigen::Vector3d A = Eigen::Vector3d::Zero(), B = Eigen::Vector3d::Zero();

    A << m_V(m_F[face][edge], 0), m_V(m_F[face][edge], 1), m_V(m_F[face][edge], 2);
    B << m_V(m_F[face][next], 0), m_V(m_F[face][next], 1), m_V(m_F[face][next], 2);
    
    Eigen::Vector3d D = B - A;

    // 
    if (normalize) 
    {
        D.normalize();
    }

    // 
    if (fixZeros) 
    {
        utils::fixZeros(D, threshold);
    }

    return D;
}

const std::vector<size_t> & VF::face(size_t index) const
{
    assert(index >= 0 && index < m_addFaceIndex);

    return m_F[index];
}

bool VF::facesHaveEvenNumberOfSides() const
{
    size_t nFaces = m_F.size(), i = 0;

    for (i = 0; i < nFaces; i += 1) 
    {
        if (!HasEvenNumberOfSides(i)) 
        {
            return false;
        }
    }

    return true;
}

void VF::flip()
{
    size_t fIdx, n, N, vIdx, t, idx;

    for (fIdx = 0; fIdx < m_addFaceIndex; fIdx += 1) 
    {
        N = m_F[fIdx].size();
        n = m_F[fIdx].size() / 2;

        for (vIdx = 0; vIdx < n; vIdx += 1) 
        {
            t = m_F[fIdx][vIdx];
            idx = N - vIdx - 1;
            m_F[fIdx][vIdx] = m_F[fIdx][idx];
            m_F[fIdx][idx] = t;
        }
    }
}

DynamicAttributes & VF::faceAttributes(size_t index)
{
    assert(index >= 0 && index < m_addFaceIndex);

    return m_facesAttributes[index];
}

Eigen::MatrixXd VF::faceVertices(size_t face) const
{
    assert(face >= 0 && face < m_F.size());

    size_t nVertices = m_F[face].size(), i = 0;
    
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(nVertices, 4);

    for (i = 0; i < nVertices; i += 1) 
    {
        M.row(i) << m_V.row(m_F[face][i]);
    }

    return M;
}

bool VF::HasEvenNumberOfSides(size_t face) const
{
    assert(face >= 0 && face < m_addFaceIndex);

    return m_F[face].size() % 2 == 0;
}

bool VF::IsComplete() const
{
    return m_addVertexIndex == m_V.rows() && m_addFaceIndex == m_F.size();
}

bool VF::IsCoplanar(size_t face, const toolkit::Plane & plane, double threshold) const
{
    assert(face >= 0 && face < m_addFaceIndex);

    size_t nVertices = m_F[face].size(), i = 0;

    Eigen::Vector3d V = Eigen::Vector3d::Zero();

    for (i = 0; i < nVertices; i += 1) 
    {
        V << m_V(m_F[face][i], 0), m_V(m_F[face][i], 1), m_V(m_F[face][i], 2);

        if (!plane.IsPointInPlane(V, threshold))
        {
            return false;
        }
    }

    return true;
}

bool VF::IsPlanar(size_t face, double threshold) const
{
    assert(face >= 0 && face < m_addFaceIndex);

    return IsCoplanar(face, Plane(face), threshold);
}

bool VF::IsPointIn(size_t face, const Eigen::Vector3d & P, double threshold) const
{
    assert(face >= 0 && face < m_addFaceIndex);

    // Test 1: Check if point P lies in the same plane of the face, if not then return false
    if (!Plane(face).IsPointInPlane(P, threshold)) 
    {
        return false;
    }

    size_t nVertices = m_F[face].size() - 1, i = 0, j = 0;

    Eigen::Vector3d 
        V0 = Eigen::Vector3d::Zero(), 
        Vi = Eigen::Vector3d::Zero(), 
        Vj = Eigen::Vector3d::Zero();

    // Keep the coordinates of the first vertex of the face
    V0 << m_V(m_F[face][0], 0), m_V(m_F[face][0], 1), m_V(m_F[face][0], 2);

    // Traverse through the edges of the face. Here we define a triangle between vertex v0 and the 
    // end points of the current edge
    for (i = 1; i < nVertices; i += 1) 
    {
        j = i + 1;

        Vi << m_V(m_F[face][i], 0), m_V(m_F[face][i], 1), m_V(m_F[face][i], 2);
        Vj << m_V(m_F[face][j], 0), m_V(m_F[face][j], 1), m_V(m_F[face][j], 2);

        // Test 2: Check if P lies within the triangle defined by points v0, v1 and v2. If that's 
        // the case then return true
        if (utils::isPointInTriangle(P, V0, Vi, Vj, threshold)) 
        {
            return true;
        }
    }

    // 
    toolkit::LineSegment edge;

    // Traverse through the edges of the face and check if the point lies at any of them
    for (i = 0; i < nVertices - 1; i += 1) 
    {
        j = (i == nVertices - 1) ? 0 : i + 1;

        Vi << m_V(m_F[face][i], 0), m_V(m_F[face][i], 1), m_V(m_F[face][i], 2);
        Vj << m_V(m_F[face][j], 0), m_V(m_F[face][j], 1), m_V(m_F[face][j], 2);

        edge.Set(Vi, Vj);

        // Test 3: Check if P lies within the current edge. If that's the case then return true
        if (edge.IsPointIn(P, threshold)) 
        {
            return true;
        }
    }
    
    // Return false since P didn't pass any the tests
    return false;
}

bool VF::IsPointInSurface(const Eigen::Vector3d & P, double threshold) const
{
    size_t fIdx = 0;

    for (fIdx = 0; fIdx < m_addFaceIndex; fIdx += 1)
    {
        if (!IsPointIn(fIdx, P, threshold))
        {
            return false;
        }
    }

    return true;
}

Eigen::Vector3d VF::Normal(size_t face, bool normalize, bool fixZeros, double threshold) const
{
    assert(face >= 0 && face < m_addFaceIndex);

    Eigen::Vector3d 
        N = Eigen::Vector3d::Zero(), 
        V0 = Eigen::Vector3d::Zero(), 
        Vi = Eigen::Vector3d::Zero(), 
        Vj = Eigen::Vector3d::Zero();

    V0 << m_V(m_F[face][0], 0), m_V(m_F[face][0], 1), m_V(m_F[face][0], 2);

    size_t i = 0, nVertices = m_F[face].size() - 1, j = 0;

    for (i = 1; i < nVertices; i += 1) 
    {
        j = i + 1;
        
        Vi << m_V(m_F[face][i], 0), m_V(m_F[face][i], 1), m_V(m_F[face][i], 2);
        Vj << m_V(m_F[face][j], 0), m_V(m_F[face][j], 1), m_V(m_F[face][j], 2);

        N += (Vj - Vi).cross(V0 - Vi);
    }

    if (nVertices > 0) 
    {
        N /= (double)(nVertices - 1);
    }

    if (normalize) 
    {
        N.normalize();
    }

    if (fixZeros) 
    {
        utils::fixZeros(N, threshold);
    }

    return N;
}

void VF::NormalizeVertices()
{
    Eigen::Vector3d V = Eigen::Vector3d::Zero();

    Eigen::Index i = 0;

    for (i = 0; i < m_addVertexIndex; i += 1) 
    {
        V << m_V(i, 0), m_V(i, 1), m_V(i, 2);
        V.normalize();

        m_V.row(i) << V.x(), V.y(), V.z(), 1;
    }
}

toolkit::Plane VF::Plane(size_t face, bool normalize, bool fixZeros, double threshold) const
{
    assert(face >= 0 && face < m_addFaceIndex);

    return toolkit::Plane(
        Normal(face, normalize, fixZeros, threshold), 
        centroid(face, fixZeros, threshold));
}

int VF::PointLocation(size_t face, const Eigen::Vector3d & point, double threshold) const
{
    assert(face >= 0 && face < m_addFaceIndex);

    return Plane(face, false, true, threshold).PointLocation(point, threshold);
}

void VF::Reset()
{
    // 
    m_F.erase(m_F.begin(), m_F.end());
    m_F.clear();

    m_attributes.Clear();
    m_facesAttributes.erase(m_facesAttributes.begin(), m_facesAttributes.end());
}

void VF::Rotate(const Eigen::Vector3d & K, double angle)
{
    double sin_a = sin(angle);
    double cos_a = cos(angle);
    double _cos_a = 1.0 - cos_a;

    Eigen::Vector3d Kn = K.normalized(), V = Eigen::Vector3d::Zero(), R = Eigen::Vector3d::Zero();

    Eigen::Index i = 0;

    for (i = 0; i < m_addVertexIndex; i += 1) 
    {
        V << m_V(i, 0), m_V(i, 1), m_V(i, 2);
        R << (V * cos_a) + (Kn.cross(V) * sin_a) + (Kn * Kn.dot(V) * _cos_a);

        m_V.row(i) << R.x(), R.y(), R.z(), 1;
    }
}

void VF::Scale(double s, const Eigen::Vector3d & C)
{
    double _s = 1.0 - s;

    Eigen::Matrix4d S;

    S <<    s,          0,          0,          0,
            0,          s,          0,          0,
            0,          0,          s,          0,
            C.x() * _s, C.y() * _s, C.z() * _s, 1;

    m_V *= S;
}

void VF::Set(const VF & vf)
{
    m_addVertexIndex = 0;
    m_addFaceIndex = 0;
    m_F.erase(m_F.begin(), m_F.end());
    m_F.clear();

    size_t nVertices = vf.countVertices(), nFaces = vf.countFaces(), i;

    m_V.resize(nVertices, 4);
    m_F.resize(nFaces);

    for (i = 0; i < nVertices; i += 1) 
    {
        addVertex(vf.Vertex(i));
    }

    for (i = 0; i < nFaces; i += 1) 
    {
        addFace(vf.face(i));
    }
}

void VF::Translate(const Eigen::Vector3d & D)
{
    Eigen::Matrix4d T;

    T << 1,     0,     0,     0,
         0,     1,     0,     0,
         0,     0,     1,     0,
         D.x(), D.y(), D.z(), 1;

    m_V *= T;
}

VF VF::TriangulateFacesByVertices() const
{
    // The new number of vertices of the triangulates geometry is the current number of vertices 
    // plus the number of faces (one new vertex per face)
    size_t nNewVertices = (size_t)m_addVertexIndex + m_addFaceIndex;

    size_t nNewFaces = 0, fIdx, newIdx = 0, nFaceVertices = 0, vIdx, jIdx;

    for (fIdx = 0; fIdx < m_addFaceIndex; fIdx += 1) 
    {
        nNewFaces += m_F[fIdx].size();
    }

    // Initialize the triangulated geometry
    VF triangulated(nNewVertices, nNewFaces);

    // Copy the original vertices into the triangulated geometry
    for (vIdx = 0; vIdx < (size_t)m_addVertexIndex; vIdx += 1)
    {
        triangulated.addVertex(Vertex(vIdx));
    }

    // Define the new vertices and faces of the triangulated geometry
    for (fIdx = 0; fIdx < m_addFaceIndex; fIdx += 1) 
    {
        newIdx = triangulated.addVertex(centroid(fIdx));

        nFaceVertices = m_F[fIdx].size();

        for (vIdx = 0; vIdx < nFaceVertices; vIdx += 1) 
        {
            jIdx = (vIdx == nFaceVertices - 1) ? 0 : vIdx + 1;

            triangulated.addFace(m_F[fIdx][vIdx], m_F[fIdx][jIdx], newIdx);
        }
    }

    // Return the triangulated geometry
    return triangulated;
}

/*double VF::Volume() const
{
    // Initialize the volume of the polyhedron
    double volume = 0.0;

    size_t nIndices = F.size();

    // Traverse through the faces
    for (size_t i = 0; i , nIndices; i += 1)
    {
        

        // Get the pointer to the current face
        std::shared_ptr<Face> face = *it;

        // Keep the reference to the start vertex of the incident half edge of the face. It will be
        // used for triangulating the face
        const Eigen::Vector3d & v0 = face->halfedge->start->Coords();

        // Get the reference to the next half edge to the incident half edge of the face
        std::shared_ptr<Halfedge> halfedge = face->halfedge->next;

        // Traverse through the half edges of the face and define the points of the triangles using
        // v0 and the end points of the half edges. Stop at the previous half edge to the incident 
        // half edge of the face
        do
        {
            // Get the references to the other two points of the current triangle
            const Eigen::Vector3d & v1 = halfedge->start->Coords();
            const Eigen::Vector3d & v2 = halfedge->twin->start->Coords();

            // Calculate the cross product between the sides of the triangle
            Eigen::Vector3d N = (v2 - v1).cross(v0 - v1);

            // Calculate the area of the current triangle (actually it is twice the area)
            double twiceArea = abs(N.norm());

            // Normalize the normal vector of the current triangle
            N.normalize();

            // Calculate the expression of the volume from the current face and accumulate it into 
            // the volume value
            volume += (v1.dot(N) * twiceArea);

            // Move to the next half edge of the face
            halfedge = halfedge->next;

        } while (halfedge != face->halfedge->previous);
    }

    // Divide the volume by 6 and return it
    return volume / 6.0;
}*/

Eigen::Vector3d VF::Vertex(size_t index) const
{
    Eigen::Index i = index;

    assert(i >= 0 && i < m_addVertexIndex);

    return Eigen::Vector3d(m_V(i, 0), m_V(i, 1), m_V(i, 2));
}

void VF::Vertex(size_t index, const Eigen::Vector3d& P)
{
    Eigen::Index i = index;

    assert(i < m_addVertexIndex);

    m_V.row(i) << P.x(), P.y(), P.z(), 1;
}

void VF::Vertex(size_t index, double x, double y, double z)
{
    Eigen::Index i = index;

    assert(i < m_addVertexIndex);

    m_V.row(i) << x, y, z, 1;
}

double VF::Volume() const
{
    double volume = 0, area = 0;

    size_t nFaces = m_F.size(), fIdx = 0, nIndices = 0, i = 0, j = 0;

    Eigen::Vector3d 
        C = Eigen::Vector3d::Zero(), 
        N = Eigen::Vector3d::Zero(), 
        X = Eigen::Vector3d::Zero(), 
        V0 = Eigen::Vector3d::Zero(), 
        Vi = Eigen::Vector3d::Zero(), 
        Vj = Eigen::Vector3d::Zero();

    for (fIdx = 0; fIdx < nFaces; fIdx += 1) 
    {
        // The number of vertices of the current face minus one
        nIndices = m_F[fIdx].size() - 1;

        // Keep the information of the first vertex of the face
        V0 << m_V(m_F[fIdx][0], 0), m_V(m_F[fIdx][0], 1), m_V(m_F[fIdx][0], 2);

        // Initialize the variables for the centroid, normal and area of the current face
        C << V0;
        N << 0, 0, 0;
        area = 0;

        // Traverse through the vertices of the current face, update the centroid, normal and area
        // values accordingly
        for (i = 1; i < nIndices; i += 1) 
        {
            // Get the information of the vertices from the current edge
            j = i + 1;
            Vi << m_V(m_F[fIdx][i], 0), m_V(m_F[fIdx][i], 1), m_V(m_F[fIdx][i], 2);
            Vj << m_V(m_F[fIdx][j], 0), m_V(m_F[fIdx][j], 1), m_V(m_F[fIdx][j], 2);

            // Update the centroid, normal and area values
            C += Vi;
            X << (Vj - Vi).cross(V0 - Vi);
            N += X;
            area += abs(X.norm());
        }

        // Add the last vertex and calculate the centroid (divide by the number of vertices of the 
        // face)
        C += Vj;
        C /= (double)(nIndices + 1);

        // Calculate the normal vector (divide by the number of triangles of the face) and 
        // normalize it
        N /= (double)(nIndices - 1);
        N.normalize();

        // Update the volume value (do not divide by 2, it'll be done at the end)
        volume += C.dot(N) * area;
    }

    return volume / 6.0;
}

void VF::Write() const
{
    for (Eigen::Index i = 0; i < m_addVertexIndex; i += 1) 
    {
        std::cout << "V: (" << m_V(i, 0) << ", " << m_V(i, 1) << ", " << m_V(i, 2) << ")" << std::endl;
    }

    for (auto it = m_F.begin(); it != m_F.end(); ++it)
    {
        const std::vector<size_t> & indices = *it;

        std::cout << "F:";

        for (auto jt = indices.begin(); jt != indices.end(); ++jt)
        {
            std::cout << ' ' << std::to_string(*jt);
        }

        std::cout << std::endl;
    }
}

void VF::Write(std::stringstream & ss) const
{
    for (Eigen::Index i = 0; i < m_addVertexIndex; i += 1)
    {
        ss << "v " << m_V(i, 0) << ' ' << m_V(i, 1) << ' ' << m_V(i, 2) << std::endl;
    }

    for (auto it = m_F.begin(); it != m_F.end(); ++it)
    {
        const std::vector<size_t> & indices = *it;

        ss << 'f';

        for (auto i = indices.begin(); i != indices.end(); ++i)
        {
            ss << ' ' << ((*i) + 1);
        }

        ss << std::endl;
    }
}

void VF::WriteGeogebraJs(
    std::stringstream & ss, 
    const std::string & prefix, 
    int r, 
    int g, 
    int b, 
    int a, 
    double threshold) const
{
    // Get the number of vertices of the geometry
    /*size_t nVertices = vertices.size();

    // Traverse through the vertices of the geometry
    for (size_t i = 0; i < nVertices; i += 1)
    {
        // Set the index of the vertex
        vertices[i]->Attributes().Set<size_t>(ATTRIB_INDEX, i);

        // Write the commands that generate the current vertex in GeoGebra
        vertices[i]->WriteGeogebraJs(ss, prefix, r, g, b, threshold);
    }

    // Get the number of faces of the geometry
    size_t nFaces = faces.size();

    // Traverse through the faces and write their content
    for (size_t i = 0; i < nFaces; i += 1)
    {
        // Set the index of the face
        faces[i]->Attributes().Set<size_t>(ATTRIB_INDEX, i);

        // If the current face is planar then write its indices using a single polygon. Otherwise, 
        // write its indices as triangles
        if (faces[i]->IsPlanar(threshold))
        {
            faces[i]->WritePlanarGeogebraJs(ss, prefix);
        }
        else
        {
            faces[i]->WriteTriangularGeogebraJs(ss, prefix);
        }
    }

    // remove the index attributes
    RemoveVerticesAttribute(ATTRIB_INDEX);
    RemoveFacesAttribute(ATTRIB_INDEX);*/
}

void VF::WriteObjFile(const std::string & filepath, double threshold) const
{
    Eigen::Vector3d P = Eigen::Vector3d::Zero();

    std::ofstream file(filepath, std::ios::out);
    
    for (Eigen::Index i = 0; i < m_addVertexIndex; i += 1) 
    {
        P << m_V(i, 0), m_V(i, 1), m_V(i, 2);

        utils::fixZeros(P, threshold);

        file << "v " << P.x() << ' ' << P.y() << ' ' << P.z() << std::endl;
    }

    for (auto it = m_F.begin(); it != m_F.end(); ++it) 
    {
        const std::vector<size_t> & indices = *it;

        file << 'f';

        for (auto i = indices.begin(); i != indices.end(); ++i) 
        {
            file << ' ' << ((*i) + 1);
        }

        file << std::endl;
    }

    file.close();
}
