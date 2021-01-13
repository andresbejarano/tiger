#include "dcel.h"
#include "utils.h"
#include "Eigen/Geometry"
#include <iostream>
#include <list>
#include <math.h>
#include <fstream>


dcel::Vertex::Vertex() : 
    m_coords(0.0, 0.0, 0.0), 
    halfedge(nullptr)
{
}

dcel::Vertex::Vertex(const Eigen::Vector3d & C) : 
    m_coords(C), 
    halfedge(nullptr)
{
}

dcel::Vertex::Vertex(double x, double y, double z) :
    m_coords(x, y, z),
    halfedge(nullptr)
{
}

dcel::Vertex::~Vertex()
{
    Clear();
}

void dcel::Vertex::CheckConsistency() const
{
    // Vertex must have a half edge
    assert(halfedge);

    // *this must be the same as halfedge->start
    assert(shared_from_this() == halfedge->start);
}

void dcel::Vertex::Clear()
{
    m_attributes.Clear();
    m_coords << 0, 0, 0;
    halfedge = nullptr;
}

const Eigen::Vector3d & dcel::Vertex::Coords() const
{
    return m_coords;
}

void dcel::Vertex::Coords(const Eigen::Vector3d & V, bool fixZeros, double threshold)
{
    m_coords << V;

    // Fix the values close to zero if indicated
    if (fixZeros)
    {
        utils::fixZeros(m_coords, threshold);
    }
}

void dcel::Vertex::Coords(double x, double y, double z, bool fixZeros, double threshold)
{
    m_coords << x, y, z;

    // Fix the values close to zero if indicated
    if (fixZeros)
    {
        utils::fixZeros(m_coords, threshold);
    }
}

std::string dcel::Vertex::CoordsToString(
    const std::string & sep, 
    const std::string & opening, 
    const std::string & closing, 
    double threshold) const
{
    // Get a copy of the coordinates of the vertex and fix its zeros
    Eigen::Vector3d C(m_coords);
    utils::fixZeros(C, threshold);

    // Initialize a String Stream and define its content
    std::stringstream ss;
    ss << opening << C.x() << sep << C.y() << sep << C.z() << closing;

    // Convert the String Stream to String and return it
    return ss.str();
}

void dcel::Vertex::FixZeros(double threshold)
{
    utils::fixZeros(m_coords, threshold);
}

bool dcel::Vertex::IsSurrounded() const
{
	// Initialize the surrounded indicator of the vertex. Assume it is surrounded
	bool surrounded = true;

	// Get the pointer to the incident half edge of the vertex
	std::shared_ptr<Halfedge> currentHalfedge = halfedge;

	// Traverse through the half edges around the vertex
	do 
	{
		// If the current half edge does not have an incident face then indicate the vertex is no 
		// longer surrounded
		if (!currentHalfedge->face) 
		{
			surrounded = false;
		}

		// Move to the next half edge around the vertex
		currentHalfedge = currentHalfedge->twin->next;

	} while (surrounded && currentHalfedge != halfedge);

	// Return the surrounded indicator
	return surrounded;
}

size_t dcel::Vertex::NumberOfIncidentFaces() const
{
    // Initialize the counter for the incident faces
    size_t faceCount = 0;

    // Get the incident half edge of the vertex
    std::shared_ptr<Halfedge> currentHalfedge = halfedge;

    // Traverse through the half edges around the vertex
    do
    {
        // If the current half edge has an incident face then increase the counter by one
        if (currentHalfedge->face)
        {
            faceCount += 1;
        }

        // Move to the next half edge around the vertex
        currentHalfedge = currentHalfedge->twin->next;

    } while (currentHalfedge != halfedge);

    // Return the count of incident faces
    return faceCount;
}

void dcel::Vertex::WriteGeogebraJs(
    std::stringstream & ss, 
    const std::string prefix, 
    int r, 
    int g, 
    int b, 
    double threshold) const
{
    // Get the index of the vertex
    size_t index;
    assert(m_attributes.Get<size_t>(ATTRIB_INDEX, index));

    // Write the command that generates the current vertex in GeoGebra
    ss << "ggbApplet.evalCommand(\"" << prefix << "V" << index << " = Point(" << 
        CoordsToString(", ", "{", "}", threshold) << ")\");" << std::endl;

    // 
    ss << "ggbApplet.setColor(\"" << prefix << "V" << index << "\", " << r << ", " << g << ", " << 
        b << ");" << std::endl;
}

dcel::Face::Face() : halfedge(nullptr)
{
}

dcel::Face::~Face()
{
    Clear();
}

double dcel::Face::Area() const
{
    // Initialize the area of the face
    double area = 0.0;

    // Get the next half edge from the incident half edge of the face
    std::shared_ptr<Halfedge> currentHalfedge = halfedge->next;

    // Get the reference to the coordinates of the start vertex from the incident half edge of the 
    // face
    Eigen::Vector3d AB = Eigen::Vector3d::Zero(), AC = Eigen::Vector3d::Zero();
    const Eigen::Vector3d & A = halfedge->start->Coords();

    // Traverse through the half edges of the face
    do
    {
        // Get the vectors from A to the end points of the half edge
        AB << currentHalfedge->start->Coords() - A;
        AC << currentHalfedge->twin->start->Coords() - A;

        // Calculate the area of triangle between A and the end points of the current half edge
        area += (AB.cross(AC).norm());

        // Move to the next half edge
        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != halfedge->previous);

    // Return the area of the face
    return area / 2.0;
}

Eigen::Vector3d dcel::Face::Centroid(bool fixZeros, double threshold) const
{
    // 
    Eigen::Vector3d C = Eigen::Vector3d::Zero();

    // Initialize the counter for the number of sides of the face
    size_t nSides = 0;

    // Get the incident half edge
    std::shared_ptr<Halfedge> currentHalfedge = halfedge;

    // Traverse through the half edges
    do
    {
        // Add the values of the start vertex coordinates
        C += currentHalfedge->start->Coords();

        // Increment the number of faces by one
        nSides += 1;

        // Go to the next half edge
        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != halfedge);

    // Divide the point by the number of sides of the face
    if (nSides > 0)
    {
        C /= double(nSides);
    }

    // Fix the values close to zero if indicated
    if (fixZeros)
    {
        utils::fixZeros(C, threshold);
    }

    return C;
}

void dcel::Face::CheckConsistency() const
{
    // The face must have a half edge
    assert(halfedge);

    // Get the pointer to the incident half edge of the face
    std::shared_ptr<Halfedge> currentHalfedge = halfedge;

    // Traverse through the half edges of the face
    do
    {
        // Each half edge must be incident to this face
        assert(shared_from_this() == currentHalfedge->face);

        // Move to the next half edge
        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != halfedge);
}

void dcel::Face::Clear()
{
    m_attributes.Clear();
    halfedge = nullptr;
}

size_t dcel::Face::CountNeighbors() const
{
    size_t nNeighbors = 0;

    std::shared_ptr<Halfedge> currentHalfedge = halfedge;

    do 
    {
        if (currentHalfedge->twin->face) 
        {
            nNeighbors += 1;
        }

        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != halfedge);

    return nNeighbors;
}

size_t dcel::Face::CountSides() const
{
    // Initialize the counter of the number of sides of the face
    size_t nSides = 0;

    // Get the incident half edge
    std::shared_ptr<Halfedge> currentHalfedge = halfedge;

    // Traverse through the half edges
    do
    {
        // Increment the number of sides by one
        nSides += 1;

        // Go to the next half edge
        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != halfedge);

    // Return the number of sides
    return nSides;
}

size_t dcel::Face::CountTriangles() const
{
    return CountSides() - 2;
}

size_t dcel::Face::CountVertices() const
{
    return CountSides();
}

VF dcel::Face::Flip() const
{
    // 
    size_t nVertices = CountVertices();

    // 
    VF flipped(nVertices, 1);

    // 
    std::vector<size_t> indices(nVertices);

    // 
    std::shared_ptr<Halfedge> currentHalfedge = halfedge;

    size_t i = 0;

    // Traverse through the half edges of the face by
    do 
    {
        // 
        flipped.addVertex(currentHalfedge->start->Coords());

        // 
        indices.push_back(i++);

        // Move to the previous half edge
        currentHalfedge = currentHalfedge->previous;
    
    } while (currentHalfedge != halfedge);

    // 
    flipped.addFace(indices);

    return flipped;
}

bool dcel::Face::HasAllNeighbors() const
{
    // Get the incident half edge
    std::shared_ptr<Halfedge> currentHalfedge = halfedge;

    // Traverse through the half edges
    do
    {
        // If the twin half edge has no incident face then return false
        if (!currentHalfedge->twin->face)
        {
            return false;
        }

        // Go to the next half edge
        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != halfedge);

    // Indicate all neighboring faces exist
    return true;
}

bool dcel::Face::HasEvenNumberOfSides() const
{
    // Get the number of sides of the face and check if it is an even number
    return (CountSides() % 2 == 0);
}

bool dcel::Face::IsAtBoundary() const
{
    // If the face does not have at least one of its neighbors then it is at a boundary
    return !HasAllNeighbors();
}

bool dcel::Face::IsCoplanar(const toolkit::Plane & plane, double threshold) const
{
    // Get the pointer to the incident half edge of the face
    std::shared_ptr<dcel::Halfedge> currentHalfedge = halfedge;

    // Traverse through the half edges of the face and check if their incident start vertices are 
    // coplanar with the given plane
    do
    {
        // Check if the coordinates of the start vertex of the current half edge lie in the given 
        // plane. If not then return false
        if (!plane.IsPointInPlane(currentHalfedge->start->Coords(), threshold))
        {
            return false;
        }

        // Move to the next half edge
        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != halfedge);

    // Return true since all vertices of the face were found to be coplanar with the plane
    return true;
}

bool dcel::Face::IsPlanar(double threshold) const
{
    // Get the plane of the face
    toolkit::Plane plane = Plane(true, true, threshold);

    // Get the pointer to the incident half edge of the face
    std::shared_ptr<Halfedge> currentHalfedge = halfedge;

    // Traverse through the half edges of the face
    do 
    {
        // If the coordinates of the start vertex from the current half edge are not in the plane
        // of the face then return false. This could happen (although it is expected not to happen)
        if (!plane.IsPointInPlane(currentHalfedge->start->Coords(), threshold)) 
        {
            return false;
        }

        // Move to the next half edge
        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != halfedge);

    // Return true since all face vertices are coplanar
    return true;
}

bool dcel::Face::IsPointIn(const Eigen::Vector3d & P, double threshold) const
{
    // Test 1: Check if point P lies in the same plane of the face, if not then return false
    if (!Plane().IsPointInPlane(P, threshold))
    {
        return false;
    }

    // Get the next half edge to the incident half edge of the face
    std::shared_ptr<Halfedge> currentHalfedge = halfedge->next;

    // Get the coordinates of the start vertex from the incident half edge of the face
    const Eigen::Vector3d & v0 = halfedge->start->Coords();

    // Traverse through the half edges of the face. Finish at the previous half edge to the incident
    // half edge of the face. Here we define a triangle between vertex v0 and the end points of the 
    // current half edge
    do
    {
        // Get the references to the coordinates of the end points of the current half edge
        const Eigen::Vector3d & v1 = currentHalfedge->start->Coords();
        const Eigen::Vector3d & v2 = currentHalfedge->twin->start->Coords();

        // Test 2: Check if P lies within the triangle defined by points v0, v1 and v2. If that's 
        // the case then return true
        if (utils::isPointInTriangle(P, v0, v1, v2, threshold))
        {
            return true;
        }

        // Move to the next half edge
        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != halfedge->previous);

    // Get the incident half edge of the face
    currentHalfedge = halfedge;

    // Traverse through the half edges of the face
    do
    {
        // Test 3: Check if P lies within the current half edge. If that's the case then return 
        // true
        if (currentHalfedge->IsPointIn(P, threshold))
        {
            return true;
        }

        // Move to the next half edge
        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != halfedge);

    // Return false since P didn't pass any the tests
    return false;
}

Eigen::Vector3d dcel::Face::Normal(bool normalize, bool fixZeros, double threshold) const
{
    Eigen::Vector3d N = Eigen::Vector3d::Zero(), C = Eigen::Vector3d::Zero();

    // Initialize the counter for the number of sides of the face
    size_t nSides = 0;

    // Get the next half edge to the incident half edge
    std::shared_ptr<Halfedge> current = halfedge->next;

    // Get the coordinates of the start vertex of the incident half edge
    const Eigen::Vector3d & P = halfedge->start->Coords();

    // Traverse through the half edges
    do
    {
        // Build the two vectors for the current triangle section of the face and calculate its 
        // normal vector
        C << (current->start->Coords() - P).cross(current->twin->start->Coords() - P);

        // Add the cross product vector to the normal vector
        N += C;

        // Increment the number of sides by one
        nSides += 1;

        // Go to the next half edge
        current = current->next;

    } while (current->twin->start != halfedge->start);

    // Divide the normal vector by the number of sides
    if (nSides > 0)
    {
        N /= double(nSides);
    }

    // Normalize the normal vector if indicated
    if (normalize)
    {
        N.normalize();
    }

    // Fiz the zero values of the normal vector if indicated
    if (fixZeros)
    {
        utils::fixZeros(N, threshold);
    }

    // Return the normal vector of the face
    return N;
}

template <typename T>
void dcel::Face::SetVerticesAttribute(const std::string & name, T value)
{
    // Get the incident half edge of the face
    std::shared_ptr<Halfedge> currentHalfedge = halfedge;

    // Traverse through the half edges of the face
    do
    {
        // Set the start vertex of the half edge with the given attribute name and value
        currentHalfedge->start->Attributes().Set<T>(name, value);

        // Move to the next half edge
        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != halfedge);
}

void dcel::Face::Write() const
{
    std::shared_ptr<Halfedge> currentHalfedge = halfedge;

    do 
    {
        utils::write(currentHalfedge->start->Coords(), true);

        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != halfedge);
}

void dcel::Face::WritePlanarGeogebraJs(std::stringstream & ss, const std::string & prefix) const
{
    // Get the index of the face
    size_t index = 0;
    assert(m_attributes.Get<size_t>(ATTRIB_INDEX, index));

    // Start the command for generating the polygon representing the face
    ss << "ggbApplet.evalCommand(\"" << prefix << "F" << index << " = Polygon(";
    
    // Get the pointer to the incident half edge of the face
    std::shared_ptr<Halfedge> currentHalfedge = halfedge;

    // 
    bool firstVertex = true;

    // Traverse through the half edges of the face
    do 
    {
        // Get the index of the start vertex from the current half edge
        assert(currentHalfedge->start->Attributes().Get<size_t>(ATTRIB_INDEX, index));

        // 
        if (firstVertex)
        {
            ss << prefix << "V" << index;

            firstVertex = false;
        }
        else
        {
            ss << ", " << prefix << "V" << index;
        }

        // Move to the next half edge
        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != halfedge);

    // Close the command
    ss << ")\");" << std::endl;
}

void dcel::Face::WriteTriangularGeogebraJs(
    std::stringstream & ss, 
    const std::string & prefix) const
{
    // Get the index of the face
    size_t faceIdx;
    assert(m_attributes.Get<size_t>(ATTRIB_INDEX, faceIdx));

    // Get the index of the start vertex from the incident half edge of the face
    size_t v0Idx;
    assert(halfedge->start->Attributes().Get<size_t>(ATTRIB_INDEX, v0Idx));

    // Get the pointer to the next half edge from the incident half edge of the face
    std::shared_ptr<Halfedge> currentHalfedge = halfedge->next;

    // Keep track of the triangle indices of the face
    size_t triangleIdx = 0;

    // Traverse through the half edges of the face. Finish at the previous half edge from the 
    // incident half edge of thef ace
    do 
    {
        // Get the indices of the end points from the current half edge
        size_t v1Idx, v2Idx;
        assert(currentHalfedge->start->Attributes().Get<size_t>(ATTRIB_INDEX, v1Idx));
        assert(currentHalfedge->twin->start->Attributes().Get<size_t>(ATTRIB_INDEX, v2Idx));

        // Write the command for generating the current triangle of the face
        ss << "ggbApplet.evalCommand(\"" << prefix << "F" << faceIdx << "T" << triangleIdx << 
            " = Polygon(" << prefix << "V" << v0Idx << ", " << prefix << "V" << v1Idx << ", " << 
            prefix << "V" << v2Idx << ")\");" << std::endl;

        // Update the triangle index
        triangleIdx += 1;

        // Move to the next half edge
        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != halfedge->previous);
}


toolkit::Plane dcel::Face::Plane(bool normalize, bool fixZeros, double threshold) const
{
    return toolkit::Plane(Normal(normalize, fixZeros, threshold), Centroid(fixZeros, threshold));
}

int dcel::Face::PointLocation(const Eigen::Vector3d & P, double threshold) const
{
    return Plane(false, true, threshold).PointLocation(P, threshold);
}

VF dcel::Face::vf() const
{
    // Get the number of vertices of the face. It is the same number of sides of the face
    size_t nVertices = CountSides(), index = 0;

    // Initialize the vertex coordinates and vertex indices of the face. Then, resize its vectors
    VF vf(nVertices, 1);

    // Initialize the vector for the vertex indices of the face. Then, get its reference and resize
    // it
    std::vector<size_t> indices(nVertices);

    // Ge the incident half edge of the face
    std::shared_ptr<Halfedge> currentHalfedge = halfedge;

    // Traverse through the half edges of the face
    do
    {
        // Make a copy of the coordinates from the start vertex of the half edge
        vf.addVertex(currentHalfedge->start->Coords());

        // Set the index attribute for the start vertex of the current half edge, then increase the 
        // index value
        indices[index++] = index;

        // Move to the next half edge
        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != halfedge);

    vf.addFace(indices);

    // Return the object with the vertices and vertex indices
    return vf;
}

dcel::Halfedge::Halfedge() :
    start(nullptr),
    previous(nullptr),
    next(nullptr),
    twin(nullptr),
    face(nullptr)
{
}

dcel::Halfedge::~Halfedge()
{
    Clear();
}

void dcel::Halfedge::CheckConsistency() const
{
    // Must have a pointer to a start vertex
    assert(start);

    // Must have a pointer to a twin half edge
    assert(twin);

    // this must be the twin's tein
    assert(shared_from_this() == twin->twin);

    // Start vertex must be twin's end vertex
    assert(start == twin->twin->start);

    // This must be prevous' next half edge
    assert(shared_from_this() == previous->next);

    // This must be next's previous half edge
    assert(shared_from_this() == next->previous);

    // Incident face must be the same as previous' face
    assert(face == previous->face);

    // Incident face must be the same as next's face
    assert(face == next->face);
}

Eigen::Vector3d dcel::Halfedge::Direction(bool normalize, bool fixZeros, double threshold) const
{
    // Calculate the direction vector of the half edge
    Eigen::Vector3d D = twin->start->Coords() - start->Coords();

    // Normalize the direction vector if indicated
    if (normalize)
    {
        D.normalize();
    }

    // Fix the zero values of the direction vector if indicated
    if (fixZeros)
    {
        utils::fixZeros(D, threshold);
    }

    // Return the direction vector
    return D;
}

void dcel::Halfedge::Clear()
{
    m_attributes.Clear();
    start = nullptr;
    previous = nullptr;
    next = nullptr;
    twin = nullptr;
    face = nullptr;
}

toolkit::Ray dcel::Halfedge::Ray(bool normalize) const
{
    return toolkit::Ray(
        start->Coords(), 
        twin->start->Coords() - start->Coords(), 
        normalize);
}

bool dcel::Halfedge::Intersect(
    std::shared_ptr<Halfedge> H,
    Eigen::Vector3d & P,
    double & s,
    double & t,
    double threshold) const
{
    // 
    Eigen::Vector3d t1 = start->Coords() - H->start->Coords();
    Eigen::Vector3d t2 = H->Direction();

    // 
    if (abs(t2.x()) < threshold && abs(t2.y()) < threshold && abs(t2.z()) < threshold)
    {
        return false;
    }

    // 
    Eigen::Vector3d t3 = Direction();

    //
    if (abs(t3.x()) < threshold && abs(t3.y()) < threshold && abs(t3.z()) < threshold)
    {
        return false;
    }

    //
    double d1343 = (t1.x() * t2.x()) + (t1.y() * t2.y()) + (t1.z() * t2.z());
    double d4321 = (t2.x() * t3.x()) + (t2.y() * t3.y()) + (t2.z() * t3.z());
    double d1321 = (t1.x() * t3.x()) + (t1.y() * t3.y()) + (t1.z() * t3.z());
    double d4343 = (t2.x() * t2.x()) + (t2.y() * t2.y()) + (t2.z() * t2.z());
    double d2121 = (t3.x() * t3.x()) + (t3.y() * t3.y()) + (t3.z() * t3.z());

    //
    double denom = (d2121 * d4343) - (d4321 * d4321);

    // 
    if (abs(denom) < threshold)
    {
        return false;
    }

    //
    double numer = (d1343 * d4321) - (d1321 * d4343);

    // Calculate the parameters for points Pa and Pb
    s = numer / denom;
    t = (d1343 + (d4321 * s)) / d4343;

    // Calculate points Pa and Pb
    Eigen::Vector3d pa = (t3 * s) + start->Coords();
    Eigen::Vector3d pb = (t2 * t) + H->start->Coords();

    // Calculate the midpoint between points Pa and Pb
    P << (pa + pb) / 2.0;

    // Indicate there is a result and set M as the intersection point between the half edges
    return true;
}

double dcel::Halfedge::Length() const
{
    return Direction().norm();
}

Eigen::Vector3d dcel::Halfedge::Lerp(double t) const
{
    return start->Coords() + ((twin->start->Coords() - start->Coords()) * t);
}

toolkit::LineSegment dcel::Halfedge::LineSegment() const
{
    return toolkit::LineSegment(start->Coords(), twin->start->Coords());
}

Eigen::Vector3d dcel::Halfedge::Midpoint(bool fixZeros, double threshold) const
{
    // Calculate the midpoint of the half edge
    Eigen::Vector3d M = (start->Coords() + twin->start->Coords()) / 2.0;

    // Fix the zero values if indicated
    if (fixZeros)
    {
        utils::fixZeros(M, threshold);
    }

    // Return the midpoint of the half edge
    return M;
}

Eigen::Vector3d dcel::Halfedge::Normal(bool normalize, bool fixZeros, double threshold) const
{
    // Get the normal vector of the face incident to the half edge. If no face then get a 0 vector
    Eigen::Vector3d N1 = (face) ? face->Normal() : Eigen::Vector3d(0.0, 0.0, 0.0);

    // Get the normal vector of the face incident to the twin half edge. If no face then get a 0 
    // vector
    Eigen::Vector3d N2 = (twin->face) ? twin->face->Normal() : Eigen::Vector3d(0.0, 0.0, 0.0);

    // Calculate the average between both normal vectors
    Eigen::Vector3d N = (N1 + N2) / 2.0;

    // Normalize the normal vector if indicated
    if (normalize)
    {
        N.normalize();
    }

    // Fix the zero values of the normal vector if indicated
    if (fixZeros)
    {
        utils::fixZeros(N, threshold);
    }

    // Return the normal vector
    return N;
}

bool dcel::Halfedge::IsPointIn(const Eigen::Vector3d & P, double threshold) const
{
    // Get the vector from the start vertex of the half edge and point P
    Eigen::Vector3d AP = P - start->Coords();

    // Get the direction vector of the half edge
    Eigen::Vector3d AB = Direction();

    // Get the cross product between AP and the direction vector of the half edge
    Eigen::Vector3d test = AP.cross(AB);

    // Fix the zeros of the test vector
    utils::fixZeros(test, threshold);

    // If the test result is not the zero vector then P does not lie along the line defined by the 
    // end points of the half edge
    if (!(test(0) == 0.0 && test(1) == 0.0 && test(2) == 0.0))
    {
        return false;
    }

    // Get the dot product between AP and AB, and AB with itself
    double dotAP = AP.dot(AB);
    double dotAB = AB.dot(AB);

    // Fix zeros
    utils::fixZero(dotAP, threshold);
    utils::fixZero(dotAB, threshold);

    // If dotAP is equal to 0 then P matches with the start vertex of the half edge. If it is equal 
    // to dotAB then P matches with the end vertex of the half edge. If it is greater than zero but
    // less than dotAB then P lies in the line segment between the end points of the half edge. In 
    // such cases return true. Otherwise, return false
    return (dotAP == 0.0 || dotAP == dotAB || (dotAP > 0.0 && dotAP < dotAB));
}

dcel::DCEL::DCEL() 
{
}

dcel::DCEL::DCEL(const VF & vf)
{
    Set(vf);
}

dcel::DCEL::DCEL(const DCEL & dcel)
{
    Set(dcel);
}

dcel::DCEL::DCEL(const std::string & filename)
{
    LoadDcelFile(filename);
}

dcel::DCEL::~DCEL()
{
    Clear();
}

bool dcel::DCEL::AreFacesEvenSided() const
{
    // Traverse through the faces of the DCEL and count their number of sides
    for (auto itF = m_faces.begin(); itF != m_faces.end(); ++itF)
    {
        // If the face has an odd number of sides then return false
        if ((*itF)->CountSides() % 2 == 1)
        {
            return false;
        }
    }

    // Return true since all faces have an even number of sides
    return true;
}

void dcel::DCEL::AxisAlignedBoundingBox(Eigen::Vector3d & min, Eigen::Vector3d & max) const
{
    // Initialize both min and max points using the coordinates of the first vertex of the geometry
    min << m_vertices[0]->Coords();
    max << min;

    size_t nVertices = m_vertices.size();

    Eigen::Vector3d P = Eigen::Vector3d::Zero();

    // Traverse through the vertices of the geometry and update the corner points respectively
    for (size_t i = 1; i < nVertices; i += 1)
    {
        P << m_vertices[i]->Coords();

        if (P.x() < min.x())
        {
            min(0) = P.x();
        }

        if (P.y() < min.y())
        {
            min(1) = P.y();
        }

        if (P.z() < min.z())
        {
            min(2) = P.z();
        }

        if (P.x() > max.x())
        {
            max(0) = P.x();
        }

        if (P.y() > max.y())
        {
            max(1) = P.y();
        }

        if (P.z() > max.z())
        {
            max(2) = P.z();
        }
    }
}

Eigen::Vector3d dcel::DCEL::Centroid(bool fixZeros, double threshold) const
{
    // Initialize the variable for calculating the volume of the polyhedron
    double V = 0;

    // Initialize the object for storing the centroid coordinates
    Eigen::Vector3d C(0.0, 0.0, 0.0);

    // Traverse through the faces of the geometry
    for (auto it = m_faces.begin(); it != m_faces.end(); ++it)
    {
        // Get the pointer to the current face
        std::shared_ptr<Face> face = *it;

        // Get the reference to the coordinates of the start vertex of the incident half edge of 
        // the current face
        const Eigen::Vector3d & v0 = face->halfedge->start->Coords();

        // Get the next half edge to the incident half edge of the current face
        std::shared_ptr<Halfedge> currentHalfedge = face->halfedge->next;

        // Traverse through the half edges of the polyhedron. Stop at the previous half edge to the 
        // incident half edge of the current face. By doing this we traverse the current face by 
        // triangles
        do
        {
            // Get the references to the coordinates of the start and end vertices of the current 
            // half edge
            const Eigen::Vector3d & v1 = currentHalfedge->start->Coords();
            const Eigen::Vector3d & v2 = currentHalfedge->twin->start->Coords();

            // Calculate the normal vector of the triangle formed by v0, v1 and v2
            Eigen::Vector3d N = (v1 - v0).cross(v2 - v0);

            // Update the volume of the polyhedron
            V += (v0.dot(N) / 6.0);

            // Update the centroid coordinates
            C(0) += N.x() * (pow(v0.x() + v1.x(), 2) + pow(v1.x() + v2.x(), 2) + pow(v2.x() + v0.x(), 2));
            C(1) += N.y() * (pow(v0.y() + v1.y(), 2) + pow(v1.y() + v2.y(), 2) + pow(v2.y() + v0.y(), 2));
            C(2) += N.z() * (pow(v0.z() + v1.z(), 2) + pow(v1.z() + v2.z(), 2) + pow(v2.z() + v0.z(), 2));

            // Move to the next half edge in the face
            currentHalfedge = currentHalfedge->next;

        } while (currentHalfedge != face->halfedge->previous);
    }

    // Update the centroid coordinates
    C *= (1.0 / (24.0 * 2.0 * V));

    // Fix the zero values of the centroid if indicated
    if (fixZeros)
    {
        utils::fixZeros(C, threshold);
    }

    return C;
}

void dcel::DCEL::CheckConsistency() const
{
    // Traverse through the vertices of the geometry and check their consistency
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        (*it)->CheckConsistency();
    }

    // Traverse through the faces of the geometry and check their consistency
    for (auto it = m_faces.begin(); it != m_faces.end(); ++it)
    {
        (*it)->CheckConsistency();
    }

    // Traverse through the half edges of the geometry and check their consistency
    for (auto it = m_halfedges.begin(); it != m_halfedges.end(); ++it)
    {
        (*it)->CheckConsistency();
    }
}

void dcel::DCEL::Clear()
{
    m_attributes.Clear();

    size_t nV = m_vertices.size(), nF = m_faces.size(), nH = m_halfedges.size(), i = 0;

    for (i = 0; i < nV; i += 1)
    {
        m_vertices[i]->Clear();
        m_vertices[i] = nullptr;
    }

    m_vertices.clear();

    for (i = 0; i < nF; i += 1)
    {
        m_faces[i]->Clear();
        m_faces[i] = nullptr;
    }

    m_faces.clear();

    for (i = 0; i < nH; i += 1)
    {
        m_halfedges[i]->Clear();
        m_halfedges[i] = nullptr;
    }

    m_halfedges.clear();
}

/*void dcel::DCEL::Clip(const toolkit::Plane & plane, VF & pvf, VF & nvf, double threshold) const
{
    // Reset the vertex coordinates and vertex indices for both geometries
    pvf.reset();
    nvf.reset();

    // Set the location of the vertices with respect to the clipping plane, location values are 
    // stored at each vertex in the ATTRIB_LOCATION dynamic attribute
    SetVerticesLocation(plane, threshold);

    // Traverse through the vertices of the original geometry and add them to the respective 
    // clipped geometry according to its location with respect to the clipping plane
    for (auto it = vertices.begin(); it != vertices.end(); ++it) 
    {
        // Get the pointer to the current vertex
        std::shared_ptr<Vertex> V = *it;

        // Get the location value of the vertex
        int location;
        assert(V->Attributes().Get<int>(ATTRIB_LOCATION, location));
        assert(location == 0 || location == 1 || location == -1);

        // Add the vertex in the respective geometry according to its location with respect to the
        // clipping plane. Store its index since we need it when defining the vertex indices of the
        // clipped geometries
        switch (location) 
        {
            // Vertex is located at the plane, it must be added to both positive and negative 
            // geometries
            case 0: 
            {
                size_t positiveIndex = pvf.AddVertex(V->Coords());
                size_t negativeIndex = nvf.AddVertex(V->Coords());
                V->Attributes().Set<size_t>(ATTRIB_POSITIVE_INDEX, positiveIndex);
                V->Attributes().Set<size_t>(ATTRIB_NEGATIVE_INDEX, negativeIndex);
                break;
            }

            // Vertex is located at the positive half space, it must be added to the positive 
            // geometry
            case 1: 
            {
                size_t positiveIndex = pvf.AddVertex(V->Coords());
                V->Attributes().Set<size_t>(ATTRIB_POSITIVE_INDEX, positiveIndex);
                break;
            }

            // Vertex is located at the negative half space, it must be added to the negative 
            // geometry
            case -1: 
            {
                size_t negativeIndex = nvf.AddVertex(V->Coords());
                V->Attributes().Set<size_t>(ATTRIB_NEGATIVE_INDEX, negativeIndex);
                break;
            }
        }
    }

    // Set all half edges as not visited
    SetHalfedgesAttribute(ATTRIB_VISITED, false);

    // Traverse through the half edges of the original geometry and check if the plane intersects 
    // with each one of them
    for (auto it = halfedges.begin(); it != halfedges.end(); ++it) 
    {
        // Get the pointer to the current half edge
        std::shared_ptr<Halfedge> H = *it;

        // Get the visited attribute of the half edge
        bool visited;
        assert(H->Attributes().Get<bool>(ATTRIB_VISITED, visited));

        // If the current half edge has been visited then continue with the next half edge
        if (visited) 
        {
            continue;
        }

        // Get the location of the end points of the half edge
        int startLocation, endLocation;
        assert(H->start->Attributes().Get<int>(ATTRIB_LOCATION, startLocation));
        assert(H->twin->start->Attributes().Get<int>(ATTRIB_LOCATION, endLocation));
        assert(startLocation == 0 || startLocation == 1 || startLocation == -1);
        assert(endLocation == 0 || endLocation == 1 || endLocation == -1);

        // If the endpoints lie at different half spaces then calculate the intersection between 
        // the plane and the half edge, then store it in both clipped geometries and keep the point
        // indices in the current half edge and its twin
        if (startLocation * endLocation == -1) 
        {
            // Get the line segment representin the current half edge
            const toolkit::LineSegment linesegment = H->LineSegment();

            // Calculate the intersection parameter between the plane and the line segment
            double t;
            assert(utils::PlaneLineSegmentIntersection(plane, linesegment, t, threshold));

            // Calculate the coordinates of the intersection point
            const Eigen::Vector3d intersection = linesegment.At(t);

            // Store the intersection point in both clipped geometries, keep the respective indices
            size_t positiveIndex = pvf.AddVertex(intersection);
            size_t negativeIndex = nvf.AddVertex(intersection);

            // Store the indices in the current half edge and its twin
            H->Attributes().Set<size_t>(ATTRIB_POSITIVE_INDEX, positiveIndex);
            H->Attributes().Set<size_t>(ATTRIB_NEGATIVE_INDEX, negativeIndex);
            H->twin->Attributes().Set<size_t>(ATTRIB_POSITIVE_INDEX, positiveIndex);
            H->twin->Attributes().Set<size_t>(ATTRIB_NEGATIVE_INDEX, negativeIndex);
        }

        // Set the current half edge and its twin as visited
        H->Attributes().Set<bool>(ATTRIB_VISITED, true);
        H->twin->Attributes().Set<bool>(ATTRIB_VISITED, true);
    }

    // Initialize a flag that indicates whether a face made exclusively of intersection points has
    // been defined or not
    bool intersectionFace = true;

    // Traverse through the faces of the original geometry and generate the vertex indices of the 
    // clipped geometries
    for (auto it = faces.begin(); it != faces.end(); ++it) 
    {
        // Get the pointer to the current face
        std::shared_ptr<Face> F = *it;

        // Get the pointer to the incident half edge of the current face
        std::shared_ptr<Halfedge> currentHalfedge = F->halfedge;

        // Initialize the vectors for storing the vertex indices of the clipped faces
        std::vector<size_t> positiveIndices, negativeIndices;

        // Initialize a flag that indicates whether the current face is made exclusively of 
        // intersection points or not
        bool allIntersectionPoints = true;

        // Traverse through the half edges of the current face
        do 
        {
            // Get the location of the end points of the current half edge
            int startLocation, endLocation;
            assert(currentHalfedge->start->Attributes().Get<int>(ATTRIB_LOCATION, startLocation));
            assert(currentHalfedge->twin->start->Attributes().Get<int>(ATTRIB_LOCATION, endLocation));
            assert(startLocation == 0 || startLocation == 1 || startLocation == -1);
            assert(endLocation == 0 || endLocation == 1 || endLocation == -1);

            // Add the vertex index in the respective geometry according to its location with 
            // respect to the clipping plane
            switch (startLocation) 
            {
                // Vertex is located at the plane, it must be added to both positive and negative 
                // geometries
                case 0: 
                {
                    size_t positiveIndex, negativeIndex;
                    assert(currentHalfedge->start->Attributes().Get<size_t>(ATTRIB_POSITIVE_INDEX, positiveIndex));
                    assert(currentHalfedge->start->Attributes().Get<size_t>(ATTRIB_NEGATIVE_INDEX, negativeIndex));
                    positiveIndices.push_back(positiveIndex);
                    negativeIndices.push_back(negativeIndex);
                    break;
                }

                // Vertex is located at the positive half space, it must be added to the positive 
                // geometry
                case 1: 
                {
                    size_t positiveIndex;
                    assert(currentHalfedge->start->Attributes().Get<size_t>(ATTRIB_POSITIVE_INDEX, positiveIndex));
                    positiveIndices.push_back(positiveIndex);

                    // The vertex is not at the intersection plane, so the face is not made 
                    // exclusively of intersection points
                    allIntersectionPoints = false;

                    break;
                }

                // Vertex is located at the negative half space, it must be added to the negative 
                // geometry
                case -1: 
                {
                    size_t negativeIndex;
                    assert(currentHalfedge->start->Attributes().Get<size_t>(ATTRIB_NEGATIVE_INDEX, negativeIndex));
                    negativeIndices.push_back(negativeIndex);

                    // The vertex is not at the intersection plane, so the face is not made 
                    // exclusively of intersection points
                    allIntersectionPoints = false;

                    break;
                }
            }

            // If both end points lie at different half spaces then get the vertex index of the
            // intersection point to both positive and negative geometries
            if (startLocation * endLocation == -1) 
            {
                size_t positiveIndex, negativeIndex;
                assert(currentHalfedge->Attributes().Get<size_t>(ATTRIB_POSITIVE_INDEX, positiveIndex));
                assert(currentHalfedge->Attributes().Get<size_t>(ATTRIB_NEGATIVE_INDEX, negativeIndex));
                positiveIndices.push_back(positiveIndex);
                negativeIndices.push_back(negativeIndex);
            }

            // Move to the next half edge in the face
            currentHalfedge = currentHalfedge->next;

        } while (currentHalfedge != F->halfedge);

        // Store the vertex indices respecitvely
        if (positiveIndices.size() > 0) 
        {
            pvf.F.push_back(positiveIndices);
        }

        if (negativeIndices.size() > 0) 
        {
            nvf.F.push_back(negativeIndices);
        }

        // Update the intersection face flag
        intersectionFace = intersectionFace && allIntersectionPoints;
    }

    // If no face of the clipped geometries is made exclusively of intersection points then we 
    // need to generate them
    if (!intersectionFace) 
    {
        // Set the half edges of the original geometry as not visited
        SetHalfedgesAttribute(ATTRIB_VISITED, false);

        // 
        std::shared_ptr<Halfedge> intersectHalfedge = nullptr;

        // 
        int startLocation, endLocation;

        // Traverse through the half edges of the original geometry and keep the pointer to the 
        // first half edge that intersects the clipping plane
        for (auto it = halfedges.begin(); it != halfedges.end(); it += 1) 
        {
            // Get the pointer to the current half edge
            std::shared_ptr<Halfedge> H = *it;

            // Get the visited attribute of the half edge
            bool visited;
            assert(H->Attributes().Get<bool>(ATTRIB_VISITED, visited));

            // If the current half edge has been visited then continue with the next one
            if (visited) 
            {
                continue;
            }

            // Get the location of the end points of the half edge
            assert(H->start->Attributes().Get<int>(ATTRIB_LOCATION, startLocation));
            assert(H->twin->start->Attributes().Get<int>(ATTRIB_LOCATION, endLocation));
            assert(startLocation == 0 || startLocation == 1 || startLocation == -1);
            assert(endLocation == 0 || endLocation == 1 || endLocation == -1);

            // If the current half edge has both end points at different half spaces then exit the 
            // for loop since we don't need to keep searching
            if (startLocation * endLocation == -1) 
            {
                intersectHalfedge = H;
                break;
            }
            else 
            {
                H->Attributes().Set<bool>(ATTRIB_VISITED, true);
                H->twin->Attributes().Set<bool>(ATTRIB_VISITED, true);
            }
        }

        // assert(H);

        // 
        //if (startLocation * endLocation == -1) 
        if (intersectHalfedge)
        {
            // Be sure we start at a half edge with intersection point
            assert(intersectHalfedge->start->Attributes().Get<int>(ATTRIB_LOCATION, startLocation));
            assert(intersectHalfedge->twin->start->Attributes().Get<int>(ATTRIB_LOCATION, endLocation));
            assert(startLocation == 0 || startLocation == 1 || startLocation == -1);
            assert(endLocation == 0 || endLocation == 1 || endLocation == -1);
            assert(startLocation * endLocation == -1);

            // Initialize the vectors for storing the vertex indices of the clipped faces
            std::vector<size_t> positiveIndices, negativeIndices;

            // 
            size_t positiveIndex, negativeIndex;
            assert(intersectHalfedge->Attributes().Get<size_t>(ATTRIB_POSITIVE_INDEX, positiveIndex));
            assert(intersectHalfedge->Attributes().Get<size_t>(ATTRIB_NEGATIVE_INDEX, negativeIndex));
            positiveIndices.push_back(positiveIndex);
            negativeIndices.push_back(negativeIndex);
            
            // Make a copy of the coordinates of the intersection point. We need it to define the 
            // orientation of the final faces for the clipped geometries
            Eigen::Vector3d lastP(pvf.V[positiveIndex]);

            // Get the pointer to the twin half edge of the intersected half edge
            std::shared_ptr<Halfedge> currentHalfedge = intersectHalfedge->twin->next;

            // 
            do
            {
                // Get the location of the end points of the half edge
                assert(currentHalfedge->start->Attributes().Get<int>(ATTRIB_LOCATION, startLocation));
                assert(currentHalfedge->twin->start->Attributes().Get<int>(ATTRIB_LOCATION, endLocation));
                assert(startLocation == 0 || startLocation == 1 || startLocation == -1);
                assert(endLocation == 0 || endLocation == 1 || endLocation == -1);

                // If the current half edge has both end points at different half spaces then add 
                // its vertex index to the last faces
                if (startLocation * endLocation == -1)
                {
                    // 
                    assert(currentHalfedge->Attributes().Get<size_t>(ATTRIB_POSITIVE_INDEX, positiveIndex));
                    assert(currentHalfedge->Attributes().Get<size_t>(ATTRIB_NEGATIVE_INDEX, negativeIndex));

                    // Get the reference to the intersection point
                    const Eigen::Vector3d & currentP = pvf.V[positiveIndex];

                    // Check if the last intersection point and the current one are making a left 
                    // turn with respect to the clipping plane
                    int isLeft = utils::IsLeftTurn(lastP, currentP, plane, threshold);
                    assert(isLeft != 0);

                    // 
                    if (isLeft == 1)
                    {
                        positiveIndices.insert(positiveIndices.begin(), positiveIndex);
                        negativeIndices.push_back(negativeIndex);
                    }
                    else
                    {
                        positiveIndices.push_back(positiveIndex);
                        negativeIndices.insert(negativeIndices.begin(), negativeIndex);
                    }

                    // Update the last intersection point
                    lastP << currentP;

                    // Move to the next half edge of the twin half edge
                    currentHalfedge = currentHalfedge->twin->next;
                }
                else
                {
                    // Move to the next half edge
                    currentHalfedge = currentHalfedge->next;
                }

            } while (currentHalfedge != intersectHalfedge);

            // Store the vertex indices respecitvely
            pvf.F.push_back(positiveIndices);
            nvf.F.push_back(negativeIndices);
        }
    }

    // Remove the dynamic attributes from the original geometry
    RemoveVerticesAttribute(ATTRIB_LOCATION);
    RemoveVerticesAttribute(ATTRIB_POSITIVE_INDEX);
    RemoveVerticesAttribute(ATTRIB_NEGATIVE_INDEX);
    RemoveHalfedgesAttribute(ATTRIB_VISITED);
    RemoveHalfedgesAttribute(ATTRIB_POSITIVE_INDEX);
    RemoveHalfedgesAttribute(ATTRIB_NEGATIVE_INDEX);
}*/

VF dcel::DCEL::Clip(const toolkit::Plane & plane, const double threshold) const
{
    // Initialize the lists to store the vertices and faces of the clipped geometry
    std::list<Eigen::Vector3d> clippedVertices;
    std::list<std::vector<size_t>> clippedFaces;

    //
    ClippedElements(plane, clippedVertices, clippedFaces, threshold);

    VF vf(clippedVertices, clippedFaces);

    //
    clippedVertices.clear();
    clippedFaces.clear();

    return vf;
}

void dcel::DCEL::ClippedFaces(
    std::list<std::vector<size_t>> & clippedFaces, 
    bool & faceAtClippingPlane, 
    double threshold) const
{
    clippedFaces.clear();

    faceAtClippingPlane = true;

    bool allIntersectionPoints = true;

    int startLocation = 0, endLocation = 0;

    size_t index = 0;

    // Traverse through the faces of the original geometry and generate the faces of the clipped 
    // geometry
    for (auto it = m_faces.begin(); it != m_faces.end(); ++it)
    {
        // Get the pointer to the current face
        std::shared_ptr<Face> F = *it;

        // Get the pointer to the incident half edge of the current face
        std::shared_ptr<Halfedge> currentHalfedge = F->halfedge;

        // Initialize the vectors for storing the vertex indices of the clipped faces
        std::vector<size_t> indices;

        // Initialize a flag that indicates whether the current face is made exclusively of 
        // intersection points or not
        allIntersectionPoints = true;

        // Traverse through the half edges of the current face
        do
        {
            // Get the location of the end points of the current half edge
            assert(currentHalfedge->start->Attributes().Get<int>(ATTRIB_LOCATION, startLocation));
            assert(currentHalfedge->twin->start->Attributes().Get<int>(ATTRIB_LOCATION, endLocation));
            assert(startLocation == 0 || startLocation == 1 || startLocation == -1);
            assert(endLocation == 0 || endLocation == 1 || endLocation == -1);

            // Add the vertex index in the respective geometry according to its location with 
            // respect to the clipping plane
            if (startLocation == 0 || startLocation == -1)
            {
                assert(currentHalfedge->start->Attributes().Get<size_t>(ATTRIB_INDEX, index));
                indices.push_back(index);

                // The vertex is not at the intersection plane, so the face is not made exclusively
                // of intersection points
                if (startLocation == -1)
                {
                    allIntersectionPoints = false;
                }
            }

            // If both end points lie at different half spaces then get the vertex index of the
            // intersection point to both positive and negative geometries
            if (startLocation * endLocation == -1)
            {
                assert(currentHalfedge->Attributes().Get<size_t>(ATTRIB_INDEX, index));
                indices.push_back(index);
            }

            // Move to the next half edge in the face
            currentHalfedge = currentHalfedge->next;

        } while (currentHalfedge != F->halfedge);

        // Store the vertex indices respecitvely
        if (indices.size() > 0)
        {
            clippedFaces.push_back(indices);
        }

        // Update the intersection face flag
        faceAtClippingPlane = faceAtClippingPlane && allIntersectionPoints;
    }
}

void dcel::DCEL::ClippedElements(
    const toolkit::Plane & plane, 
    std::list<Eigen::Vector3d> & clippedVertices, 
    std::list<std::vector<size_t>> & clippedFaces, 
    double threshold) const
{
    // Initialize a flag that indicates there are intersection points between the geometry and the 
    // plane
    size_t pointsInPlane = 0;

    // Populate the vector with the clipped vertices that lie at the positive half space of the 
    // clipping plane
    ClippedVertices(plane, clippedVertices, pointsInPlane, threshold);

    // Initialize a flag that indicates whether a face made exclusively of intersection points has
    // been defined or not
    bool faceAtClippingPlane = false, visited = false;

    // 
    ClippedFaces(clippedFaces, faceAtClippingPlane, threshold);

    int startLocation, endLocation;

    // Check if there are more than two clipped vertices that lie at the clipping plane and there 
    // is no face made of such vertices. If so then we need to define such face
    if (pointsInPlane > 2 && !faceAtClippingPlane)
    {
        // Set the half edges of the original geometry as not visited
        SetHalfedgesAttribute(ATTRIB_VISITED, false);

        // 
        std::shared_ptr<Halfedge> intersectHalfedge = nullptr;

        // Traverse through the half edges of the original geometry and keep the pointer to the 
        // first half edge that intersects the clipping plane
        for (auto it = m_halfedges.begin(); it != m_halfedges.end(); it += 1)
        {
            // Get the pointer to the current half edge
            std::shared_ptr<Halfedge> H = *it;

            // Get the visited attribute of the half edge
            assert(H->Attributes().Get<bool>(ATTRIB_VISITED, visited));

            // If the current half edge has been visited then continue with the next one
            if (visited)
            {
                continue;
            }

            // Get the location of the end points of the half edge
            assert(H->start->Attributes().Get<int>(ATTRIB_LOCATION, startLocation));
            assert(H->twin->start->Attributes().Get<int>(ATTRIB_LOCATION, endLocation));
            assert(startLocation == 0 || startLocation == 1 || startLocation == -1);
            assert(endLocation == 0 || endLocation == 1 || endLocation == -1);

            // If the current half edge has both end points at different half spaces then exit the 
            // for loop since we don't need to keep searching
            if (startLocation * endLocation == -1)
            {
                // 
                double dot = H->Direction().dot(plane.Normal());

                // 
                if (abs(dot) > threshold)
                {
                    intersectHalfedge = (dot > 0) ? H : H->twin;
                    break;
                }
            }

            // 
            H->Attributes().Set<bool>(ATTRIB_VISITED, true);
            H->twin->Attributes().Set<bool>(ATTRIB_VISITED, true);
        }

        // Pay attention to this part!!!
        // THERE MUST EXIST A HALF EDGE HERE
        if (intersectHalfedge) 
        {
            // Be sure we start at a half edge with intersection point
            assert(intersectHalfedge->start->Attributes().Get<int>(ATTRIB_LOCATION, startLocation));
            assert(intersectHalfedge->twin->start->Attributes().Get<int>(ATTRIB_LOCATION, endLocation));
            assert(startLocation == 0 || startLocation == 1 || startLocation == -1);
            assert(endLocation == 0 || endLocation == 1 || endLocation == -1);
            assert(startLocation * endLocation == -1);

            // Initialize the vectors for storing the vertex indices of the clipped faces
            std::vector<size_t> indices;

            // 
            size_t index;
            assert(intersectHalfedge->Attributes().Get<size_t>(ATTRIB_INDEX, index));
            indices.push_back(index);

            // Make a copy of the coordinates of the intersection point. We need it to define the 
            // orientation of the final faces for the clipped geometries
            //Eigen::Vector3d lastP(vf.V[index]);

            // Get the pointer to the twin half edge of the intersected half edge
            std::shared_ptr<Halfedge> currentHalfedge = intersectHalfedge->twin->next;

            // 
            do
            {
                // Get the location of the end points of the half edge
                assert(currentHalfedge->start->Attributes().Get<int>(ATTRIB_LOCATION, startLocation));
                assert(currentHalfedge->twin->start->Attributes().Get<int>(ATTRIB_LOCATION, endLocation));
                assert(startLocation == 0 || startLocation == 1 || startLocation == -1);
                assert(endLocation == 0 || endLocation == 1 || endLocation == -1);

                // If the current half edge has both end points at different half spaces then add 
                // its vertex index to the last faces
                if (startLocation * endLocation == -1)
                {
                    // 
                    assert(currentHalfedge->Attributes().Get<size_t>(ATTRIB_INDEX, index));

                    // Get the reference to the intersection point
                    //const Eigen::Vector3d & currentP = vf.V[index];

                    // Check if the last intersection point and the current one are making a left 
                    // turn with respect to the clipping plane
                    //int isLeft = utils::IsLeftTurn(lastP, currentP, plane, threshold);
                    //assert(isLeft != 0);

                    // 
                    //if (isLeft == 1)
                    //{
                    //    indices.insert(indices.begin(), index);
                    //}
                    //else
                    //{
                    //    indices.push_back(index);
                    //}

                    // Update the last intersection point
                    //lastP << currentP;

                    indices.push_back(index);

                    // Move to the next half edge of the twin half edge
                    currentHalfedge = currentHalfedge->twin->next;
                }
                else
                {
                    // Move to the next half edge
                    currentHalfedge = currentHalfedge->next;
                }

            } while (currentHalfedge != intersectHalfedge);

            // Store the vertex indices respecitvely
            clippedFaces.push_back(indices);
        }
    }

    // Remove the dynamic attributes from the original geometry
    RemoveVerticesAttribute(ATTRIB_LOCATION);
    RemoveVerticesAttribute(ATTRIB_INDEX);
    RemoveHalfedgesAttribute(ATTRIB_VISITED);
    RemoveHalfedgesAttribute(ATTRIB_INDEX);
}

void dcel::DCEL::ClippedVertices(
    const toolkit::Plane & plane, 
    std::list<Eigen::Vector3d> & clippedVertices, 
    size_t & pointsInPlane, 
    double threshold) const
{
    clippedVertices.clear();

    size_t vIdx = 0;

    pointsInPlane = 0;

    int location = 0, startLocation = 0, endLocation = 0;

    bool visited = false;

    double t = 0;

    std::shared_ptr<Vertex> V = nullptr;
    std::shared_ptr<Halfedge> H = nullptr;

    toolkit::LineSegment linesegment;

    // Set the location of the vertices with respect to the clipping plane, location values are 
    // stored at each vertex in the ATTRIB_LOCATION dynamic attribute
    SetVerticesLocation(plane, threshold);

    // Traverse through the vertices of the original geometry and add them to the clipped geometry 
    // according to its location with respect to the clipping plane
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        // Get the pointer to the current vertex
        V = *it;

        // Get the location value of the vertex
        assert(V->Attributes().Get<int>(ATTRIB_LOCATION, location));
        assert(location == 0 || location == 1 || location == -1);

        // Add the vertex in the respective geometry according to its location with respect to the
        // clipping plane. Store its index since we need it when defining the vertex indices of the
        // clipped geometries
        if (location == 0 || location == -1)
        {
            // Insert the coordinates of the vertex in the list of clipped vertices
            clippedVertices.push_back(V->Coords());

            // Store the index of the clipped vertex in the original vertex
            V->Attributes().Set<size_t>(ATTRIB_INDEX, vIdx++);

            // If the current vertex is at the clipping plane then indicate the plane touches a 
            // vertex of the original geometry
            if (location == 0)
            {
                pointsInPlane += 1;
            }
        }
    }

    SetHalfedgesAttribute(ATTRIB_VISITED, false);

    // Traverse through the half edges of the original geometry and check if the plane intersects 
    // with each one of them
    for (auto it = m_halfedges.begin(); it != m_halfedges.end(); ++it)
    {
        // Get the pointer to the current half edge
        H = *it;

        // Get the visited attribute of the half edge
        assert(H->Attributes().Get<bool>(ATTRIB_VISITED, visited));

        // If the current half edge has been visited then continue with the next half edge
        if (visited)
        {
            continue;
        }

        // Get the location of the end points of the half edge
        assert(H->start->Attributes().Get<int>(ATTRIB_LOCATION, startLocation));
        assert(H->twin->start->Attributes().Get<int>(ATTRIB_LOCATION, endLocation));
        assert(startLocation == 0 || startLocation == 1 || startLocation == -1);
        assert(endLocation == 0 || endLocation == 1 || endLocation == -1);

        // If the endpoints lie at different half spaces then calculate the intersection between 
        // the plane and the half edge, then store it in both clipped geometries and keep the point
        // indices in the current half edge and its twin
        if (startLocation * endLocation == -1)
        {
            // Get the line segment representing the current half edge
            linesegment = H->LineSegment();

            // Calculate the intersection parameter between the plane and the line segment
            assert(utils::planeLineSegmentIntersection(plane, linesegment, t, threshold));
            
            // Calculate the intersection point and store it in the list of clipped vertices
            clippedVertices.push_back(linesegment.At(t));

            // Indicate there is an intersection point between the geometry and the plane
            pointsInPlane += 1;

            // Store the indices in the current half edge and its twin
            H->Attributes().Set<size_t>(ATTRIB_INDEX, vIdx);
            H->twin->Attributes().Set<size_t>(ATTRIB_INDEX, vIdx++);
        }

        // Set the current half edge and its twin as visited
        H->Attributes().Set<bool>(ATTRIB_VISITED, true);
        H->twin->Attributes().Set<bool>(ATTRIB_VISITED, true);
    }
}

dcel::DCEL dcel::DCEL::Clone() const
{
    return DCEL(vf());
}

void dcel::DCEL::CloseLoops()
{
    // Traverse through the half edges
    for (auto itH = m_halfedges.begin(); itH != m_halfedges.end(); ++itH)
    {
        // If the current half edge has an incident face then move to the next half edge
        if ((*itH)->face)
        {
            continue;
        }

        // 
        std::vector<std::shared_ptr<Halfedge>>::iterator jtH = itH;

        // Traverse through the remaining half edges
        for (jtH = ++jtH; jtH != m_halfedges.end(); ++jtH)
        {
            // If the current half edge has an incident face then move to the next half edge
            if ((*jtH)->face)
            {
                continue;
            }

            // If the end vertex of the i-th half edge is the same as the start vertex of the j-th 
            // half edge and the i-th half edge does not have a next and the j-th half edge does 
            // not have a previous then link them together
            if ((*itH)->twin->start == (*jtH)->start && !(*itH)->next && !(*jtH)->previous)
            {
                (*itH)->next = (*jtH);
                (*jtH)->previous = (*itH);
            }

            // If the start vertex of the i-th half edge is the same as the end vertex of the j-th 
            // half edge and the i-th half edge does not have a previous and the j-th half edge 
            // does not have a next then link them together
            else if ((*itH)->start == (*jtH)->twin->start && !(*itH)->previous && !(*jtH)->next)
            {
                (*itH)->previous = (*jtH);
                (*jtH)->next = (*itH);
            }
        }
    }
}

dcel::DCEL dcel::DCEL::Copy() const
{
    return DCEL(vf());
}

VF dcel::DCEL::Dual() const
{
    // Initialize the lists to store the vertices and faces of the dual geometry
    std::list<Eigen::Vector3d> dualVertices;
    std::list<std::vector<size_t>> dualFaces;

    // 
    DualElements(dualVertices, dualFaces);

    // 
    VF vf(dualVertices, dualFaces);

    // 
    dualVertices.clear();
    dualFaces.clear();

    return vf;
	
}

void dcel::DCEL::DualElements(
    std::list<Eigen::Vector3d> & dualVertices, 
    std::list<std::vector<size_t>> & dualFaces) const
{
    // 
    dualVertices.clear();
    dualFaces.clear();

    // 
    size_t vIdx = 0;

    // Set the faces of the geometry as not visited. This use this attribute for keeping track of 
    // the faces with a dual vertex in the dual geometry
    SetFacesAttribute(ATTRIB_VISITED, false);

    // Traverse through the vertices of the geometry
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        // Get the pointer to the current vertex
        std::shared_ptr<Vertex> vertex = *it;

        // If the current vertex is not surrounded then continue with the next one
        if (!vertex->IsSurrounded())
        {
            continue;
        }

        // 
        std::vector<size_t> indices;

        // Get the pointer to the incident half edge of the vertex
        std::shared_ptr<Halfedge> currentHalfedge = vertex->halfedge;

        // Traverse through the half edges around the vertex and define the geometry of its 
        // respective dual face
        do
        {
            // Get the pointer to the incident face to the current half edge
            std::shared_ptr<Face> face = currentHalfedge->face;

            // Get the visited attribute value
            bool isVisited = false;
            assert(face->Attributes().Get<bool>(ATTRIB_VISITED, isVisited));

            // Initialize the index of the dual vertex of the face
            size_t dualIndex = 0;

            // If the current face is not visited then store the coordinates of its respective dual
            // vertex and store the respective info. Otherwise, get the index of the dual vertex
            if (!isVisited)
            {
                // Store the coordinates of the barycenter from the current face in the dual 
                // geometry. Then, indicate the face has been visited and store its dual vertex 
                // index
                //dualIndex = dual.AddVertex(face->Barycenter());
                dualVertices.push_back(face->Centroid());

                // Indicate the face has been visited (in other words, there exist its respective 
                // dual vertex coordinates in the dual geometry. Then, store the index of the dual 
                // vertex
                face->Attributes().Set<bool>(ATTRIB_VISITED, true);
                face->Attributes().Set<size_t>(ATTRIB_DUAL_INDEX, vIdx);

                dualIndex = vIdx++;
            }
            else
            {
                // Get the index of the dual vertex of the face
                assert(face->Attributes().Get<size_t>(ATTRIB_DUAL_INDEX, dualIndex));
            }

            // Insert the index of the dual vertex into the indices of the dual face (insert it at 
            // the begining of the indices vector)
            indices.insert(indices.begin(), dualIndex);

            // Move to the next half edge around the vertex
            currentHalfedge = currentHalfedge->twin->next;

        } while (currentHalfedge != vertex->halfedge);

        // 
        dualFaces.push_back(indices);
    }

    // Remove the visited and dual index dynamic attributes from the faces of the geometry
    RemoveFacesAttribute(ATTRIB_VISITED);
    RemoveFacesAttribute(ATTRIB_DUAL_INDEX);
}

const std::vector<std::shared_ptr<dcel::Face>> & dcel::DCEL::Faces() const
{
    return m_faces;
}

const std::vector<std::shared_ptr<dcel::Halfedge>>& dcel::DCEL::Halfedges() const
{
    return m_halfedges;
}

/*void dcel::DCEL::Flip() 
{
    std::shared_ptr<Halfedge> t = nullptr;

    for (auto it = m_halfedges.begin(); it != m_halfedges.end(); ++it) 
    {
        t = (*it)->previous;
        (*it)->previous = (*it)->next;
        (*it)->next = t;
    }
}*/

bool dcel::DCEL::Intersects(const DCEL & dcel, double threshold) const
{
    //assert();
    return false;
}

bool dcel::DCEL::Intersects(const toolkit::Plane & plane, double threshold) const
{
    // Set the location of the vertices with respect of the plane
    SetVerticesLocation(plane, threshold);

    // Set the half edges of the geometry as not visited
    SetHalfedgesAttribute(ATTRIB_VISITED, false);

    // Initialize the intersection indicator. It will be true if at least one vertex is in the 
    // plane, or the end points of a half edge lie in different half spaces defined by the plane
    bool intersects = false;

    // Traverse through the half edges and check if any of them has been split. If so, stop 
    // searching since the plane intersects the geometry
    for (auto it = m_halfedges.begin(); it != m_halfedges.end(); ++it)
    {
        // Get the pointer to the current half edge
        std::shared_ptr<Halfedge> H = *it;

        // Get the visited value of the half edge
        bool visited;
        H->Attributes().Get<bool>(ATTRIB_VISITED, visited);

        // If the half edge has been visited then continue with the next one
        if (visited) 
        {
            continue;
        }

        // Get the location value from the end points of the half edge
        int locA, locB;
        assert(H->start->Attributes().Get<int>(ATTRIB_LOCATION, locA));
        assert(H->twin->start->Attributes().Get<int>(ATTRIB_LOCATION, locB));

        // If at least one of the vertices is at the plane, or they are at different half spaces 
        // then the plane intersects the geometry
        if (locA * locB <= 0) 
        {
            intersects = true;
            break;
        }

        // Set the current half edge and its twin as visited
        H->Attributes().Set<bool>(ATTRIB_VISITED, true);
        H->twin->Attributes().Set<bool>(ATTRIB_VISITED, true);
    }

    // Remove the dynamic attributes
    RemoveVerticesAttribute(ATTRIB_LOCATION);
    RemoveHalfedgesAttribute(ATTRIB_VISITED);

    // Return the intersection indicator
    return intersects;
}

bool dcel::DCEL::IsPointIn(const Eigen::Vector3d & point, double threshold) const
{
    // Traverse through the faces and check the location of the point with respect of them
    for (auto it = m_faces.begin(); it != m_faces.end(); ++it)
    {
        // Get the pointer to the current face
        std::shared_ptr<Face> F = *it;

        // Get the location of the point with respect of the plane of the face
        int location = F->PointLocation(point, threshold);

        // If the point is in front of the plane then it is not inside or in the geometry, then 
        // return false
        if (location > 0) 
        {
            return false;
        }

        // If the point is in the plane then check if it is in the face (including its edges and 
        // vertices)
        if (location == 0) 
        {
            return F->IsPointIn(point, threshold);
        }
    }

    // Since the point is at the back of all faces then return true
    return true;
}

VF dcel::DCEL::LineSegments() const
{
    // Label all half edges as not visited
    SetHalfedgesAttribute(ATTRIB_VISITED, false);

    // Initialize the vertex coordinates and vertex indices for the line segments representing the
    // edges of the geometry
    VF vf(m_vertices.size(), m_halfedges.size() / 2);

    // Initialize the vertex index for the line segments
    size_t i = 0;

    // Traverse through the half edges of the geometry and define their geometric information
    for (auto it = m_halfedges.begin(); it != m_halfedges.end(); ++it)
    {
        // Get the pointer to the current half edge
        std::shared_ptr<Halfedge> halfedge = *it;

        // Get the visited value of the current half edge
        bool visited;
        assert(halfedge->Attributes().Get<bool>(ATTRIB_VISITED, visited));

        // If the current half edge has been visited then continue to the next half edge
        if (visited)
        {
            continue;
        }

        // Get the references to the coordinates of the end points from the current halfedge
        //const Eigen::Vector3d & v0 = halfedge->start->Coords();
        //const Eigen::Vector3d & v1 = halfedge->twin->start->Coords();

        // Insert the coordinates of the end points
        vf.addVertex(halfedge->start->Coords());
        vf.addVertex(halfedge->twin->start->Coords());
        //vf.V.emplace_back(v0);
        //vf.V.emplace_back(v1);

        // Define the indices for the current line segment
        vf.addFace(i, i + 1);
        //vf.F.emplace_back();
        //std::vector<size_t> & indices = vf.F[vf.F.size()];
        //indices.push_back(i);
        //indices.push_back(i + 1);

        // Update the vertex indices for the next line segment
        i += 2;

        // Label the current half edge and its twin as visited
        halfedge->Attributes().Set<bool>(ATTRIB_VISITED, true);
        halfedge->twin->Attributes().Set<bool>(ATTRIB_VISITED, true);
    }

    // Remove the visited dynamic attribute from all half edges
    RemoveHalfedgesAttribute(ATTRIB_VISITED);

    return vf;
}

void dcel::DCEL::LoadDcelFile(const std::string & filename)
{
    Clear();

    std::ifstream file(filename);

    // First line has the number of vertices, faces and half edges of the geometry
    size_t nV, nF, nH;
    file >> nV >> nF >> nH;
    assert(nV > 0 && nF > 0 && nH > 0);

    m_vertices.resize(nV);
    m_faces.resize(nF);
    m_halfedges.resize(nH);

    char type;
    size_t idx, sidx, pidx, nidx, tidx, fidx, hidx;

    // Load the information of the vertices
    for (size_t i = 0; i < nV; i += 1) 
    {
        double x, y, z;
        
        // Line format: v idx X Y Z hidx
        file >> type >> idx >> x >> y >> z >> hidx;
        assert(type == 'v');
        assert(idx >= 0 && idx < nV);
        assert(hidx >= 0 && hidx < nH);
        assert(m_vertices[idx] == nullptr);
        
        if (m_halfedges[hidx] == nullptr)
        {
            m_halfedges[hidx] = std::make_shared<Halfedge>();
        }

        m_vertices[idx] = std::make_shared<Vertex>(x, y, z);
        m_vertices[idx]->halfedge = m_halfedges[hidx];
    }

    // Load the information of the faces
    for (size_t i = 0; i < nF; i += 1) 
    {
        // Line format: f idx hidx
        file >> type >> idx >> hidx;
        assert(type == 'f');
        assert(idx >= 0 && idx < nF);
        assert(hidx >= 0 && hidx < nH);
        assert(m_faces[idx] == nullptr);

        if (m_halfedges[hidx] == nullptr)
        {
            m_halfedges[hidx] = std::make_shared<Halfedge>();
        }

        m_faces[idx] = std::make_shared<Face>();
        m_faces[idx]->halfedge = m_halfedges[hidx];
    }

    // Load the information of the half edges
    for (size_t i = 0; i < nH; i += 1) 
    {
        // Line format: h idx sidx pidx nidx tidx fidx
        file >> type >> idx >> sidx >> pidx >> nidx >> tidx >> fidx;
        assert(type == 'h');
        assert(idx >= 0 && idx < nH);
        assert(sidx >= 0 && sidx < nV);
        assert(pidx >= 0 && pidx < nH);
        assert(nidx >= 0 && nidx < nH);
        assert(tidx >= 0 && tidx < nH);
        assert(fidx >= 0 && fidx < nF);
        assert(m_vertices[sidx] && m_faces[fidx]);

        if (m_halfedges[idx] == nullptr)
        {
            m_halfedges[idx] = std::make_shared<Halfedge>();
        }

        if (m_halfedges[pidx] == nullptr)
        {
            m_halfedges[pidx] = std::make_shared<Halfedge>();
        }

        if (m_halfedges[nidx] == nullptr)
        {
            m_halfedges[nidx] = std::make_shared<Halfedge>();
        }

        if (m_halfedges[tidx] == nullptr)
        {
            m_halfedges[tidx] = std::make_shared<Halfedge>();
        }

        m_halfedges[idx]->start = m_vertices[sidx];
        m_halfedges[idx]->previous = m_halfedges[pidx];
        m_halfedges[idx]->next = m_halfedges[nidx];
        m_halfedges[idx]->twin = m_halfedges[tidx];
        m_halfedges[idx]->face = m_faces[fidx];
    }

    file.close();

    CheckConsistency();
}

void dcel::DCEL::NormalizeVertices(double radius, const Eigen::Vector3d & C)
{
    // Traverse through the vertices of the DCEL, normalize their coordinates and scale to the 
    // given radius
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        // Make a copy of the normalized coordinates of the current vertex. Then, translate it with
        // respect of the given center and multiply by the given radius
        Eigen::Vector3d P = ((*it)->Coords().normalized() + C) * radius;

        // Normalize the coordinates and multiply by the given radius
        //P.normalize();
        //P += C;
        //P *= radius;

        // Set the new coordinates of the current vertex
        (*it)->Coords(P);
    }
}

size_t dcel::DCEL::NumberOfInternalEdges() const
{
    // Initialize the count of internal edges
    size_t nInternal = 0;

    // Set the half edges of the geometry as not visited
    SetHalfedgesAttribute(ATTRIB_VISITED, false);

    // Traverse through the half edges of the geometry
    for (auto it = m_halfedges.begin(); it != m_halfedges.end(); ++it)
    {
        // Get the pointer to the current half edge
        std::shared_ptr<Halfedge> halfedge = *it;

        // Get the visited value of the half edge
        bool visited;
        assert(halfedge->Attributes().Get<bool>(ATTRIB_VISITED, visited));

        // If the half edge has been visited then continue to the next one
        if (visited)
        {
            continue;
        }

        // If both half edge and its twin have an incident face then update the count of internal 
        // edges
        if (halfedge->face && halfedge->twin->face)
        {
            nInternal += 1;
        }

        // Set both half edge and its twin as visited
        halfedge->Attributes().Set<bool>(ATTRIB_VISITED, true);
        halfedge->twin->Attributes().Set<bool>(ATTRIB_VISITED, true);
    }

    // Remove the visited dynamic attribute from the half edges
    RemoveHalfedgesAttribute(ATTRIB_VISITED);

    // Return the number of internal edges
    return nInternal;
}

size_t dcel::DCEL::NumberOfTriangles() const
{
    // Initialize the counter for the number of triangles required for representing the geometry
    size_t nTriangles = 0;

    // Traverse through the faces of the geometry and count the number of triangles required for 
    // each face in the geometry
    for (auto it = m_faces.begin(); it != m_faces.end(); ++it)
    {
        nTriangles += (*it)->CountTriangles();
    }

    // Return the number of triangles
    return nTriangles;
}

size_t dcel::DCEL::NumberOfVertices() const
{
    return m_vertices.size();
}

void dcel::DCEL::QuadrangulateFaces()
{
    // Split half edges into two half edges. This function sets the ATTRIB_ORIGINAL attribute on 
    // the original vertices of the DCEL
    SplitHalfedgesByMidpoint();

    // Get the number of original faces in the geometric domain
    size_t nOriginalFaces = m_faces.size();

    // Traverse through the original faces and subdivide them into quadrilaterals
    for (size_t i = 0; i < nOriginalFaces; i += 1)
    {
        SubdivideFaceIntoQuadrilaterals(m_faces[i]);
    }

    // Remove the ATTRIB_ORIGINAL and ATTRIB_VISITED attributes from the vertices in the DCEL
    RemoveVerticesAttribute(ATTRIB_ORIGINAL);
    RemoveVerticesAttribute(ATTRIB_VISITED);
}

void dcel::DCEL::RemoveFacesAttribute(const std::string & name) const
{
    for (auto it = m_faces.begin(); it != m_faces.end(); ++it)
    {
        (*it)->Attributes().Erase(name);
    }
}

void dcel::DCEL::RemoveHalfedgesAttribute(const std::string & name) const
{
    for (auto it = m_halfedges.begin(); it != m_halfedges.end(); ++it)
    {
        (*it)->Attributes().Erase(name);
    }
}

void dcel::DCEL::RemoveVerticesAttribute(const std::string & name) const
{
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        (*it)->Attributes().Erase(name);
    }
}

void dcel::DCEL::Rotate(const Eigen::Vector3d & K, const double angle)
{
    // Normalize the axis vector
    const Eigen::Vector3d nK = K.normalized();

    // Traverse through the vertices and rotate them
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        std::shared_ptr<Vertex> V = *it;
        Eigen::Vector3d R = utils::axisAngleRotation(V->Coords(), nK, angle);
        V->Coords(R);
    }
}

void dcel::DCEL::Scale(double factor, const Eigen::Vector3d & C)
{
    // Traverse through the vertices vector and scale the vertices
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        // Scale the coordinates with respect of the scale reference point and the scale factor
        const Eigen::Vector3d P = C + (((*it)->Coords() - C) * factor);

        // Set the new coordinates for the current vertex
        (*it)->Coords(P);
    }
}

void dcel::DCEL::Set(const VF & vf)
{
    // Clear the DCEL
    Clear();

    // 
    size_t nVertices = vf.countVertices(), i = 0, j = 0, nFaces = vf.countFaces(), nIndices = 0, 
        prevIndex = 0, currIndex = 0, nextIndex = 0, nextIdx = 0;
    
    // Resize the vertices vector using the number of vertices from the given geometry
    m_vertices.resize(nVertices, nullptr);

    // Traverse through the array of vertices and generate the DCEL vertices
    for (i = 0; i < nVertices; i += 1)
    {
        m_vertices[i] = std::make_shared<Vertex>(vf.Vertex(i));
    }

    // Initialize an object to be used as a map for the half edges to be found by their start and 
    // end vertices indexes
    std::unordered_map<std::string, std::shared_ptr<Halfedge>> vToH;

    // 
    m_faces.resize(nFaces, nullptr);

    m_halfedges.resize(vf.countEdges() * 2, nullptr);

    // Traverse through the vertex indices and generate the faces and half edges of the DCEL
    for (i = 0; i < nFaces; i += 1)
    {
        // Initialize the current DCEL face
        m_faces[i] = std::make_shared<Face>();

        // Get the reference to the current indices vector
        const std::vector<size_t> & indices = vf.face(i);

        // Get the number of vertices for the current face
        nIndices = indices.size();

        // Traverse through the vertex indices of the current face
        for (j = 0; j < nIndices; j += 1)
        {
            // Get the indexes for the previous, current and next vertices of the face
            prevIndex = (j == 0) ? indices[nIndices - 1] : indices[j - 1];
            currIndex = indices[j];
            nextIndex = (j == nIndices - 1) ? indices[0] : indices[j + 1];

            // Get the respective DCEL vertices
            std::shared_ptr<Vertex> prevVertex = m_vertices[prevIndex];
            std::shared_ptr<Vertex> currVertex = m_vertices[currIndex];
            std::shared_ptr<Vertex> nextVertex = m_vertices[nextIndex];

            // Declare the variables for the previous and next half edges
            std::shared_ptr<Halfedge> previous = nullptr;
            std::shared_ptr<Halfedge> next = nullptr;

            // Define the key for the previous half edge
            std::string prevKey = std::to_string(prevIndex) + '-' + std::to_string(currIndex);

            // If the previous half edge does not exist then generate it and its twin. Otherwise, 
            // get it directly
            if (vToH.find(prevKey) == vToH.end())
            {
                // Generate the previous half edge and set its start vertex, end vertex and incident
                // face
                previous = std::make_shared<Halfedge>();
                previous->start = prevVertex;
                previous->face = m_faces[i];

                // Generate the twin half edge and set its start vertex, end vertex and twin half 
                // edge
                std::shared_ptr<Halfedge> twin = std::make_shared<Halfedge>();
                twin->start = currVertex;
                twin->twin = previous;

                // Set the twin to the previous half edge
                previous->twin = twin;

                // If the previous vertex doesn't have an incident half edge then use the previous 
                // half edge
                if (!prevVertex->halfedge)
                {
                    prevVertex->halfedge = previous;
                }

                // If the current vertex doesn't have an incident half edge then use the twin half 
                // edge
                if (!currVertex->halfedge)
                {
                    currVertex->halfedge = twin;
                }

                // If the current face doesn't have an incident half edge then use the previous half
                // edge
                if (!m_faces[i]->halfedge)
                {
                    m_faces[i]->halfedge = previous;
                }

                // Define the key for the twin half edge
                std::string twinKey = std::to_string(currIndex) + '-' + std::to_string(prevIndex);

                // Insert the previous and twin half edges into the half edges map
                vToH.insert(std::make_pair(prevKey, previous));
                vToH.insert(std::make_pair(twinKey, twin));

                // Insert the previous and twin half edge into the half edges array
                m_halfedges[nextIdx++] = previous;
                m_halfedges[nextIdx++] = twin;
            }
            else
            {
                // Get the previous half edge
                previous = vToH.find(prevKey)->second;

                // If the previous half edge doesn't have an incident face then use the current one
                if (!previous->face)
                {
                    previous->face = m_faces[i];
                }

                // If the current face doesn't have an incideht half edge then use the previous 
                // half edge
                if (!m_faces[i]->halfedge)
                {
                    m_faces[i]->halfedge = previous;
                }
            }

            // Define the key for the next half edge
            std::string nextKey = std::to_string(currIndex) + '-' + std::to_string(nextIndex);

            // If the next half edge does not exist then generate it and its twin. Otherwise, get 
            // it directly
            if (vToH.find(nextKey) == vToH.end())
            {
                // Generate the next half edge and set its start vertex, end vertex and incident 
                // face
                next = std::make_shared<Halfedge>();
                next->start = currVertex;
                next->face = m_faces[i];

                // Generate the twin half edge and set its start vertex, end vertex and twin half 
                // edge
                std::shared_ptr<Halfedge> twin = std::make_shared<Halfedge>();
                twin->start = nextVertex;
                twin->twin = next;

                // Set the twin to the next half edge
                next->twin = twin;

                // If the current vertex doesn't have an incident half edge then use the next half 
                // edge
                if (!currVertex->halfedge)
                {
                    currVertex->halfedge = next;
                }

                // If the next vertex doesn't have an incident half edge then use the twin half 
                // edge
                if (!nextVertex->halfedge)
                {
                    nextVertex->halfedge = twin;
                }

                // If the current face doesn't have an incident half edge then use the next half 
                // edge
                if (!m_faces[i]->halfedge)
                {
                    m_faces[i]->halfedge = next;
                }

                // Define the key for the twin half edge
                std::string twinKey = std::to_string(nextIndex) + '-' + std::to_string(currIndex);

                // Insert the next and twin half edges into the half edges map
                vToH.insert(std::make_pair(nextKey, next));
                vToH.insert(std::make_pair(twinKey, twin));

                // Insert the next and twin half edge into the half edges array
                m_halfedges[nextIdx++] = next;
                m_halfedges[nextIdx++] = twin;
            }
            else
            {
                // Get the next half edge
                next = vToH.find(nextKey)->second;

                // If the next half edge doesn't have an incident face then use the current one
                if (!next->face)
                {
                    next->face = m_faces[i];
                }

                // If the current face doesn't have an incideht half edge then use the next half 
                // edge
                if (!m_faces[i]->halfedge)
                {
                    m_faces[i]->halfedge = next;
                }
            }

            // Link the previous and next half edges
            previous->next = next;
            next->previous = previous;
        }
    }

    // Close the remaining loops in the DCEL
    CloseLoops();

    // Erase and clear the map between the keys and the half edges
    vToH.erase(vToH.begin(), vToH.end());
    vToH.clear();
}

void dcel::DCEL::Set(const DCEL & dcel)
{
    Set(dcel.vf());
}

void dcel::DCEL::SetFacesIndex() const
{
    size_t nFaces = m_faces.size(), i = 0;

    for (i = 0; i < nFaces; i += 1)
    {
        m_faces[i]->Attributes().Set<size_t>(ATTRIB_INDEX, i);
    }
}

template <typename T>
void dcel::DCEL::SetFacesAttribute(const std::string & name, T value) const
{
    // Traverse through the faces and set the attribute
    for (auto it = m_faces.begin(); it != m_faces.end(); ++it)
    {
        (*it)->Attributes().Set<T>(name, value);
    }
}

void dcel::DCEL::SetHalfedgesIndex() const
{
    // Get the number of half edges in the geometry
    size_t nHalfedges = m_halfedges.size();

    // Traverse through the half edges and set the ATTRIB_INDEX dynamic attribute
    for (size_t i = 0; i < nHalfedges; i += 1)
    {
        m_halfedges[i]->Attributes().Set<size_t>(ATTRIB_INDEX, i);
    }
}

template <typename T>
void dcel::DCEL::SetHalfedgesAttribute(const std::string & name, T value) const
{
    // Traverse through the half edges and set the attribute
    for (auto it = m_halfedges.begin(); it != m_halfedges.end(); ++it)
    {
        (*it)->Attributes().Set<T>(name, value);
    }
}

template <typename T>
void dcel::DCEL::SetVerticesAttribute(const std::string & name, T value) const
{
    // Traverse through the vertices and set the attribute
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        (*it)->Attributes().Set<T>(name, value);
    }
}

void dcel::DCEL::SetVerticesIndex() const
{
    // Get the number of vertices in the geometry
    size_t nVertices = m_vertices.size();

    // Traverse through the vertices and set the ATTRIB_INDEX dynamic attribute
    for (size_t i = 0; i < nVertices; i += 1)
    {
        m_vertices[i]->Attributes().Set<size_t>(ATTRIB_INDEX, i);
    }
}

void dcel::DCEL::SetVerticesLocation(const toolkit::Plane & plane, double threshold) const
{
    // Traverse through the vertices of the geometry and set their location with respect of the 
    // plane
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        // Get the pointer to the current vertex
        std::shared_ptr<Vertex> V = *it;

        // Get the location of the point with respect of the plane and store it in the location 
        // dynamic attribute
        int location = plane.PointLocation(V->Coords(), threshold);
        V->Attributes().Set<int>(ATTRIB_LOCATION, location);
    }
}

Eigen::Vector3d dcel::DCEL::SimpleCentroid(bool fixZeros, double threshold) const
{
    // 
    Eigen::Vector3d C(0, 0, 0);

    // 
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        C += (*it)->Coords();
    }

    // 
    C /= (double)m_vertices.size();

    // 
    if (fixZeros) 
    {
        utils::fixZeros(C, threshold);
    }

    return C;
}

void dcel::DCEL::SplitHalfedgeByMidpoint(std::shared_ptr<Halfedge> halfedge)
{
    // Get the coordinates of the half edge's mid point
    Eigen::Vector3d coordinates = halfedge->Midpoint();

    // Generate the DCEL midpoint vertex, indicate it is not an original vertex
    std::shared_ptr<Vertex> midpoint = std::make_shared<Vertex>(coordinates);
    midpoint->Attributes().Set<bool>(ATTRIB_ORIGINAL, false);

    // Push the new vertex into the vertices array
    m_vertices.push_back(midpoint);

    // Keep the reference to the twin half edge (we need it for later)
    std::shared_ptr<Halfedge> twinHalfedge = halfedge->twin;

    // Generate a new DCEL half edge from the midpoint to the current half edge's end point
    std::shared_ptr<Halfedge> newHalfedge = std::make_shared<Halfedge>();
    newHalfedge->start = midpoint;
    newHalfedge->previous = halfedge;
    newHalfedge->next = halfedge->next;
    newHalfedge->face = halfedge->face;
    newHalfedge->twin = nullptr;
    newHalfedge->Attributes().Set<bool>(ATTRIB_SPLIT, true);

    // Make the new half edge incident to the new midpoint vertex
    midpoint->halfedge = newHalfedge;

    // Push the new half edge into the half edges array
    m_halfedges.push_back(newHalfedge);

    // Fix the previous reference of the current half edge's next
    halfedge->next->previous = newHalfedge;

    // Edit the current half edge to go from the start to the midpoint
    halfedge->next = newHalfedge;
    halfedge->twin = nullptr;
    halfedge->Attributes().Set<bool>(ATTRIB_SPLIT, true);

    // Generate a new DCEL half edge from midpoint to twin's end point
    std::shared_ptr<Halfedge> newTwinHalfedge = std::make_shared<Halfedge>();
    newTwinHalfedge->start = midpoint;
    newTwinHalfedge->previous = twinHalfedge;
    newTwinHalfedge->next = twinHalfedge->next;
    newTwinHalfedge->face = twinHalfedge->face;
    newTwinHalfedge->twin = halfedge;
    newTwinHalfedge->Attributes().Set<bool>(ATTRIB_SPLIT, true);

    // Push the new twin half edge into the half edges array
    m_halfedges.push_back(newTwinHalfedge);

    // Set the current half edge twin to the new twin half edge
    halfedge->twin = newTwinHalfedge;

    // Fix the previous reference of the twin's next
    twinHalfedge->next->previous = newTwinHalfedge;

    // Edit the current twin half edge to go from the start to the midpoint
    twinHalfedge->next = newTwinHalfedge;
    twinHalfedge->twin = newHalfedge;
    twinHalfedge->Attributes().Set<bool>(ATTRIB_SPLIT, true);

    // Set the new half edge's twin to the twin half edge
    newHalfedge->twin = twinHalfedge;
}

void dcel::DCEL::SplitHalfedgesByMidpoint()
{
    // Set all vertices in the DCEL as original
    SetVerticesAttribute<bool>(ATTRIB_ORIGINAL, true);

    // Set all half edges in the DCEL with the ATTRIB_SPLIT attribute with value false
    SetHalfedgesAttribute<bool>(ATTRIB_SPLIT, false);

    // Set all half edges in the DCEL as original
    SetHalfedgesAttribute<bool>(ATTRIB_ORIGINAL, true);

    // Get the number of original half edges in the geometric domain
    size_t nOriginalHalfedges = m_halfedges.size();

    // Traverse through the original half edges
    for (size_t i = 0; i < nOriginalHalfedges; i += 1)
    {
        // Get the pointer to the current half edge
        std::shared_ptr<Halfedge> halfedge = m_halfedges[i];

        // Check if the current half edge has been split
        bool isSplit;
        assert(halfedge->Attributes().Get<bool>(ATTRIB_SPLIT, isSplit));

        // If the current half edge has been split then continue to the next half edge
        if (isSplit)
        {
            continue;
        }

        // If the current half edge is not one of the originals then continue to the next half edge
        if (!halfedge->Attributes().Has(ATTRIB_ORIGINAL))
        {
            continue;
        }

        // Split the current half edge by its midpoint
        SplitHalfedgeByMidpoint(halfedge);
    }

    // Remove the split and original dynamic attributes from the half edges
    RemoveHalfedgesAttribute(ATTRIB_SPLIT);
    RemoveHalfedgesAttribute(ATTRIB_ORIGINAL);

    // Check the consistency of the geometry
    CheckConsistency();
}

void dcel::DCEL::SubdivideFaceByMidpoints(std::shared_ptr<Face> face)
{
	// 
	assert(face);

	// Get the original dynamic attribute of the start vertex from the incident half edge of the 
	// face
	bool isOriginal = false;
	assert(face->halfedge->start->Attributes().Get<bool>(ATTRIB_ORIGINAL, isOriginal));

	// Get the pointer to a half edge whose start vertex is original
	std::shared_ptr<Halfedge> currentHalfedge = (isOriginal) ?
		face->halfedge : face->halfedge->next;

	// Traverse through the half edges of the geometry. Define the new faces by triangulating every 
	// two non-original vertices
	do
	{
		// Get the pointer to the half edge for the next iteration. We need to get this information now
		// since the incidency of the current half edge will change
		std::shared_ptr<Halfedge> nextIterationHalfedge = currentHalfedge->next->next;

		// Initialize a new face. Then, insert it into the faces vector of the geometry
		std::shared_ptr<Face> newFace = std::make_shared<Face>();
        m_faces.push_back(newFace);

		// Initialize the new half edges. Then, insert them into the half edges vector of the geometry
		std::shared_ptr<Halfedge> newHalfedge = std::make_shared<Halfedge>();
		std::shared_ptr<Halfedge> newHalfedgeTwin = std::make_shared<Halfedge>();

		// Set the start vertices of the new half edges
		newHalfedge->start = currentHalfedge->twin->start;
		newHalfedgeTwin->start = currentHalfedge->previous->start;

		// Set both half edges as twins
		newHalfedge->twin = newHalfedgeTwin;
		newHalfedgeTwin->twin = newHalfedge;

		// Fix the incidency associated to the original face
		newHalfedgeTwin->next = currentHalfedge->next;
		currentHalfedge->next->previous = newHalfedgeTwin;
		newHalfedgeTwin->previous = currentHalfedge->previous->previous;
		currentHalfedge->previous->previous->next = newHalfedgeTwin;
		face->halfedge = newHalfedgeTwin;
		newHalfedgeTwin->face = face;

		// Fix the incidency associated to the new face
		newHalfedge->previous = currentHalfedge;
		currentHalfedge->next = newHalfedge;
		newHalfedge->next = currentHalfedge->previous;
		currentHalfedge->previous->previous = newHalfedge;
		currentHalfedge->face = newFace;
		currentHalfedge->previous->face = newFace;
		newHalfedge->face = newFace;
		newFace->halfedge = currentHalfedge;

		// Get the pointer to the half edge for the next iteration
		currentHalfedge = nextIterationHalfedge;

		// Check if the start vertex of the current half edge is original
		assert(currentHalfedge->start->Attributes().Get<bool>(ATTRIB_ORIGINAL, isOriginal));

	} while (isOriginal);
}

void dcel::DCEL::SubdivideFaceIntoQuadrilaterals(std::shared_ptr<Face> face)
{
    // Set the ATTRIB_VISITED attribute on all incident vertices to the face as false
    face->SetVerticesAttribute<bool>(ATTRIB_VISITED, false);

    // Get the reference to the initial half edge. It is required for the initial half edge to have 
    // a non-original vertex
    std::shared_ptr<Halfedge> currentHalfedge;
    currentHalfedge = !face->halfedge->start->Attributes().Has(ATTRIB_ORIGINAL) ?
        face->halfedge : face->halfedge->next;

    // Get the coordinates of the face's barycenter and generate the new vertex
    Eigen::Vector3d coordinates;
    assert(face->Attributes().Get<Eigen::Vector3d>(ATTRIB_CENTER, coordinates));

    // Generate a new DCELVertex with the coordinates of the center point of the face. Indicate the
    // vertex is not original and it is not visited
    std::shared_ptr<Vertex> center = std::make_shared<Vertex>(coordinates);
    center->Attributes().Set<bool>(ATTRIB_ORIGINAL, false);
    center->Attributes().Set<bool>(ATTRIB_VISITED, false);

    // Push the new vertex into the vertices array of the DCEL
    m_vertices.push_back(center);

    // NOTE: The first sub-face is always processed in a different way since it is the only time we 
    // will generate two half edges simultaneously. That's why it is not generated in the loop way 
    // below.

    // Generate a new DCELFace and push it into the faces array of the DCEL
    std::shared_ptr<Face> newFace = std::make_shared<Face>();
    m_faces.push_back(newFace);

    // Generate a new DCELHalfedge from the center to the current half edge's start vertex. 
    // Initialize the respective twin half edge as well
    std::shared_ptr<Halfedge> newHalfedge = std::make_shared<Halfedge>();
    std::shared_ptr<Halfedge> newHalfedgeTwin = std::make_shared<Halfedge>();

    // Generate a new DCELHalfedge from the current half edge's next's end vertex to the center. 
    // Initialize the respective twin half edge as well
    std::shared_ptr<Halfedge> newPreviousHalfedge = std::make_shared<Halfedge>();
    std::shared_ptr<Halfedge> newPreviousHalfedgeTwin = std::make_shared<Halfedge>();

    // Push the new half edges into the DCEL half edges array
    m_halfedges.push_back(newHalfedge);
    m_halfedges.push_back(newHalfedgeTwin);
    m_halfedges.push_back(newPreviousHalfedge);
    m_halfedges.push_back(newPreviousHalfedgeTwin);

    // Define the attribute values of the new half edges
    newHalfedge->start = center;
    newHalfedge->previous = newPreviousHalfedge;
    newHalfedge->next = currentHalfedge;
    newHalfedge->face = newFace;
    newHalfedge->twin = newHalfedgeTwin;

    // Define the attribute values of the new half edge's twin
    newHalfedgeTwin->start = currentHalfedge->start;
    newHalfedgeTwin->previous = currentHalfedge->previous;
    newHalfedgeTwin->next = newPreviousHalfedgeTwin;
    newHalfedgeTwin->face = face;
    newHalfedgeTwin->twin = newHalfedge;

    // Define the attribute values of the previous new half edge
    newPreviousHalfedge->start = currentHalfedge->next->twin->start;
    newPreviousHalfedge->previous = currentHalfedge->next;
    newPreviousHalfedge->next = newHalfedge;
    newPreviousHalfedge->face = newFace;
    newPreviousHalfedge->twin = newPreviousHalfedgeTwin;

    // Define the attribute values of the previous new half edge's twin
    newPreviousHalfedgeTwin->start = center;
    newPreviousHalfedgeTwin->previous = newHalfedgeTwin;
    newPreviousHalfedgeTwin->next = currentHalfedge->next->next;
    newPreviousHalfedgeTwin->face = face;
    newPreviousHalfedgeTwin->twin = newPreviousHalfedge;

    // Fix the half edges affected in the original face
    currentHalfedge->previous->next = newHalfedgeTwin;
    currentHalfedge->next->next->previous = newPreviousHalfedgeTwin;

    // Fix the half edges affected in the new face
    currentHalfedge->face = newFace;
    currentHalfedge->previous = newHalfedge;
    currentHalfedge->next->face = newFace;
    currentHalfedge->next->next = newPreviousHalfedge;

    // Make the initial new half edge as the incident one for the new face and the center
    newFace->halfedge = newHalfedge;
    center->halfedge = newHalfedge;

    // It is possible the incident half edge for the original face was the current half edge or its 
    // next one. Let's play safe and make it to be the new half edge's twin's previous
    face->halfedge = newHalfedgeTwin->previous;

    // Set the respective vertices with the ATTRIB_VISITED attribute as true
    newHalfedge->twin->start->Attributes().Set<bool>(ATTRIB_VISITED, true);
    center->Attributes().Set<bool>(ATTRIB_VISITED, true);
    newPreviousHalfedge->start->Attributes().Set<bool>(ATTRIB_VISITED, true);

    // Not actually visited, but we are marking all of the new face vertices as visited though
    newHalfedge->next->twin->start->Attributes().Set<bool>(ATTRIB_VISITED, true);

    // Now locate the current half edge to the next half edge in the original face
    currentHalfedge = newPreviousHalfedgeTwin->next;

    // Get the reference to the current half edge's previous previous's half edge. We need it for 
    // setting up correctly the twin half edge. Also, no matter which is the current half edge now, 
    // the previous previous's half edge will always be the same half edge
    std::shared_ptr<Halfedge> previousPreviousHalfedge = newHalfedgeTwin;

    // Check if the end vertex of the next half edge of the current half edge has been visited
    bool isVisited;
    assert(currentHalfedge->next->twin->start->Attributes().Get<bool>(ATTRIB_VISITED, isVisited));

    // Repeat while the end vertex of the next half edge of the current half edge has not been 
    // visited
    while (!isVisited)
    {
        // Generate a new DCELFace and push it into the faces array of the DCEL
        newFace = std::make_shared<Face>();
        m_faces.push_back(newFace);

        // Generate the new half edge and its twin. Push both into the DCEL half edges array
        newHalfedge = std::make_shared<Halfedge>();
        newHalfedgeTwin = std::make_shared<Halfedge>();
        m_halfedges.push_back(newHalfedge);
        m_halfedges.push_back(newHalfedgeTwin);

        // Define the attributes of the new half edge
        newHalfedge->start = currentHalfedge->next->twin->start;
        newHalfedge->previous = currentHalfedge->next;
        newHalfedge->next = currentHalfedge->previous;
        newHalfedge->face = newFace;
        newHalfedge->twin = newHalfedgeTwin;

        // Define the attributes of the new half edge's twin
        newHalfedgeTwin->start = center;
        newHalfedgeTwin->previous = previousPreviousHalfedge;
        newHalfedgeTwin->next = currentHalfedge->next->next;
        newHalfedgeTwin->face = face;
        newHalfedgeTwin->twin = newHalfedge;

        // Fix the half edges affected in the original face
        previousPreviousHalfedge->next = newHalfedgeTwin;
        currentHalfedge->next->next->previous = newHalfedgeTwin;

        // Fix the half edges affected in the new face
        currentHalfedge->face = newFace;
        currentHalfedge->previous->face = newFace;
        currentHalfedge->next->face = newFace;
        currentHalfedge->previous->previous = newHalfedge;
        currentHalfedge->next->next = newHalfedge;

        // Make the new half edge as the incident for the new face
        newFace->halfedge = newHalfedge;

        // Set the ATTRIB_VISITED attributes of the respective vertices as visited. The center and 
        // the current half edge's start are already visited
        currentHalfedge->twin->start->Attributes().Set<bool>(ATTRIB_VISITED, true);
        newHalfedge->start->Attributes().Set<bool>(ATTRIB_VISITED, true);

        // Move to the next respective half edge in the original face
        currentHalfedge = newHalfedgeTwin->next;

        // Check if the end vertex of the next half edge of the current half edge has been visited
        assert(currentHalfedge->next->twin->start->Attributes().Get<bool>(ATTRIB_VISITED, isVisited));
    }

    // Check the consistency of the geometry
    CheckConsistency();
}

void dcel::DCEL::SubdivideFacesByMidpoints()
{
	// Split the half edges of the geometry by their midpoint
	SplitHalfedgesByMidpoint();

	// Get the number of original faces from the geometry
	size_t nOriginalFaces = m_faces.size();

	// Traverse through the original faces of the geometry and subdivide them by their edge 
	// midpoints
	for (size_t faceIdx = 0; faceIdx < nOriginalFaces; faceIdx += 1)
	{
		SubdivideFaceByMidpoints(m_faces[faceIdx]);
	}

	// Remove the ATTRIB_ORIGINAL and ATTRIB_VISITED attributes from the vertices in the DCEL
	RemoveVerticesAttribute(ATTRIB_ORIGINAL);
	RemoveVerticesAttribute(ATTRIB_VISITED);
}

void dcel::DCEL::Translate(const Eigen::Vector3d & V)
{
    // 
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        std::shared_ptr<Vertex> vertex = *it;
        vertex->Coords(vertex->Coords() + V);
    }
}

void dcel::DCEL::TriangulateFacesByVertices()
{
	// Get the number of original faces in the geometry
	size_t nFaces = m_faces.size();

	// Traverse through the faces of the geometry and triangulate them
	for (size_t i = 0; i < nFaces; i += 1)
	{
		// Get the pointer to the current face
		std::shared_ptr<dcel::Face> face = m_faces[i];

		// Get the coordinates of the center associated to the face
		Eigen::Vector3d C;
		assert(face->Attributes().Get<Eigen::Vector3d>(ATTRIB_CENTER, C));

		// Triangulate the current face using the coordinates of the associated center
		TriangulateFaceByVertices(face, C);
	}
}

void dcel::DCEL::TriangulateFaceByVertices(std::shared_ptr<Face> face, const Eigen::Vector3d & C)
{
    // 
    assert(face);

    // Get the pointer to the incident half edge of the face
    std::shared_ptr<Halfedge> halfedge = face->halfedge;

    // Get the pointer to the previous half edge of the incident half edge of the face
    std::shared_ptr<Halfedge> lastHalfedge = halfedge->previous;

    // Generate a new vertex using the given coordinates. Then, add it to the vertices vector
    std::shared_ptr<Vertex> center = std::make_shared<Vertex>(C);
    m_vertices.push_back(center);

    // Generate a new face, this is the first triangle from the given face. Then, add it to the 
    // faces vector
    std::shared_ptr<Face> newFace = std::make_shared<Face>();
    m_faces.push_back(newFace);

    // Generate four new half edges. Then, add them to the halfedges vector
    std::shared_ptr<Halfedge> nextHalfedge = std::make_shared<Halfedge>();
    std::shared_ptr<Halfedge> prevHalfedge = std::make_shared<Halfedge>();
    std::shared_ptr<Halfedge> nextHalfedgeTwin = std::make_shared<Halfedge>();
    std::shared_ptr<Halfedge> prevHalfedgeTwin = std::make_shared<Halfedge>();
    m_halfedges.push_back(nextHalfedge);
    m_halfedges.push_back(prevHalfedge);
    m_halfedges.push_back(nextHalfedgeTwin);
    m_halfedges.push_back(prevHalfedgeTwin);

    // Define the start and end vertices of the new half edges
    nextHalfedge->start = halfedge->twin->start;
    prevHalfedge->start = center;
    nextHalfedgeTwin->start = center;
    prevHalfedgeTwin->start = halfedge->start;

    // Set the previous and next half edges
    lastHalfedge->next = prevHalfedgeTwin;
    prevHalfedgeTwin->next = nextHalfedgeTwin;
    nextHalfedgeTwin->next = halfedge->next;
    halfedge->next->previous = nextHalfedgeTwin;
    nextHalfedgeTwin->previous = prevHalfedgeTwin;
    prevHalfedgeTwin->previous = lastHalfedge;
    nextHalfedge->next = prevHalfedge;
    prevHalfedge->previous = nextHalfedge;
    halfedge->next = nextHalfedge;
    nextHalfedge->previous = halfedge;
    halfedge->previous = prevHalfedge;
    prevHalfedge->next = halfedge;

    // Set the twin half edges
    nextHalfedge->twin = nextHalfedgeTwin;
    nextHalfedgeTwin->twin = nextHalfedge;
    prevHalfedge->twin = prevHalfedgeTwin;
    prevHalfedgeTwin->twin = prevHalfedge;

    // Set the incident face of the half edges
    halfedge->face = newFace;
    nextHalfedge->face = newFace;
    prevHalfedge->face = newFace;
    prevHalfedgeTwin->face = face;
    nextHalfedgeTwin->face = face;

    // Set the incident half edge on the faces
    face->halfedge = lastHalfedge;
    newFace->halfedge = halfedge;

    // Set the incident half edge to the center vertex
    center->halfedge = prevHalfedge;

    // Move to the next half edge as it was originally in the given face (this was the next half 
    // edge of the current half edge)
    halfedge = nextHalfedgeTwin->next;

    // Traverse through the half edges of the original face
    do
    {
        // Generate a new face. Then, add it to the faces vector
        newFace = std::make_shared<Face>();
        m_faces.push_back(newFace);

        // Generate two new half edges. Then, add them to the halfedges vector
        nextHalfedge = std::make_shared<Halfedge>();
        nextHalfedgeTwin = std::make_shared<Halfedge>();
        m_halfedges.push_back(nextHalfedge);
        m_halfedges.push_back(nextHalfedgeTwin);

        // Define the start and end vertices of the new half edges
        nextHalfedge->start = halfedge->twin->start;
        nextHalfedgeTwin->start = center;

        // Set the previous and next half edges
        lastHalfedge->next->next = nextHalfedgeTwin;
        nextHalfedgeTwin->next = halfedge->next;
        halfedge->next->previous = nextHalfedgeTwin;
        nextHalfedgeTwin->previous = lastHalfedge->next;
        halfedge->next = nextHalfedge;
        nextHalfedge->previous = halfedge;
        nextHalfedge->next = halfedge->previous;
        halfedge->previous->previous = nextHalfedge;

        // Set the twin half edges
        nextHalfedge->twin = nextHalfedgeTwin;
        nextHalfedgeTwin->twin = nextHalfedge;

        // Set the incident face of the half edges
        halfedge->face = newFace;
        nextHalfedge->face = newFace;
        halfedge->previous->face = newFace;
        nextHalfedgeTwin->face = face;

        // Set the incident half edge of the new face
        newFace->halfedge = halfedge;

        // Move to the next half edge as it was originally in the given face (this was the next 
        // half edge of the current half edge)
        halfedge = nextHalfedgeTwin->next;

    } while (halfedge != lastHalfedge);

    // Check the consistency of the geometry
    CheckConsistency();
}

const std::vector<std::shared_ptr<dcel::Vertex>>& dcel::DCEL::Vertices() const
{
    return m_vertices;
}

VF dcel::DCEL::vf() const
{
    // Initialize the vertex coordinates and vertex indices
    VF vf(m_vertices.size(), m_faces.size());

    // Initialize the variable for keep track of the vertex indices
    size_t i = 0;

    // Traverse through the vertices of the DCEL and set their respective index
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        // Get the pointer to the current vertex
        std::shared_ptr<Vertex> vertex = *it;

        // Insert the coordinates of the current vertex into the vertex coordinates vector
        vf.addVertex(vertex->Coords());

        // Set the index of the vertex as its attribute
        vertex->Attributes().Set<size_t>(ATTRIB_INDEX, i++);

        // Update the vertex index
        //i += 1;
    }

    // 
    std::vector<size_t> indices;

    // Traverse through the faces of the DCEL
    for (auto it = m_faces.begin(); it != m_faces.end(); ++it)
    {
        // Get the pointer to the current face
        std::shared_ptr<Face> face = *it;

        // Get the number of vertices of the face
        size_t nVertices = face->CountVertices();

        // Reset the indices vector
        indices.clear();
        indices.resize(nVertices);

        // 
        size_t idx = 0;

        // Initialize the array for the vertex indices of the current face. Then, get its reference
        //vf.F.emplace_back();
        //std::vector<size_t> & vertices = vf.F[vf.F.size() - 1];

        // Get the incident half edge to the face
        std::shared_ptr<Halfedge> currentHalfedge = face->halfedge;

        // Traverse through the half edges of the face
        do
        {
            // Insert the index of the start vertex into the vertices array
            assert(currentHalfedge->start->Attributes().Get<size_t>(ATTRIB_INDEX, i));
            //indices.push_back(i);
            indices[idx++] = i;

            // Move to the next half edge in the face
            currentHalfedge = currentHalfedge->next;

        } while (currentHalfedge != face->halfedge);

        // 
        vf.addFace(indices);
    }

    // Traverse through the vertices of the DCEL and remove their respective index
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        (*it)->Attributes().Erase(ATTRIB_INDEX);
    }

    // Return the vertex coordinates and vertex indices
    return vf;
}

double dcel::DCEL::Volume() const
{
    // Initialize the volume of the polyhedron
    double volume = 0.0;

    // Traverse through the faces
    for (auto it = m_faces.begin(); it != m_faces.end(); ++it)
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
}

void dcel::DCEL::Write() const
{
    // Set the vertex index on the vertices of the geometry
    SetVerticesIndex();

    // Write the label indicating the start of the geometry content
    std::cout << "geometry " << m_vertices.size() << " " << m_faces.size() << std::endl;

    // Traverse through the vertices of the geometry and write their content
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        // Get the pointer to the current vertex
        std::shared_ptr<Vertex> vertex = *it;

        // Get the index of the vertex
        size_t index;
        assert(vertex->Attributes().Get<size_t>(ATTRIB_INDEX, index));

        // Get the reference to the vertex coordinates
        const Eigen::Vector3d & C = vertex->Coords();

        // Write the information of the vertex
        std::cout << "v " << index << " (" << C.x() << ", " << C.y() << ", " << C.z() << ")" << std::endl;
    }

    // Traverse through the faces of the geometry and write their content
    for (auto it = m_faces.begin(); it != m_faces.end(); ++it)
    {
        // Get the pointer to the current face
        std::shared_ptr<Face> face = *it;

        // Get the pointer to the incident half edge of the face
        std::shared_ptr<Halfedge> halfedge = face->halfedge;

        // Start the line with the vertex indices of the current face
        std::cout << "f";

        // Traverse through the half edges of the face
        do
        {
            // Get the index of the start vertex of the half edge
            size_t index;
            assert(halfedge->start->Attributes().Get<size_t>(ATTRIB_INDEX, index));

            // Write the index
            std::cout << " " << index;

            // Move to the next half edge
            halfedge = halfedge->next;

        } while (halfedge != face->halfedge);

        // Finish the line with the vertex indices of the current face
        std::cout << std::endl;
    }

    // Write the label indicating the end of the geometry content
    std::cout << "/geometry" << std::endl;

    // Remove the vertex index attribute from the vertices of the geometry
    RemoveVerticesAttribute(ATTRIB_INDEX);
}

void dcel::DCEL::Write(std::ofstream & file) const
{
    // Set the vertex index on the vertices of the geometry
    SetVerticesIndex();

    // Write the label indicating the start of the geometry content
    file << "geometry " << m_vertices.size() << " " << m_faces.size() << std::endl;

    // Traverse through the vertices of the geometry and write their content
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        // Get the pointer to the current vertex
        std::shared_ptr<Vertex> vertex = *it;

        // Get the index of the vertex
        size_t index;
        assert(vertex->Attributes().Get<size_t>(ATTRIB_INDEX, index));

        // Get the reference to the vertex coordinates
        const Eigen::Vector3d & C = vertex->Coords();

        // Write the information of the vertex
        file << "v " << index << " " << C.x() << " " << C.y() << " " << C.z() << std::endl;
    }

    // Traverse through the faces of the geometry and write their content
    for (auto it = m_faces.begin(); it != m_faces.end(); ++it)
    {
        // Get the pointer to the current face
        std::shared_ptr<Face> face = *it;

        // Get the pointer to the incident half edge of the face
        std::shared_ptr<Halfedge> halfedge = face->halfedge;

        // Start the line with the vertex indices of the current face
        file << "f";

        // Traverse through the half edges of the face
        do
        {
            // Get the index of the start vertex of the half edge
            size_t index;
            assert(halfedge->start->Attributes().Get<size_t>(ATTRIB_INDEX, index));

            // Write the index
            file << " " << index;

            // Move to the next half edge
            halfedge = halfedge->next;

        } while (halfedge != face->halfedge);

        // Finish the line with the vertex indices of the current face
        file << std::endl;
    }

    // Write the label indicating the end of the geometry content
    file << "/geometry" << std::endl;

    // Remove the vertex index attribute from the vertices of the geometry
    RemoveVerticesAttribute(ATTRIB_INDEX);

}

void dcel::DCEL::Write(std::stringstream & ss) const
{
    // Set the vertex index on the vertices of the geometry
    SetVerticesIndex();

    // Write the label indicating the start of the geometry content
    ss << "geometry " << m_vertices.size() << " " << m_faces.size() << std::endl;

    // Traverse through the vertices of the geometry and write their content
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        // Get the pointer to the current vertex
        std::shared_ptr<Vertex> vertex = *it;

        // Get the index of the vertex
        size_t index;
        assert(vertex->Attributes().Get<size_t>(ATTRIB_INDEX, index));

        // Get the reference to the vertex coordinates
        const Eigen::Vector3d & C = vertex->Coords();

        // Write the information of the vertex
        ss << "v " << index << " " << C.x() << " " << C.y() << " " << C.z() << std::endl;
    }

    // Traverse through the faces of the geometry and write their content
    for (auto it = m_faces.begin(); it != m_faces.end(); ++it)
    {
        // Get the pointer to the current face
        std::shared_ptr<Face> face = *it;

        // Get the pointer to the incident half edge of the face
        std::shared_ptr<Halfedge> halfedge = face->halfedge;

        // Start the line with the vertex indices of the current face
        ss << "f";

        // Traverse through the half edges of the face
        do
        {
            // Get the index of the start vertex of the half edge
            size_t index;
            assert(halfedge->start->Attributes().Get<size_t>(ATTRIB_INDEX, index));

            // Write the index
            ss << " " << index;

            // Move to the next half edge
            halfedge = halfedge->next;

        } while (halfedge != face->halfedge);

        // Finish the line with the vertex indices of the current face
        ss << std::endl;
    }

    // Write the label indicating the end of the geometry content
    ss << "/geometry" << std::endl;

    // Remove the vertex index attribute from the vertices of the geometry
    RemoveVerticesAttribute(ATTRIB_INDEX);
}

void dcel::DCEL::WriteDcelFile(const std::string & filename) const
{
    SetVerticesIndex();
    SetFacesIndex();
    SetHalfedgesIndex();

    std::ofstream file(filename, std::ios::out);

    // First line is the number of vertices, faces and half edges
    file << m_vertices.size() << ' ' << m_faces.size() << ' ' << m_halfedges.size() << std::endl;

    // Write the information of the vertices
    for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
    {
        std::shared_ptr<Vertex> V = *it;
        const Eigen::Vector3d & C = V->Coords();

        // The index of the vertex
        size_t idx;
        assert(V->Attributes().Get(ATTRIB_INDEX, idx));

        // The index of the incident half edge
        size_t hidx;
        assert(V->halfedge->Attributes().Get(ATTRIB_INDEX, hidx));

        // Line format: v idx X Y Z hidx
        file << 'v ' << idx << ' ' << C.x() << ' ' << C.y() << ' ' << C.z() << ' ' << hidx << std::endl;
    }

    // Write the information of the faces
    for (auto it = m_faces.begin(); it != m_faces.end(); ++it)
    {
        std::shared_ptr<Face> F = *it;

        // The index of the face
        size_t idx;
        assert(F->Attributes().Get(ATTRIB_INDEX, idx));

        // The index of the incident half edge
        size_t hidx;
        assert(F->halfedge->Attributes().Get(ATTRIB_INDEX, hidx));

        // Line format: f idx hidx
        file << 'f ' << idx << ' ' << hidx << std::endl;
    }

    // Write the information of the half edges
    for (auto it = m_halfedges.begin(); it != m_halfedges.end(); ++it)
    {
        std::shared_ptr<Halfedge> H = *it;

        // The index of the half edge
        size_t idx;
        assert(H->Attributes().Get(ATTRIB_INDEX, idx));

        // The indices of the incident elements
        size_t sidx, pidx, nidx, tidx, fidx;
        assert(H->start->Attributes().Get(ATTRIB_INDEX, sidx));
        assert(H->previous->Attributes().Get(ATTRIB_INDEX, pidx));
        assert(H->next->Attributes().Get(ATTRIB_INDEX, nidx));
        assert(H->twin->Attributes().Get(ATTRIB_INDEX, tidx));
        assert(H->face->Attributes().Get(ATTRIB_INDEX, fidx));

        // Line format: h idx sidx pidx nidx tidx fidx
        file << 'h ' << idx << ' ' << sidx << ' ' << pidx << ' ' << nidx << ' ' << tidx << ' ' << fidx << std::endl;
    }
    
    file.close();

    RemoveVerticesAttribute(ATTRIB_INDEX);
    RemoveFacesAttribute(ATTRIB_INDEX);
    RemoveHalfedgesAttribute(ATTRIB_INDEX);
}

void dcel::DCEL::WriteGeogebraJs(
    std::stringstream & ss, 
    const std::string & prefix, 
    int r, 
    int g, 
    int b, 
    int a, 
    double threshold) const
{
    // Get the number of vertices of the geometry
    size_t nVertices = m_vertices.size();

    // Traverse through the vertices of the geometry
    for (size_t i = 0; i < nVertices; i += 1) 
    {
        // Set the index of the vertex
        m_vertices[i]->Attributes().Set<size_t>(ATTRIB_INDEX, i);

        // Write the commands that generate the current vertex in GeoGebra
        m_vertices[i]->WriteGeogebraJs(ss, prefix, r, g, b, threshold);
    }

    // Get the number of faces of the geometry
    size_t nFaces = m_faces.size();

    // Traverse through the faces and write their content
    for (size_t i = 0; i < nFaces; i += 1) 
    {
        // Set the index of the face
        m_faces[i]->Attributes().Set<size_t>(ATTRIB_INDEX, i);

        // If the current face is planar then write its indices using a single polygon. Otherwise, 
        // write its indices as triangles
        if (m_faces[i]->IsPlanar(threshold))
        {
            m_faces[i]->WritePlanarGeogebraJs(ss, prefix);
        }
        else 
        {
            m_faces[i]->WriteTriangularGeogebraJs(ss, prefix);
        }
    }

    // remove the index attributes
    RemoveVerticesAttribute(ATTRIB_INDEX);
    RemoveFacesAttribute(ATTRIB_INDEX);
}
