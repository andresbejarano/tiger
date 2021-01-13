#include "assembly.h"
#include "utils.h"
#include <iostream>
#include <fstream>

Assembly::Assembly() 
{
}

Assembly::Assembly(const std::vector<VF> & blocks)
{
    size_t nblocks = blocks.size(), i = 0;

    m_blocks.resize(nblocks, nullptr);

    for (i = 0; i < nblocks; i += 1)
    {
        m_blocks[i] = std::make_shared<Block>(blocks[i], i);
    }
}

Assembly::~Assembly() 
{
    Clear();
}

void Assembly::AdaptiveTileOffsetClipping(
    const std::shared_ptr<Tessellation> tessellation, 
    const std::string & topFunction, 
    const std::string & bottomFunction, 
    std::vector<VF> & clipped, 
    double threshold) const
{
    assert(tessellation);

    size_t nblocks = m_blocks.size(), i = 0, nLengths = 0;

    clipped.clear();
    clipped.resize(nblocks);

    std::vector<double> lengths;

    double topLength = 0, bottomLength = 0;

    // Traverse through the blocks and find the requested lenghts from the neighborhood of the 
    // respective tile associated to each block
    for (i = 0; i < nblocks; i += 1)
    {
        assert(m_blocks[i]->FaceIndex() == i);

        // Get the top central lenghts from the neighborhood
        nLengths = NeighborhoodCentralLengths(tessellation, i, lengths, 1);
        topLength = utils::vectorFunction(lengths, topFunction);

        // Get the bottom central lengths from the neighborhood
        nLengths = NeighborhoodCentralLengths(tessellation, i, lengths, -1);
        bottomLength = utils::vectorFunction(lengths, bottomFunction);

        // Use the top and bottom central lenghts and clip the block with the plane offset method
        clipped[i] = m_blocks[i]->PlaneOffsetClipping(
            tessellation->Geometry().Plane(i), topLength, bottomLength, threshold);
    }
}

bool Assembly::AreBlocksModified() const
{
    std::shared_ptr<Block> block = nullptr;

    bool isModified = false;

    // Traverse through the blocks and check if they have the ATTRIB_MODIFIED dynamic attribute
    for (auto it = m_blocks.begin(); it != m_blocks.end(); ++it)
    {
        // Get the pointer to the current block
        block = *it;

        // If the block doesn't have the ATTRIB_MODIFIED dynamic attribute then continue to the 
        // next block
        if (!block->Attributes().Has(ATTRIB_MODIFIED))
        {
            continue;
        }

        // Get the value of the ATTRIB_MODIFIED dynamic attribute
        assert(block->Attributes().Get<bool>(ATTRIB_MODIFIED, isModified));

        // If the block is modified then return true
        if (isModified)
        {
            return true;
        }
    }

    // Since no block has been modified then return false
    return false;
}

size_t Assembly::CalculateInterfacePolygons(
    std::shared_ptr<Tessellation> tessellation,
    std::vector<std::shared_ptr<VF>> & interfaces,
    double threshold) const
{
    // There must be a tessellation and a non-empty assembly of blocks
    assert(tessellation);
    assert(tessellation->IsDCELLoaded());
    assert(m_blocks.size() > 0);

    //tessellation->LoadDCEL();

    // Get the reference to the geometry of the tessellation
    std::shared_ptr<dcel::DCEL> geometry = tessellation->DCEL();

    geometry->SetFacesIndex();

    // Get the number of internal edges of the tessellation. This is the number of interface 
    // polygons between the blocks of the assembly
    size_t nInterfaces = geometry->NumberOfInternalEdges(), nExistingInterfacePolygons = 0, idx = 0, p1Index = 0, p2Index = 0;

    // Clear the interface polygons vector and resize it
    interfaces.clear();
    interfaces.resize(nInterfaces, nullptr);

    // Set all half edges of the tessellation as not visited
    geometry->SetHalfedgesAttribute<bool>(ATTRIB_VISITED, false);

    std::shared_ptr<dcel::Halfedge> halfedge = nullptr;

    bool visited = false;

    int dir = 0;

    std::shared_ptr<Block> P1 = nullptr, P2 = nullptr;

    // Traverse through the half edges of the tessellation
    for (auto it = geometry->Halfedges().begin(); it != geometry->Halfedges().end(); ++it)
    {
        P1 = nullptr;
        P2 = nullptr;

        // Get the pointer to the current half edge
        halfedge = *it;

        // Check whether the current half edge has been visited
        assert(halfedge->Attributes().Get<bool>(ATTRIB_VISITED, visited));

        // If the current half edge doesn't have a face, or the twin half edge doesn't have a face,
        // or the current half edge has been visited then continue to the next half edge
        if (!halfedge->face || !halfedge->twin->face || visited)
        {
            continue;
        }

        // Get the direction value of the current half edge
        assert(halfedge->Attributes().Get<int>(ATTRIB_DIRECTION_VALUE, dir));

        // Get the indices to the adjacent blocks. block 1 will be the block located at the 
        // negative direction, block 2 will be the block at the positive direction
        if (dir < 0)
        {
            assert(halfedge->face->Attributes().Get<size_t>(ATTRIB_INDEX, p1Index));
            assert(halfedge->twin->face->Attributes().Get<size_t>(ATTRIB_INDEX, p2Index));
        }
        else
        {
            assert(halfedge->twin->face->Attributes().Get<size_t>(ATTRIB_INDEX, p1Index));
            assert(halfedge->face->Attributes().Get<size_t>(ATTRIB_INDEX, p2Index));
        }

        // Get the pointers to the respective incident blocks
        P1 = m_blocks[p1Index];
        P2 = m_blocks[p2Index];

        std::shared_ptr<VF> currentInterface = std::make_shared<VF>(P1->InterfacePolygon(P2, threshold));

        if (currentInterface->countVertices() == 0 || currentInterface->countFaces() == 0)
        {
            // The blocks don't have a common interface. It could happen when the resultant blocks
            // are not aligned. In such case we disable the respective blocks
            P1->Disable();
            P2->Disable();

            currentInterface->clear();
            currentInterface = nullptr;
        }
        else 
        {
            interfaces[idx] = currentInterface;

            interfaces[idx]->Attributes().Set<size_t>(ATTRIB_POSITIVE_BLOCK_INDEX, p1Index);
            interfaces[idx]->Attributes().Set<size_t>(ATTRIB_NEGATIVE_BLOCK_INDEX, p2Index);

            // Update the number of existing interface polygons
            nExistingInterfacePolygons += 1;

            currentInterface = nullptr;
        }

        //std::cout << "Finished block " << idx << std::endl;

        // Update the index for the next interface polygon
        idx += 1;

        // Set the half edge and its twin as visited
        halfedge->Attributes().Set<bool>(ATTRIB_VISITED, true);
        halfedge->twin->Attributes().Set<bool>(ATTRIB_VISITED, true);
    }

    // Remove the visited dynamic attribute from the half edges of the tessellation
    geometry->RemoveFacesAttribute(ATTRIB_INDEX);
    geometry->RemoveHalfedgesAttribute(ATTRIB_VISITED);

    // 
    return nExistingInterfacePolygons;
}

bool Assembly::CalculateOverlapping(
    std::shared_ptr<Tessellation> tessellation, 
    std::vector<std::shared_ptr<VF>> & overlapping, 
    double threshold) const
{
    //TODO
    return false;
}

void Assembly::Centroids(std::vector<Eigen::Vector3d> & C, bool fixZeros, double threshold) const
{
    size_t nblocks = m_blocks.size(), i;

    C.clear();
    C.resize(nblocks);

    for (i = 0; i < nblocks; i += 1)
    {
        C[i] << m_blocks[i]->Geometry().centroid(fixZeros, threshold);
    }
}

void Assembly::Clear() 
{
    size_t nblocks = m_blocks.size(), i = 0;

    for (i = 0; i < nblocks; i += 1) 
    {
        m_blocks[i]->Clear();
        m_blocks[i] = nullptr;
    }

    m_blocks.clear();
}

size_t Assembly::CountBlocks() const
{
    return m_blocks.size();
}

size_t Assembly::CountEdges() const
{
    size_t nEdges = 0;

    for (auto it = m_blocks.begin(); it != m_blocks.end(); ++it) 
    {
        nEdges += (*it)->Geometry().countEdges();
    }

    return nEdges;
}

size_t Assembly::CountTriangles() const
{
    size_t nTriangles = 0;

    for (auto it = m_blocks.begin(); it != m_blocks.end(); ++it)
    {
        nTriangles += (*it)->Geometry().countTriangles();
    }

    return nTriangles;
}

size_t Assembly::CountVertices() const
{
    size_t nVertices = 0;

    for (auto it = m_blocks.begin(); it != m_blocks.end(); ++it)
    {
        nVertices += (*it)->Geometry().countVertices();
    }

    return nVertices;
}

void Assembly::DisableIntersectedBlocks(const toolkit::Plane & plane)
{
    size_t nBlocks = m_blocks.size();

    for (size_t i = 0; i < nBlocks; i += 1) 
    {
        m_blocks[i]->DisableIfIntersects(plane);
    }
}

void Assembly::EnableAll() 
{
    for (auto it = m_blocks.begin(); it != m_blocks.end(); ++it) 
    {
        (*it)->Enable();
    }
}

void Assembly::TileOffsetClipping(
    const std::shared_ptr<Tessellation> domain, 
    double extrados, 
    double intrados, 
    std::vector<VF> & clipped, 
    double threshold) const
{
    assert(domain);

    size_t nblocks = m_blocks.size(), i = 0;

    clipped.clear();
    clipped.resize(nblocks);

    for (i = 0; i < nblocks; i += 1)
    {
        clipped[i] = m_blocks[i]->PlaneOffsetClipping(
            domain->Geometry().Plane(i), extrados, intrados, threshold);
    }
}

const std::vector<std::shared_ptr<Block>> & Assembly::GetBlocks() const
{
    return m_blocks;
}

bool Assembly::IsEmpty() const 
{
    return m_blocks.empty();
}

/*void Assembly::LinkBlocksToTiles(std::shared_ptr<Tessellation> tessellation)
{
    // There must exist a tessellation
    assert(tessellation);
    assert(tessellation->IsDCELLoaded());

    //tessellation->LoadDCEL();

    std::shared_ptr<dcel::DCEL> geometry = tessellation->DCEL();

    // Set the index of the faces. Indices are stored in the ATTRIB_INDEX dynamic attribute
    //domain->GetGeometry().SetFacesIndex();
    geometry->SetFacesIndex();

    // Get the reference to the vector with the pointers to the faces of the tessellation
    //std::vector<std::shared_ptr<dcel::Face>> faces = domain->GetGeometry().faces;
    //const std::vector<std::shared_ptr<dcel::Face>> & faces = geometry->Faces();

    std::shared_ptr<dcel::Face> face = nullptr;

    size_t faceIndex = 0, blockIndex = 0;

    // Traverse through the faces of the tessellation
    for (auto it = geometry->Faces().begin(); it != geometry->Faces().end(); ++it)
    {
        // Get the pointer to the current face
        face = *it;

        // Get the index of the current face
        assert(face->Attributes().Get<size_t>(ATTRIB_INDEX, faceIndex));

        // Get the index of the block associated to the current face
        face->Attributes().Get<size_t>(ATTRIB_BLOCK_INDEX, blockIndex);

        // Set the index of the associated face to the current block
        m_blocks[blockIndex]->Attributes().Set<size_t>(ATTRIB_FACE_INDEX, faceIndex);
    }
}*/

std::vector<std::vector<double>> Assembly::Loads(double density) const
{
    size_t nblocks = m_blocks.size(), i = 0;

    std::vector<std::vector<double>> loads(nblocks);

    for (i = 0; i < nblocks; i += 1)
    {
        loads[i] = m_blocks[i]->Loads(density);
    }

    return loads;
}

size_t Assembly::NeighborhoodCentralLengths(
    const std::shared_ptr<Tessellation> tessellation,
    size_t index, 
    std::vector<double> & lengths, 
    int direction, 
    double threshold) const
{
    assert(index >= 0 && index < m_blocks.size());

    // Set the index attribute at the faces of the geometry
    tessellation->DCEL()->SetFacesIndex();

    assert(m_blocks[index]->FaceIndex() == index);

    // Get the pointer to the face associated to the given block index
    std::shared_ptr<dcel::Face> face = tessellation->DCEL()->Faces()[index];

    lengths.clear();
    lengths.resize(face->CountNeighbors());

    std::shared_ptr<dcel::Halfedge> halfedge = face->halfedge;

    size_t idx = 0, blockIndex = 0;

    do 
    {
        if (halfedge->twin->face)
        {
            assert(halfedge->twin->face->Attributes().Get<size_t>(ATTRIB_INDEX, blockIndex));

            assert(m_blocks[blockIndex]->FaceIndex() == blockIndex);

            lengths[idx++] = m_blocks[blockIndex]->CentralLength(tessellation, direction, threshold);
        }

        halfedge = halfedge->next;

    } while (halfedge != face->halfedge);

    tessellation->DCEL()->RemoveFacesAttribute(ATTRIB_INDEX);

    return lengths.size();
}

void Assembly::RemoveBlocksAttribute(const std::string & name) 
{
    for (auto it = m_blocks.begin(); it != m_blocks.end(); ++it) 
    {
        (*it)->Attributes().Erase(name);
    }
}

/*void Assembly::Set(std::vector<std::shared_ptr<block>> & blocks)
{
    // Clear the blocks and resize the vector
    Clear();

    // Get the number of given blocks
    size_t nblocks = blocks.size();

    // Resize the vector for storing the pointers of the blocks
    m_blocks.resize(nblocks);

    // Traverse through the vector with the given blocks and copy the pointers
    for (size_t i = 0; i < nblocks; i += 1)
    {
        m_blocks[i] = blocks[i];
    }
}*/

/*void Assembly::Set(const std::vector<VF> & vfs) 
{
    Clear();

    size_t nblocks = vfs.size(), i = 0;

    m_blocks.resize(nblocks, nullptr);

    for (i = 0; i < nblocks; i += 1) 
    {
        m_blocks[i] = std::make_shared<Block>(vfs[i]);
    }
}*/

/*void Assembly::SetBlockIndex() 
{
    size_t nblocks = m_blocks.size(), i = 0;

    for (i = 0; i < nblocks; i += 1)
    {
        m_blocks[i]->Attributes().Set<size_t>(ATTRIB_INDEX, i);
    }
}*/

template <typename T>
void Assembly::SetBlocksAttribute(const std::string & name, T value)
{
    for (auto it = m_blocks.begin(); it != m_blocks.end(); ++it)
    {
        (*it).Attributes().Set<T>(name, value);
    }
}

void Assembly::ToggleBlocksByBoundary(std::shared_ptr<Tessellation> tessellation)
{
    assert(tessellation);
    assert(tessellation->IsDCELLoaded());

    std::shared_ptr<Block> block = nullptr;

    std::shared_ptr<dcel::Face> face = nullptr;

    // Traverse through the blocks of the assembly
    for (auto it = m_blocks.begin(); it != m_blocks.end(); ++it)
    {
        // Get the pointer to the current block
        block = *it;

        // Get the pointer to the associated face in the tessellation
        face = tessellation->DCEL()->Faces()[block->FaceIndex()];

        // Set the block as enabled or disabled. A block is enabled if its associated face in the 
        block->SetEnabled(face->HasAllNeighbors());
    }
}

void Assembly::Write() const
{
    size_t nblocks = m_blocks.size(), i = 0;

    std::cout << "assembly " << m_blocks.size() << std::endl;

    for (i = 0; i < nblocks; i += 1)
    {
        std::cout << i << std::endl;
        m_blocks[i]->Write();
    }

    std::cout << "/assembly" << std::endl;
}

void Assembly::Write(std::stringstream & ss) const 
{
    size_t nblocks = m_blocks.size(), i = 0;

    ss << "assembly " << m_blocks.size() << std::endl;

    for (i = 0; i < nblocks; i += 1) 
    {
        ss << i << std::endl;
        m_blocks[i]->Write(ss);
    }

    ss << "/assembly" << std::endl;
}

void Assembly::WriteGeogebraJs(
    std::stringstream && ss, 
    const std::string & prefix, 
    int r, 
    int g, 
    int b, 
    int a, 
    double threshold) const
{
    // TODO (WHAT?!!! I forgot)

    std::stringstream blockPrefix;

    size_t nblocks = m_blocks.size(), i = 0;

    for (i = 0; i < nblocks; i += 1) 
    {
        blockPrefix.str(std::string());
        blockPrefix.clear();
    }
}

void Assembly::WriteObjFile(const std::string & filepath, double threshold) const
{
    std::stringstream verticesSS, facesSS;

    std::shared_ptr<Block> block = nullptr;

    size_t nVertices = 0, nFaces = 0, idx = 0, nWrittenVertices = 0;

    Eigen::Vector3d P = Eigen::Vector3d::Zero();
    
    // Traverse through the blocks and write their geometry in the streamstring objects
    for (auto it = m_blocks.begin(); it != m_blocks.end(); ++it) 
    {
        block = *it;

        if (!block) 
        {
            continue;
        }

        nVertices = block->Geometry().countVertices();
        
        // Traverse through the vertices of the block and write them
        for (idx = 0; idx < nVertices; idx += 1) 
        {
            P << block->Geometry().Vertex(idx);

            verticesSS << "v " << P.x() << " " << P.y() << " " << P.z() << std::endl;
        }

        nFaces = block->Geometry().countFaces();

        // Traverse through the faces of the block and write them
        for (idx = 0; idx < nFaces; idx += 1) 
        {
            facesSS << "f";

            const std::vector<size_t> & indices = block->Geometry().face(idx);

            for (auto fIt = indices.begin(); fIt != indices.end(); ++fIt) 
            {
                facesSS << " " << (*fIt + nWrittenVertices + 1);
            }

            facesSS << std::endl;
        }

        nWrittenVertices += nVertices;
    }

    // Write the file
    std::ofstream file(filepath, std::ios::out);

    file << verticesSS.str() << facesSS.str();

    file.close();
}
