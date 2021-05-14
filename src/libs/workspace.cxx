#include <tiger/workspace/workspace.h>
#include <tiger/workspace/message.h>
#include <tiger/utils.h>
#include <tiger/geometries.h>
#include <tiger/toolkit/algorithms.h>

#include <queue>
#include <Eigen/Geometry>
#include <iostream>
#include <fstream>

Workspace::Workspace() : 
    m_assembly(nullptr), 
    m_directory(""), 
    m_equilibriumResults(nullptr), 
    m_interfaces(nullptr), 
    m_tessellation(nullptr)
{
}

Workspace::Workspace(const std::string & directory) :
    m_assembly(nullptr),
    m_directory(directory), 
    m_equilibriumResults(nullptr), 
    m_interfaces(nullptr), 
    m_tessellation(nullptr)
{
}

Workspace::~Workspace()
{
    // Clear the content of the workspace
    Clear();
}

bool Workspace::AdaptiveTileOffsetClipping(
    const std::string & topFunction, 
    const std::string & bottomFunction, 
    std::vector<VF> & blocks, 
    double threshold) const
{
    // There must exist a tessellation and blocks in the assembly
    assert(m_tessellation);
    assert(m_assembly);
    assert(m_assembly->CountBlocks() > 0);

    // 
    m_assembly->AdaptiveTileOffsetClipping(m_tessellation, topFunction, bottomFunction, blocks, threshold);

    // 
    return blocks.size() > 0;
}

bool Workspace::CheckAbaqusCodeRequirements(Message& msg) const 
{
    // If no blocks then report it
    if (!m_assembly || m_assembly->CountBlocks() == 0)
    {
        msg.set(Message::MSG_NO_BLOCKS);
        return false;
    }

    // TODO: Check physic attributes

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckAssemblyCodeRequirements(Message& msg) const 
{
    // If no blocks then report it
    if (!m_assembly || m_assembly->CountBlocks() == 0)
    {
        msg.set(Message::MSG_NO_BLOCKS);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckCentroidToOriginRequirements(Message& msg) const
{
    // If no tessellation then report it
    if (!m_tessellation)
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckClippingRequirements(Message& msg) const
{
    // If no tessellation then report it
    if (!m_tessellation)
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    // If no blocks then report it
    if (!m_assembly || m_assembly->CountBlocks() == 0)
    {
        msg.set(Message::MSG_NO_BLOCKS);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckDualTessellationRequirements(Message& msg) const
{
    // If no tessellation then report it
    if (!m_tessellation)
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckEquilibriumAnalysisRequirements(Message& msg) const 
{
    // If no blocks then report it
    if (!m_assembly || m_assembly->CountBlocks() == 0)
    {
        msg.set(Message::MSG_NO_BLOCKS);
        return false;
    }

    // If no interface polygons then report it
    if (!m_interfaces || m_interfaces->CountInterfaces() == 0)
    {
        msg.set(Message::MSG_NO_INTERFACES);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckFlipTessellationRequirements(Message& msg) const
{
    // If no tessellation then report it
    if (!m_tessellation)
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckFixBottomBlocksRequirements(Message& msg) const
{
    // If no tessellation then report it
    if (!m_tessellation)
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    // If no blocks then report it
    if (!m_assembly || m_assembly->CountBlocks() == 0)
    {
        msg.set(Message::MSG_NO_BLOCKS);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckHeightBisectionRequirements(Message& msg) const 
{
    // If no tessellation then report it
    if (!m_tessellation)
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    // If not all tiles have an even number of sides then report it
    if (!m_tessellation->CheckTilesAreEvenSided())
    {
        msg.set(Message::MSG_NO_EVEN_SIDED_TILES);
        return false;
    }

    // If not all tiles have the ATTRIB_CENTER attribute then report ot
    if (!m_tessellation->DCEL() && !m_tessellation->CheckTileCenters())
    {
        msg.set(Message::MSG_NO_TILE_CENTERS);
        return false;
    }

    // If not all half edges have both ATTRIB_DIRECTION_VALUE and ATTRIB_DIRECTION_VECTOR 
    // attributes then report it
    if (!m_tessellation->CheckHalfedgeDirections())
    {
        msg.set(Message::MSG_NO_EDGE_DIRECTIONS);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckInterfacePolygonsRequirements(Message& msg) const
{
    // If no tessellation then report it
    if (!m_tessellation)
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    // If no blocks then report it
    if (!m_assembly || m_assembly->CountBlocks() == 0)
    {
        msg.set(Message::MSG_NO_BLOCKS);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckRotateTessellationRequirements(Message& msg) const
{
    // If no tessellation then report it
    if (!m_tessellation)
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckScaleBlocksRequirements(Message& msg) const
{
    // If no blocks then report it
    if (!m_assembly || m_assembly->CountBlocks() == 0)
    {
        msg.set(Message::MSG_NO_BLOCKS);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckScaleTessellationRequirements(Message& msg) const 
{
    // If no tessellation then report it
    if (!m_tessellation)
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckSetEdgeDirectionsRequirements(Message& msg) const
{
    // If no tessellation then report it
    if (!m_tessellation)
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    if (!m_tessellation->IsDCELLoaded())
    {
        m_tessellation->LoadDCEL();
    }
    
    // Get the reference to the geometry of the tessellation
    std::shared_ptr<dcel::DCEL> geometry = m_tessellation->DCEL();

    std::shared_ptr<dcel::Face> face = nullptr;

    // Traverse through the tiles of the tessellation and check they have the center dynamic 
    // attribute
    for (auto it = geometry->Faces().begin(); it != geometry->Faces().end(); ++it)
    {
        // Get the pointer to the current face
        face = *it;

        // If the face does not have an even number of sides then report it
        if (!face->HasEvenNumberOfSides())
        {
            msg.set(Message::MSG_NO_EVEN_SIDED_TILES);
            return false;
        }
    }

    geometry = nullptr;

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckSetFaceCentersRequirements(Message& msg) const 
{
    // If no tessellation then report it
    if (!m_tessellation)
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckTessellationVertexNormalizationRequirements(Message& msg) const
{
    // If no tessellation then report it
    if (!m_tessellation)
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckTileSubdivisionRequirements(Message& msg) const
{
    // If no tessellation then report it
    if (!m_tessellation)
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckTiltingAngleRequirements(Message& msg) const 
{
    // If no tessellation then report it
    if (!m_tessellation)
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    // If not all tiles have an even number of sides then report it
    if (!m_tessellation->CheckTilesAreEvenSided())
    {
        msg.set(Message::MSG_NO_EVEN_SIDED_TILES);
        return false;
    }

    // If not all tiles have the ATTRIB_CENTER attribute then report ot
    if (!m_tessellation->CheckTileCenters())
    {
        msg.set(Message::MSG_NO_TILE_CENTERS);
        return false;
    }

    // If not all half edges have both ATTRIB_DIRECTION_VALUE and ATTRIB_DIRECTION_VECTOR 
    // attributes then report it
    if (!m_tessellation->CheckHalfedgeDirections())
    {
        msg.set(Message::MSG_NO_EDGE_DIRECTIONS);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckTruncateBlocksRequirements(Message& msg) const 
{
    // If no tessellation then report it
    if (!m_tessellation)
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    // If no blocks then report it
    if (!m_assembly || m_assembly->CountBlocks() == 0)
    {
        msg.set(Message::MSG_NO_BLOCKS);
        return false;
    }

    // If a block has been modified then report it
    if (m_assembly->AreBlocksModified()) 
    {
        msg.set(Message::MSG_BLOCKS_MODIFIED);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

bool Workspace::CheckVertexNormalizationRequirements(Message& msg) const 
{
    // If no tessellation then report it
    if (!m_tessellation) 
    {
        msg.set(Message::MSG_NO_TESSELLATION);
        return false;
    }

    // Return a positive result since all requirements are satisfied
    msg.set(Message::MSG_OK);
    return true;
}

void Workspace::Clear()
{
    // Erase the tessellation from the workspace (this also erases all subsequent content from
    // the workspace as well)
    EraseTessellation();
}

void Workspace::EraseAssembly()
{
    // Erase the interface polygons from the workspace
    EraseInterfaces();

    if (m_assembly)
    {
        m_assembly->Clear();
        m_assembly = nullptr;
    }
}

void Workspace::EraseEdgeDirections() 
{
    // Erase the assembly from the workspace (interface polygons are removed as well)
    EraseAssembly();

    // If there is a tessellation then remove the direction values and direction vectors
    if (m_tessellation && m_tessellation->IsDCELLoaded()) 
    {
        std::shared_ptr<dcel::DCEL> geometry = m_tessellation->DCEL();

        std::shared_ptr<dcel::Halfedge> halfedge = nullptr;

        // Traverse through the half edges of the tessellation and remove the dynamic 
        // attributes related to the edge directions
        for (auto it = geometry->Halfedges().begin(); it != geometry->Halfedges().end(); ++it)
        {
            halfedge = *it;

            halfedge->Attributes().Erase(ATTRIB_DIRECTION_VALUE);
            halfedge->Attributes().Erase(ATTRIB_DIRECTION_VECTOR);
        }

        halfedge = nullptr;
        geometry = nullptr;
    }
}

void Workspace::EraseEquilibriumResults()
{
    if (m_equilibriumResults) 
    {
        m_equilibriumResults->Reset();
        m_equilibriumResults = nullptr;
    }
}

void Workspace::EraseTessellation() 
{
    // Erase the face centers from the workspace (edge directions, blocks and interface polygons 
    // are removed as well)
    EraseTileCenters();

    if (m_tessellation) 
    {
        m_tessellation->Clear();
        m_tessellation = nullptr;
    }
}

void Workspace::EraseInterfaces()
{
    EraseEquilibriumResults();

    if (m_interfaces) 
    {
        m_interfaces->Clear();
        m_interfaces = nullptr;
    }
}

void Workspace::EraseTileCenters()
{
    // Erase the edge directions from the workspace (blocks and interface 
    // polygons are removed as well)
    EraseEdgeDirections();

    // If there is a tessellation then erase the face centers from the 
    // tessellation
    if (m_tessellation && m_tessellation->IsDCELLoaded())
    {
        std::shared_ptr<dcel::DCEL> geometry = m_tessellation->DCEL();

        // Traverse through the tiles of the tessellation and remove the 
        // dynamic attributes related to the face centers
        for (auto it = geometry->Faces().begin(); it != geometry->Faces().end(); ++it)
        {
            (*it)->Attributes().Erase(ATTRIB_CENTER);
        }

        geometry = nullptr;
    }
}

void Workspace::GenerateEquilibriumSystem(
    const std::list<std::tuple<size_t, size_t>> & R,
    std::vector<Eigen::Vector3d> & C,
    std::vector<std::vector<double>> & W,
    double density, 
    double gravity, 
    std::vector<std::shared_ptr<VF>> & I, 
    std::list<std::tuple<size_t, size_t>> & nR) const
{
    // Clear the given vectors
    C.clear();
    W.clear();
    I.clear();
    nR.clear();

    // Define the map for relating the original block indices and their new 
    // indices
    std::map<size_t, size_t> blockIndices;

    // Define the map for relating the original interface polygon indices and 
    // their new indices
    std::map<size_t, size_t> interfaceIndices;

    // Traverse through the given relations
    for (auto it = R.begin(); it != R.end(); ++it) 
    {
        // Get the original index of the block
        size_t origblockIndex = std::get<0>(*it);

        // Initialize the new index of the block
        size_t newblockIndex;

        // Search for the original block index in the block indices map
        auto findblockIndex = blockIndices.find(origblockIndex);

        // If the original block index is not in the block indices map then 
        // define the new block index and insert both original and new block 
        // indices into the block indices map
        if (findblockIndex == blockIndices.end()) 
        {
            // Define the new block index
            newblockIndex = blockIndices.size();

            // Store both original and new block indices into the block indices
            // map
            blockIndices.insert(std::make_pair(origblockIndex, newblockIndex));

            // Get the centroid of the block and insert it in the centroids vector
            C.push_back(m_assembly->GetBlocks()[origblockIndex]->Geometry().centroid());

            // Get the loads of the block and insert it in the weights vector
            //W.push_back(m_assembly->GetBlocks()[origblockIndex]->Loads(density, gravity));
            const std::array<double, 6>& blockLoads = m_assembly->GetBlocks()[origblockIndex]->GetLoads();
            W.push_back({ blockLoads[0], blockLoads[1], blockLoads[2], blockLoads[3], blockLoads[4], blockLoads[5] });

            // Get the pointer to the current block, indicate it is used in the
            // equilibrium analysis and set its new index (as it will be used 
            // during the equilibrium analysis)
            std::shared_ptr<Block> block = m_assembly->GetBlocks()[origblockIndex];
            block->Attributes().Set<bool>(ATTRIB_IN_EQUILIBRIUM, true);
            block->Attributes().Set<size_t>(ATTRIB_EQUILIBRIUM_INDEX, newblockIndex);
        }
        else 
        {
            // Get the new block index
            newblockIndex = findblockIndex->second;
        }

        // Get the original index of the interface
        size_t origInterfaceIndex = std::get<1>(*it);

        // Initialize the new index of the interface
        size_t newInterfaceIndex;

        // Search for the original interface index in the interface indices map
        auto findInterfaceIndex = interfaceIndices.find(origInterfaceIndex);

        // If the original interface index is not in the interface indices map then define the new 
        // interface index and insert both original and new interface indices into the interface 
        // indices map
        if (findInterfaceIndex == interfaceIndices.end())
        {
            // Define the new interface index
            newInterfaceIndex = interfaceIndices.size();

            // Store both original and new interface indices into the interface indices map
            interfaceIndices.insert(std::make_pair(origInterfaceIndex, newInterfaceIndex));

            // Get the geometry of the interface and insert it in the interfaces vector
            I.push_back(m_interfaces->GetInterfaces()[origInterfaceIndex]);

            // Get the pointer to the interface polygon, indicate it is used in the equilibrium 
            // analysis and set its new index (as it will be used during the equilibrium analysis)
            auto intf = m_interfaces->GetInterfaces()[origInterfaceIndex];
            intf->Attributes().Set<bool>(ATTRIB_IN_EQUILIBRIUM, true);
            intf->Attributes().Set<size_t>(ATTRIB_EQUILIBRIUM_INDEX, newInterfaceIndex);
        }
        else
        {
            // Get the new interface index
            newInterfaceIndex = findInterfaceIndex->second;
        }

        // Define a new relations tuple using the new block and interface indices
        nR.push_back(std::make_tuple(newblockIndex, newInterfaceIndex));
    }
}

std::shared_ptr<Assembly> Workspace::GetAssembly() const 
{
    return m_assembly;
}

bool Workspace::GetFaceRotatedPlanes(
    std::shared_ptr<dcel::Face> face, 
    std::vector<toolkit::Plane> & planes) const
{
    // Clear the vector where the planes will be stored
    planes.clear();

    // Get the pointer to the incident half edge of the face
    std::shared_ptr<dcel::Halfedge> currentHalfedge = face->halfedge;

    // Traverse through the half edges of the face
    do
    {
        // If the current half edge doesn't have a rotated vector then return null
        if (!currentHalfedge->Attributes().Has(ATTRIB_ROTATED_VECTOR))
        {
            // Clear any content in the planes vector and return it empty
            planes.clear();
            return false;
        }

        // Get the reference to the rotated vector of the current half edge
        Eigen::Vector3d R;
        assert(currentHalfedge->Attributes().Get<Eigen::Vector3d>(ATTRIB_ROTATED_VECTOR, R));

        // Get the reference to the midpoint of the current half edge
        Eigen::Vector3d M;
        assert(currentHalfedge->Attributes().Get<Eigen::Vector3d>(ATTRIB_MIDPOINT, M));

        // Define the plane object by setting its components and push it into the planes array
        planes.emplace_back(R, M);

        // Move to the next half edge
        currentHalfedge = currentHalfedge->next;

    } while (currentHalfedge != face->halfedge);

    // Return the array with the planes of the face
    return true;
}

std::shared_ptr<Tessellation> Workspace::GetTessellation() const
{
	return m_tessellation;
}

bool Workspace::HasEquilibriumResults() const
{
    return m_equilibriumResults != nullptr;
}

bool Workspace::HasInterfacePolygons() const
{
    return m_interfaces != nullptr;
}

std::shared_ptr<InterfacePolygons> Workspace::GetInterfacePolygons() const
{
    return m_interfaces;
}

bool Workspace::GetBlockInterfaceRelations(std::list<std::tuple<size_t, size_t>> & PI) const 
{
    // Clear the block-interface relations list
    PI.clear();

    // Get the number of interface polygons
    size_t nInterfaces = m_interfaces->CountInterfaces(), p1Index = 0, p2Index = 0;

    // Traverse through the interface polygons
    for (size_t i = 0; i < nInterfaces; i += 1)
    {
        // Get the pointer to the current interface polygon
        std::shared_ptr<VF> currentInterface = m_interfaces->GetInterfaces()[i];

        if (!currentInterface)
        {
            continue;
        }

        // Get the index of the block located at the positive side of the interface polygon
        assert(currentInterface->Attributes().Get<size_t>(ATTRIB_POSITIVE_BLOCK_INDEX, p1Index));

        // If the block is enabled then store the tuple between the block and the current interface
        // polygon
        if (m_assembly->GetBlocks()[p1Index]->IsEnabled())
        {
            PI.push_back(std::make_tuple(p1Index, i));
        }

        // Get the index of the block located at the negative side of the interface polygon
        assert(currentInterface->Attributes().Get<size_t>(ATTRIB_NEGATIVE_BLOCK_INDEX, p2Index));

        // If the block is enabled then store the tuple between the block and the current interface
        // polygon
        if (m_assembly->GetBlocks()[p2Index]->IsEnabled())
        {
            PI.push_back(std::make_tuple(p2Index, i));
        }

        currentInterface = nullptr;
    }

    // Indicate if at least one tuple was generated
    return PI.size() > 0;
}

bool Workspace::GetBlocksFromRotatedVectors(std::vector<VF> & blocks) const 
{
    // There must exist a tessellation
    assert(m_tessellation);
    assert(m_tessellation->IsDCELLoaded());

    // Get the number of tiles of the geometry
    size_t nTiles = m_tessellation->DCEL()->Faces().size();

    // Clear the given vector and resize it
    blocks.clear();
    blocks.resize(nTiles);

    // Traverse through the tiles of the tessellation and generate the blocks
    for (size_t i = 0; i < nTiles; i += 1)
    {
        // 
        blocks[i].Set(GetBlock(m_tessellation->DCEL()->Faces()[i]));
    }

    // Indicate that at least one block was generated
    return blocks.size() > 0;
}

std::shared_ptr<EquilibriumAnalysis::Result> Workspace::GetEquilibriumResults() const
{
    return m_equilibriumResults;
}

bool Workspace::HasTessellation() const 
{
    return (m_tessellation != nullptr);
}

bool Workspace::HasRotatedVectors() const 
{
    if (!m_tessellation)
    {
        return false;
    }

    if (!m_tessellation->IsDCELLoaded()) 
    {
        return false;
    }

    // Get the reference to the geometry of the tessellation
    std::shared_ptr<dcel::DCEL> geometry = m_tessellation->DCEL();

    // Traverse through the half edges of the tessellation
    for (auto it = geometry->Halfedges().begin(); it != geometry->Halfedges().end(); ++it)
    {
        // If the current half edge does not have the rotated vector then return false
        if (!(*it)->Attributes().Has(ATTRIB_ROTATED_VECTOR))
        {
            return false;
        }
    }

    geometry = nullptr;

    // Return true since all half edges have the rotated vector dynamic attribute
    return true;
}

bool Workspace::HeightBisectionMethod(double height, bool boundary, std::vector<VF> & blocks) const 
{
    // There must exist a tessellation
    assert(m_tessellation);

    // Calculate the top and bottom points of the tiles in the tessellation using the given 
    // height value
    m_tessellation->SetTileTopBottomSectionPoints(height);

    // Calculate the rotated vectors using the top and bottom section points of each face in the 
    // tessellation, then set them as attributes of the half edges
    m_tessellation->SetEdgeRotatedVectorsUsingTopBottomSectionPoints(boundary);

    // Generate the blocks from the intersection of the planes defined by the rotated vectors in 
    // the half edges. Then return ...
    return GetBlocksFromRotatedVectors(blocks);
}

VF Workspace::GetBlock(const std::shared_ptr<dcel::Face> face) const
{
    // Initialize the lists to store the blocks and tiles of the block
    std::list<Eigen::Vector3d> blockVertices;
    std::list<std::vector<size_t>> blocktiles;

    // 
    GetBlockElements(face, blockVertices, blocktiles);

    // 
    VF vf(blockVertices, blocktiles);

    // 
    blockVertices.clear();
    blocktiles.clear();

    return vf;
}

void Workspace::GetBlockElements(
    const std::shared_ptr<dcel::Face> face, 
    std::list<Eigen::Vector3d> & blockVertices, 
    std::list<std::vector<size_t>> & blocktiles) const
{
    // 
    blockVertices.clear();
    blocktiles.clear();

    // Initialize the vector where the rotated planes incident to the half edges of the face will 
    // be stored
    std::vector<toolkit::Plane> planes;

    // Get the rotated planes from the current face. If no rotated planes then move to the next
    // face in the tessellation
    if (!GetFaceRotatedPlanes(face, planes))
    {
        return;
    }

    // Initialize the vertex coordinates and vertex indices that describe the geometry of the
    // block from the current face
    //VF vf;

    std::vector<Eigen::Vector3d> points;

    // Get the intersection points from the planes. These are the vertices of the block. If no 
    // intersection points then move to the next face
    if (!algorithms::getIntersectionPointsFromPlanes(planes, points))
    {
        return;
    }

    // Store the vertices of the block in the block vertices list
    for (auto it = points.begin(); it != points.end(); ++it) 
    {
        blockVertices.push_back(*it);
    }

    // Get the direction value of the incident half edge of the face. This direction value is 
    // associated to the first plane in the planes array and therefore to the first point in 
    // the intersection points array
    int dir;
    assert(face->halfedge->Attributes().Get<int>(ATTRIB_DIRECTION_VALUE, dir));

    // Get the number of intersection points
    size_t nPoints = points.size();

    // Traverse through the intersection points and generate the tiles of the triangular strip 
    // of the block
    for (size_t i = 0; i < nPoints; i += 1)
    {
        // Calculate the indexes of the previous, current and next points
        size_t prevIndex = (i == 0) ? nPoints - 1 : i - 1;
        size_t currIndex = i;
        size_t nextIndex = (i == nPoints - 1) ? 0 : i + 1;

        // Initialize the vector for storing the vertex indices for the current face of the 
        // block. Then, get its reference
        std::vector<size_t> indices;

        // If the direction is positive it means the intersection point is above the face. 
        // Otherwise, it is below the face. This is crucial for the proper order of the vertex 
        // indexes of the face
        if (dir > 0)
        {
            indices.emplace_back(currIndex);
            indices.emplace_back(prevIndex);
            indices.emplace_back(nextIndex);
        }
        else
        {
            indices.emplace_back(currIndex);
            indices.emplace_back(nextIndex);
            indices.emplace_back(prevIndex);
        }

        // Toggle the direction for the next vertex
        dir *= -1;

        // 
        blocktiles.push_back(indices);
    }

    // If the block has more than four vertices then the top and bottom tiles have to be 
    // defined. Such tiles are polygons with half number of sides with respect to the face
    if (nPoints > 4)
    {
        std::vector<size_t> topIndices;

        // Define the vertex indices for the top face of the block
        for (size_t j = (dir > 0) ? 0 : 1; j < nPoints; j += 2)
        {
            //topFace.emplace_back(j);
            topIndices.emplace_back(j);
        }

        blocktiles.push_back(topIndices);

        // Initialize the vector for the vertex indices of the bottom face of the block. Then, 
        // get its reference
        std::vector<size_t> bottomIndices;

        // Define the vertex indices for the bottom face of the block
        for (size_t j = (dir > 0) ? nPoints - 1 : nPoints - 2; j >= 0; j -= 2)
        {
            //bottomFace.emplace_back(j);
            bottomIndices.emplace_back(j);

            // Break the cycle if j is less or equal to 1 (size_t variables do not handle 
            // negative values
            if (j <= 1)
            {
                break;
            }
        }

        blocktiles.push_back(bottomIndices);
    }
}

bool Workspace::RunEquilibriumAnalysis(
    double density,
    double friction, 
    const std::string lengthUnit, 
    bool verbose, 
    bool files, 
    double cWeight,
    double tWeight,
    double uWeight,
    double vWeight) 
{
    // Get the relations between the blocks and the interface polygons. Exit the function if there 
    // are no relations
    std::list<std::tuple<size_t, size_t>> PI;
    if (!GetBlockInterfaceRelations(PI)) 
    {
        return false;
    }

    // Initialize the vector for storing the centroids of the blocks
    std::vector<Eigen::Vector3d> centroids;

    // Initialize the vector for storing the loads on the blocks
    std::vector<std::vector<double>> loads;

    // Initialize the vector for storing the pointers to the interface polygons
    std::vector<std::shared_ptr<VF>> interfaces;

    // Initialize the list for storing the 
    std::list<std::tuple<size_t, size_t>> nR;

    // Determine the value for gravity based on the given length unit
    double gravity = ((lengthUnit.compare("m") == 0) ? -9.8 : -980);
    std::cout << "Gravity = " << gravity << std::endl;

    // 
    GenerateEquilibriumSystem(PI, centroids, loads, density, gravity, interfaces, nR);

    // Generate the files if indicated
    if (files) 
    {
        WriteCentroidsFile(centroids);
        //writeblocksFile();
        WriteInterfacesFile(interfaces);
        WriteRelationsFile(nR);
        WriteGeogebraJsFile(centroids, interfaces);
    }

    // Define the prefix for the files associated to the current model
    std::stringstream ss;
    ss << "model_" << density << "_" << friction;
    std::string prefix = ss.str();

    // Initialize a new Equilibrium Analysis object
    EraseEquilibriumResults();
    m_equilibriumResults = std::make_shared<EquilibriumAnalysis::Result>();
    
    // Run the equilibrium analysis
    return EquilibriumAnalysis::Run(
        centroids, 
        loads, 
        interfaces, 
        nR, 
        friction, 
        *m_equilibriumResults,
        verbose, 
        files, 
        prefix, 
        cWeight, 
        tWeight, 
        uWeight, 
        vWeight);
}

bool Workspace::SaveAsTicFile(const std::string filename, const std::string directory) const 
{
    // Build the file path
    std::string filepath = directory + "/" + filename + ".tic";

    // Open the file
    std::fstream file;
    file.open(filepath.c_str(), std::fstream::in | std::fstream::trunc);

    // Initialize the string stream where the content will be written
    std::stringstream ss;

    // Write the content of the workspace in the string stream
    Write(ss);

    // Convert the content to string and write it to the file
    file << ss.str();

    // Close the file
    file.close();

    // Indicate everything went fine
    return true;
}

void Workspace::ScaleTessellation(double scale)
{
    assert(m_tessellation);
    
    EraseTileCenters();
    
    m_tessellation->Geometry().Scale(scale);
}

void Workspace::SetAssembly(std::vector<VF> & blocks, bool disableBoundary)
{
    assert(blocks.size() == m_tessellation->DCEL()->Faces().size());

    // Remove the assembly from the workspace. This removes interface polygons as well
    EraseAssembly();

    // Set the blocks to the assembly
    m_assembly = std::make_shared<Assembly>(blocks);

    // Disable the blocks at the boundary if indicated
    if (disableBoundary) 
    {
        assert(m_tessellation);
        m_assembly->ToggleBlocksByBoundary(m_tessellation);
    }
}

std::vector<size_t> Workspace::SelectedBlocks(const std::vector<std::string>& strings) const 
{
    size_t nBlocks = m_assembly->CountBlocks();

    std::vector<size_t> selected;

    for (auto it = strings.begin(); it != strings.end(); ++it) 
    {
        // If the current value is "all" then reset the selected vector and 
        // populate it with all block indices
        if ((*it).compare("all") == 0) 
        {
            selected = std::vector<size_t>(nBlocks);

            for (size_t i = 0; i < nBlocks; i += 1) 
            {
                selected[i] = i;
            }

            return selected;
        }

        // Split the current string into tokens. The number of resultant tokens
        // must be either 1 or 2
        std::vector<std::string> tokens = utils::split(*it, '-');
        assert(tokens.size() == 1 || tokens.size() == 2);

        //
        if (tokens.size() == 1) 
        {
            // Get the index. It must be a value less than the number of blocks
            // of the assembly.
            size_t index = std::stoi(tokens[0]);
            assert(index < nBlocks);

            selected.push_back(std::stoi(tokens[0]));
        }
        else 
        {
            // Get the indices of the indicated range. Both indices must be 
            // less than the number of blocks of the assembly. Also, the first 
            // index must be less than the second index.
            size_t first = std::stoi(tokens[0]);
            size_t second = std::stoi(tokens[1]);
            assert(first < nBlocks);
            assert(second < nBlocks);
            assert(first < second);

            for (size_t i = first; i <= second; i += 1) 
            {
                selected.push_back(i);
            }
        }
    }

    // Verify at least one block was selected
    assert(selected.size() > 0);

    return selected;
}

void Workspace::SetAssemblyUsingAdaptiveTileOffsetClippedBlocks(
    const std::string & topFunction, 
    const std::string & bottomFunction, 
    double threshold)
{
    // There must be a tessellation and blocks in the assembly
    assert(m_tessellation);
    assert(m_assembly);
    assert(m_assembly->CountBlocks() > 0);

    // Clip the blocks using the Adaptive Tile Offset Clipping method
    std::vector<VF> blocks;
    assert(AdaptiveTileOffsetClipping(topFunction, bottomFunction, blocks, threshold));

    // Set the blocks to the assembly
    SetAssembly(blocks);
}

void Workspace::SetAssemblyUsingHeightBisectionMethod(double height, bool boundary) 
{
    // There must be a tessellation
    assert(m_tessellation);

    // Generate the blocks using the Height-Bisection method
    std::vector<VF> blocks;
    assert(HeightBisectionMethod(height, boundary, blocks));

    // Set the blocks to the assembly
    SetAssembly(blocks);
}

void Workspace::SetAssemblyUsingTileOffsetClippedBlocks(
    double extrados, 
    double intrados, 
    double threshold)
{
    // There must be a tessellation and blocks in the assembly
    assert(m_tessellation);
    assert(m_assembly);
    assert(m_assembly->CountBlocks() > 0);

    // Clip the blocks using the Tile Offset Clipping method
    std::vector<VF> blocks;
    assert(TileOffsetClipping(extrados, intrados, blocks, threshold));

    // Set the blocks to the assembly
    SetAssembly(blocks);
}

void Workspace::SetAssemblyUsingTiltingAngleMethod(double angle)
{
    // There must be a tessellation
    assert(m_tessellation);

    // Generate the blocks using the Tilting Angle method
    std::vector<VF> blocks;
    assert(TiltingAngleMethod(angle, blocks));

    // Set the blocks to the assembly
    SetAssembly(blocks);
}

void Workspace::SetTessellation(const VF & vf)
{
    // Clear the content of the workspace (a new tessellation invalidates anything in both 
    // workspace and scene
    Clear();

    // Set the vertex coordinates and vertex indices to the tessellation
    m_tessellation = std::make_shared<Tessellation>(vf);
}

void Workspace::SetTessellationEdgeDirections(int dir) 
{
    // Tehere must exists a tessellation
    assert(m_tessellation);

    // Erase the direction values and direction vectors from the both scene and edges of the 
    // tessellation
    EraseEdgeDirections();

    // Set the direction values, direction vectors, and midpoints of the edges of the tessellation
    m_tessellation->SetEdgesDirection(dir);
}

void Workspace::SetTessellationTileCenters()
{
    // There must be a tessellation
    assert(m_tessellation);

    // Erase the face centers from the workspace
    EraseTileCenters();

    // Set the specified center to the tiles of the tessellation
    m_tessellation->SetTilesCenter();
}

void Workspace::SetTessellationScale(double scale) 
{
    // There must be a tessellation
    assert(m_tessellation);

    // Scale the tessellation using the given scale factor
    m_tessellation->Geometry().Scale(scale);

    // Get the vertex coordinates and vertex indices of the scaled tessellation
    VF vf = m_tessellation->Geometry();
    
    // Set the new geometry as the new tessellation
    SetTessellation(vf);
}

void Workspace::SetInterfacePolygons(std::vector<std::shared_ptr<VF>> & interfaces)
{
    // Assert there are interface polygons
    assert(interfaces.size() > 0);

    // Remove the interface polygons from the workspace
    EraseInterfaces();

    // 
    m_interfaces = std::make_shared<InterfacePolygons>();
    m_interfaces->Set(interfaces);
}

void Workspace::SubdivideTessellationTiles(Tessellation::TILE_SUBDIVISION_TYPE type)
{
    assert(m_tessellation);

    switch (type) 
    {
		case Tessellation::TILE_SUBDIVISION_TYPE::QUADRILATERALS:
		{
			m_tessellation->QuadrangulateTiles();
			break;
		}

		case Tessellation::TILE_SUBDIVISION_TYPE::MIDPOINTS: 
		{
			m_tessellation->SubdivideTilesByMidpoints();
			break;
		}

        case Tessellation::TILE_SUBDIVISION_TYPE::TRIANGLES:
        {
			m_tessellation->TriangulateTilesByVertices();
            break;
        }
    }

    VF vf = m_tessellation->Geometry().clone();

    SetTessellation(vf);
}

bool Workspace::TileOffsetClipping(
    double extrados, 
    double intrados, 
    std::vector<VF> & blocks, 
    double threshold) const
{
    // There must exist a tessellation and blocks in the assembly
    assert(m_tessellation);
    assert(m_assembly);
    assert(m_assembly->CountBlocks() > 0);

    // 
    m_assembly->TileOffsetClipping(m_tessellation, extrados, intrados, blocks, threshold);

    // 
    return blocks.size() > 0;
}

bool Workspace::TiltingAngleMethod(double angle, std::vector<VF> & blocks) const 
{
    // There must exists a tessellation
    assert(m_tessellation);
    
    // Calculate the rotated vectors using the given angle value and set them as attributes of the 
    // half edges in the DCEL
    m_tessellation->SetEdgeRotatedVectorsUsingAngle(angle);

    // Generate the blocks from the intersection of the planes defined by the rotated vectors in 
    // the half edges
    return GetBlocksFromRotatedVectors(blocks);
}

void Workspace::Write() const
{
    // If there is a tessellation then write its content
    if (m_tessellation)
    {
        m_tessellation->Write();
    }

    // If there is an assembly then write its content
    if (m_assembly)
    {
        m_assembly->Write();
    }

    // If there are interfaces then write their content
    if (m_interfaces)
    {
        m_interfaces->Write();
    }
}

void Workspace::Write(std::stringstream & ss) const 
{
    // If there is a tessellation then write its content
    if (m_tessellation) 
    {
        m_tessellation->Write(ss);
    }

    // If there is an assembly then write its content
    if (m_assembly) 
    {
        m_assembly->Write(ss);
    }

    // If there are interfaces then write their content
    if (m_interfaces) 
    {
        m_interfaces->Write(ss);
    }
}

void Workspace::WriteAllAssemblyGeogebraJsFile(const std::string filename, int r, int g, int b, double threshold) const
{
    // There must exist an assembly
    assert(m_assembly);

    // Initialize the string stream for storing the commands that generate the geometry of the 
    // assembly
    std::stringstream ss;

    // Generate and write the JavaScript commands that generate the geometry of the tessellation in
    // GeoGebra
    m_tessellation->Geometry().WriteGeogebraJs(ss, "geomdom", r, g, b, 1, threshold);

    // Open the file, write its content and close it
    std::ofstream file(filename, std::ios::out);
    file << ss.str() << std::endl;
    file.close();
}

void Workspace::WriteCentroidsFile(
    const std::vector<Eigen::Vector3d> & C,
    const std::string filename) const 
{
    // Get the number of centroids
    size_t nCentroids = C.size();

    // Open the file
    std::ofstream file(filename, std::ios::out);

    // Traverse through the centroids and write them in the file
    for (size_t i = 0; i < nCentroids; i += 1) 
    {
        file << i << " (" << C[i].x() << ", " << C[i].y() << ", " << C[i].z() << ")" << std::endl;
    }

    // Close the file
    file.close();
}

void Workspace::WriteGeogebraJSGeometryFile(const std::string & filename, double threshold) const
{
    // Initialize the string stream for writing the commands that generate the geometries
    std::stringstream ss;

    // If there is a tessellation then write its content
    if (m_tessellation) 
    {
        if (!m_tessellation->IsDCELLoaded()) 
        {
            m_tessellation->LoadDCEL();
        }

        // Get the reference to the geometry of the tessellation
        std::shared_ptr<dcel::DCEL> geometry = m_tessellation->DCEL();

        // Write the GepGebra commands that generate the geometry of the tessellation
        geometry->WriteGeogebraJs(ss, "geomdom", 128, 128, 128, 1, threshold);

        geometry = nullptr;
    }

    // If there are blocks then write their content
    if (m_assembly) 
    {
        // Initialize the string stream for defining the prefixes of the blocks
        std::stringstream blockPrefix;

        // Get the number of blocks of the assembly
        size_t nblocks = m_assembly->CountBlocks();

        // Traverse through the blocks and write the GeoGebra commands that generate them
        for (size_t i = 0; i < nblocks; i += 1) 
        {
            // Clear the string stream for the block prefixes
            blockPrefix.str(std::string());
            blockPrefix.clear();

            // Define the prefix for the current block
            blockPrefix << "block" << i;

            // Write the GeoGebra commands that generate the geometry of the current block
            m_assembly->GetBlocks()[i]->Geometry().WriteGeogebraJs(
                ss, blockPrefix.str(), 255, 0, 0, 1, threshold);
        }
    }

    // If there are interface polygons then write their content
    if (m_interfaces) 
    {
        // Initialize the string stream for defining the prefixes of the interface polygons
        std::stringstream interfacePrefix;

        // Get the number of interface polygons between the blocks
        size_t ninterfaces = m_interfaces->CountInterfaces();

        // Traverse through the interface polygons and write the GeoGebra commands that generate 
        // them
        for (size_t i = 0; i < ninterfaces; i += 1)
        {
            // Clear the string stream for the interface prfixes
            interfacePrefix.str(std::string());
            interfacePrefix.clear();

            // Define the prefix for the current interface polygon
            interfacePrefix << "interface" << i;

            // Write the GeoGebra commands that generate the geometry of the current interface 
            // polygon
            m_interfaces->GetInterfaces()[i]->WriteGeogebraJs(
                ss, interfacePrefix.str(), 255, 255, 0, 1, threshold);
        }
    }

    // Open the file
    std::ofstream file(filename, std::ios::out);

    // Write the geometries into the file. Then, close it
    file << ss.str() << std::endl;
    file.close();
}

void Workspace::WriteGeogebraJsFile(
    const std::vector<Eigen::Vector3d> & C, 
    const std::vector<std::shared_ptr<VF>> & I, 
    const std::string filename) const 
{
    /*assert(C.size() > 0);
    assert(I.size() > 0);

    // Open the file
    std::ofstream file(filename, std::ios::out);

    // Get the number of centroids
    size_t nCentroids = C.size();

    // Traverse through the centroids and write their content
    for (size_t centrIdx = 0; centrIdx < nCentroids; centrIdx += 1)
    {
        // Get the reference to the coordinates of the current centroid
        Eigen::Vector3d c = C[centrIdx];
        utils::FixZeros(c);

        // Write the command that generates the current centroid, then set its color and size
        file << "ggbApplet.evalCommand(\"c" << centrIdx << " = Point({" << c.x() << ", " << 
            c.y() << ", " << c.z() << "})\");" << std::endl;
        file << "ggbApplet.setColor(\"c" << centrIdx << "\", 0, 0, 255);" << std::endl;
        file << "ggbApplet.setPointSize(\"c" << centrIdx << "\", 9);" << std::endl;
    }

    // Get the number of interfaces
    size_t ninterfaces = I.size();

    // Traverse through the interfaces and write their content in the file
    for (size_t intfIdx = 0; intfIdx < ninterfaces; intfIdx += 1)
    {
        // Get the reference to the geometry of the current interface
        const dcel::DCEL & geometry = I[intfIdx]->GetGeometry();

        // Get the number of vertices of the interface
        size_t nVertices = geometry.vertices.size();

        // Traverse through the vertices of the interface and write their content
        for (size_t vertexIdx = 0; vertexIdx < nVertices; vertexIdx += 1)
        {
            // Get the reference to the coordinates of the current vertex
            Eigen::Vector3d c = geometry.vertices[vertexIdx]->Coords();
            utils::FixZeros(c);

            // Write the current vertex
            file << "ggbApplet.evalCommand(\"i" << intfIdx << "v" << vertexIdx << " = Point({" << 
                c.x() << ", " << c.y() << ", " << c.z() << "})\");" << std::endl;

            // Set the id to the vertex
            geometry.vertices[vertexIdx]->Attributes().Set<size_t>(ATTRIB_INDEX, vertexIdx);
        }

        // Get the number of tiles of the geometry (there should be one though)
        size_t ntiles = geometry.tiles.size();

        // Traverse through the tiles of the interface and write their content
        for (size_t faceIdx = 0; faceIdx < ntiles; faceIdx += 1)
        {
            // Get the pointer to the incident half edge of the face
            std::shared_ptr<dcel::Halfedge> halfedge = geometry.tiles[faceIdx]->halfedge;

            // Start the line for the current interface face
            file << "ggbApplet.evalCommand(\"i" << intfIdx << " = Polygon(";

            bool firstVertex = true;

            // Traverse through the half edges of the face
            do
            {
                // Get the index of the start vertex of the half edge
                size_t vertexIdx;
                assert(halfedge->start->Attributes().Get<size_t>(ATTRIB_INDEX, vertexIdx));

                if (firstVertex)
                {
                    file << "i" << intfIdx << "v" << vertexIdx;

                    firstVertex = false;
                }
                else
                {
                    file << ", i" << intfIdx << "v" << vertexIdx;
                }

                // Move to the next half edge
                halfedge = halfedge->next;

            } while (halfedge != geometry.tiles[faceIdx]->halfedge);

            // Close the line for the current interface face
            file << ")\");" << std::endl;

            // Write the centroid of the interface polygon
            file << "ggbApplet.evalCommand(\"ci" << intfIdx << " = Centroid(i" << intfIdx << 
                ")\");" << std::endl;
            file << "ggbApplet.setColor(\"ci" << intfIdx << "\", 255, 0, 0);" << std::endl;
            file << "ggbApplet.setPointSize(\"ci" << intfIdx << "\", 9);" << std::endl;

            // Write the normal vector of the interface polygon
            Eigen::Vector3d N = geometry.tiles[faceIdx]->Normal(true, true);
            file << "ggbApplet.evalCommand(\"Ni" << intfIdx << " = UnitVector(Segment(ci" << 
                intfIdx << ", ci" << intfIdx << " + (" << N.x() << ", " << N.y() << ", " << 
                N.z() << ")))\");" << std::endl;
            file << "ggbApplet.setColor(\"Ni" << intfIdx << "\", 0, 255, 255);" << std::endl;

            // Write the tangential U vector of the interface polygon
            Eigen::Vector3d U = geometry.tiles[faceIdx]->halfedge->Direction(true, true);
            file << "ggbApplet.evalCommand(\"Ui" << intfIdx << " = UnitVector(Segment(ci" << 
                intfIdx << ", ci" << intfIdx << " + (" << U.x() << ", " << U.y() << ", " << 
                U.z() << ")))\");" << std::endl;
            file << "ggbApplet.setColor(\"Ui" << intfIdx << "\", 255, 0, 255);" << std::endl;

            Eigen::Vector3d V = N.cross(U).normalized();
            utils::FixZeros(V);
            file << "ggbApplet.evalCommand(\"Vi" << intfIdx << " = UnitVector(Segment(ci" <<
                intfIdx << ", ci" << intfIdx << " + (" << V.x() << ", " << V.y() << ", " <<
                V.z() << ")))\");" << std::endl;
            file << "ggbApplet.setColor(\"Vi" << intfIdx << "\", 255, 255, 0);" << std::endl;
        }

        // Remove the indices from the vertices of the geometry
        geometry.RemoveVerticesAttribute(ATTRIB_INDEX);
    }

    // Close the file
    file.close();*/
}

void Workspace::WriteTessellationGeogebraJsFile(
    const std::string & filename, 
    int r, 
    int g, 
    int b, 
    double threshold) const
{
    // There must exist a tessellation
    assert(m_tessellation);

    // Initialize the string stream for storing the commands that generate the geometry of the 
    // tessellation
    std::stringstream ss;

    // Generate and write the JavaScript commands that generate the geometry of the tessellation in
    // GeoGebra
    m_tessellation->Geometry().WriteGeogebraJs(ss, "geomdom", r, g, b, 1, threshold);

    // Open the file, write its content and close it
    std::ofstream file(filename, std::ios::out);
    file << ss.str() << std::endl;
    file.close();
}

void Workspace::WriteTessellationObjFile(const std::string & filename, double threshold) const
{
    // There must exist a tessellation
    assert(m_tessellation);

    // Generate and write the OBJ file with the geometry of the tessellation
    m_tessellation->Geometry().WriteObjFile(filename, threshold);
}

void Workspace::WriteInterfacesFile(
    const std::vector<std::shared_ptr<VF>> & I, 
    const std::string filename) const
{
    // Open the file
    std::ofstream file(filename, std::ios::out);

    // Initialize the stringstream for writing the content of the interface polygons
    std::stringstream ss;

    // Traverse through the interface polygons and write their geometry
    for (auto it = I.begin(); it != I.end(); ++it) 
    {
        (*it)->Write(ss);
    }

    // Write the content in the file
    file << ss.str() << std::endl;

    // Close the file
    file.close();
}

void Workspace::WriteBlocksFile(const std::string filename) const 
{
    // Open the file
    std::ofstream file(filename, std::ios::out);

    // Initialize the stringstream for writing the content of the blocks
    std::stringstream ss;

    // Write the content of the blocks
    m_assembly->Write(ss);

    // Write the content in the file
    file << ss.str() << std::endl;

    // Close the file
    file.close();
}

void Workspace::WriteRelationsFile(
    const std::list<std::tuple<size_t, size_t>> & R, 
    const std::string filename) const 
{
    // Get the number of relations
    size_t nRelations = R.size();

    // Open the file
    std::ofstream file(filename, std::ios::out);

    // 
    size_t idx = 0;

    // Traverse through the relations and write them in the file
    for (auto it = R.begin(); it != R.end(); ++it)
    {
        // 
        file << idx << " P " << std::get<0>(*it) << " I " << std::get<1>(*it) << std::endl;

        // 
        idx += 1;
    }

    // Close the file
    file.close();
}
