#include <tiger/equilibriumanalysis.h>
#include <Eigen/Geometry>

#include <gurobi_c++.h>

#include <fstream>

//
// @param const char * s1
// @param const char * s2
//
std::string name(const char * s1, const char* s2)
{
    std::stringstream ss;
    ss << s1 << s2;
    return ss.str();
}

//
// @param const char * s
// @param size_t i
// @param const char* type
//
std::string name(const char * s, size_t i, const char* type)
{
    std::stringstream ss;
    ss << s << i << type;
    return ss.str();
}

//
// @param const char * s1
// @param size_t i1
// @param const char * s2
// @param size_t i2
// @param const char * type
//
std::string name(const char * s1, size_t i1, const char * s2, size_t i2, const char * type)
{
    std::stringstream ss;
    ss << s1 << i1 << s2 << i2 << type;
    return ss.str();
}

//
// Counts the number of vertices from the given interface polygons.
// @param const std::vector<std::shared_ptr<VF>> & I The reference to the 
// vector with the geometries of the interface polygons.
// @return size_t The number of vertices.
//
size_t countInterfacePolygonVertices(const std::vector<std::shared_ptr<VF>> & I) 
{
    size_t count = 0;

    // Traverse through the interface polygons and count their vertices
    for (auto it = I.begin(); it != I.end(); ++it) 
    {
        count += (*it)->countVertices();
    }

    return count;
}

bool EquilibriumAnalysis::Run(
    const std::vector<Eigen::Vector3d> & C, 
    const std::vector<std::vector<double>> & W, 
    const std::vector<std::shared_ptr<VF>> & I, 
    const std::list<std::tuple<size_t, size_t>> & BI, 
    double friction, 
    Result & result, 
    bool verbose, 
    bool files, 
    std::string fileprefix, 
    double cWeight, 
    double tWeight, 
    double uWeight, 
    double vWeight, 
    int logToConsole)
{
    // Reset the results object
    result.Reset();

    // Get the number of blocks
    size_t nBlocks = C.size();

    // Get the number of interfaces
    size_t nInterfaces = I.size();

    // Initialize the total number of interface polygon vertices
    size_t nVertices = countInterfacePolygonVertices(I);

    // Define the name of the output file
    std::string outputFilename = name(fileprefix.c_str(), ".out");

    if (verbose)
    {
        std::cout << std::endl;
        std::cout << "+---------------------------------+" << std::endl;
        std::cout << "| Equilibrium Analysis Parameters |" << std::endl;
        std::cout << "+---------------------------------+" << std::endl;
        std::cout << "#Blocks = " << nBlocks << std::endl;
        std::cout << "#Interfaces = " << nInterfaces << std::endl;
        std::cout << "#Vertices = " << nVertices << std::endl;
        std::cout << "#Variables = " << (nVertices * 4) << std::endl;
        std::cout << "#Contacts = " << BI.size() << std::endl;
        std::cout << "Friction Coeff. = " << friction << std::endl;
        std::cout << "Compression Weight = " << cWeight << std::endl;
        std::cout << "Tension Weight = " << tWeight << std::endl;
        std::cout << "U-tangential Weight = " << uWeight << std::endl;
        std::cout << "V-tangential Weight = " << vWeight << std::endl;
        std::cout << std::endl;
    }

    if (files)
    {
        // Open the output file
        std::ofstream file(outputFilename.c_str(), std::ofstream::out | std::ofstream::app);

        file << std::endl;
        file << "+---------------------------------+"           << std::endl;
        file << "| Equilibrium Analysis Parameters |"           << std::endl;
        file << "+---------------------------------+"           << std::endl;
        file << "#blocks = "                << nBlocks          << std::endl;
        file << "#Interfaces = "            << nInterfaces      << std::endl;
        file << "#Vertices = "              << nVertices        << std::endl;
        file << "#Variables = "             << (nVertices * 4)  << std::endl;
        file << "#Contacts = "              << BI.size()        << std::endl;
        file << "Friction Coeff. = "        << friction         << std::endl;
        file << "Compression Weight = "     << cWeight          << std::endl;
        file << "Tension Weight = "         << tWeight          << std::endl;
        file << "U-tangential Weight = "    << uWeight          << std::endl;
        file << "V-tangential Weight = "    << vWeight          << std::endl;
        file << std::endl;

        // Close the output file
        file.close();
    }

    // Define the quadratic program and try to run it
    try
    {
        // Initialize the Gurobi environment
        GRBEnv env;

        // Initialize the Gurobi model in the environment
        GRBModel model(env);
        //model.set(GRB_STR_PAR_LOGFILE, "");
        //model.set(GRB_INT_PAR_LOGTOCONSOLE, logToConsole);
        //model.set(GRB_INT_PAR_OUTPUTFLAG, 1);

        // 
        // +-------------------------------+
        // | Define the decision variables |
        // +-------------------------------+
        //

        // Determine the number of decision variables. There is a force vector 
        // per interface vertex. Each force is decomposed into 4 components: 
        // nc (compression, axial), nt (tension, axial), u (friction, 
        // tangential) and v (friction, tangential). The quadratic program 
        // needs to find the magnitudes of such force components. Then, there 
        // are nVertices * 4 decision variables
        size_t nVariables = nVertices * 4;

        // Initialize a vector for storing the decision variables
        std::vector<GRBVar> vars(nVariables);

        // Initialize the vector for storing the start index of the decision 
        // variables associated to the vertices of each interface polygon. We 
        // need to keep record of the start indices since not all interface 
        // polygons have the same number of vertices
        std::vector<size_t> startVertexIndex(nInterfaces);

        // Initialize the variable for tracking the indices of the vertices as 
        // they are placed in the decision variables vector
        size_t index = 0;

        // Traverse through the interface polygons and define the decision 
        // variables per vertex
        for (size_t k = 0; k < nInterfaces; k += 1)
        {
            // Store the start index of the decision variables associated to 
            // the vertices of the current interface polygon
            startVertexIndex[k] = index;

            // Traverse through the vertices of the interface polygon and 
            // define their respective decision variables
            for (size_t v = 0; v < I[k]->countVertices(); v += 1)
            {
                // The magnitude of the compression component of the force. It 
                // must be positive
                vars[index + 0] = model.addVar(
                    0.0, 
                    GRB_INFINITY, 
                    0, 
                    GRB_CONTINUOUS, 
                    name("k", k, "_v", v, "_c"));

                // The magnitude of the tension component of the force. It must
                // be positive
                vars[index + 1] = model.addVar(
                    0.0, 
                    GRB_INFINITY, 
                    0, 
                    GRB_CONTINUOUS, 
                    name("k", k, "_v", v, "_t"));

                // The magnitude of the u-tangential component of the force. It
                // is free
                vars[index + 2] = model.addVar(
                    -GRB_INFINITY, 
                    GRB_INFINITY, 
                    0, 
                    GRB_CONTINUOUS, 
                    name("k", k, "_v", v, "_u"));

                // The magnitude of the v-tangential component of the force. It
                // is free
                vars[index + 3] = model.addVar(
                    -GRB_INFINITY, 
                    GRB_INFINITY, 
                    0, 
                    GRB_CONTINUOUS, 
                    name("k", k, "_v", v, "_v"));

                // Update the variable index for the next vertex
                index += 4;
            }
        }


        //
        // +-------------------------------+
        // | Define the objective function |
        // +-------------------------------+
        //

        // Initialize the quadratic expression for the objective function
        GRBQuadExpr objective;

        // Traverse through the decision variables associated to the tension 
        // magnitudes of the forces
        for (size_t i = 0; i < nVariables; i += 4)
        {
            // Set the expression in the objective function related to the 
            // force vector of the current vertex. It is the sum of the squares
            // of the force components associated to the vertices of the 
            // interface polygons. Each term in the objective function is 
            // multiplied by the respective weight of the components
            objective +=
                (cWeight * vars[i + 0] * vars[i + 0]) +
                (tWeight * vars[i + 1] * vars[i + 1]) +
                (uWeight * vars[i + 2] * vars[i + 2]) +
                (vWeight * vars[i + 3] * vars[i + 3]);
        }

        // Set the objective function to the model. Indicate it is to be 
        // minimized
        model.setObjective(objective, GRB_MINIMIZE);


        // 
        // +-----------------------------------------------------------------+
        // | Define the equality constraints (net force and net momentum for |
        // | each block)                                                     |
        // +-----------------------------------------------------------------+
        //

        // Initialize the vector for storing the equality constraints. These 
        // are the net force and net momentum equations for each block. Their 
        // respective sum must be equal to 0. There are six equality 
        // constraints per block
        std::vector<GRBLinExpr> blockEqualityConstraints(nBlocks * 6);

        // Traverse through the tuples representing the interaction between 
        // blocks and interfaces and define the respective net force and net 
        // momentum constraints
        for (auto itBI = BI.begin(); itBI != BI.end(); ++itBI)
        {
            // Get the indices of the block and interface polygon respectively
            size_t blockIndex = std::get<0>(*itBI);
            size_t intfIndex = std::get<1>(*itBI);

            // Calculate the reference frame of the interface polygon. It is 
            // made of the normalized normal and tangential vectors of the 
            // interface
            Eigen::Vector3d N = I[intfIndex]->Normal(0, true);
            Eigen::Vector3d U = I[intfIndex]->direction(0, 0, true);
            Eigen::Vector3d V = N.cross(U).normalized();

            // Determine the direction of the reference vectors of the 
            // interface polygon with respect to the centroid of the block. If 
            // the centroid is facing the opposite direction of the face then 
            // invert the direction of the vectors
            if (I[intfIndex]->PointLocation(0, C[blockIndex]) < 0)
            {
                N *= -1.0;
                U *= -1.0;
                V *= -1.0;
            }

            // Get the start index for the equations associated to the block
            size_t blockEquationIndex = blockIndex * 6;

            // Get the start index of the decision variables associated to the 
            // vertices of the current interface polygon
            index = startVertexIndex[intfIndex];

            // Traverse through the vertices of the interface polygon and 
            // update the respective equality constraints 
            for (size_t v = 0; v < I[intfIndex]->countVertices(); v += 1)
            {
                // Update the net force constraint for the current block along 
                // the X axis
                blockEqualityConstraints[blockEquationIndex + 0] +=
                    (N.x() * vars[index + 0]) - (N.x() * vars[index + 1]) +
                    (U.x() * vars[index + 2]) + (V.x() * vars[index + 3]);

                // Update the net force constraint for the current block along 
                // the Y axis
                blockEqualityConstraints[blockEquationIndex + 1] +=
                    (N.y() * vars[index + 0]) - (N.y() * vars[index + 1]) +
                    (U.y() * vars[index + 2]) + (V.y() * vars[index + 3]);

                // Update the net force constraint for the current block along 
                // the Z axis
                blockEqualityConstraints[blockEquationIndex + 2] +=
                    (N.z() * vars[index + 0]) - (N.z() * vars[index + 1]) +
                    (U.z() * vars[index + 2]) + (V.z() * vars[index + 3]);

                // Calculate the relative position vector of the current vertex
                // of the interface polygon with respect to the centroid of the
                // current block
                Eigen::Vector3d Vij = C[blockIndex] - I[intfIndex]->Vertex(v);

                // Calculate the cross product between the reference vectors of
                // the current interface polygon and the relative position 
                // vector
                Eigen::Vector3d N_x_Vij = N.cross(Vij);
                Eigen::Vector3d U_x_Vij = U.cross(Vij);
                Eigen::Vector3d V_x_Vij = V.cross(Vij);

                // Update the net momentum constraint for the current block 
                // along the X axis
                blockEqualityConstraints[blockEquationIndex + 3] +=
                    (N_x_Vij.x() * vars[index + 0]) - (N_x_Vij.x() * vars[index + 1]) +
                    (U_x_Vij.x() * vars[index + 2]) + (V_x_Vij.x() * vars[index + 3]);

                // Update the net momentum constraint for the current block 
                // along the Y axis
                blockEqualityConstraints[blockEquationIndex + 4] +=
                    (N_x_Vij.y() * vars[index + 0]) - (N_x_Vij.y() * vars[index + 1]) +
                    (U_x_Vij.y() * vars[index + 2]) + (V_x_Vij.y() * vars[index + 3]);

                // Update the net momentum constraint for the current block 
                // along the Z axis
                blockEqualityConstraints[blockEquationIndex + 5] +=
                    (N_x_Vij.z() * vars[index + 0]) - (N_x_Vij.z() * vars[index + 1]) +
                    (U_x_Vij.z() * vars[index + 2]) + (V_x_Vij.z() * vars[index + 3]);

                // Update the variable index for the next vertex
                index += 4;
            }
        }

        // Traverse through the block equality constraints. Complete them and 
        // add them to the model
        for (size_t j = 0; j < nBlocks; j += 1)
        {
            // Get the start index for the equations associated to the block
            size_t blockEquationIndex = j * 6;

            // Complete the net force constraints of the block by adding its 
            // force loads. Then, add the net force constraints of the current 
            // block to the model
            model.addConstr(
                blockEqualityConstraints[blockEquationIndex + 0] + W[j][0],
                GRB_EQUAL,
                0.0,
                name("constr_block_", j, "_Fx"));

            model.addConstr(
                blockEqualityConstraints[blockEquationIndex + 1] + W[j][1],
                GRB_EQUAL,
                0.0,
                name("constr_block_", j, "_Fy"));

            model.addConstr(
                blockEqualityConstraints[blockEquationIndex + 2] + W[j][2],
                GRB_EQUAL,
                0.0,
                name("constr_block_", j, "_Fz"));

            // Complete the net momentum constraints of the block by adding its
            // momentum loads. Then, add the net momentum constraints of the 
            // current block to the model
            model.addConstr(
                blockEqualityConstraints[blockEquationIndex + 3] + W[j][3],
                GRB_EQUAL,
                0.0,
                name("constr_block_", j, "_Mx"));

            model.addConstr(
                blockEqualityConstraints[blockEquationIndex + 4] + W[j][4],
                GRB_EQUAL,
                0.0,
                name("constr_block_", j, "_My"));

            model.addConstr(
                blockEqualityConstraints[blockEquationIndex + 5] + W[j][5],
                GRB_EQUAL,
                0.0,
                name("constr_block_", j, "_Mz"));
        }


        //
        // +-------------------------------------------------------+
        // | Define the inequality constraints (friction pyramids) |
        // +-------------------------------------------------------+
        //

        // Calculate a conservative value for the friction cone. It helps to 
        // define a friction pyramid that fits within the actual friction cone
        double conservativeFriction = friction / sqrt(2.0);

        // Traverse through the interface polygons and define the friction 
        // coefficients for the tangential components associated to the forces 
        // at the interface vertices
        for (size_t k = 0; k < nInterfaces; k += 1)
        {
            // Get the start index of the decision variables associated to the 
            // vertices of the current interface polygon
            index = startVertexIndex[k];

            // Traverse through the vertices of the interface polygon and add 
            // its respective friction inequality constraints to the model
            for (size_t v = 0; v < I[k]->countVertices(); v += 1)
            {
                // Define the four faces of the friction pyramid. The 
                // compression component of the normal force and the 
                // conservative friction value defines the pyramid
                model.addConstr(
                    vars[index + 2] - (conservativeFriction * vars[index + 0]),
                    GRB_LESS_EQUAL,
                    0.0,
                    name("constr_interface_", k, "_v", v, "_friction_1"));

                model.addConstr(
                    -vars[index + 2] - (conservativeFriction * vars[index + 0]),
                    GRB_LESS_EQUAL,
                    0.0,
                    name("constr_interface_", k, "_v", v, "_friction_2"));

                model.addConstr(
                    vars[index + 3] - (conservativeFriction * vars[index + 0]),
                    GRB_LESS_EQUAL,
                    0.0,
                    name("constr_interface_", k, "_v", v, "_friction_3"));

                model.addConstr(
                    -vars[index + 3] - (conservativeFriction * vars[index + 0]),
                    GRB_LESS_EQUAL,
                    0.0,
                    name("constr_interface_", k, "_v", v, "_friction_4"));

                // Update the variable index for the next vertex
                index += 4;
            }
        }


        //
        // +-------------------------+
        // | Run the quadratic model |
        // +-------------------------+
        //

        // If indicated, write the files with the content of the model
        if (files)
        {
            // Write the files representing the model
            model.write(name(fileprefix.c_str(), ".mps"));
            model.write(name(fileprefix.c_str(), ".rew"));
            model.write(name(fileprefix.c_str(), ".lp"));
            model.write(name(fileprefix.c_str(), ".rlp"));
        }

        // Run the model (fingers crossed!)
        model.optimize();

        // If the model reached an optimal solution then show it; otherwise, 
        // say what happened
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        {
            // Indicate there is a solution to the Quadratic program. Then, 
            // store the optimal energy value
            result.isSolution = true;
            result.energy = model.get(GRB_DoubleAttr_ObjVal);

            // Initialize the index for traversing the decision variables
            index = 0;

            // Traverse through the interface polygons
            for (size_t k = 0; k < nInterfaces; k += 1)
            {

                // Traverse through the vertices of the interface polygon and 
                // define their respective decision variables
                for (size_t v = 0; v < I[k]->countVertices(); v += 1)
                {
                    // Define the key for storing the current force. The key is
                    // defined by the indices of the current interface polygon 
                    // and vertex
                    std::tuple<size_t, size_t> key = std::make_tuple(k, v);

                    // Insert a force into the forces map. Then, get its 
                    // reference
                    result.forces.insert(std::make_pair(key, Force()));
                    Force & force = result.forces.at(key);

                    // Initialize and define the force object for the current 
                    // vertex
                    force.Set(
                        k, 
                        v, 
                        vars[index + 0].get(GRB_DoubleAttr_X),
                        vars[index + 1].get(GRB_DoubleAttr_X),
                        vars[index + 2].get(GRB_DoubleAttr_X),
                        vars[index + 3].get(GRB_DoubleAttr_X));

                    // Update the minimum and maximum force magnitude values
                    result.maxCompression = std::max(result.maxCompression, force.compression);
                    result.maxTension     = std::max(result.maxTension,     force.tension);
                    result.maxUTangential = std::max(result.maxUTangential, force.uTangential);
                    result.maxVTangential = std::max(result.maxVTangential, force.vTangential);
                    result.minCompression = std::min(result.minCompression, force.compression);
                    result.minTension     = std::min(result.minTension,     force.tension);
                    result.minUTangential = std::min(result.minUTangential, force.uTangential);
                    result.minVTangential = std::min(result.minVTangential, force.vTangential);

                    // Update the force component energies
                    result.compressionEnergy += (force.compression * force.compression);
                    result.tensionEnergy     += (force.tension     * force.tension);
                    result.uTangentialEnergy += (force.uTangential * force.uTangential);
                    result.vTangentialEnergy += (force.vTangential * force.vTangential);
                    
                    // Update the variable index for the next vertex
                    index += 4;
                }
            }

            // Return true indicating everything went OK
            return true;
        }
        else
        {
            // Check the infeasibility of the model. These lines don't solve 
            // the issue but show what variable/constraint might be causing the
            // solution to be infeasible
            model.computeIIS();
            model.write(name(fileprefix.c_str(), ".ilp"));

            // 
            result.isSolution = false;

            // Return false since the model could not be solved
            return false;
        }
    }
    catch (GRBException e)
    {
        // Gurobi failed! Write the error code
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        return false;
    }
    catch (...)
    {
        // Something else failed! Sorry :(
        std::cout << "Exception during optimization" << std::endl;
        return false;
    }
}

EquilibriumAnalysis::Force::Force() : 
    compression(0.0), 
    interfaceIndex(0), 
    tension(0.0), 
    uTangential(0.0), 
    vTangential(0.0), 
    vertexIndex(0)
{
}

void EquilibriumAnalysis::Force::Normalize(double value)
{
    assert(value != 0.0);

    compression /= value;
    tension /= value;
    uTangential /= value;
    vTangential /= value;
}

void EquilibriumAnalysis::Force::Set(
    size_t intfIdx, 
    size_t vIdx, 
    double C, 
    double T, 
    double UT, 
    double VT)
{
    interfaceIndex = intfIdx;
    vertexIndex = vIdx;
    compression = C;
    tension = T;
    uTangential = UT;
    vTangential = VT;
}

void EquilibriumAnalysis::Force::Write() const
{
    std::cout << std::endl;
    std::cout << "Interface " << interfaceIndex << ", Vertex " << vertexIndex << std::endl;
    std::cout << "Compression  = " << compression << std::endl;
    std::cout << "Tension      = " << tension     << std::endl;
    std::cout << "U Tangential = " << uTangential << std::endl;
    std::cout << "V Tangential = " << vTangential << std::endl;
}

void EquilibriumAnalysis::Result::GetMinMaxForces(Force::TYPE type, double & min, double & max) const
{
    switch (type) 
    {
        case Force::TYPE::COMPRESSION: 
        {
            min = minCompression;
            max = maxCompression;
            break;
        }

        case Force::TYPE::TENSION:
        {
            min = minTension;
            max = maxTension;
            break;
        }

        case Force::TYPE::UTANGENTIAL:
        {
            min = minUTangential;
            max = maxUTangential;
            break;
        }

        case Force::TYPE::VTANGENTIAL:
        {
            min = minVTangential;
            max = maxVTangential;
            break;
        }
    }
}

void EquilibriumAnalysis::Result::GetCTMinMaxForces(double & min, double & max) const
{
    min = std::min(minCompression, minTension);
    max = std::max(maxCompression, maxTension);
}

void EquilibriumAnalysis::Result::Normalize(double value)
{
    assert(value != 0.0);

    for (auto it = forces.begin(); it != forces.end(); ++it) 
    {
        it->second.Normalize(value);
    }

    UpdateMinMaxValues();
}

EquilibriumAnalysis::Result::Result() :
    compressionEnergy(0.0),
    energy(0.0),
    forces(), 
    isSolution(false), 
    maxCompression(0.0), 
    maxTension(0.0), 
    maxUTangential(0.0), 
    maxVTangential(0.0), 
    minCompression(0.0), 
    minTension(0.0), 
    minUTangential(0.0), 
    minVTangential(0.0), 
    tensionEnergy(0.0), 
    uTangentialEnergy(0.0), 
    vTangentialEnergy(0.0)
{
}

void EquilibriumAnalysis::Result::Reset()
{
    // Indicate there is no solution
    isSolution = false;

    // Clear the forces map
    forces.clear();

    // Reset the energy and component values
    compressionEnergy = 0;
    energy            = 0;
    maxCompression    = 0;
    maxTension        = 0;
    maxUTangential    = 0;
    maxVTangential    = 0;
    minCompression    = 0;
    minTension        = 0;
    minUTangential    = 0;
    minVTangential    = 0;
    tensionEnergy     = 0;
    uTangentialEnergy = 0;
    vTangentialEnergy = 0;
}

void EquilibriumAnalysis::Result::UpdateMinMaxValues()
{
    auto it = forces.begin();
    assert(it != forces.end());

    minCompression = it->second.compression;
    maxCompression = it->second.compression;
    minTension = it->second.tension;
    maxTension = it->second.tension;
    minUTangential = it->second.uTangential;
    maxUTangential = it->second.uTangential;
    minVTangential = it->second.vTangential;
    maxVTangential = it->second.vTangential;

    ++it;

    for (it; it != forces.end(); ++it) 
    {
        if (it->second.compression < minCompression) 
        {
            minCompression = it->second.compression;
        }

        if (it->second.compression > maxCompression) 
        {
            maxCompression = it->second.compression;
        }

        if (it->second.tension < minTension) 
        {
            minTension = it->second.tension;
        }

        if (it->second.tension > maxTension) 
        {
            maxTension = it->second.tension;
        }

        if (it->second.uTangential < minUTangential) 
        {
            minUTangential = it->second.uTangential;
        }

        if (it->second.uTangential > maxUTangential) 
        {
            maxUTangential = it->second.uTangential;
        }

        if (it->second.vTangential < minVTangential) 
        {
            minVTangential = it->second.vTangential;
        }

        if (it->second.vTangential > maxVTangential) 
        {
            maxVTangential = it->second.vTangential;
        }
    }
}

void EquilibriumAnalysis::Result::Write() const
{
    WriteForceComponents();
    WriteForceRanges();
}

void EquilibriumAnalysis::Result::WriteForceComponents() const
{
    std::cout << "+------------------+" << std::endl;
    std::cout << "| Force Components |" << std::endl;
    std::cout << "+------------------+" << std::endl;

    //size_t k, v;

    // Traverse through the forces and write their content
    for (auto it = forces.begin(); it != forces.end(); ++it)
    {
        // Get the indices of the interface polygon and vertex respectively
        //k = std::get<0>(it->first);
        //v = std::get<1>(it->first);

        // Get the reference to the current force object
        //const Force & force = it->second;

        // Write the force components
        it->second.Write();
    }
}

void EquilibriumAnalysis::Result::WriteForceRanges() const
{
    std::cout << "+-----------------------+" << std::endl;
    std::cout << "| Energy & Force Bounds |" << std::endl;
    std::cout << "+-----------------------+" << std::endl;
    std::cout << "Energy = " << energy << std::endl;
    std::cout << "Compressions = [" << minCompression << ", " << maxCompression << "]" << std::endl;
    std::cout << "Tensions = [" << minTension << ", " << maxTension << "]" << std::endl;
    std::cout << "uTangentials = [" << minUTangential << ", " << maxUTangential << "]" << std::endl;
    std::cout << "vTangentials = [" << minVTangential << ", " << maxVTangential << "]" << std::endl;
}
