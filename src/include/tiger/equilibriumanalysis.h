#ifndef _EQUILIBRIUM_ANALYSIS_H_
#define _EQUILIBRIUM_ANALYSIS_H_

#pragma once

#include <tiger/ds/vf.h>
#include <map>

//
// The namespace for storing the classes and methods related to the Equilibrium
// Analysis.
//
namespace EquilibriumAnalysis
{

    //
    // The class for storing the information of a force after running the 
    // equilibrium analysis. Each force object stores its respective component 
    // magnitudes and the indices of the interface polygon and respective 
    // vertex.
    //
    class Force 
    {

    public:

        enum TYPE { COMPRESSION = 1, TENSION, UTANGENTIAL, VTANGENTIAL };

        // The magnitude of the compresion component
        double compression;

        // The index of the interface polygon
        size_t interfaceIndex;

        // The magnitude of the tension component
        double tension;

        // The magnitude of the U tangential component
        double uTangential;

        // The magnitude of the V tangential component
        double vTangential;

        // The index of the vertex in the interface polygon
        size_t vertexIndex;

        //
        // Constructor of the class.
        //
        Force();

        //
        // @param double value The normalizing value.
        //
        void Normalize(double value);

        //
        // Sets the content of the force.
        // @param size_t intfIdx The index of the interface polygon.
        // @param size_t vIdx The index of the vertex in the interface polygon.
        // @param double C The magnitude of the compression component.
        // @param double T The magnitude of the tension component.
        // @param double UT The magnitude of the U tangential component.
        // @param double VT The magnitude of the V tangential component.
        //
        void Set(size_t intfIdx, size_t vIdx, double C, double T, double UT, double VT);

        //
        // Writes the content of the force.
        //
        void Write() const;
    };

    //
    // The class representing the results of an equilibrium analysis. Results 
    // store the magnitudes of the forces along with some metadata for 
    // visualization purposes.
    //
    class Result 
    {

    public:

        // The compression component of the optimal energy value. It is the sum
        // of squared compression components from all forces
        double compressionEnergy;

        // The optimal energy value
        double energy;

        // The map for storing the forces at the vertices of the interface 
        // polygons. The key of each force is defined by the indices of the 
        // respective interface polygon and vertex
        std::map<std::tuple<size_t, size_t>, Force> forces;

        // Indicates if the equilibrium analysis found a solution of the 
        // quadratic program
        bool isSolution;

        // The maximum compression magnitude among all forces
        double maxCompression;

        // The maximum tension magnitude among all forces
        double maxTension;

        // The maximum U tangential magnitude among all forces
        double maxUTangential;

        // The maximum V tangential magnitude among all forces
        double maxVTangential;

        // The minimum compression magnitude among all forces
        double minCompression;

        // The minimum tension magnitude ammong all forces
        double minTension;

        // The minimum U tangential magnitude among all forces
        double minUTangential;

        // The minimum V tangential magnitude among all forces
        double minVTangential;

        // The tension component of the optimal energy value. It is the sum of 
        // squared tension components from all forces
        double tensionEnergy;

        // The U tangential component of the optimal energy value. It is the 
        // sum of squared U tangential components from all forces
        double uTangentialEnergy;

        // The V tangential component of the optimal energy value. It is the 
        // sum of squared V tangential components from all forces
        double vTangentialEnergy;

    public:

        //
        // Constructor of the class.
        //
        Result();

        //
        // @param Force::TYPE type
        // @param double & min
        // @param double & max
        //
        void GetMinMaxForces(Force::TYPE type, double & min, double & max) const;

        //
        // @param double & min
        // @param double & max
        //
        void GetCTMinMaxForces(double & min, double & max) const;

        //
        //
        void Normalize(double value);

        //
        // Resets the content of the results.
        //
        void Reset();

        //
        //
        void UpdateMinMaxValues();

        //
        //
        void Write() const;

        //
        //
        void WriteForceComponents() const;

        //
        //
        void WriteForceRanges() const;

    };

    //
    // Runs the equilibrium analysis. 
    // @param const std::vector<Eigen::Vector3d> & C The reference to the 
    // vector with the centroid of the blocks.
    // @param const std::vector<std::vector<double>> & W The reference to the 
    // vector with the loads that apply to the blocks.
    // @param const std::list<dcel::DCEL> & I The reference to the list with 
    // the geometry of the interface polygons between blocks.
    // @param const std::vector<std::tuple<size_t, size_t>> & BI The reference 
    // to the vector with the tuples indicating how blocks and interface 
    // polygons interact. First index references the block, second index 
    // references the interface polygon
    // @param double friction The friction coefficient.
    // @param double & energy The optimal value.
    // @param bool verbose Indicates whether to print results as they are 
    // calculated.
    // @param bool files Indicates whether to write the files with the content 
    // of the model.
    // @param std::string fileprefix
    // @param double cWeight
    // @param double tWeight
    // @param double uWeight
    // @param double vWeight
    // @return bool
    //
    bool Run(
        const std::vector<Eigen::Vector3d> & C,
        const std::vector<std::vector<double>> & W,
        const std::vector<std::shared_ptr<VF>> & I,
        const std::list<std::tuple<size_t, size_t>> & BI,
        double friction,
        Result & result, 
        bool verbose = false,
        bool files = false,
        std::string fileprefix = "model",
        double cWeight = 1.0e5,
        double tWeight = 1.0e5,
        double uWeight = 1.0e3,
        double vWeight = 1.0e3, 
        int logToConsole = 1);

}

#endif
