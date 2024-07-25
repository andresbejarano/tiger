#include <tiger/utils.h>
#include <Eigen/Geometry>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>

Eigen::Vector3d utils::axisAngleRotation(
    const Eigen::Vector3d & V,
    const Eigen::Vector3d & K,
    double angle)
{
    double sin_a = sin(angle), cos_a = cos(angle), _cos_a = 1.0 - cos_a;

    Eigen::Vector3d 
        nK = K.normalized(), 
        R = (V * cos_a) + (nK.cross(V) * sin_a) + (nK * (nK.dot(V)) * _cos_a);

    return R;
}

Eigen::Vector3d utils::centroid(
    const std::list<Eigen::Vector3d> & points, 
    bool fixZeros, 
    double threshold)
{
    if (points.size() == 0) 
    {
        return Eigen::Vector3d::Zero();
    }

    Eigen::Vector3d C = Eigen::Vector3d::Zero();

    for (auto it = points.begin(); it != points.end(); ++it) 
    {
        C += *it;
    }

    return C / (double)points.size();
}

Eigen::Vector3d utils::centroid(
    const std::vector<Eigen::Vector3d> & points, 
    bool fixZeros, 
    double threshold)
{
    if (points.size() == 0)
    {
        return Eigen::Vector3d::Zero();
    }

    Eigen::Vector3d C = Eigen::Vector3d::Zero();

    for (auto it = points.begin(); it != points.end(); ++it)
    {
        C += *it;
    }

    return C / (double)points.size();
}

double utils::det(
    double a, double b, double c,
    double d, double e, double f,
    double g, double h, double i)
{
    return (a * ((e * i) - (f * h))) - (b * ((d * i) - (f * g))) + (c * ((d * h) - (e * g)));
}

bool utils::findPath(
    Eigen::MatrixXd & adj,
    std::vector<bool> & visited,
    size_t i,
    std::vector<size_t> & path)
{
    // If point i has been visited then the loop is closed. Otherwise, keep 
    // visiting
    if (visited[i])
    {
        // Traverse through the visited vector and search for a non visited 
        // index
        for (auto it = visited.begin(); it != visited.end(); ++it)
        {
            // If the current index is not visited then return false
            if (!(*it))
            {
                return false;
            }
        }

        // Return true since all indices are visited
        return true;
    }
    else
    {
        // Set the current point as visited
        visited[i] = true;

        // Get the number of points
        size_t nPoints = visited.size(), j = 0;

        // Traverse through the possible directions from i
        for (j = 0; j < nPoints; j += 1)
        {
            // Avoid invalid directions
            if (adj(i, j) == -1.0)
            {
                continue;
            }

            // Check if we can find a path by moving from i to j
            if (findPath(adj, visited, j, path))
            {
                path.insert(path.begin(), i);
                return true;
            }
        }

        // No path was possible from i
        visited[i] = false;
        return false;
    }
}

void utils::fixZero(double & val, double threshold)
{
    if (abs(val) <= threshold)
    {
        val = 0.0;
    }
}

void utils::fixZeros(Eigen::Vector3d & V, double threshold)
{
    if (abs(V.x()) <= threshold)
    {
        V(0) = 0.0;
    }

    if (abs(V.y()) <= threshold)
    {
        V(1) = 0.0;
    }

    if (abs(V.z()) <= threshold)
    {
        V(2) = 0.0;
    }
}

void utils::fixZeros(std::vector<Eigen::Vector3d> & V, double threshold)
{
    // Traverse through the 3D vectors and fiz their zero values
    for (auto it = V.begin(); it != V.end(); ++it)
    {
        fixZeros(*it, threshold);
    }
}

bool utils::getUniquePoints(
    std::vector<Eigen::Vector3d> & points,
    std::vector<Eigen::Vector3d> & unique,
    double threshold)
{
    // Clear the vector where unique points will be stored
    unique.clear();

    bool isUnique = false;

    // Traverse through the array of points
    for (auto i = points.begin(); i != points.end(); ++i)
    {
        // Get the reference to the current point
        Eigen::Vector3d & V = *i;

        // Assume the current point Vi is unique
        isUnique = true;

        // Traverse through the array of unique points
        for (auto j = unique.begin(); j != unique.end(); ++j)
        {
            // Get the reference to the current unique point
            Eigen::Vector3d & U = *j;

            // If the coordinate values from Vi and Uj are equal then indicate 
            // that Vi is no longer unique
            if (abs(V.x() - U.x()) <= threshold &&
                abs(V.y() - U.y()) <= threshold &&
                abs(V.z() - U.z()) <= threshold)
            {
                isUnique = false;
            }
        }

        // If Vi is unique then add its clone to the array of unique points
        if (isUnique)
        {
            unique.emplace_back(V);
        }
    }

    // Indicate whether there is at least one unique point
    return unique.size() > 0;
}

bool utils::isPointInTriangle(
    const Eigen::Vector3d & P,
    const Eigen::Vector3d & v0,
    const Eigen::Vector3d & v1,
    const Eigen::Vector3d & v2,
    double threshold)
{
    return sameSide(P, v0, v1, v2, threshold) &&
        sameSide(P, v1, v0, v2, threshold) &&
        sameSide(P, v2, v0, v1, threshold);
}

void utils::print(const char * s)
{
    std::cout << s << std::endl;
}

void utils::print(const char * s, double v)
{
    std::cout << s << " = " << v << std::endl;
}

void utils::print(const char * s, const Eigen::Vector3d & V)
{
    std::cout << s << " = (" << V.x() << ", " << V.y() << ", " << V.z() << ")" << std::endl;
}

void utils::print(const char * s, size_t i, const Eigen::Vector3d & V)
{
    std::cout << s << i << " = (" << V.x() << ", " << V.y() << ", " << V.z() << ")" << std::endl;
}

void utils::print(const char * s1, size_t i1, const char * s2, size_t i2)
{
    std::cout << s1 << i1 << s2 << i2 << std::endl;
}

void utils::print(const char * s1, size_t i1, const char * s2, size_t i2, const Eigen::Vector3d & V)
{
    std::cout << s1 << i1 << s2 << i2 << " = (" << V.x() << ", " << V.y() << ", " << V.z() << ")" << std::endl;
}

void utils::print(const std::vector<Eigen::Vector3d> & V)
{
    for (auto it = V.begin(); it != V.end(); ++it)
    {
        std::cout << "[" << *it << "]" << std::endl;
    }
}

bool utils::sameSide(
    const Eigen::Vector3d & P,
    const Eigen::Vector3d & Q,
    const Eigen::Vector3d & A,
    const Eigen::Vector3d & B,
    double threshold)
{
    // 
    Eigen::Vector3d 
        dir = (B - A).normalized(), 
        cp1 = dir.cross(P - A).normalized(), 
        cp2 = dir.cross(Q - A).normalized();

    // 
    double test = cp1.dot(cp2);
    fixZero(test, threshold);

    // If the dot product is greater or equal to zero then both points P and Q 
    // are at the same side of vector AB
    return (test >= 0.0);
}

std::vector<std::string> utils::split(const std::string& text, char sep) 
{
    assert(text.size() > 0);

    std::vector<std::string> tokens;

    size_t n = text.size();
    size_t i = 0;

    for (size_t j = 0; j < n; j += 1) 
    {
        if (text.at(j) == sep) 
        {
            tokens.push_back(text.substr(i, j - i));
            i = j + 1;
        }
    }

    tokens.push_back(text.substr(i, n - i));

    // Verify there was at least one token
    assert(tokens.size() > 0);

    return tokens;
}

std::string utils::toString(const Eigen::Vector3d& V, double threshold) 
{
    std::stringstream ss;

    ss << "(" 
        << (V.x() <= threshold ? 0.0 : V.x()) << ", " 
        << (V.y() <= threshold ? 0.0 : V.y()) << ", " 
        << (V.z() <= threshold ? 0.0 : V.z()) << ")";

    return ss.str();
}

double utils::vectorFunction(std::vector<double> & values, const std::string & function)
{
    if (function == "AVG") 
    {
        return std::accumulate(values.begin(), values.end(), 0.0) / (double)values.size();

    }
    else if (function == "MIN") 
    {
        return *std::min_element(values.begin(), values.end());
    }

    // No valid function name given
    assert(false);
    return 0.0;
}

void utils::write(const Eigen::Vector3d & V, bool newline)
{
    std::cout << "(" << V.x() << ", " << V.y() << ", " << V.z() << ")";

    if (newline) 
    {
        std::cout << std::endl;
    }
}

void utils::write(const std::array<double, 6>& V) 
{
    std::cout << "[ " << V[0];

    for (size_t i = 1; i < 6; i += 1) 
    {
        std::cout << ",\t" << V[i];
    }

    std::cout << " ]" << std::endl;
}

void utils::write(const std::list<Eigen::Vector3d> & points)
{
    for (auto it = points.begin(); it != points.end(); ++it) 
    {
        const Eigen::Vector3d & v = *it;

        std::cout << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")" << std::endl;
    }
}

void utils::write(const std::list<std::vector<size_t>> & L)
{
    for (auto itL = L.begin(); itL != L.end(); ++itL) 
    {
        const std::vector<size_t> & V = *itL;

        size_t nV = V.size();

        std::cout << V[0];

        for (size_t i = 1; i < nV; i += 1) 
        {
            std::cout << " " << V[i];
        }

        std::cout << std::endl;
    }
}

void utils::write(const std::vector<Eigen::Vector3d> & V)
{
    for (auto it = V.begin(); it != V.end(); ++it)
    {
        const Eigen::Vector3d & v = *it;

        std::cout << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")" << std::endl;
    }
}

void utils::write(const std::vector<double> & V) 
{
    size_t nV = V.size();

    std::cout << V[0];

    for (size_t i = 1; i < nV; i += 1)
    {
        std::cout << " " << V[i];
    }

    std::cout << std::endl;
}

void utils::write(const std::vector<size_t>& V) 
{
    size_t nV = V.size();

    std::cout << V[0];

    for (size_t i = 1; i < nV; i += 1)
    {
        std::cout << " " << V[i];
    }

    std::cout << std::endl;
}

void utils::write(const std::vector<std::string>& V) 
{
    size_t nV = V.size();

    std::cout << V[0];

    for (size_t i = 1; i < nV; i += 1) 
    {
        std::cout << " " << V[i];
    }

    std::cout << std::endl;
}
