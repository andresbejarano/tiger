#include "geometries.h"
#include "utils.h"
#include <Eigen/Geometry>

const double TWO_PI = 2.0 * acos(-1.0);

/*VF geometries::BendedSquare(
    double radius, 
    double length, 
    size_t rs, 
    size_t ls, 
    int direction, 
    const Eigen::Vector3d & U, 
    const Eigen::Vector3d & V)
{
    VF vf((rs + 1) * (ls + 1), rs * ls);

    double sign = direction >= 0 ? 1.0 : -1.0;

    double halfLength = length / 2.0;
    double lengthStep = length / (double)ls;
    double arcAngle = halfLength / radius;
    double angleStep = (2.0 * arcAngle) / (double)rs;

    Eigen::Vector3d W = U.cross(V).normalized(), C = sign * radius * W;

    double angle = -arcAngle, sin_a = 0, cos_a = 0;

    for (size_t r = 0; r < rs; r += 1) 
    {
        sin_a = sin(angle);
        cos_a = cos(angle);



        for (size_t l = 0; l < ls; l += 1) 
        {
        }

        angle += angleStep;
    }
    
    return vf;
}*/

VF geometries::BendedSquare(double radius, double length, size_t rs, size_t ls, int bending)
{
    // Initialize the vertex coordinates and vertex indices structure
    VF vf((rs + 1) * (ls + 1), rs * ls);

    // Calculate the required constants
    double halfLength = length / 2.0;
    double lengthStep = length / double(ls);
    double arcAngle = halfLength / radius;
    double arcStep = (2.0 * arcAngle) / double(rs);

    // Calculate the center of the circumference and the center point (at the origin)
    Eigen::Vector3d H(0.0, double(bending) * radius, 0.0), C = Eigen::Vector3d::Zero();

    // Define the vertices
    for (size_t i = 0; i <= ls; i += 1)
    {
        for (size_t j = 0; j <= rs; j += 1)
        {
            // Calculate the sine and cosine of the current angle
            double sin_a = sin(-arcAngle + (double(j) * arcStep));
            double cos_a = cos(-arcAngle + (double(j) * arcStep));

            // Calculate the coordinates of the current vertex
            double x = H.x() + ((C.x() - H.x()) * cos_a) - ((C.y() - H.y()) * sin_a);
            double y = H.y() + ((C.x() - H.x()) * sin_a) + ((C.y() - H.y()) * cos_a);
            double z = halfLength - (double(i) * lengthStep);

            // Generate the vertex objext and push it into V
            vf.addVertex(x, y, z);
        }
    }

    // Define the vertex indices for the faces
    for (size_t i = 0; i < ls; i += 1)
    {
        for (size_t j = 0; j < rs; j += 1)
        {
            // Calculate the vertex indices that define the current face
            size_t v0 = j + (i * (rs + 1));
            size_t v1 = v0 + 1;
            size_t v2 = (j + 1) + ((i + 1) * (rs + 1));
            size_t v3 = v2 - 1;

            vf.addFace(v0, v1, v2, v3);
        }
    }

    // Return the object with the vertices and faces
    return vf;
}

VF geometries::Cone(
    const Eigen::Vector3d & A,
    const Eigen::Vector3d & B,
    double radius,
    size_t ls,
    size_t rs,
    bool open) 
{
    // Initialize the vertex coordinates and vertex indices structure
    /*VF vf;

    // Get a random point
    double xRand = floor(Random01() * 10.0) + 1.0;
    double yRand = floor(Random01() * 10.0) + 1.0;
    double zRand = floor(Random01() * 10.0) + 1.0;
    Eigen::Vector3d P(xRand, yRand, zRand);

    // Define the axis vector of the cone
    Eigen::Vector3d AB = B - A;

    // Calculate a temproal vector from A to the random point
    Eigen::Vector3d temp = P - A;

    // Calculate the R basis vector
    Eigen::Vector3d R = AB.cross(temp);
    R.normalize();

    // Calculate the S basis vector
    Eigen::Vector3d S = AB.cross(R);
    S.normalize();

    // Calculate the angle of the radial segments
    double as = 2.0 * PI / (double)rs;

    // Insert point A into the vertices array. This point correspond to the singular point of the 
    // cone
    vf.addVertex(A);

    // Traverse through the axis vector of the cone
    for (size_t t = 1; t <= ls; t += 1)
    {
        // Calculate the parameter value for the current cone level
        double u = (double)t / (double)ls;

        // Calculate the value of the first term of the parametric equation
        Eigen::Vector3d term1 = A + (u * AB);

        // Traverse through the radial segments of the current cone level
        for (size_t a = 0; a < rs; a += 1)
        {
            // Calculate the angle of the current radial segment. Then, get its sine and cosine
            double v = (double)a * (double)as;
            double sin_v = sin(v);
            double cos_v = cos(v);

            // Calculate the second term of the parametric equation
            Eigen::Vector3d term2 = (u * radius * sin_v) * R;

            // Calculate the third term of the parametric equation
            Eigen::Vector3d term3 = (u * radius * cos_v) * S;

            // Calculate the current point of the cone and add it to the vertices array
            Eigen::Vector3d Q = term1 + term2 + term3;
            vf.addVertex(Q);
        }
    }

    // Traverse through the indices of the first level of the cone and generate the faces that are 
    // incident to the singularity point
    for (size_t a = 1; a <= rs; a += 1)
    {
        // Determine the next vertex index in the current face
        size_t next = (a < rs) ? a + 1 : 1;

        // Generate the vertex indexes array and push it into F
        vf.F.emplace_back();
        std::vector<size_t> & indices = vf.F[vf.F.size() - 1];
        indices.push_back(0);
        indices.push_back(a);
        indices.push_back(next);
    }

    // Define the vertex indexes for the faces
    for (size_t t = 1; t < ls; t += 1)
    {
        // 
        size_t index = ((t - 1) * rs) + 1;

        // 
        size_t next = (t * rs) + 1;

        // 
        for (size_t a = 0; a < rs; a += 1)
        {
            // 
            size_t v0 = index + a;
            size_t v1 = next + a;
            size_t v2;
            size_t v3;

            // Calculate the vertex indexes that define the current face. The last face of each ring
            // uses the first and last indexes in the current ring and the next one
            if (a < rs - 1)
            {
                v2 = v1 + 1;
                v3 = v0 + 1;
            }
            else
            {
                v2 = next;
                v3 = index;
            }

            // Generate the vertex indexes array and push it into F
            vf.F.emplace_back();
            std::vector<size_t> & indices = vf.F[vf.F.size() - 1];
            indices.push_back(v0);
            indices.push_back(v1);
            indices.push_back(v2);
            indices.push_back(v3);
        }
    }

    // 
    if (!open)
    {
        // 
        vf.F.emplace_back();
        std::vector<size_t> & indices = vf.F[vf.F.size() - 1];

        // 
        for (size_t a = 0; a < rs; a += 1)
        {
            indices.push_back(vf.V.size() - a - 1);
        }
    }

    // Return the vertex coordinates and vertex indices
    return vf;*/
    return Square(1);
}

VF geometries::Cyclide(double a, double b, double c, double d, size_t Rs, size_t rs)
{
    //const Eigen::Vector3d Z = X.cross(Y);

    size_t n = Rs * rs;

    // Initialize the vertex coordinates and vertex indices structure
    VF vf(n, n);

    // Calculate the angle steps
    double uStep = TWO_PI / (double)Rs;
    double vStep = TWO_PI / (double)rs;

    double u, sin_u, cos_u, v, sin_v, cos_v, denom, x, y, z;

    // Define the vertices of the torus
    for (size_t l = 0; l < Rs; l += 1)
    {
        u = (double)l * uStep;
        sin_u = sin(u);
        cos_u = cos(u);

        //Eigen::Vector3d T1 = C + (X * (R * cos_u)) + (Y * (R * sin_u));
        //Eigen::Vector3d T2 = (X * -sin_u) + (Y * cos_u);
        //Eigen::Vector3d T3 = T2.cross(Z).normalized();

        for (size_t m = 0; m < rs; m += 1)
        {
            v = (double)m * vStep;
            sin_v = sin(v);
            cos_v = cos(v);

            //Eigen::Vector3d Q = T1 + (T3 * (r * cos_v)) + (Z * (r * sin_v));

            denom = a - (c * cos_u * cos_v);

            x = ((d * (c - (a * cos_u * cos_v))) + (b * b * cos_u)) / denom;
            y = (b * sin_u * (a - (d * cos_v))) / denom;
            z = (b * sin_v * ((c * cos_u) - d)) / denom;

            //Eigen::Vector3d Q(x, y, z);
            vf.addVertex(x, y, z);
        }
    }

    // Declare the variables for the vertex indices for the current face
    size_t v0, v1, v2, v3;

    // Define the vertex indices for the faces
    for (size_t l = 0; l < Rs; l += 1)
    {
        for (size_t m = 0; m < rs; m += 1)
        {
            // Determine the vertex indices of the current face according to the current vertex 
            // indices values
            if (l < Rs - 1)
            {
                if (m < rs - 1)
                {
                    // This is the common case: All squares not at the end of the ring or at the 
                    // last ring. This case happens mostly of the times
                    v0 = (l * rs) + m;
                    v1 = ((l + 1) * rs) + m;
                    v2 = v1 + 1;
                    v3 = v0 + 1;
                }
                else
                {
                    // This is the last square in a minor ring (not the last one). This case happens
                    // once per minor ring
                    v0 = (l * rs) + m;
                    v1 = ((l + 1) * rs) + m;
                    v2 = v0 + 1;
                    v3 = v0 - rs + 1;
                }
            }
            else
            {
                if (m < rs - 1)
                {
                    // This is a square (not the last) in the last minor ring. It uses the first 
                    // vertices generated for the torus. This case happens mostly in the last minor 
                    // ring
                    v0 = (l * rs) + m;
                    v1 = m;
                    v2 = v1 + 1;
                    v3 = v0 + 1;
                }
                else
                {
                    // This is the very last square in the last minor ring. This case happens once
                    v0 = (l * rs) + m;
                    v1 = m;
                    v2 = 0;
                    v3 = v0 - rs + 1;
                }
            }

            vf.addFace(v0, v1, v2, v3);
        }
    }

    return vf;
}

VF geometries::Cylinder(
    const Eigen::Vector3d & C, 
    const Eigen::Vector3d & K, 
    double radius, 
    double length, 
    size_t rs, 
    size_t ls, 
    bool open)
{
    // Initialize the vertex coordinates and vertex indices structure
    VF vf(rs * (ls + 1), rs * ls);

    // Calculate the required constants
    double halfLength = length / 2.0;
    double angle = 2.0 * utils::PI / double(rs);

    // 
    Eigen::Vector3d temp = K.normalized() * halfLength;

    // Define the end points of the cylinder
    Eigen::Vector3d A = C - temp;
    //Eigen::Vector3d B = C + temp;

    // Get the vector AB
    //Eigen::Vector3d AB = B - A;
    Eigen::Vector3d AB = 2.0 * temp;

    // Calculate a point somewhere but along the AB line segment (TO DO in next versions)
    double xRand = floor(Random01() * 10.0) + 1.0;
    double yRand = floor(Random01() * 10.0) + 1.0;
    double zRand = floor(Random01() * 10.0) + 1.0;
    Eigen::Vector3d P(xRand, yRand, zRand);

    // Calculate the basis vector R
    temp = P - C;
    Eigen::Vector3d R = K.cross(temp).normalized();

    // Calculate the basis vector S
    Eigen::Vector3d S = K.cross(R).normalized();

    size_t l = 0, r = 0;

    double t = 0.0, a = 0.0, sin_a = 0.0, cos_a = 0.0;

    Eigen::Vector3d 
        term1 = Eigen::Vector3d::Zero(), 
        r_R_sin_a = Eigen::Vector3d::Zero(), 
        r_R_cos_a = Eigen::Vector3d::Zero();

    // Define the vertices of the cylinder
    for (l = 0; l <= ls; l += 1)
    {
        // Calculate the first term of the parametric equation (it is the same for all vertices at 
        // this length step)
        t = double(l) / double(ls);
        term1 << (AB * t) + A;

        for (r = 0; r < rs; r += 1)
        {
            // Calculate angle a and its sine and cosine values
            a = double(r) * angle;
            sin_a = sin(a);
            cos_a = cos(a);

            // Calculate the second and third term of the parametric equation
            r_R_cos_a << R * (radius * cos_a);
            r_R_sin_a << S * (radius * sin_a);

            // Calculate the vertex and push it into V
            vf.addVertex(term1 + r_R_cos_a + r_R_sin_a);
        }
    }

    size_t v0 = 0, v1 = 0, v2 = 0, v3 = 0;

    // Define the vertex indices for the faces
    for (l = 0; l < ls; l += 1)
    {
        for (r = 0; r < rs; r += 1)
        {
            // Calculate the vertex indices that define the current face. The last face of each 
            // ring uses the first and last indices in the current ring and the next one
            if (r == rs - 1)
            {
                v0 = (l * rs) + r;
                v1 = v0 - (rs - 1);
                v2 = v0 + 1;
                v3 = v2 + (rs - 1);
            }
            else
            {
                v0 = (l * rs) + r;
                v1 = v0 + 1;
                v2 = ((l + 1) * rs) + (r + 1);
                v3 = v2 - 1;
            }

            vf.addFace(v0, v1, v2, v3);
        }
    }

    return vf;
}

VF geometries::Cylinder(
    const Eigen::Vector3d & A,
    const Eigen::Vector3d & B,
    double radius,
    size_t rs,
    size_t ls, 
    bool open) 
{
    // Calculate the center point between the end points
    Eigen::Vector3d C = (A + B) / 2.0;

    // Calculate the direction vector between the end points
    Eigen::Vector3d V = B - A;

    // Get the magnitude of the direction vector
    double length = V.norm();

    // Return the vertex coordinates and vertex indices of a cylinder centered at C and along 
    // vector V
    return Cylinder(C, V, radius, length, rs, ls, open);
}

VF geometries::Dodecahedron(double radius)
{
    double c = (sqrt(5.0) - 1.0) / 2.0;
    double c1 = 1.0 + c;
    double c2 = 1.0 - c * c;
    
    VF vf(20, 12);

    vf.addVertex(1.0, 1.0, 1.0);
    vf.addVertex(1.0, 1.0, -1.0);
    vf.addVertex(1.0, -1.0, 1.0);
    vf.addVertex(1.0, -1.0, -1.0);
    vf.addVertex(-1.0, 1.0, 1.0);
    vf.addVertex(-1.0, 1.0, -1.0);
    vf.addVertex(-1.0, -1.0, 1.0);
    vf.addVertex(-1.0, -1.0, -1.0);
    vf.addVertex(0.0, c1, c2);
    vf.addVertex(0.0, c1, -c2);
    vf.addVertex(0.0, -c1, c2);
    vf.addVertex(0.0, -c1, -c2);
    vf.addVertex(c1, c2, 0.0);
    vf.addVertex(c1, -c2, 0.0);
    vf.addVertex(-c1, c2, 0.0);
    vf.addVertex(-c1, -c2, 0.0);
    vf.addVertex(c2, 0.0, c1);
    vf.addVertex(c2, 0.0, -c1);
    vf.addVertex(-c2, 0.0, c1);
    vf.addVertex(-c2, 0.0, -c1);

    vf.addFace({ 2, 16, 18, 6, 10 });
    vf.addFace({ 6, 18, 4, 14, 15 });
    vf.addFace({ 4, 18, 16, 0, 8 });
    vf.addFace({ 0, 16, 2, 13, 12 });
    vf.addFace({ 1, 9, 8, 0, 12 });
    vf.addFace({ 5, 14, 4, 8, 9 });
    vf.addFace({ 7, 15, 14, 5, 19 });
    vf.addFace({ 11, 10, 6, 15, 7 });
    vf.addFace({ 3, 13, 2, 10, 11 });
    vf.addFace({ 17, 1, 12, 13, 3 });
    vf.addFace({ 17, 3, 11, 7, 19 });
    vf.addFace({ 5, 9, 1, 17, 19 });

    vf.NormalizeVertices();
    vf.Scale(radius);
    
    return vf;
}

VF geometries::EllipticParaboloid(
    double a, 
    double b, 
    double width, 
    double height, 
    size_t ws, 
    size_t hs, 
    const Eigen::Vector3d& X, 
    const Eigen::Vector3d& Y)
{
    return Grid(
        width, 
        height, 
        ws, 
        hs, 
        [&a, &b](double& x, double& y) -> double 
        {
            return ((x * x) / (a * a)) + ((y * y) / (b * b)); 
        }, 
        X, 
        Y);
}

VF geometries::EquilateralTriangle(
    double length, 
    const Eigen::Vector3d & U, 
    const Eigen::Vector3d & V)
{
    VF vf(3, 1);

    double halfLength = length / 2.0;
    double height = (length * sqrt(3.0)) / 2.0;

    vf.addVertex((halfLength * U) - ((height / 3.0) * V));
    vf.addVertex(((2.0 * height) / 3.0) * V);
    vf.addVertex((-halfLength * U) - ((height / 3.0) * V));

    vf.addFace(0, 1, 2);

    return vf;
}

VF geometries::Hexahedron(double radius)
{
    VF vf(8, 6);

    vf.addVertex(1.0, 1.0, 1.0);
    vf.addVertex(1.0, 1.0, -1.0);
    vf.addVertex(1.0, -1.0, 1.0);
    vf.addVertex(1.0, -1.0, -1.0);
    vf.addVertex(-1.0, 1.0, 1.0);
    vf.addVertex(-1.0, 1.0, -1.0);
    vf.addVertex(-1.0, -1.0, 1.0);
    vf.addVertex(-1.0, -1.0, -1.0);

    vf.addFace({ 0, 4, 6, 2 });
    vf.addFace({ 1, 3, 7, 5 });
    vf.addFace({ 0, 2, 3, 1 });
    vf.addFace({ 6, 4, 5, 7 });
    vf.addFace({ 2, 6, 7, 3 });
    vf.addFace({ 4, 0, 1, 5 });

    vf.NormalizeVertices();
    vf.Scale(radius);
    
    return vf;
}

VF geometries::Grid(
    double width, 
    double height, 
    size_t ws, 
    size_t hs, 
    const std::function<double(double& x, double& y)>& Z, 
    const Eigen::Vector3d& X, 
    const Eigen::Vector3d& Y)
{
    // Initialize the vertex coordinates and vertex indices structure
    VF vf((ws + 1) * (hs + 1), ws * hs);

    // Calculate the half lengths of the square
    double halfWidth = width / 2.0;
    double halfHeight = height / 2.0;

    // Calculate the width and height step values
    double widthStep = width / (double)ws;
    double heightStep = height / (double)hs;

    Eigen::Vector3d P = Eigen::Vector3d::Zero(), N = X.cross(Y).normalized();

    double w, h, n;

    size_t i, j;

    // Define the vertex coordinates of the geometry
    for (i = 0; i <= hs; i += 1)
    {
        for (j = 0; j <= ws; j += 1)
        {
            w = -halfWidth + ((double)j * widthStep);
            h = -halfHeight + ((double)i * heightStep);
            n = Z(w, h);

            P << (w * X) + (h * Y) + (n * N);

            vf.addVertex(P);
        }
    }

    size_t v0, v1, v2, v3;

    // Define the vertex indices for the faces
    for (i = 0; i < hs; i += 1)
    {
        for (j = 0; j < ws; j += 1)
        {
            // Calculate the vertex indices that define the current face
            v0 = j + (i * (ws + 1));
            v1 = v0 + 1;
            v2 = (j + 1) + ((i + 1) * (ws + 1));
            v3 = v2 - 1;

            vf.addFace(v0, v1, v2, v3);
        }
    }

    // Return the vertex coordinates and vertex indices of the geometry
    return vf;
}

VF geometries::HyperbolicParaboloid(
    double a, 
    double b, 
    double width, 
    double height, 
    size_t ws, 
    size_t hs, 
    const Eigen::Vector3d& X, 
    const Eigen::Vector3d& Y)
{
    return Grid(
        width,
        height,
        ws,
        hs,
        [&a, &b](double& x, double& y) -> double
        {
            return ((y * y) / (b * b)) - ((x * x) / (a * a));
        },
        X,
            Y);
}

VF geometries::Icosahedron(double radius)
{
    double c = (1.0 + sqrt(5.0)) / 2.0;

    VF vf(12, 20);
    
    vf.addVertex(0.0, 1.0, c);
    vf.addVertex(0.0, 1.0, -c);
    vf.addVertex(0.0, -1.0, c);
    vf.addVertex(0.0, -1.0, -c);
    vf.addVertex(1.0, c, 0.0);
    vf.addVertex(1.0, -c, 0.0);
    vf.addVertex(-1.0, c, 0.0);
    vf.addVertex(-1.0, -c, 0.0);
    vf.addVertex(c, 0.0, 1.0);
    vf.addVertex(c, 0.0, -1.0);
    vf.addVertex(-c, 0.0, 1.0);
    vf.addVertex(-c, 0.0, -1.0);

    vf.addFace(2, 0, 10);
    vf.addFace(2, 8, 0);
    vf.addFace(8, 4, 0);
    vf.addFace(4, 6, 0);
    vf.addFace(6, 10, 0);
    vf.addFace(10, 7, 2);
    vf.addFace(7, 5, 2);
    vf.addFace(5, 8, 2);
    vf.addFace(9, 8, 5);
    vf.addFace(9, 4, 8);
    vf.addFace(1, 4, 9);
    vf.addFace(1, 6, 4);
    vf.addFace(11, 6, 1);
    vf.addFace(11, 10, 6);
    vf.addFace(11, 7, 10);
    vf.addFace(3, 5, 7);
    vf.addFace(9, 5, 3);
    vf.addFace(3, 7, 11);
    vf.addFace(1, 3, 11);
    vf.addFace(3, 1, 9);

    vf.NormalizeVertices();
    vf.Scale(radius);
    
    return vf;
}

VF geometries::Octahedron(double radius) 
{
    VF vf(6, 8);

    vf.addVertex(1.0, 0.0, 0.0);
    vf.addVertex(-1.0, 0.0, 0.0);
    vf.addVertex(0.0, 1.0, 0.0);
    vf.addVertex(0.0, -1.0, 0.0);
    vf.addVertex(0.0, 0.0, 1.0);
    vf.addVertex(0.0, 0.0, -1.0);
    
    vf.addFace(0, 2, 4);
    vf.addFace(2, 1, 4);
    vf.addFace(1, 3, 4);
    vf.addFace(3, 0, 4);
    vf.addFace(5, 0, 3);
    vf.addFace(5, 2, 0);
    vf.addFace(5, 1, 2);
    vf.addFace(5, 3, 1);

    vf.NormalizeVertices();
    vf.Scale(radius);
    
    return vf;
}

VF geometries::Polygon(
    size_t sides, 
    double length, 
    const Eigen::Vector3d & U, 
    const Eigen::Vector3d & V)
{
    VF vf(sides, 1);

    double angle = 0.0;
    double angleStep = 2.0 * utils::PI / (double)sides;

    // Calculate the radius of the circumference that contains the n-sided regular polygon with the
    // given side length
    double radius = length / (2.0 * cos(utils::PI / (double)sides));

    Eigen::Vector3d P = Eigen::Vector3d::Zero();

    std::vector<size_t> indices(sides, 0);

    for (size_t i = 0; i < sides; i += 1) 
    {
        P << radius * ((cos(angle) * U) + (sin(angle) * V));

        vf.addVertex(P);

        indices[i] = i;

        angle += angleStep;
    }

    vf.addFace(indices);

    return vf;
}

double geometries::Random01()
{
    return ((double)rand() / (RAND_MAX)) + 1;
}

VF geometries::Rectangle(
    double width, 
    double height, 
    const Eigen::Vector3d & W, 
    const Eigen::Vector3d & H)
{
    VF vf(4, 1);

    Eigen::Vector3d _W = (width / 2.0) * W, _H = (height / 2.0) * H;

    vf.addVertex( _W + _H);
    vf.addVertex(-_W + _H);
    vf.addVertex(-_W - _H);
    vf.addVertex( _W - _H);

    vf.addFace(0, 1, 2, 3);

    return vf;
}

VF geometries::Saddle(
    double width, 
    double height, 
    size_t ws, 
    size_t hs, 
    const Eigen::Vector3d& X, 
    const Eigen::Vector3d& Y)
{
    return Grid(
        width, 
        height, 
        ws, 
        hs, 
        [](double& x, double& y) -> double 
        {
            return x * y; 
        }, 
        X, Y);
}

VF geometries::Square(
    double length, 
    const Eigen::Vector3d & U, 
    const Eigen::Vector3d & V)
{
    VF vf(4, 1);

    Eigen::Vector3d _U = (length / 2.0) * U, _V = (length / 2.0) * V;

    vf.addVertex( _U + _V);
    vf.addVertex(-_U + _V);
    vf.addVertex(-_U - _V);
    vf.addVertex( _U - _V);

    vf.addFace(0, 1, 2, 3);

    return vf;
}

VF geometries::SquareGrid(
    double width, 
    double height, 
    size_t ws, 
    size_t hs, 
    const Eigen::Vector3d& X, 
    const Eigen::Vector3d& Y)
{
    return Grid(
        width, 
        height, 
        ws, 
        hs, 
        [](double& x, double& y) -> double 
        {
            return 0; 
        }, 
        X, 
        Y);
}

VF geometries::Tetrahedron(double radius)
{
    double c = 1.0 / sqrt(2.0);

    VF vf(4, 4);

    vf.addVertex(1.0,  0.0, -c);
    vf.addVertex(-1.0,  0.0, -c);
    vf.addVertex(0.0,  1.0,  c);
    vf.addVertex(0.0, -1.0,  c);
    
    vf.addFace(0, 1, 2);
    vf.addFace(0, 2, 3);
    vf.addFace(1, 0, 3);
    vf.addFace(1, 3, 2);

    vf.NormalizeVertices();
    vf.Scale(radius);
    
    return vf;
}

VF geometries::Torus(
    double R,
    double r,
    size_t Rs,
    size_t rs, 
    const Eigen::Vector3d& C,
    const Eigen::Vector3d& X,
    const Eigen::Vector3d& Y)
{
    const Eigen::Vector3d Z = X.cross(Y);

    size_t n = Rs * rs;

    // Initialize the vertex coordinates and vertex indices structure
    VF vf(n, n);

    // Calculate the angle steps
    double uStep = TWO_PI / (double)Rs;
    double vStep = TWO_PI / (double)rs;

    double u, sin_u, cos_u, v, sin_v, cos_v;

    // Define the vertices of the torus
    for (size_t l = 0; l < Rs; l += 1) 
    {
        u = (double)l * uStep;
        sin_u = sin(u);
        cos_u = cos(u);

        Eigen::Vector3d T1 = C + (X * (R * cos_u)) + (Y * (R * sin_u));
        Eigen::Vector3d T2 = (X * -sin_u) + (Y * cos_u);
        Eigen::Vector3d T3 = T2.cross(Z).normalized();

        for (size_t m = 0; m < rs; m += 1) 
        {
            v = (double)m * vStep;
            sin_v = sin(v);
            cos_v = cos(v);

            Eigen::Vector3d Q = T1 + (T3 * (r * cos_v)) + (Z * (r * sin_v));

            vf.addVertex(Q);
        }
    }

    // Declare the variables for the vertex indices for the current face
    size_t v0, v1, v2, v3;

    // Define the vertex indices for the faces
    for (size_t l = 0; l < Rs; l += 1)
    {
        for (size_t m = 0; m < rs; m += 1)
        {
            // Determine the vertex indices of the current face according to the current vertex 
            // indices values
            if (l < Rs - 1)
            {
                if (m < rs - 1)
                {
                    // This is the common case: All squares not at the end of the ring or at the 
                    // last ring. This case happens mostly of the times
                    v0 = (l * rs) + m;
                    v1 = ((l + 1) * rs) + m;
                    v2 = v1 + 1;
                    v3 = v0 + 1;
                }
                else
                {
                    // This is the last square in a minor ring (not the last one). This case happens
                    // once per minor ring
                    v0 = (l * rs) + m;
                    v1 = ((l + 1) * rs) + m;
                    v2 = v0 + 1;
                    v3 = v0 - rs + 1;
                }
            }
            else
            {
                if (m < rs - 1)
                {
                    // This is a square (not the last) in the last minor ring. It uses the first 
                    // vertices generated for the torus. This case happens mostly in the last minor 
                    // ring
                    v0 = (l * rs) + m;
                    v1 = m;
                    v2 = v1 + 1;
                    v3 = v0 + 1;
                }
                else
                {
                    // This is the very last square in the last minor ring. This case happens once
                    v0 = (l * rs) + m;
                    v1 = m;
                    v2 = 0;
                    v3 = v0 - rs + 1;
                }
            }

            vf.addFace(v0, v1, v2, v3);
        }
    }

    return vf;
}

VF geometries::TruncatedCone(
    const Eigen::Vector3d & C, 
    const Eigen::Vector3d & K, 
    double br, 
    double tr, 
    double length, 
    size_t rs, 
    size_t ls)
{
    // Initialize the vertex coordinates and vertex indices structure
    VF vf(rs * (ls + 1), rs * ls);

    // Calculate the required constants
    double halfLength = length / 2.0;
    double angle = 2.0 * utils::PI / (double)rs;

    // Define the end points of the cylinder
    Eigen::Vector3d temp = K.normalized() * halfLength;
    Eigen::Vector3d A = C - temp;
    Eigen::Vector3d B = C + temp;

    // Get the vector AB
    Eigen::Vector3d AB = B - A;

    // Calculate a point somewhere but along the AB line segment (TO DO in next versions)
    double xRand = floor(Random01() * 10.0) + 1.0;
    double yRand = floor(Random01() * 10.0) + 1.0;
    double zRand = floor(Random01() * 10.0) + 1.0;
    Eigen::Vector3d P(xRand, yRand, zRand);

    // Calculate the basis vector R
    temp = P - C;
    Eigen::Vector3d R = K.cross(temp).normalized();

    // Calculate the basis vector S
    Eigen::Vector3d S = K.cross(R).normalized();

    size_t l, r;
    double t, radius, a, sin_a, cos_a;

    // Define the vertices of the cylinder
    for (l = 0; l <= ls; l += 1)
    {
        // Calculate the first term of the parametric equation (it is the same for all vertices at 
        // this length step)
        t = (double)l / (double)ls;
        Eigen::Vector3d term1 = (AB * t) + A;

        // Calculate the radius for the current section of the truncated cone
        radius = br + (t * (tr - br));

        for (r = 0; r < rs; r += 1)
        {
            // Calculate angle a and its sine and cosine values
            a = (double)r * angle;
            sin_a = sin(a);
            cos_a = cos(a);

            // Calculate the second and third term of the parametric equation
            Eigen::Vector3d r_R_cos_a = R * (radius * cos_a);
            Eigen::Vector3d r_R_sin_a = S * (radius * sin_a);

            // Calculate the vertex and push it into V
            Eigen::Vector3d P = term1 + r_R_cos_a + r_R_sin_a;
            vf.addVertex(P);
        }
    }

    // Declare the variables for the vertex indices for the current face
    size_t v0, v1, v2, v3;

    // Define the vertex indices for the faces
    for (l = 0; l < ls; l += 1)
    {
        for (r = 0; r < rs; r += 1)
        {
            // Calculate the vertex indices that define the current face. The last face of each ring
            // uses the first and last indices in the current ring and the next one
            if (r == rs - 1)
            {
                v0 = (l * rs) + r;
                v1 = v0 - (rs - 1);
                v2 = v0 + 1;
                v3 = v2 + (rs - 1);
            }
            else
            {
                v0 = (l * rs) + r;
                v1 = v0 + 1;
                v2 = ((l + 1) * rs) + (r + 1);
                v3 = v2 - 1;
            }

            vf.addFace(v0, v1, v2, v3);
        }
    }
    
    return vf;
}
