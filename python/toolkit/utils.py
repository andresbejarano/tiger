# -*- coding: utf-8 -*-

import math
from toolkit.pvector import PVector
from toolkit.plane import Plane
from toolkit.linesegment import LineSegment
from toolkit.ray import Ray

PI = math.acos(-1.0)
HALF_PI = PI / 2.0
TWO_PI = PI * 2.0


# Returns the determinant of the 2x2 matrix specified as:
# |a, b|
# |c, d|
def det2(a:int|float, b:int|float, c:int|float, d:int|float) -> float:
    return (a * d) - (c * b)


# Returns the determinant of the 3x3 matrix specified as:
# |a, b, c|
# |d, e, f|
# |g, h, i|
def det3(a:int|float, b:int|float, c:int|float, 
        d:int|float, e:int|float, f:int|float, 
        g:int|float, h:int|float, i:int|float) -> float:
    t1 = (a * ((e * i) - (f * h)))
    t2 = (b * ((d * i) - (f * g)))
    t3 = (c * ((d * h) - (e * g)))
    return t1 - t2 + t3


# Populates a list with unique points.
def getuniquepoints(L:list, e:int|float=1e-8) -> tuple[bool, list]:
    assert L is not None, "L is missing"
    assert isinstance(e, (int, float)), "e not a real number"
    unique = []
    for V in L:
        isunique = True
        for U in unique:
            diff = (V - U).abs().fixzeros(e)
            if diff.iszerovector():
                isunique = False
                break
        if isunique:
            unique.append(V.copy())
    return len(unique) > 0, unique


# Calculates N * ((Q - P) x (A - P)), where * is dot product and x is cross 
# product. This value indicates if a vector from P to Q goes in CCW direction
# with respect to the normal vector N and point A (the resultant vector and N 
# must point at the same half space)
def isleftturn(P:PVector, Q:PVector, plane:Plane, e:int|float=1e-8, r:int=10) -> int:
    assert P is not None, "Missing P"
    assert isinstance(P, PVector), "P not a PVector"
    assert Q is not None, "Missing Q"
    assert isinstance(Q, PVector), "Q not a PVector"
    assert plane is not None, "Missing plane"
    assert isinstance(plane, Plane), "plane not a Plane"
    assert isinstance(e, (int, float)), "e not a real number"
    assert isinstance(r, int), "r not an integer number"
    dot = round((Q - P).cross(plane.point - P).dot(plane.normal), r)
    return 0 if abs(dot) <= e else (1 if dot > 0 else -1)


# Checks if a point lies within the triangle defined by its vertices.
# https:// stackoverflow.com/questions/995445/determine-if-a-3d-point-is-within-a-triangle
def ispointintriangle(P:PVector, v0:PVector, v1:PVector, v2:PVector, e:int|float=1e-8, r:int=10) -> bool:
    assert P is not None, "Missing P"
    assert v0 is not None, "Missing v0"
    assert v1 is not None, "Missing v1"
    assert v2 is not None, "Missing v2"
    assert isinstance(P, PVector), "P not a PVector"
    assert isinstance(v0, PVector), "v0 not a PVector"
    assert isinstance(v1, PVector), "v1 not a PVector"
    assert isinstance(v2, PVector), "v2 not a PVector"
    assert isinstance(e, (int, float)), "e not a real number"
    assert isinstance(r, int), "r not an integer number"
    return sameside(P, v0, v1, v2, e, r) and sameside(P, v1, v0, v2, e, r) and sameside(P, v2, v0, v1, e, r)


# Calculates the intersection between two line segments. If both line segments 
# are skewed then it calculates the midpoint between the closest points from 
# both line segments. The algorithm is adapted from the solution by Paul Bourke
# in "The shortest line between two lines in 3D" found in 
# http://paulbourke.net/geometry/pointlineplane/
# L1: The reference to a line segment.
# L2: The reference to a line segment.
# e: The threshold for values close to zero.
# Returns:
# bool: Indicates whether there is an intersection point or not.
# P: The reference to store the intersection point.
# s: The parameter of the intersection point along L1.
# t: The parameter of the intersection point along L2.
def linesegmentsintersection(L1:LineSegment, L2:LineSegment, e:int|float=1e-8) -> tuple[bool, PVector, float, float]:
    assert L1 is not None, "Missing L1"
    assert isinstance(L1, LineSegment), "L1 not a LineSegment"
    assert L2 is not None, "Missing L2"
    assert isinstance(L2, LineSegment), "L2 not a LineSegment"
    assert isinstance(e, (int, float)), "e not a real number"
    
    t1 = L1.A - L2.A
    t2 = L2.direction()
    if abs(t2.x) <= e and abs(t2.y) <= e and abs(t2.z) <= e:
        return False, None, None, None
    
    t3 = L1.direction()
    if abs(t3.x) <= e and abs(t3.y) <= e and abs(t3.z) <= e:
        return False, None, None, None
    
    d1343 = (t1.x * t2.x) + (t1.y * t2.y) + (t1.z * t2.z)
    d4321 = (t2.x * t3.x) + (t2.y * t3.y) + (t2.z * t3.z)
    d1321 = (t1.x * t3.x) + (t1.y * t3.y) + (t1.z * t3.z)
    d4343 = (t2.x * t2.x) + (t2.y * t2.y) + (t2.z * t2.z)
    d2121 = (t3.x * t3.x) + (t3.y * t3.y) + (t3.z * t3.z)
    
    denom = (d2121 * d4343) - (d4321 * d4321)
    
    if abs(denom) <= e:
        return False, None, None, None
    
    numer = (d1343 * d4321) - (d1321 * d4343)
    
    # Calculate the parameters for points Pa and Pb
    s = numer / denom
    t = (d1343 + (d4321 * s)) / d4343
    
    # Calculate points Pa and Pb
    pa = (t3 * s) + L1.A
    pb = (t2 * t) + L2.A
    
    # Calculate the midpoint between points Pa and Pb
    P = (pa + pb) / 2.0
    
    return True, P, s, t


def planelinesegmentinserction(plane:Plane, linesegment:LineSegment, e:int|float=1e-8, r:int=10) -> tuple[bool, float]:
    assert plane is not None, "Missing plane"
    assert isinstance(plane, Plane), "plane not a Plane"
    assert linesegment is not None, "Missing linesegment"
    assert isinstance(linesegment, LineSegment), "linesegment not a LineSegment"
    assert isinstance(e, (int, float)), "e not a real number"
    assert isinstance(r, int), "r not an integer number"
    ray = Ray(linesegment.A, linesegment.direction())
    return planerayintersection(plane, ray, e, r)


def planerayintersection(plane:Plane, ray:Ray, e:int|float=1e-8, r:int=10) -> tuple[bool, float]:
    assert plane is not None, "Missing plane"
    assert isinstance(plane, Plane), "plane not a Plane"
    assert ray is not None, "Missing ray"
    assert isinstance(ray, Ray), "ray not a Ray"
    assert isinstance(e, (int, float)), "e not a real number"
    assert isinstance(r, int), "r not an integer number"
    
    # Calculate the dot product between the direction vector of the ray and the
    # normal vector of the plane
    denom = round(ray.direction.dot(plane.normal), r)
    
    # If the absolute value of the dot product is zero then return false since
    # there is no intersection between the plane and the ray
    if abs(denom) <= e:
        return False, None
    
    num = round((plane.point - ray.point).dot(plane.normal), r)
    
    # Calculate the parameter value along the ray for the intersection point 
    # between the plane and the ray. Then return true since there is an 
    # intersection point
    t = num / denom
    
    return True, t


# Calculates the centroid (arithmetic mean) of a list of points.
def pointscentroid(L:list|set|tuple) -> PVector:
    assert L is not None, "P is missing"
    assert isinstance(L, (list, set, tuple)), "L not a list, set, or tuple"
    n = len(L)
    S = PVector()
    if n == 0:
        return S
    for l in L:
        assert isinstance(l, PVector), "l not a PVector"
        S += l
    return S / n


# Calculates the intersection between two rays. The calculation relies on the
# line segments intersection method. So, each ray is transformed into a line
# segment so such a method is used. If both rays are skewed then it calculates 
# the midpoint between the closest points from both rays.
# Returns:
# bool: Indicates whether there is an intersection point or not.
# P: The reference to store the intersection point.
# s: The parameter of the intersection point along R1.
# t: The parameter of the intersection point along R2.
def raysintersection(R1:Ray, R2:Ray, e:int|float=1e-8) -> tuple[bool, PVector, float, float]:
    assert R1 is not None, "Missing R1"
    assert isinstance(R1, Ray), "R1 not a Ray"
    assert R2 is not None, "Missing R2"
    assert isinstance(R2, Ray), "R2 not a Ray"
    assert isinstance(e, (int, float)), "e not a real number"
    L1 = LineSegment(R1.point, R1.at(1))
    L2 = LineSegment(R2.point, R2.at(1))
    return linesegmentsintersection(L1, L2, e)


def sameside(P:PVector, Q:PVector, A:PVector, B:PVector, e:int|float=1e-8, r:int=10) -> float:
    assert P is not None, "Missing P"
    assert Q is not None, "Missing Q"
    assert A is not None, "Missing A"
    assert B is not None, "Missing B"
    assert isinstance(P, PVector), "P not a PVector"
    assert isinstance(Q, PVector), "Q not a PVector"
    assert isinstance(A, PVector), "A not a PVector"
    assert isinstance(B, PVector), "B not a PVector"
    assert isinstance(e, (int, float)), "e not a real number"
    assert isinstance(r, int), "r not an integer number"
    D = (B - A).normalize()
    cp1 = D.cross(P - A).normalize()
    cp2 = D.cross(Q - A).normalize()
    test = round(cp1.dot(cp2), r)
    return 0 if test <= e else test


def threeplanesintersection(A:Plane, B:Plane, C:Plane, e:int|float=1e-8, r:int=10) -> tuple[bool, PVector]:
    assert A is not None, "Missing A"
    assert isinstance(A, Plane), "A not a Plane"
    assert B is not None, "Missing B"
    assert isinstance(B, Plane), "B not a Plane"
    assert C is not None, "Missing C"
    assert isinstance(C, Plane), "C not a Plane"
    assert isinstance(e, (int, float)), "e not a real number"
    assert isinstance(r, int), "r not an integer number"
    
    Aa, Ab, Ac, Ad = A.equationcoefficients()
    Ba, Bb, Bc, Bd = B.equationcoefficients()
    Ca, Cb, Cc, Cd = C.equationcoefficients()
    
    # Get the determinant of the system
    D = det3(Aa, Ab, Ac, Ba, Bb, Bc, Ca, Cb, Cc)
    
    # If the determinant is zero return false since there is no intersection 
    # point between the planes
    if abs(D) <= e:
        return False, None
    
    # Use the Cramer's rule to obtain the coordinates of the intersection point
    D1 = det3(Ad, Ab, Ac, Bd, Bb, Bc, Cd, Cb, Cc)
    D2 = det3(Aa, Ad, Ac, Ba, Bd, Bc, Ca, Cd, Cc)
    D3 = det3(Aa, Ab, Ad, Ba, Bb, Bd, Ca, Cb, Cd)
    
    # Set the coordinates of the intersection point and return true
    return True, PVector(D1 / D, D2 / D, D3 / D)


def twoplanesintersection(P1:Plane, P2:Plane, e:int|float=1e-8) -> tuple[bool, Ray]:
    assert P1 is not None, "Missing P1"
    assert isinstance(P1, Plane), "P1 is not a Plane"
    assert P2 is not None, "Missing P2"
    assert isinstance(P2, Plane), "P2 is not a Plane"
    assert isinstance(e, (int, float)), "e not a real number"
    N = P1.normal.cross(P2.normal)
    a1, b1, c1, d1 = P1.equationcoefficients()
    a2, b2, c2, d2 = P2.equationcoefficients()
    a3, b3, c3 = N.x, N.y, N.z
    d1 *= -1
    d2 *= -1
    DET = det3(a1, b1, c1, a2, b2, c2, a3, b3, c3)
    
    # If the determinant is 0 then the planes are parallel and there is no
    # intersection between them
    if abs(DET) <= e:
        return False, None
    
    x = (d2 * det2(b1, c1, b3, c3) - d1 * det2(b2, c2, b3, c3))/ DET
    y = (d2 * det2(a3, c3, a1, c1) - d1 * det2(a3, c3, a2, c2))/ DET
    z = (d2 * det2(a1, b1, a3, b3) - d1 * det2(a2, b2, a3, b3))/ DET
    return True, Ray(PVector(x, y, z), N)


def todegrees(rads:int|float) -> float:
    assert isinstance(rads, (int, float)), "rads not a real number"
    return 180.0 * rads / PI


def toradians(degs:int|float) -> float:
    assert isinstance(degs, (int, float)), "degs not a real number"
    return PI * degs / 180.0
