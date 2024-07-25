
from ds.vf import VF
import geometries.planar as planar
import math
import toolkit.utils as utils
from workspace.workspace import Workspace
from toolkit.ray import Ray
from toolkit.linesegment import LineSegment
from toolkit.plane import Plane


def evolve(f:VF, angles:list, height:int|float) -> dict:
    assert f is not None, "Missing f"
    assert isinstance(f, VF), "f not a VF"
    assert angles is not None, "Missing angles"
    assert isinstance(angles, list), "angles not a list"
    assert isinstance(height, (int, float)), "height not a real number"
    
    # Initialize a dictionary to store the resultant points of the evolution.
    # The keys will be the resultant points. The dictionary will preserve the
    # insertion order and will prevent storing replicated points
    points = {}
    return points


S = planar.square()
print(S.normal(0))




sidelength = 1
height = sidelength / (2 ** 0.5)
angle = utils.HALF_PI - math.atan(1 / (2 ** 0.5))
S = planar.square()

W = Workspace()
W.setdomain(S)


W.setedgemidpoints()
W.setedgerotationangles(angle)
W.setedgevectorsfromnormals()
W.setedgedirections()
W.rotateedgevectors()
W.setedgeplanes()
W.setfaceheights(height, height)


face = W.domain.F[0]
h = face.halfedge


# Calculate the rays representing the line segments incident to each vertex of
# the face
linesegments = {}
while True:
    N = W.halfedgeplanes[h.previous].normal.cross(W.halfedgeplanes[h].normal)
    ray = Ray(h.start.coords, N.normalize())
    linesegments[h.start] = LineSegment(h.start.coords, ray.at(1))
    h = h.next
    if h == face.halfedge:
        break

# Find the intersection points
points = {}
centroid = face.centroid()
normal = face.normal().normalize()
topplane = Plane(normal, centroid + normal * W.faceheights[face]['top'])
h = face.halfedge
while True:
    
    # Calculate the intersection point between the line segments incident to 
    # the endpoints of the current half edge
    L1 = linesegments[h.start]
    L2 = linesegments[h.twin.start]
    intersect, Q, s, t = utils.linesegmentsintersection(L1, L2)
    
    # If there is an intersection between the line segments then check where it
    # lies with respect to the top plane. Otherwise, calculate the intersection
    # points between the first line segment and the top plane (the second line
    # segment will be handled by the next iteration of the while loop).
    if intersect:
        
        where = topplane.pointlocation(Q)
        if where <= 0:
            points[Q] = None
        else:
            intersect, t = utils.planelinesegmentinserction(topplane, L1)
            assert intersect, "Whoops"
            points[L1.at(t).fixzeros()] = None
    
    else:
        intersect, t = utils.planelinesegmentinserction(topplane, L1)
        assert intersect, "Whoops"
        points[L1.at(t).fixzeros()] = None
    
    h = h.next
    if h == face.halfedge:
        break

for p in points:
    print(p)