# -*- coding: utf-8 -*-

# The following code shows the generation of a Mecke block using the Evolution 
# Steps Method.

from ds.vf import VF
import toolkit.utils as utils
from toolkit.plane import Plane
import workspace.workspace as workspace
import vtk
import apps.vtk.vtkworkspace as vtkworkspace
import math


# Set the angles for each segment
radius = 1
Aa = 0.2
Ab = utils.HALF_PI - Aa

# Create the irregular octagon
S = VF()
angle = 0
for i in range(4):
    S.addvertex(radius * math.cos(angle), radius * math.sin(angle))
    angle += Aa
    S.addvertex(radius * math.cos(angle), radius * math.sin(angle))
    angle += Ab

S.addface([i for i in range(8)])

# Set up the workspace
W = workspace.Workspace()
W.setdomain(S)



# Set the edge directions properly
W.halfedgedirections = {}
h = W.domain.F[0].halfedge
W.halfedgedirections[h.previous] = 1
W.halfedgedirections[h] = 1
W.halfedgedirections[h.next] = 1
h = h.next.next
W.halfedgedirections[h] = -1
h = h.next.next
W.halfedgedirections[h.previous] = 1
W.halfedgedirections[h] = 1
W.halfedgedirections[h.next] = 1
h = h.next.next
W.halfedgedirections[h] = -1

W.setedgemidpoints()
W.setedgevectorsfromnormals()
W.setedgerotationangles(utils.toradians(45))
W.rotateedgevectors()
W.setedgeplanes()

W.setfaceheights(4, 4)


face = W.domain.F[0]
normal = face.normal().normalize()
centroid = face.centroid()
topplane = Plane(normal, centroid + normal * W.faceheights[face]['top'])
toppoints = {}

# Traverse through the half edges of the face. For each iteration, consider the
# ray that passes through h.start. Such a ray is the intersection between the
# planes incident to h.previous and current h, respectively.
h = face.halfedge
while True:
    
    # Find the ray that intersects the planes incident to the previous half 
    # edge and the current one. This ray must pass through h.start
    isintersect, ray = utils.twoplanesintersection(W.halfedgeplanes[h.previous], W.halfedgeplanes[h])
    assert isintersect, "Whoops #1"
    
    # Find the intersection between the ray and the top plane
    isintersect, t = utils.planerayintersection(topplane, ray)
    assert isintersect, "Whoops #2"
    
    # Use a dictionary to store the resultant points in their calculated order 
    # and without duplications
    toppoints[ray.at(t).fixzeros()] = None
    
    # Move to the next half edge of the face. Break the cycle if we return to 
    # the first half edge
    h = h.next
    if h == face.halfedge:
        break


bottomplane = Plane(normal, centroid - normal * W.faceheights[face]['bottom'])
bottompoints = {}

while True:
    
    # Find the ray that intersects the planes incident to the current half edge
    # and the previous one. This ray must pass through h.start
    isintersect, ray = utils.twoplanesintersection(W.halfedgeplanes[h], W.halfedgeplanes[h.previous])
    assert isintersect, "Whoops #1"
    
    # Find the intersection between the ray and the bottom plane
    isintersect, t = utils.planerayintersection(bottomplane, ray)
    assert isintersect, "Whoops #2"
    
    # Use a dictionary to store the resultant points in their calculated order 
    # and without duplications
    bottompoints[ray.at(t).fixzeros()] = None
    
    h = h.previous
    if h == face.halfedge:
        break


points = list(toppoints) + list(bottompoints)
for p in points:
    print(p)

if False:
    width = 800
    height = 600
    VW = vtkworkspace.VTKWorkspace(width, height)
    VW.loaddomainactors(W.domain)
    
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetSize(width, height)
    renderWindow.AddRenderer(VW.renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()
    renderWindow.Render()
    renderWindowInteractor.Start()
