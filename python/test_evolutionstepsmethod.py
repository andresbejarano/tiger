# -*- coding: utf-8 -*-

# The following code shows the TIC generation process based on the Evolution 
# Steps Method.

import geometries.planar as planar
import toolkit.utils as utils
from toolkit.plane import Plane
import workspace.workspace as workspace
import vtk
import apps.vtk.vtkworkspace as vtkworkspace

# Pick a geometric domain
#S = planar.squaredgrid(2, 2, 2, 2)
S = planar.square()

# Set up the workspace
W = workspace.Workspace()
W.setdomain(S)

W.setedgemidpoints()
#print()
#print("Half edge midpoints")
#for h in W.halfedgemidpoints:
#    print(W.halfedgemidpoints[h])

#W.setedgedirections()
#print()
#print("Half edge directions")
#for h in W.halfedgedirections:
#    print(W.halfedgedirections[h])


face = W.domain.F[0]



W.halfedgedirections = {}
h = face.halfedge
while True:
    W.setedgedirection(h, -1)
    h = h.next
    if h == face.halfedge:
        break
    
W.setedgevectorsfromnormals()

W.setedgerotationangles(utils.toradians(60))

W.rotateedgevectors()

W.setedgeplanes()





W.setfaceheights(0.5, 0.5)
#print()
#print("face heights")
#for f in W.faceheights:
#    print(f"top:{W.faceheights[f]['top']}")
#    print(f"bottom:{W.faceheights[f]['bottom']}")



normal = face.normal().normalize()
topplane = Plane(normal, face.centroid() + normal * W.faceheights[face]['top'])
points = []

h = face.halfedge
while True:
    
    # Find the ray that intersects the planes incident to the previous half 
    # edge and the current one. This ray must pass through h.start
    isintersect, ray = utils.twoplanesintersection(W.halfedgeplanes[h.previous], W.halfedgeplanes[h])
    assert isintersect, "Whoops #1"
    
    # Find the intersection between the ray and the top plane
    isintersect, t = utils.planerayintersection(topplane, ray)
    assert isintersect, "Whoops #2"
    
    point = ray.at(t)
    points.append(point)
    
    
    h = h.next
    if h == face.halfedge:
        break





#W.calculateblocks()

# Calculate the block vertices
#print()
#print("Block vertices")
#for face in W.domain.F:
#    for v in W.assembly[face][0]:
#        t = v.clone().fixzeros()
#        print(f"({t.x}, {t.y}, {t.z})")
#    for f in W.assembly[face][1]:
#        print(f)

if False:
    width = 800
    height = 600
    VW = vtkworkspace.VTKWorkspace(width, height)
    VW.loadassemblyactors(W.assembly)
    
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetSize(width, height)
    renderWindow.AddRenderer(VW.renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()
    renderWindow.Render()
    renderWindowInteractor.Start()
