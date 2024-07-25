# -*- coding: utf-8 -*-

# The following code shows the generation of a Truncated Tetrahedron using the 
# Evolution Steps Method.

from ds.vf import VF
import geometries.planar as planar
import toolkit.utils as utils
from toolkit.plane import Plane
import workspace.workspace as workspace
import vtk
import apps.vtk.vtkworkspace as vtkworkspace
import math


sidelength = 1
S1 = planar.square(sidelength)

# Set up the workspace
W1 = workspace.Workspace()
W1.setdomain(S1)

W1.setedgemidpoints()
W1.setedgedirections(-1)
W1.setedgevectorsfromnormals()
W1.setedgerotationangles(utils.HALF_PI - math.atan(1 / (2 ** 0.5)))
W1.rotateedgevectors()
W1.setedgeplanes()


# The top height is a truncated tetrahedron, the bottom height is for the 
# regular tetrahedron
W1.setfaceheights(0.5, sidelength / (2 ** 0.5))





face = W1.domain.F[0]
normal = face.normal().normalize()
centroid = face.centroid()

topplane = Plane(normal, centroid + normal * W1.faceheights[face]['top'])
firsttoppoints = {}

# Traverse through the half edges of the face. For each iteration, consider the
# ray that passes through h.start. Such a ray is the intersection between the
# planes incident to h.previous and current h, respectively.
h = face.halfedge
while True:
    
    # Find the ray that intersects the planes incident to the previous half 
    # edge and the current one. This ray must pass through h.start
    isintersect, ray = utils.twoplanesintersection(W1.halfedgeplanes[h.previous], W1.halfedgeplanes[h])
    assert isintersect, "Whoops #1"
    
    # Find the intersection between the ray and the top plane
    isintersect, t = utils.planerayintersection(topplane, ray)
    assert isintersect, "Whoops #2"
    
    # Use a dictionary to store the resultant points in their calculated order 
    # and without duplications
    firsttoppoints[ray.at(t).fixzeros()] = None
    
    # Move to the next half edge of the face. Break the cycle if we return to 
    # the first half edge
    h = h.next
    if h == face.halfedge:
        break


bottomplane = Plane(normal, centroid - normal * W1.faceheights[face]['bottom'])
bottompoints = {}

while True:
    
    # Find the ray that intersects the planes incident to the current half edge
    # and the previous one. This ray must pass through h.start
    isintersect, ray = utils.twoplanesintersection(W1.halfedgeplanes[h], W1.halfedgeplanes[h.previous])
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




# Second positive evolution step
S2 = VF()
for p in list(firsttoppoints):
    S2.addvertexfrompvector(p)
S2.addface([0, 1, 2, 3])



W2 = workspace.Workspace()
W2.setdomain(S2)
W2.setedgemidpoints()

W2.setedgevectorsfromnormals()
W2.setedgerotationangles(utils.HALF_PI - math.atan(1 / (2 ** 0.5)))
W2.rotateedgevectors()
W2.setedgeplanes()


# The top height is a truncated tetrahedron, the bottom height is for the 
# regular tetrahedron
W2.setfaceheights(0.5, sidelength / (2 ** 0.5))









# Add the block to the assembly by associating its geometry to the face it 
# comes from
W.assembly = {}
W.assembly[face] = [list(firsttoppoints) + list(bottompoints), blockfaces]


















windowwidth = 800
windowheight = 600
VW = vtkworkspace.VTKWorkspace(windowwidth, windowheight)
VW.loadassemblyactors(W.assembly)

renderWindow = vtk.vtkRenderWindow()
renderWindow.SetSize(windowwidth, windowheight)
renderWindow.AddRenderer(VW.renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)
renderWindowInteractor.Initialize()
renderWindow.Render()
renderWindowInteractor.Start()
