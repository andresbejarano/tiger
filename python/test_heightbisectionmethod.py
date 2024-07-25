# -*- coding: utf-8 -*-

# The following code shows the TIC generation process based on the Height 
# Bisection Method.

import geometries.planar as planar
import workspace.workspace as workspace
import vtk
import apps.vtk.vtkworkspace as vtkworkspace

# Pick a geometric domain
S = planar.squaredgrid(2, 2, 2, 2)

# Set up the workspace
W = workspace.Workspace()
W.setdomain(S)

W.setedgemidpoints()
#print()
#print("Half edge midpoints")
#for h in W.halfedgemidpoints:
#    print(W.halfedgemidpoints[h])

W.setedgedirections()
#print()
#print("Half edge directions")
#for h in W.halfedgedirections:
#    print(W.halfedgedirections[h])

W.setfaceheights(2, 2)
#print()
#print("face heights")
#for f in W.faceheights:
#    print(f"top:{W.faceheights[f]['top']}")
#    print(f"bottom:{W.faceheights[f]['bottom']}")

W.setheightbisectionelements()
#print()
#print("Half edge angles")
#for h in W.halfedgerotationangles:
#    print(W.halfedgerotationangles[h])
#print()
#print("Half edge rotated vectors")
#for h in W.halfedgevectors:
#    P = W.halfedgevectors[h].clone().fixzeros()
#    print(f"({P.x}, {P.y}, {P.z})")

W.setedgeplanes()

W.calculateblocks()

# Calculate the block vertices
#print()
#print("Block vertices")
#for face in W.domain.F:
#    for v in W.assembly[face][0]:
#        t = v.clone().fixzeros()
#        print(f"({t.x}, {t.y}, {t.z})")
#    for f in W.assembly[face][1]:
#        print(f)

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
