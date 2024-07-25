# -*- coding: utf-8 -*-

# The following code shows the TIC generation process based on the Tilting 
# Angle Method.

import geometries.planar as planar
import toolkit.utils as utils
import workspace.workspace as workspace
import vtk
import apps.vtk.vtkworkspace as vtkworkspace

# Pick a geometric domain
S = planar.squaredgrid(3, 3, 3, 3)

# Set up the workspace
W = workspace.Workspace()
W.setdomain(S)

W.setedgemidpoints()
#print()
#print("Half edge midpoints")
#for h in W.halfedgemidpoints:
#    print(W.halfedgemidpoints[h])

W.setedgedirections(-1)
#print()
#print("Half edge directions")
#for h in W.halfedgedirections:
#    print(W.halfedgedirections[h])

W.setedgevectorsfromnormals()
#print()
#print("Half edge vectors")
#for h in W.halfedgevectors:
#    print(W.halfedgevectors[h])
    
W.setedgerotationangles(utils.toradians(60))
#print()
#print("Half edge angles")
#for h in W.halfedgerotationangles:
#    print(W.halfedgerotationangles[h])

W.rotateedgevectors()
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

W.setvertexrays()



width = 800
height = 600
VW = vtkworkspace.VTKWorkspace(width, height)
VW.xyplaneactor.SetVisibility(False)
VW.axesactor.SetVisibility(False)
VW.renderer.SetBackground(1, 1, 1)

VW.loaddomainactors(W.domain)
VW.domaindcelactor.GetProperty().SetColor(1, 0, 0)
VW.domaindcelactor.GetProperty().SetDiffuseColor(1, 0, 0)
VW.domaindcelactor.GetProperty().SetSpecularColor(1, 0, 0)
VW.domaindcelactor.GetProperty().SetSpecular(0.51)
VW.domaindcelactor.GetProperty().SetDiffuse(0.7)
VW.domaindcelactor.GetProperty().SetAmbient(0.7)
VW.domaindcelactor.GetProperty().SetSpecularPower(30.0)
VW.domaindcelactor.GetProperty().SetOpacity(1.0)

VW.loadassemblyactors(W.assembly)
VW.assemblyactor.GetProperty().SetColor(0.87, 0.87, 0)
VW.assemblyactor.GetProperty().SetDiffuseColor(0.87, 0.87, 0)
VW.assemblyactor.GetProperty().SetSpecularColor(0.87, 0.87, 0)
VW.assemblyactor.GetProperty().SetSpecular(0.51)
VW.assemblyactor.GetProperty().SetDiffuse(0.7)
VW.assemblyactor.GetProperty().SetAmbient(0.7)
VW.assemblyactor.GetProperty().SetSpecularPower(30.0)
VW.assemblyactor.GetProperty().SetOpacity(0.5)

VW.assemblyedgesactor.GetProperty().SetColor(0.87, 0.87, 0)
VW.assemblyedgesactor.GetProperty().SetDiffuseColor(0.87, 0.87, 0)
VW.assemblyedgesactor.GetProperty().SetSpecularColor(0.87, 0.87, 0)
VW.assemblyedgesactor.GetProperty().SetSpecular(0.51)
VW.assemblyedgesactor.GetProperty().SetDiffuse(0.7)
VW.assemblyedgesactor.GetProperty().SetAmbient(0.7)
VW.assemblyedgesactor.GetProperty().SetSpecularPower(30.0)
VW.assemblyedgesactor.GetProperty().SetOpacity(1)
VW.assemblyedgesactor.GetProperty().SetLineWidth(3)

# Traverse through the faces and create the evolution rectangles
t = 0.7
for f in W.domain.F:
    points = []
    h = f.halfedge
    while True:
        
        R = W.vertexrays[f][h.start]
        points.append(R.at(t))
        
        h = h.next
        if h == f.halfedge:
            break
    
    # Create a vtkPoints object and insert the points
    vtkpoints = vtk.vtkPoints()
    for P in points:
        vtkpoints.InsertNextPoint(P.x, P.y, P.z)

    # Create the polygon connectivity
    polygon = vtk.vtkPolygon()
    polygon.GetPointIds().SetNumberOfIds(len(points))
    for i in range(len(points)):
        polygon.GetPointIds().SetId(i, i)

    # Create a cell array to store the polygon
    polygons = vtk.vtkCellArray()
    polygons.InsertNextCell(polygon)

    # Create a polydata object
    polygonPolyData = vtk.vtkPolyData()
    polygonPolyData.SetPoints(vtkpoints)
    polygonPolyData.SetPolys(polygons)

    # Create a mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polygonPolyData)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0, 0, 1)
    actor.GetProperty().SetDiffuseColor(0, 0, 1)
    actor.GetProperty().SetSpecularColor(0, 0, 1)
    actor.GetProperty().SetSpecular(0.51)
    actor.GetProperty().SetDiffuse(0.7)
    actor.GetProperty().SetAmbient(0.7)
    actor.GetProperty().SetSpecularPower(30.0)
    actor.GetProperty().SetOpacity(1)
    VW.renderer.AddActor(actor)
    
    


renderWindow = vtk.vtkRenderWindow()
renderWindow.SetSize(width, height)
renderWindow.AddRenderer(VW.renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)
renderWindowInteractor.Initialize()
renderWindow.Render()
renderWindowInteractor.Start()
