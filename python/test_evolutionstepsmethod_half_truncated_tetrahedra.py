# -*- coding: utf-8 -*-

import geometries.planar as planar
import geometries.closed as closed
import geometries.platonicsolids as platonic
from ds.vf import VF
from ds.dcel import DCEL
from workspace.workspace import Workspace
import toolkit.utils as utils
from toolkit.plane import Plane
from toolkit.ray import Ray
from toolkit.pvector import PVector
import math
import vtk
import random
from evolve import evolve


# Set the attributes for the first evolution step
angle = utils.HALF_PI - math.atan(1 / 2 ** 0.5)



#S = planar.square()
#height = 1 / (2 ** 0.5)


# Planar grid
#S = planar.squaredgrid(7, 7, 7, 7)
#height = 1 / (2 ** 0.5)

# Elliptic paraboloid
#S = planar.ellipticparaboloid(1, 1, 1, 1, 10, 10)
#S.scale(3)
#height = 0.2

# Barrel vault
#S = planar.barrelvault(2, 2, 5, 20, PVector(), PVector(1), PVector(z=1))
#height = 0.2

# Wave
#S = planar.wave(utils.PI, utils.PI, 10, 10)
#height = 0.2

# Saddle
#S = planar.saddle(2, 2, 10, 10)
#height = 0.13

# Torus
#S = closed.torus(3, 1, 30, 10)
#height = 0.25

# Sphere
S = DCEL(platonic.octahedron())
S = DCEL(S.facemidpointsubdivision())
S = DCEL(S.facemidpointsubdivision())
S = DCEL(S.facemidpointsubdivision())
S = S.normalizevertices(3).tovf()
height = 0.3




# For this TIC, let's use a Workspace object to store the shared elements 
# between seeds.
W = Workspace()
W.setdomain(S)




#

W.setedgedirections()
W.setedgerotationangles(angle)
#W.setfaceheights(height, height)
W.setedgemidpoints()
W.setedgevectorsfromnormals()
W.rotateedgevectors()
W.setedgeplanes()

e = 1e-8


# Traverse through the faces of the geometric domain and run the first positive
# and negative evolution steps
blockpoints = {}
for face in W.domain.F:
    
    blockpoints[face] = {'points':{}, 'edges':set(), 'faces':set()}
    
    vertexrays = {}
    h = face.halfedge
    while True:
        
        v = h.start.coords.clone().fixzeros(e)
        if v not in blockpoints[face]['points']:
            blockpoints[face]['points'][v] = len(blockpoints[face]['points'])
        
        N = W.halfedgeplanes[h.previous].normal.cross(W.halfedgeplanes[h].normal)
        vertexrays[h.start] = Ray(h.start.coords, N.normalize())
        
        h = h.next
        if h == face.halfedge:
            break
    
    # Set the face and top evolution planes
    normal = face.normal().normalize()
    centroid = face.centroid()
    faceplane = Plane(normal, centroid)
    topplane = Plane(normal, centroid + normal * (height / 2))
    
    # Proceed with the first positive evolution step for the current face
    points = {}
    h = face.halfedge
    while True:
        
        # Get the intersection (if any) between the rays incident to the end
        # points of the current half edge
        R1 = vertexrays[h.start]
        R2 = vertexrays[h.twin.start]
        result, Q, s, t = utils.raysintersection(R1, R2)
        v0idx = blockpoints[face]['points'][h.start.coords.clone().fixzeros(e)]
        
        # If there is an intersection point and it lies between the face and 
        # top planes then keep it. Otherwise, calculate the intersection of R1 
        # with such a plane (R2 will be handled in the next iteration).
        if result and faceplane.pointlocation(Q) >= 0 and topplane.pointlocation(Q) <= 0:
            
            Q.fixzeros()
            v1idx = len(blockpoints[face]['points'])
            if Q not in blockpoints[face]['points']:
                blockpoints[face]['points'][Q] = v1idx
            else:
                v1idx = blockpoints[face]['points'][Q]
            points[Q] = None
            blockpoints[face]['edges'].add((min(v0idx, v1idx), max(v0idx, v1idx)))
            
        else:
            result, t = utils.planerayintersection(topplane, R1)
            assert result, "Whoops"
            
            Q = R1.at(t).fixzeros()
            v1idx = len(blockpoints[face]['points'])
            if Q not in blockpoints[face]['points']:
                blockpoints[face]['points'][Q] = v1idx
            else:
                v1idx = blockpoints[face]['points'][Q]
            
            points[Q] = None
            blockpoints[face]['edges'].add((min(v0idx, v1idx), max(v0idx, v1idx)))
        
        # Move to the next half edge in the face. Break the while loop if we 
        # return to the first halfe dge
        h = h.next
        if h == face.halfedge:
            break
        
    P = [p for p in points]
    for i in range(len(P) - 1):
        v0idx = blockpoints[face]['points'][P[i]]
        v1idx = blockpoints[face]['points'][P[i + 1]]
        blockpoints[face]['edges'].add((min(v0idx, v1idx), max(v0idx, v1idx)))
    
    v0idx = blockpoints[face]['points'][P[0]]
    v1idx = blockpoints[face]['points'][P[len(P) - 1]]
    blockpoints[face]['edges'].add((min(v0idx, v1idx), max(v0idx, v1idx)))
    
    
    
    
    # Use the points from the first positive evolution step to run the second
    # positive evolution step. Since the second step uses different evolution
    # parameters, run the evolve function with its respective arguments.
    Seed2top = VF()
    for p in points:
        Seed2top.addvertexfrompvector(p)
    Seed2top.addface([i for i in range(len(points))])

    points = evolve(Seed2top, [-angle for _ in range(4)], height / 2, blockpoints[face]['points'], blockpoints[face]['edges'], blockpoints[face]['faces'])
    
    
    # Proceed with the first negative evolution step for the current face
    bottomplane = Plane(-normal, centroid - normal * height)
    points = {}
    h = face.halfedge
    while True:
        
        # Get the intersection (if any) between the rays incident to the end
        # points of the current half edge
        R1 = vertexrays[h.start]
        R2 = vertexrays[h.twin.start]
        result, Q, s, t = utils.raysintersection(R1, R2)
        v0idx = blockpoints[face]['points'][h.start.coords.clone().fixzeros(e)]
        
        # If there is an intersection point and it lies between the face and 
        # bottom planes then keep it. Otherwise, calculate the intersection of 
        # R1 with such a plane (R2 will be handled in the next iteration).
        if result and faceplane.pointlocation(Q) <= 0 and bottomplane.pointlocation(Q) <= 0:
            
            Q.fixzeros()
            v1idx = len(blockpoints[face]['points'])
            if Q not in blockpoints[face]['points']:
                blockpoints[face]['points'][Q] = v1idx
            else:
                v1idx = blockpoints[face]['points'][Q]
            points[Q] = None
            blockpoints[face]['edges'].add((min(v0idx, v1idx), max(v0idx, v1idx)))
            
        else:
            result, t = utils.planerayintersection(bottomplane, R1)
            assert result, "Whoops"
            
            Q = R1.at(t).fixzeros()
            v1idx = len(blockpoints[face]['points'])
            if Q not in blockpoints[face]['points']:
                blockpoints[face]['points'][Q] = v1idx
            else:
                v1idx = blockpoints[face]['points'][Q]
            
            points[Q] = None
            blockpoints[face]['edges'].add((min(v0idx, v1idx), max(v0idx, v1idx)))
        
        # Move to the previous half edge in the face. Break the while loop if 
        # we return to the first halfe dge
        h = h.previous
        if h == face.halfedge:
            break
    
    P = [p for p in points]
    for i in range(len(P) - 1):
        v0idx = blockpoints[face]['points'][P[i]]
        v1idx = blockpoints[face]['points'][P[i + 1]]
        blockpoints[face]['edges'].add((min(v0idx, v1idx), max(v0idx, v1idx)))
    
    v0idx = blockpoints[face]['points'][P[0]]
    v1idx = blockpoints[face]['points'][P[len(P) - 1]]
    blockpoints[face]['edges'].add((min(v0idx, v1idx), max(v0idx, v1idx)))
    
    
    # Use the points from the first negative evolution step to run the second
    # negative evolution step. Since the second step uses different evolution
    # parameters, run the evolve function with its respective arguments.
    #Seed2bottom = VF()
    #for p in points:
    #    Seed2bottom.addvertexfrompvector(p)
    #Seed2bottom.addface([i for i in range(len(points))])

    #points = evolve(Seed2bottom, [-angle for _ in range(4)], height / 2, blockpoints[face]['points'], blockpoints[face]['edges'])




# Create a renderer, render window, and interactor
windowwidth = 1024
windowheight = 768
renderer = vtk.vtkRenderer()
renderer.SetBackground(1, 1, 1)  # Blue background for contrast



# Set the SSAO pass
#scenediagonal = (windowwidth ** 2 + windowheight ** 2) ** 0.5
#ssaopass = vtk.vtkSSAOPass()
#ssaopass.SetRadius(0.1 * scenediagonal)
#ssaopass.SetBias(0.001 * scenediagonal)
#ssaopass.SetKernelSize(128)
#ssaopass.BlurOff()
#basicpasses = vtk.vtkRenderStepsPass()
#ssaopass.SetDelegatePass(basicpasses)
#renderer.SetPass(ssaopass)


colors = vtk.vtkNamedColors()

for face in W.domain.F:
    
    lowestZ = 0
    
    points = vtk.vtkPoints()
    for p in blockpoints[face]['points']:
        points.InsertNextPoint(p.x, p.y, p.z)
        lowestZ = min(lowestZ, p.z)
    
    # Create a polydata object
    pointPolydata = vtk.vtkPolyData()
    pointPolydata.SetPoints(points)
    
    # Create the convex hull of the points
    delaunay = vtk.vtkDelaunay3D()
    delaunay.SetInputData(pointPolydata)
    delaunay.Update()
    
    # To extract the surface (convex hull)
    surfaceFilter = vtk.vtkDataSetSurfaceFilter()
    surfaceFilter.SetInputConnection(delaunay.GetOutputPort())
    surfaceFilter.Update()
    
    # Create a mapper and actor to display the convex hull
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(surfaceFilter.GetOutputPort())
    
    facesactor = vtk.vtkActor()
    facesactor.SetMapper(mapper)
    r, g, b = random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)
    facesactor.GetProperty().SetAmbientColor(r, g, b)
    #facesactor.GetProperty().SetAmbientColor(colors.GetColor3d('SaddleBrown'))
    facesactor.GetProperty().SetDiffuseColor(colors.GetColor3d('Sienna'))
    facesactor.GetProperty().SetSpecularColor(colors.GetColor3d('White'))
    facesactor.GetProperty().SetSpecular(0.51)
    facesactor.GetProperty().SetDiffuse(0.7)
    facesactor.GetProperty().SetAmbient(0.7)
    facesactor.GetProperty().SetSpecularPower(30.0)
    facesactor.GetProperty().SetOpacity(1.0)
    renderer.AddActor(facesactor)
    
    # Populate a vtkCellArray object with the information of the edges
    linesegments = vtk.vtkCellArray()
    linesegments.SetNumberOfCells(len(blockpoints[face]['edges']))
    for edge in blockpoints[face]['edges']:
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, edge[0])
        line.GetPointIds().SetId(1, edge[1])
        linesegments.InsertNextCell(line)
    
    edgesdata = vtk.vtkPolyData()
    edgesdata.SetPoints(points)
    edgesdata.SetLines(linesegments)
    
    edgesmapper = vtk.vtkPolyDataMapper()
    edgesmapper.SetInputData(edgesdata)
    
    edgesactor = vtk.vtkActor()
    edgesactor.SetMapper(edgesmapper)
    #edgesactor.GetProperty().SetColor(0, 0, 0)
    edgesactor.GetProperty().SetLineWidth(1.5)
    edgesactor.GetProperty().SetAmbientColor(colors.GetColor3d('black'))
    edgesactor.GetProperty().SetDiffuseColor(colors.GetColor3d('black'))
    edgesactor.GetProperty().SetSpecularColor(colors.GetColor3d('black'))
    edgesactor.GetProperty().SetSpecular(0.51)
    edgesactor.GetProperty().SetDiffuse(0.7)
    edgesactor.GetProperty().SetAmbient(0.7)
    edgesactor.GetProperty().SetSpecularPower(30.0)
    edgesactor.GetProperty().SetOpacity(1.0)
    
    renderer.AddActor(edgesactor)


renderWindow = vtk.vtkRenderWindow()
renderWindow.SetSize(windowwidth, windowheight)
renderWindow.AddRenderer(renderer)

# Set up the camera
camera = vtk.vtkCamera()
camera.SetPosition(5, 5, 3)
camera.SetFocalPoint(0, 0, 0)
camera.SetViewUp(0, 0, 1)
renderer.SetActiveCamera(camera)

renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)
renderWindow.Render()
renderWindowInteractor.Start()
