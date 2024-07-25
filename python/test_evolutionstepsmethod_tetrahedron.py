# -*- coding: utf-8 -*-

import geometries.planar as planar
import toolkit.utils as utils
import math
import vtk
from evolve import evolve






length = 1
Seed1top = planar.square(length)

# Set the angle and height values
angles = [utils.HALF_PI - math.atan(1 / 2 ** 0.5)] * 4
angles[1] *= -1
angles[3] *= -1
height = length / 2 ** 0.5

# A dictionary to store all points with their found index value. This will be
# used to reconstruct the evolved block.
allpoints = {}
alledges = set()
allfaces = set()

# Calculate the points resultant from the first positive evolution step
toppoints1 = evolve(Seed1top, angles, height, allpoints, alledges, allfaces)

Seed1bottom = Seed1top.copy().flip()
for i in range(len(angles)):
    angles[i] *= -1

bottompoints1 = evolve(Seed1bottom, angles, height, allpoints, alledges, allfaces)

# Create a renderer, render window, and interactor
windowwidth = 800
windowheight = 600
renderer = vtk.vtkRenderer()
renderer.SetBackground(1, 1, 1)  # Blue background for contrast

# Set up the camera
camera = vtk.vtkCamera()
camera.SetPosition(3, 3, 3)
camera.SetFocalPoint(0, 0, 0)
camera.SetViewUp(0, 0, 1)
renderer.SetActiveCamera(camera)

# Set the SSAO pass
scenediagonal = (windowwidth ** 2 + windowheight ** 2) ** 0.5
ssaopass = vtk.vtkSSAOPass()
ssaopass.SetRadius(0.1 * scenediagonal)
ssaopass.SetBias(0.001 * scenediagonal)
ssaopass.SetKernelSize(128)
ssaopass.BlurOff()
basicpasses = vtk.vtkRenderStepsPass()
ssaopass.SetDelegatePass(basicpasses)
renderer.SetPass(ssaopass)

# Create the vtkPoints object with ther vertices of the geometry
points = vtk.vtkPoints()
for p in allpoints:
    points.InsertNextPoint(p.x, p.y, p.z)

# Create a polydata object
pointPolydata = vtk.vtkPolyData()
pointPolydata.SetPoints(points)





# Create a glyph filter to show the points as spheres
sphere = vtk.vtkSphereSource()
sphere.SetRadius(0.1)  # Set the radius of the spheres

glyph = vtk.vtkGlyph3D()
glyph.SetInputData(pointPolydata)
glyph.SetSourceConnection(sphere.GetOutputPort())
glyph.ScalingOff()  # Turn off scaling by scalar values

# Create a mapper
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(glyph.GetOutputPort())

# Create an actor
pointsactor = vtk.vtkActor()
pointsactor.SetMapper(mapper)
pointsactor.GetProperty().SetColor(1, 0, 0)  # Red color
renderer.AddActor(pointsactor)





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
facesactor.GetProperty().SetColor(1, 1, 0)
facesactor.GetProperty().SetOpacity(0.5)

renderer.AddActor(facesactor)

# Populate a vtkCellArray object with the information of the edges
linesegments = vtk.vtkCellArray()
linesegments.SetNumberOfCells(len(alledges))
for edge in alledges:
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
edgesactor.GetProperty().SetColor(0, 0, 0)
edgesactor.GetProperty().SetLineWidth(4)

renderer.AddActor(edgesactor)

# Create a render window and the windoe interactor
renderWindow = vtk.vtkRenderWindow()
renderWindow.SetSize(windowwidth, windowheight)
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

# Render and interact
renderWindow.Render()
renderWindowInteractor.Start()
