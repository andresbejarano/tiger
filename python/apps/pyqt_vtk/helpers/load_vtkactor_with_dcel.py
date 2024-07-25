# -*- coding: utf-8 -*-

#from ds.vf import VF
from ds.dcel import DCEL
import vtk

    
# Loads the geometry of a DCEL object into two vtkActor objects (one for 
# triangular faces, and one for the edges). The reason we need two objects is 
# to maintain the logical correspondence between the polygonal nature of the
# geometry faces while drawing the edges. If we use a single object for both
# purposes, the drawn edges will correspond to the triangulation of the faces.
def loadVtkActorsWithDcel(dcel):
    assert dcel is not None, "dcel missing"
    assert isinstance(dcel, DCEL), "dcel not a DCEL"
    
    # A dictionary to map a vertex coordinate to its index in the vertices
    # array of the geometry
    vToI = dict()
    
    # A set to store the (src, dst) tuples representing the edges of the 
    # geometry.
    edges = set()

    # Define the points of the geometry
    nV = dcel.numvertices()
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(nV)
    for i in range(nV):
        points.InsertPoint(i, dcel.V[i].coords.tolist())
        vToI[str(dcel.V[i].coords.tolist())] = i
    
    # Define the triangles representing the faces of the geometry
    nF = dcel.numfaces()
    triangles = vtk.vtkCellArray()
    triangles.SetNumberOfCells(dcel.numtriangles())
    for i in range(nF):
        
        # The start vertex of the incident half edge of the face is the common
        # vertex for all triangles representing the current face.
        v0 = vToI[str(dcel.F[i].halfedge.start.coords.tolist())]
        
        # Walk through the half edges of the current face and define its
        # triangles
        currentHalfedge = dcel.F[i].halfedge.next
        v1 = vToI[str(currentHalfedge.start.coords.tolist())]
        
        # Define the information of the first edge of the current face
        edges.add((min(v0, v1), max(v0, v1)))
        
        while True:
            
            v2 = vToI[str(currentHalfedge.next.start.coords.tolist())]
            
            triangle = vtk.vtkTriangle()
            triangle.GetPointIds().SetId(0, v0)
            triangle.GetPointIds().SetId(1, v1)
            triangle.GetPointIds().SetId(2, v2)
            triangles.InsertNextCell(triangle)
            
            #edges.add((min(v0, v1), max(v0, v1)))
            edges.add((min(v1, v2), max(v1, v2)))
            
            # Move to the next half edge of the face. Break the cycle if
            # necessary
            v1 = v2
            currentHalfedge = currentHalfedge.next
            
            if currentHalfedge == dcel.F[i].halfedge.previous:
                break
        
        # Finish the last edge of the current face
        edges.add((min(v0, v1), max(v0, v1)))
    
    # Define the polydata object that contains the information of the faces of
    # the geometric object to be rendered as an actor
    facesdata = vtk.vtkPolyData()
    facesdata.SetPoints(points)
    facesdata.SetPolys(triangles)
    
    # Define polydatanormals object to calculate the normals of the points
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(facesdata)
    normals.ComputePointNormalsOn()
    normals.ComputeCellNormalsOn()
    normals.Update()
    
    # Create a mapper and actor
    facesmapper = vtk.vtkPolyDataMapper()
    facesmapper.SetInputConnection(normals.GetOutputPort())
    
    dcelactor = vtk.vtkActor()
    dcelactor.SetMapper(facesmapper)
    
    # Populate a vtkCellArray object with the information of the edges
    linesegments = vtk.vtkCellArray()
    linesegments.SetNumberOfCells(len(edges))
    for edge in edges:
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
    
    return dcelactor, edgesactor


def testactor(G, width=800, height=600):
    assert G is not None, "Missing G"
    assert isinstance(G, DCEL), "G not a DCEL"
    assert isinstance(width, int), "width not an integer number"
    assert isinstance(height, int), "height not an integer number"
    assert width > 0, "width not a positive integer"
    assert height > 0, "height not a positive integer"
    #T = geom.hexahedron()
    #Td = DCEL(G)
    
    dcelactor, edgesactor = loadVtkActorsWithDcel(G)
    dcelactor.GetProperty().BackfaceCullingOn()
    dcelactor.GetProperty().SetColor(1, 1, 0)
    
    edgesactor.GetProperty().SetColor(0, 0, 0)  # Set color to black
    edgesactor.GetProperty().SetLineWidth(1)  # Set line width
    
    sceneDiagonal = (width ** 2 + height ** 2) ** 0.5
    
    # Define the SSAO pass
    ssaoPass = vtk.vtkSSAOPass()
    ssaoPass.SetRadius(0.1 * sceneDiagonal)
    ssaoPass.SetBias(0.001 * sceneDiagonal)
    ssaoPass.SetKernelSize(128)
    ssaoPass.BlurOff()
    basicPasses = vtk.vtkRenderStepsPass()
    ssaoPass.SetDelegatePass(basicPasses)
    
    # Create axes actor
    axesactor = vtk.vtkAxesActor()
    axesactor.SetTotalLength(1, 1, 1)  # Set the length of the axes
    axesactor.SetAxisLabels(1)  # Enable axis labels

    # Create a renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(dcelactor)
    renderer.AddActor(edgesactor)
    renderer.AddActor(axesactor)
    renderer.SetBackground(0.1, 0.2, 0.3)  # Set background to dark blue
    #renderer.SetBackground(1, 1, 1)  # Set background to white
    renderer.SetPass(ssaoPass)
    
    # Set up the camera
    camera = vtk.vtkCamera()
    camera.SetPosition(3, 3, 3)
    camera.SetFocalPoint(0, 0, 0)
    camera.SetViewUp(0, 0, 1)
    renderer.SetActiveCamera(camera)

    # Create a render window
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(width, height)
    renderWindow.SetWindowName("Loading DCEL into vtkActor")

    # Create a render window interactor
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # Initialize the interactor and start the rendering loop
    renderWindowInteractor.Initialize()
    renderWindow.Render()
    renderWindowInteractor.Start()