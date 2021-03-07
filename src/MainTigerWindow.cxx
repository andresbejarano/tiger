#define TINYOBJLOADER_IMPLEMENTATION

#include "MainTigerWindow.h"

#include <filesystem>
#include <set>

#include <QFileDialog>
#include <QMessageBox>
#include <QScreen>
#include <QStandardPaths>
#include <QVTKOpenGLNativeWidget.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellPicker.h>
#include <vtkColorTransferFunction.h>
#include <vtkCubeSource.h>
#include <vtkDataObjectToTable.h>
#include <vtkDoubleArray.h>
#include <vtkElevationFilter.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkLightsPass.h>
#include <vtkLine.h>
#include <vtkMatrix4x4.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOBJExporter.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProp3DCollection.h>
#include <vtkProperty.h>
#include <vtkQtTableView.h>
#include <vtkRendererDelegate.h>
#include <vtkRenderStepsPass.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkSSAOPass.h>
#include <vtkTextProperty.h>
#include <vtkTriangle.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVersion.h>
#include <vtkVRMLExporter.h>

#include "adaptivetileoffsetclippingparametersdialog.h"
#include "bendedsquareparametersdialog.h"
#include "centerdirectionparametersdialog.h"
#include "cyclideparametersdialog.h"
#include "cylinderparametersdialog.h"
#include "paraboloidparametersdialog.h"
#include "equilateraltriangleparametersdialog.h"
#include "equilibriumanalysisparametersdialog.h"
#include "facesubdivisionparametersdialog.h"
#include "heightbisectionparametersdialog.h"
#include "normalizetessellationparametersdialog.h"
#include "platonicsolidparametersdialog.h"
#include "rectangleparametersdialog.h"
#include "regularpolygonparametersdialog.h"
#include "rotatetessellationparametersdialog.h"
#include "scaletessellationparametersdialog.h"
#include "squAregridparametersdialog.h"
#include "squAreparametersdialog.h"
#include "tileoffsetclippingparametersdialog.h"
#include "tiltingangleparametersdialog.h"
#include "torusparametersdialog.h"
#include "truncatedconeparametersdialog.h"

#include "geometries.h"
#include "utils.h"
#include "tiny_obj_loader.h"

#if VTK_VERSION_NUMBER >= 89000000000ULL
#define VTK890 1
#endif

// 
const QString APP_TITLE = "TIGER";
//const QString NO_TESSELLATION_VISIBLE_TITLE = "Tessellation Not Visible";
const QString NO_INTERFACES_VISIBLE_CONTENT = "The interface polygons not visible. Do you want to view it?";
const QString NO_TESSELLATION_VISIBLE_CONTENT = "The tessellation is not visible. Do you want to view it?";

/*
Sets the UV vectors that span a given plane.
@param const QString & plane The reference to the plane indicator.
@param Eigen::Vector3d & U The reference to the U vector.
@param Eigen::Vector3d & V The reference to the V vector.
*/
void UVFromPlane(const QString & plane, Eigen::Vector3d & U, Eigen::Vector3d & V) 
{
    assert(plane == "XY" || plane == "XZ" || plane == "YZ");

    if (plane == "XY")
    {
        U << 1, 0, 0;
        V << 0, 1, 0;
    }
    else if (plane == "XZ")
    {
        U << 0, 0, 1;
        V << 1, 0, 0;
    }
    else 
    {
        U << 0, 0, 1;
        V << 0, 1, 0;
    }
}

vtkStandardNewMacro(TilePickerStyle);

MainTigerWindow::MainTigerWindow() : 
    m_currentInterfacePolygonsType(INTERFACES::PLAIN), 
    m_lastDirectory(QStandardPaths::locate(QStandardPaths::DocumentsLocation, QString(), QStandardPaths::LocateDirectory)),
    m_viewAssemblyBlocksBoundingBoxes(false), 
    m_viewAssemblyGeometry(true), 
	m_viewAxes(true), 
	m_viewEdgeDirections(true), 
	m_viewEdgeRotatedVectors(true), 
	m_viewTileCenters(true), 
	m_viewTessellationBoundingBox(false), 
    m_viewTessellationGeometry(true), 
	m_viewGrid(true), 
	m_viewInterfacePolygons(true), 
	m_viewSectionHeights(true), 
    m_vtkAssemblyBlocksBoundingBoxesActor(nullptr), 
    m_vtkAssemblyEdgesActor(nullptr), 
    m_vtkAssemblyFacesActor(nullptr), 
	m_vtkAxesActor(nullptr), 
	m_vtkCamera(nullptr), 
    m_vtkTileCentersActor(nullptr), 
    m_vtkEdgeDirectionsActor(nullptr), 
    m_vtkTessellationBoundBoxActor(nullptr), 
	m_vtkTessellationEdgesActor(nullptr), 
    m_vtkTessellationTilesActor(nullptr), 
	m_vtkGridActor(nullptr), 
    m_vtkInterfacePolygonsActor(nullptr), 
    m_vtkInterfacePolygonsBarActor(nullptr), 
	m_vtkRenderer(nullptr), 
	m_workspace(nullptr), 
    m_workspacesDir(QDir::currentPath() + "\\workspaces")

{
    // 
	this->setupUi(this);

    // Calculate 85% of the height from the primary screen. Then, set the width and height of the 
    // window
    double height = QGuiApplication::primaryScreen()->geometry().height() * 0.85;
    setFixedSize(height * (1.0 + sqrt(5.0)) / 2.0, height);

	// Update the viewing menu items
    this->actionViewAssemblyBoundingBoxes->setChecked(m_viewAssemblyBlocksBoundingBoxes);
	this->actionViewAssemblyGeometry->setChecked(m_viewAssemblyGeometry);
	this->actionViewAxes->setChecked(m_viewAxes);
	this->actionViewEdgeDirections->setChecked(m_viewEdgeDirections);
	this->actionViewTileCenters->setChecked(m_viewTileCenters);
	this->actionViewTessellationBoundingBox->setChecked(m_viewTessellationBoundingBox);
    this->actionViewTessellationGeometry->setChecked(m_viewTessellationGeometry);
	this->actionViewGrid->setChecked(m_viewGrid);

    switch (m_currentInterfacePolygonsType) 
    {
        case INTERFACES::PLAIN: 
        {
            this->actionViewPlainInterfaces->setChecked(m_viewInterfacePolygons);
            this->actionViewCompressionForces->setChecked(false);
            this->actionViewTensionForces->setChecked(false);
            this->actionViewUTangentialForces->setChecked(false);
            this->actionViewVTangentialForces->setChecked(false);
            break;
        }

        case INTERFACES::COMPRESSION:
        {
            this->actionViewPlainInterfaces->setChecked(false);
            this->actionViewCompressionForces->setChecked(m_viewInterfacePolygons);
            this->actionViewTensionForces->setChecked(false);
            this->actionViewUTangentialForces->setChecked(false);
            this->actionViewVTangentialForces->setChecked(false);
            break;
        }

        case INTERFACES::TENSION:
        {
            this->actionViewPlainInterfaces->setChecked(false);
            this->actionViewCompressionForces->setChecked(false);
            this->actionViewTensionForces->setChecked(m_viewInterfacePolygons);
            this->actionViewUTangentialForces->setChecked(false);
            this->actionViewVTangentialForces->setChecked(false);
            break;
        }

        case INTERFACES::UTANGENTIAL:
        {
            this->actionViewPlainInterfaces->setChecked(false);
            this->actionViewCompressionForces->setChecked(false);
            this->actionViewTensionForces->setChecked(false);
            this->actionViewUTangentialForces->setChecked(m_viewInterfacePolygons);
            this->actionViewVTangentialForces->setChecked(false);
            break;
        }

        case INTERFACES::VTANGENTIAL:
        {
            this->actionViewPlainInterfaces->setChecked(false);
            this->actionViewCompressionForces->setChecked(false);
            this->actionViewTensionForces->setChecked(false);
            this->actionViewUTangentialForces->setChecked(false);
            this->actionViewVTangentialForces->setChecked(m_viewInterfacePolygons);
            break;
        }
    }
	
	// 
	vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
    this->qvtkWidget->setRenderWindow(renderWindow);

    // 
    vtkNew<vtkNamedColors> colors;

	// Initialize and define the camera of the VTK scene
	m_vtkCamera = vtkCamera::New();
	m_vtkCamera->SetPosition(3, 3, 3);
	m_vtkCamera->SetFocalPoint(0, 0, 0);

    // 
    //vtkWeakPointer<vtkRenderStepsPass> renderSteps = vtkRenderStepsPass::New();

    // Initialize the Screen-Space Ambient Occlusion (SSAO) pass. It darkens some pixels to improve
    // depth perception
    //vtkWeakPointer<vtkSSAOPass> ssaoPass = vtkSSAOPass::New();
    //ssaoPass->SetDelegatePass(renderSteps);

    // Initialize the VTK renderer
	m_vtkRenderer = vtkRenderer::New();
	m_vtkRenderer->SetActiveCamera(m_vtkCamera);
	m_vtkRenderer->SetBackground(colors->GetColor3d("White").GetData());
    //m_vtkRenderer->SetPass(ssaoPass);

	// VTK/Qt wedded
    this->qvtkWidget->renderWindow()->AddRenderer(m_vtkRenderer);
    this->qvtkWidget->renderWindow()->LineSmoothingOn();
    this->qvtkWidget->renderWindow()->PolygonSmoothingOn();

	// Set up action signals and slots
	//connect(this->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));

	// Initialize the axes and grid actors, they start initialized by default
	InitAxesActor();
	InitGridActor();

    // Validate the workspaces directory, make it if it doesn't exist
    assert(ValidateWorkspacesDirectory(true));
    //std::cout << CountWorkspaces() << std::endl;

    // Start a new workspace
    m_workspace = std::make_shared<Workspace>();
}

MainTigerWindow::~MainTigerWindow() 
{
	ClearWorkspace();

    if (m_vtkAxesActor)
    {
        m_vtkRenderer->RemoveActor(m_vtkAxesActor);
        m_vtkAxesActor->Delete();
    }

    if (m_vtkGridActor)
    {
        m_vtkRenderer->RemoveActor(m_vtkGridActor);
        m_vtkGridActor->Delete();
    }

    //m_tilePicker->Delete();
}

bool MainTigerWindow::AreAxesVisible() const
{
	return m_viewAxes;
}

bool MainTigerWindow::AreEdgeDirectionsVisible() const
{
	return m_viewEdgeDirections;
}

bool MainTigerWindow::AreEdgeRotatedVectorsVisible() const
{
	return m_viewEdgeRotatedVectors;
}

bool MainTigerWindow::AreFaceCentersVisible() const
{
	return m_viewTileCenters;
}

bool MainTigerWindow::AreInterfacePolygonsVisible() const
{
	return m_viewInterfacePolygons;
}

bool MainTigerWindow::AreSectionHeightsVisible() const
{
	return m_viewSectionHeights;
}

void MainTigerWindow::AskViewInterfacePolygons()
{
    // Ask whether to show the interface polygons or not. If yes then show it
    if (QMessageBox::question(
        this,
        APP_TITLE, 
        NO_INTERFACES_VISIBLE_CONTENT,
        QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes)
    {
        ViewInterfacePolygons(m_currentInterfacePolygonsType);
    }
}

void MainTigerWindow::ClearWorkspace()
{
	if (m_workspace) 
	{
		m_workspace->Clear();
		m_workspace = nullptr;
	}

    DeleteWorkspaceActors();
}

void MainTigerWindow::AskViewTessellationGeometry()
{
    // Ask whether to show the tessellation or not. If yes then show it
    if (QMessageBox::question(
        this,
        APP_TITLE,
        NO_TESSELLATION_VISIBLE_CONTENT,
        QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes)
    {
        ViewTessellationGeometry();
    }
}

void MainTigerWindow::ClearAssemblyActors()
{
    if (m_vtkAssemblyFacesActor)
    {
        m_vtkRenderer->RemoveActor(m_vtkAssemblyFacesActor);
        m_vtkAssemblyFacesActor->Delete();
    }

    if (m_vtkAssemblyEdgesActor) 
    {
        m_vtkRenderer->RemoveActor(m_vtkAssemblyEdgesActor);
        m_vtkAssemblyEdgesActor->Delete();
    }
}

void MainTigerWindow::ClearTessellationActors()
{
    if (m_vtkTessellationBoundBoxActor)
    {
        m_vtkRenderer->RemoveActor(m_vtkTessellationBoundBoxActor);
        m_vtkTessellationBoundBoxActor->Delete();
    }

    if (m_vtkTessellationEdgesActor) 
    {
        m_vtkRenderer->RemoveActor(m_vtkTessellationEdgesActor);
        m_vtkTessellationEdgesActor->Delete();
    }

    if (m_vtkTessellationTilesActor)
    {
        m_vtkRenderer->RemoveActor(m_vtkTessellationTilesActor);
        m_vtkTessellationTilesActor->Delete();
    }
}

void MainTigerWindow::ClearInteracePolygonsActor()
{
    if (m_vtkInterfacePolygonsActor)
    {
        m_vtkRenderer->RemoveActor(m_vtkInterfacePolygonsActor);
        m_vtkInterfacePolygonsActor->Delete();
    }

    if (m_vtkInterfacePolygonsBarActor) 
    {
        m_vtkRenderer->RemoveActor2D(m_vtkInterfacePolygonsBarActor);
        m_vtkInterfacePolygonsBarActor->Delete();
    }
}

size_t MainTigerWindow::CountWorkspaces() const
{
    assert(ValidateWorkspacesDirectory());

    QDir dir(m_workspacesDir);
    dir.setFilter(QDir::Dirs | QDir::NoDotAndDotDot);
    return dir.count();
}

void MainTigerWindow::DeleteWorkspaceActors()
{
    ClearTessellationActors();

    if (m_vtkTileCentersActor) 
    {
        m_vtkRenderer->RemoveActor(m_vtkTileCentersActor);
        m_vtkTileCentersActor->Delete();
    }

    if (m_vtkEdgeDirectionsActor) 
    {
        m_vtkRenderer->RemoveActor(m_vtkEdgeDirectionsActor);
        m_vtkEdgeDirectionsActor->Delete();
    }

    ClearAssemblyActors();
    ClearInteracePolygonsActor();
}

/*vtkIdType MainTigerWindow::GetPickedTileCellId(int x, int y) const
{
    assert(x >= 0 && y >= 0);

    // Return no cell id if the tessellation is not visible or there is no tessellation actor
    if (!m_viewTessellation || !m_vtkTessellationTilesActor) 
    {
        return -1;
    }

    // Initialize and define the cell picker. Then check what it picks
    vtkSmartPointer<vtkCellPicker> cellPicker = vtkSmartPointer<vtkCellPicker>::New();
    cellPicker->SetTolerance(0.0005);
    cellPicker->Pick(x, y, 0, m_tilePicker->GetDefaultRenderer());

    // Get the pointer to the collection of actors picked by the cell picker
    //vtkProp3DCollection * pickedActors = cellPicker->GetProp3Ds();

    // Check if the cell picker intersected with the tessellation actor
    //if (pickedActors->IsItemPresent(m_vtkTessellationTilesActor)) 
    //{
    //    std::cout << "Click on the tessellation" << std::endl;
    //}

    return 1;
}*/

std::string MainTigerWindow::GetInterfaceTypeName(INTERFACES type) const
{
    switch (type) 
    {
        case INTERFACES::PLAIN:         return "Plain";         break;
        case INTERFACES::COMPRESSION:   return "Compression";   break;
        case INTERFACES::TENSION:       return "Tension";       break;
        case INTERFACES::UTANGENTIAL:   return "U-Tangential";  break;
        case INTERFACES::VTANGENTIAL:   return "V-Tangential";  break;
        default: return "Unknown";
    }
}

void MainTigerWindow::InitAssemblyActors(
    const std::shared_ptr<Tessellation> domain, 
    const std::shared_ptr<Assembly> assembly)
{
    assert(domain);
    assert(assembly);

    // Clear the assembly actors
    ClearAssemblyActors();

    size_t nAssemblyVertices = assembly->CountVertices();

    // Initialize the object to store the points (vertices) of the pieces
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(nAssemblyVertices);

    // 
    vtkSmartPointer<vtkDoubleArray> facePointColors = vtkSmartPointer<vtkDoubleArray>::New();
        //edgePointColors = vtkSmartPointer<vtkDoubleArray>::New();
    facePointColors->SetNumberOfValues(nAssemblyVertices);
    //edgePointColors->SetNumberOfValues(nAssemblyVertices);

    vtkSmartPointer<vtkCellArray> edges = vtkSmartPointer<vtkCellArray>::New();
    edges->SetNumberOfCells(assembly->CountEdges());

    std::set<std::tuple<size_t, size_t>> visitedEdges;

    // Initialize the object to store the triangles representing the faces of the pieces. Each 
    // face is triangulated with respect of one of its incident vertices
    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    triangles->SetNumberOfCells(assembly->CountTriangles());

    // Keep count of the number of points that have been inserted in the points object
    size_t insertedPoints = 0;

    // 
    Eigen::Vector3d C = Eigen::Vector3d::Zero();

    // 
    double P[3] = { 0.0, 0.0, 0.0 };

    std::shared_ptr<Block> block = nullptr;

    // Traverse through the pieces of the assembly
    for (auto it = assembly->GetBlocks().begin(); it != assembly->GetBlocks().end(); ++it) 
    {
        // Get the pointer to the current block
        block = *it;

        //size_t faceIndex;
        //assert(block->Attributes().Get<size_t>(ATTRIB_FACE_INDEX, faceIndex));

        //bool allNeighbors = domain->DCEL()->Faces()[block->FaceIndex()]->HasAllNeighbors();

        double facePointColor = block->IsEnabled() ? 1.0 : 0.2;
        ////double edgePointColor = allNeighbors ? 0.0 : 1.0;
        //double edgePointColor = 1.0;

        // Get the number of vertices of the current piece
        size_t nVertices = block->Geometry().countVertices();

        // Traverse through the vertices of the current piece, insert their coordinates into the 
        // points object
        for (size_t vIdx = 0; vIdx < nVertices; vIdx += 1) 
        {
            C << block->Geometry().Vertex(vIdx);

            P[0] = C.x();
            P[1] = C.y();
            P[2] = C.z();

            points->InsertPoint(insertedPoints + vIdx, P);

            facePointColors->SetValue(insertedPoints + vIdx, facePointColor);
            //edgePointColors->SetValue(insertedPoints + vIdx, edgePointColor);
        }

        // Get the number of faces of the current piece
        size_t nFaces = block->Geometry().countFaces();

        // Traverse through the faces of the current piece, triangulate them and insert into the 
        // triangles object
        for (size_t fIdx = 0; fIdx < nFaces; fIdx += 1) 
        {
            const std::vector<size_t> & indices = block->Geometry().face(fIdx);

            size_t nIndices = indices.size();

            for (size_t vIdx = 1; vIdx < nIndices - 1; vIdx += 1) 
            {
                vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                triangle->GetPointIds()->SetId(0, insertedPoints + indices[0]);
                triangle->GetPointIds()->SetId(1, insertedPoints + indices[vIdx]);
                triangle->GetPointIds()->SetId(2, insertedPoints + indices[vIdx + 1]);

                triangles->InsertNextCell(triangle);
            }

            size_t i, j;

            // Traverse the indices and define the edges of the current piece
            for (size_t vIdx = 0; vIdx < nIndices - 1; vIdx += 1) 
            {
                i = insertedPoints + indices[vIdx], j = insertedPoints + indices[vIdx + 1];

                std::tuple<size_t, size_t> edge = std::make_tuple(std::min(i, j), std::max(i, j));

                if (visitedEdges.find(edge) == visitedEdges.end()) 
                {
                    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                    line->GetPointIds()->SetId(0, std::get<0>(edge));
                    line->GetPointIds()->SetId(1, std::get<1>(edge));

                    edges->InsertNextCell(line);

                    visitedEdges.insert(edge);
                }
            }

            // Check the last edge connecting the first index with the last one
            i = insertedPoints + indices[0], j = insertedPoints + indices[nIndices - 1];

            std::tuple<size_t, size_t> edge = std::make_tuple(std::min(i, j), std::max(i, j));

            if (visitedEdges.find(edge) == visitedEdges.end())
            {
                vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                line->GetPointIds()->SetId(0, std::get<0>(edge));
                line->GetPointIds()->SetId(1, std::get<1>(edge));

                edges->InsertNextCell(line);

                visitedEdges.insert(edge);
            }
        }

        // Update the number of inserted points
        insertedPoints += nVertices;
    }

    vtkSmartPointer<vtkColorTransferFunction> colorFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
    colorFunction->AddRGBPoint(0.2, 0.2, 0.2, 0.2); // For the pieces at the boundary
    colorFunction->AddRGBPoint(1.0, 1.0, 1.0, 1.0); // For the internal pieces

    // 
    vtkSmartPointer<vtkPolyData> facesPolyData = vtkSmartPointer<vtkPolyData>::New();
    facesPolyData->SetPoints(points);
    facesPolyData->SetPolys(triangles);
    facesPolyData->GetPointData()->SetScalars(facePointColors);
    
    vtkSmartPointer<vtkPolyDataMapper> facesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    facesMapper->SetInputData(facesPolyData);
    facesMapper->SetLookupTable(colorFunction);

    // Initialize and define the actor that renders the faces of the blocks in the assembly
    m_vtkAssemblyFacesActor = vtkActor::New();
    m_vtkAssemblyFacesActor->SetMapper(facesMapper);
    m_vtkAssemblyFacesActor->GetProperty()->BackfaceCullingOn();
	m_vtkAssemblyFacesActor->GetProperty()->SetOpacity(1.0);
    m_vtkAssemblyFacesActor->SetVisibility(m_viewAssemblyGeometry);
    
    // 
    vtkSmartPointer<vtkPolyData> edgesPolyData = vtkSmartPointer<vtkPolyData>::New();
    edgesPolyData->SetPoints(points);
    edgesPolyData->SetLines(edges);

    vtkSmartPointer<vtkPolyDataMapper> edgesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    edgesMapper->SetInputData(edgesPolyData);

    // Initialize and define the actor that renders the edges of the blocks in the assembly
    m_vtkAssemblyEdgesActor = vtkActor::New();
    m_vtkAssemblyEdgesActor->SetMapper(edgesMapper);
    m_vtkAssemblyEdgesActor->GetProperty()->SetLineWidth(4);
    m_vtkAssemblyEdgesActor->GetProperty()->SetColor(0.3, 0.3, 0.3);
    m_vtkAssemblyEdgesActor->SetVisibility(m_viewAssemblyGeometry);

    // Add the assembly actors to the renderer
    m_vtkRenderer->AddActor(m_vtkAssemblyFacesActor);
    m_vtkRenderer->AddActor(m_vtkAssemblyEdgesActor);

    std::cout << "#Blocks = " << assembly->CountBlocks() << std::endl;
}

void MainTigerWindow::InitAxesActor()
{
	// Define the coordinates of the points for generating the axes
	const double O[3] = { 0.0, 0.0, 0.0 };
	const double X[3] = { 1.0, 0.0, 0.0 };
	const double Y[3] = { 0.0, 1.0, 0.0 };
	const double Z[3] = { 0.0, 0.0, 1.0 };

	// Initialize the object to store the points for the reference axes
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(4);
	points->SetPoint(0, O);
	points->SetPoint(1, X);
	points->SetPoint(2, Y);
	points->SetPoint(3, Z);

	// Define the line representing the X axis
	vtkSmartPointer<vtkLine> xAxis = vtkSmartPointer<vtkLine>::New();
	xAxis->GetPointIds()->SetId(0, 0);
	xAxis->GetPointIds()->SetId(1, 1);

	// Define the line representing the Y axis
	vtkSmartPointer<vtkLine> yAxis = vtkSmartPointer<vtkLine>::New();
	yAxis->GetPointIds()->SetId(0, 0);
	yAxis->GetPointIds()->SetId(1, 2);

	// Define the line representing the Z axis
	vtkSmartPointer<vtkLine> zAxis = vtkSmartPointer<vtkLine>::New();
	zAxis->GetPointIds()->SetId(0, 0);
	zAxis->GetPointIds()->SetId(1, 3);

	// Initialize the object to store the lines and insert the axis lines
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    lines->SetNumberOfCells(3);
	lines->InsertNextCell(xAxis);
	lines->InsertNextCell(yAxis);
	lines->InsertNextCell(zAxis);

	// Define the colors for the axes
	const unsigned char R[3] = { 255, 0, 0 };
	const unsigned char G[3] = { 0, 255, 0 };
	const unsigned char B[3] = { 0, 0, 255 };

	// Initialize the object to store the colors and insert them
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetNumberOfTuples(3);
	colors->SetTypedTuple(0, R);
	colors->SetTypedTuple(1, G);
	colors->SetTypedTuple(2, B);

	// Initialize and define the lines poly data (this object links points, lines and colors)
	vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
	linesPolyData->SetPoints(points);
	linesPolyData->SetLines(lines);
	linesPolyData->GetCellData()->SetScalars(colors);

	// Initialize and define the lines mapper
	vtkSmartPointer<vtkPolyDataMapper> linesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	linesMapper->SetInputData(linesPolyData);

	// Initialize and define the axes actor
	m_vtkAxesActor = vtkActor::New();
	m_vtkAxesActor->SetMapper(linesMapper);
	m_vtkAxesActor->GetProperty()->SetLineWidth(4);
    m_vtkAxesActor->SetVisibility(m_viewAxes);

    // Add the axes actor to the scene
    m_vtkRenderer->AddActor(m_vtkAxesActor);
}

vtkWeakPointer<vtkActor> MainTigerWindow::InitAxisAlignedBoundingBoxActor(const VF & vf)
{
    Eigen::Vector3d min = Eigen::Vector3d::Zero(), max = Eigen::Vector3d::Zero();

    vf.axisAlignedBoundingBox(min, max);

    // Define the points of the bounding box
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(8);
    points->SetPoint(0, min.x(), min.y(), min.z());
    points->SetPoint(1, min.x(), min.y(), max.z());
    points->SetPoint(2, max.x(), min.y(), max.z());
    points->SetPoint(3, max.x(), min.y(), min.z());
    points->SetPoint(4, min.x(), max.y(), min.z());
    points->SetPoint(5, min.x(), max.y(), max.z());
    points->SetPoint(6, max.x(), max.y(), max.z());
    points->SetPoint(7, max.x(), max.y(), min.z());

    // Define the edges of the bounding box
    vtkSmartPointer<vtkLine> line01 = vtkSmartPointer<vtkLine>::New();
    line01->GetPointIds()->SetId(0, 0);
    line01->GetPointIds()->SetId(1, 1);

    vtkSmartPointer<vtkLine> line12 = vtkSmartPointer<vtkLine>::New();
    line12->GetPointIds()->SetId(0, 1);
    line12->GetPointIds()->SetId(1, 2);

    vtkSmartPointer<vtkLine> line23 = vtkSmartPointer<vtkLine>::New();
    line23->GetPointIds()->SetId(0, 2);
    line23->GetPointIds()->SetId(1, 3);

    vtkSmartPointer<vtkLine> line30 = vtkSmartPointer<vtkLine>::New();
    line30->GetPointIds()->SetId(0, 3);
    line30->GetPointIds()->SetId(1, 0);

    vtkSmartPointer<vtkLine> line45 = vtkSmartPointer<vtkLine>::New();
    line45->GetPointIds()->SetId(0, 4);
    line45->GetPointIds()->SetId(1, 5);

    vtkSmartPointer<vtkLine> line56 = vtkSmartPointer<vtkLine>::New();
    line56->GetPointIds()->SetId(0, 5);
    line56->GetPointIds()->SetId(1, 6);

    vtkSmartPointer<vtkLine> line67 = vtkSmartPointer<vtkLine>::New();
    line67->GetPointIds()->SetId(0, 6);
    line67->GetPointIds()->SetId(1, 7);

    vtkSmartPointer<vtkLine> line74 = vtkSmartPointer<vtkLine>::New();
    line74->GetPointIds()->SetId(0, 7);
    line74->GetPointIds()->SetId(1, 4);

    vtkSmartPointer<vtkLine> line04 = vtkSmartPointer<vtkLine>::New();
    line04->GetPointIds()->SetId(0, 0);
    line04->GetPointIds()->SetId(1, 4);

    vtkSmartPointer<vtkLine> line15 = vtkSmartPointer<vtkLine>::New();
    line15->GetPointIds()->SetId(0, 1);
    line15->GetPointIds()->SetId(1, 5);

    vtkSmartPointer<vtkLine> line26 = vtkSmartPointer<vtkLine>::New();
    line26->GetPointIds()->SetId(0, 2);
    line26->GetPointIds()->SetId(1, 6);

    vtkSmartPointer<vtkLine> line37 = vtkSmartPointer<vtkLine>::New();
    line37->GetPointIds()->SetId(0, 3);
    line37->GetPointIds()->SetId(1, 7);

    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    lines->SetNumberOfCells(12);
    lines->InsertNextCell(line01);
    lines->InsertNextCell(line12);
    lines->InsertNextCell(line23);
    lines->InsertNextCell(line30);
    lines->InsertNextCell(line45);
    lines->InsertNextCell(line56);
    lines->InsertNextCell(line67);
    lines->InsertNextCell(line74);
    lines->InsertNextCell(line04);
    lines->InsertNextCell(line15);
    lines->InsertNextCell(line26);
    lines->InsertNextCell(line37);

    // Initialize and define the lines poly data (this object links points, lines and colors)
    vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
    linesPolyData->SetPoints(points);
    linesPolyData->SetLines(lines);
    //linesPolyData->GetCellData()->SetScalars(colors);

    // Initialize and define the lines mapper
    vtkSmartPointer<vtkPolyDataMapper> linesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    linesMapper->SetInputData(linesPolyData);

    // Define and initialize the actor with the information of the bounding box
    vtkWeakPointer<vtkActor> boundingBoxActor = vtkActor::New();
    boundingBoxActor->SetMapper(linesMapper);
    boundingBoxActor->GetProperty()->SetLineWidth(1);
    boundingBoxActor->GetProperty()->SetColor(1.0, 0.0, 0.0);

    return boundingBoxActor;
}

void MainTigerWindow::InitEdgeDirectionsActor(const std::shared_ptr<dcel::DCEL> dcel, double length)
{
    assert(dcel);

    // If there is an edge directions actor then remove it from the renderer (if viewing) and 
    // delete it
    if (m_vtkEdgeDirectionsActor) 
    {
        m_vtkRenderer->RemoveActor(m_vtkEdgeDirectionsActor);
        m_vtkEdgeDirectionsActor->Delete();
    }

    size_t nHalfedges = dcel->Halfedges().size();

    // Initialize the object to store the endpoints of the edge direction lines
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(nHalfedges);

    // Initialize an array to store the coordinates of the vertices and generate the points
    double P[3] = { 0.0, 0.0, 0.0 };

    // 
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    lines->SetNumberOfCells(nHalfedges / 2);

    // Set the half edges of the tessellation as not visited
    dcel->SetHalfedgesAttribute<bool>(ATTRIB_VISITED, false);

    bool visited = false;

    Eigen::Vector3d M = Eigen::Vector3d::Zero(), D = Eigen::Vector3d::Zero();

    size_t insertedPoints = 0;

    for (size_t heIdx = 0; heIdx < nHalfedges; heIdx += 1) 
    {
        // Get the pointer to the current half edge
        std::shared_ptr<dcel::Halfedge> halfedge = dcel->Halfedges()[heIdx];

        // Check if the current half edge has been visited
        assert(halfedge->Attributes().Get<bool>(ATTRIB_VISITED, visited));

        // If the current half edge has been visited then continue with the next half edge
        if (visited) 
        {
            continue;
        }

        // Get the coordinates of the Midpoint of the half edge
        M = halfedge->Midpoint();

        // Get the direction vector of the half edge
        assert(halfedge->Attributes().Get<Eigen::Vector3d>(ATTRIB_DIRECTION_VECTOR, D));
        D.normalize();
        D *= length;

        // Load the coordinates of the first end point and insert it into the points object
        P[0] = M.x();
        P[1] = M.y();
        P[2] = M.z();
        points->InsertPoint(insertedPoints + 0, P);

        // Load the coordinates of the second end point and insert it into the points object
        P[0] = M.x() + D.x();
        P[1] = M.y() + D.y();
        P[2] = M.z() + D.z();
        points->InsertPoint(insertedPoints + 1, P);

        // 
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, insertedPoints + 0);
        line->GetPointIds()->SetId(1, insertedPoints + 1);
        lines->InsertNextCell(line);

        // Update the number of points that have been inserted in the points object
        insertedPoints += 2;
    }

    // Remove the ATTRIB_VISITED dynamic attribute from the half edges of the tessellation
    dcel->RemoveHalfedgesAttribute(ATTRIB_VISITED);

    // Initialize and define the polydata
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);

    // Initialize and define the mapper
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polyData);

    // Initialize and define the edge directions actor
    m_vtkEdgeDirectionsActor = vtkActor::New();
    m_vtkEdgeDirectionsActor->SetMapper(mapper);
    m_vtkEdgeDirectionsActor->GetProperty()->SetColor(1.0, 0.0, 1.0);
    m_vtkEdgeDirectionsActor->GetProperty()->SetLineWidth(4);
    m_vtkEdgeDirectionsActor->SetVisibility(m_viewEdgeDirections);

    // Add the actor to the renderer
    m_vtkRenderer->AddActor(m_vtkEdgeDirectionsActor);
}

void MainTigerWindow::InitTileCentersActor(const std::shared_ptr<dcel::DCEL> dcel, double radius)
{
    assert(dcel);

    // If there is a face centers actor then remove it from the renderer (if viewing) and  delete 
    // it
    if (m_vtkTileCentersActor)
    {
        m_vtkRenderer->RemoveActor(m_vtkTileCentersActor);
        m_vtkTileCentersActor->Delete();
    }

    // Get the vertex coordinates and vertex indices of an octahedron with the given radius
    const VF OCT = geometries::Octahedron(radius);

    size_t nFaces = dcel->Faces().size();

    // Initialize the object to store the points (vertices) of the face centers. Face centers 
    // Are represented using octahedra, an octahedron has 6 vertices
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(nFaces * 6);

    // Initialize the object to store the triangles representing the face centers
    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    triangles->SetNumberOfCells(nFaces * 8);

    // Initialize an array to store the coordinates of the vertices and generate the points
    double P[3] = { 0.0, 0.0, 0.0 };

    Eigen::Vector3d V = Eigen::Vector3d::Zero(), C = Eigen::Vector3d::Zero();

    // Traverse through the faces of the tessellation and store their respective octahedron 
    // vertices coordinates in the points object
    for (size_t faceIdx = 0; faceIdx < nFaces; faceIdx += 1)
    {
        // Get the reference to the center of the current face
        assert(dcel->Faces()[faceIdx]->Attributes().Get<Eigen::Vector3d>(ATTRIB_CENTER, C));

        // Calculate the start index for the points of the octahedron for the current face center
        size_t octIdx = faceIdx * 6;

        // Traverse through the vertices of the octahedron. Translate its coordinates with respect 
        // to the location of the face center and add them to the points object
        for (size_t vIdx = 0; vIdx < 6; vIdx += 1)
        {
            V << OCT.Vertex(vIdx);

            //P[0] = C.x() + OCT.V[vIdx].x();
            //P[1] = C.y() + OCT.V[vIdx].y();
            //P[2] = C.z() + OCT.V[vIdx].z();

            P[0] = C.x() + V.x();
            P[1] = C.y() + V.y();
            P[2] = C.z() + V.z();

            points->InsertPoint(octIdx + vIdx, P);
        }

        // Traverse through the vertex indices of the octahedron. Define the faces and insert them 
        // in the triangles object
        for (size_t fIdx = 0; fIdx < 8; fIdx += 1)
        {
            const std::vector<size_t> & indices = OCT.face(fIdx);

            vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
            //triangle->GetPointIds()->SetId(0, octIdx + OCT.F[fIdx][0]);
            //triangle->GetPointIds()->SetId(1, octIdx + OCT.F[fIdx][1]);
            //triangle->GetPointIds()->SetId(2, octIdx + OCT.F[fIdx][2]);
            triangle->GetPointIds()->SetId(0, octIdx + indices[0]);
            triangle->GetPointIds()->SetId(1, octIdx + indices[1]);
            triangle->GetPointIds()->SetId(2, octIdx + indices[2]);

            triangles->InsertNextCell(triangle);
        }
    }

    // Initialize and define the polydata
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->SetPolys(triangles);

    // Initialize and define the mapper
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polydata);

    // Initialize and define the face centers actor
    m_vtkTileCentersActor = vtkActor::New();
    m_vtkTileCentersActor->SetMapper(mapper);
    m_vtkTileCentersActor->GetProperty()->BackfaceCullingOn();
    m_vtkTileCentersActor->GetProperty()->SetColor(1.0, 0.0, 1.0);
    m_vtkTileCentersActor->SetVisibility(m_viewTileCenters);

    // Add the actor to the renderer
    m_vtkRenderer->AddActor(m_vtkTileCentersActor);
}

/*void MainTigerWindow::InitTessellationActor(const dcel::DCEL & domain)
{
	// If there is a tessellation actor then remove it from the renderer (if viewing) and 
	// delete it
	if (m_vtkTessellationActor) 
	{
		if (m_viewTessellation) 
		{
			m_vtkRenderer->RemoveActor(m_vtkTessellationActor);
		}
        
        m_vtkTessellationActor->Delete();
	}

	// Initialize the object to store the points (vertices) of the tessellation
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(domain.vertices.size());

	// Initialize an array to store the coordinates of the vertices and generate the points
	double P[3] = { 0.0, 0.0, 0.0 };

	// Traverse through the vertices of the tessellation and store their coordinates in the 
	// points object
	for (size_t i = 0; i < domain.vertices.size(); i += 1) 
	{
		// Get the reference to the coordinates of the current vertex
		const Eigen::Vector3d & coords = domain.vertices[i]->Coords();

		// Set the coordinates of the point, then insert it into the points object
		P[0] = coords.x();
		P[1] = coords.y();
		P[2] = coords.z();
		points->InsertPoint(i, P);

		// Set the index dynamic attribute of the current vertex
		domain.vertices[i]->Attributes().Set<size_t>(ATTRIB_INDEX, i);
	}

    // Initialize the object to store the lines representing the wireframe of the geometric
    // domain
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

    // Initialize the visited dynamic attribute of the half edges of the tessellation as false
    domain.SetHalfedgesAttribute(ATTRIB_VISITED, false);

    // Traverse through the half edges of the tessellation and add their information to the 
    // lines object
    for (auto hIt = domain.halfedges.begin(); hIt != domain.halfedges.end(); ++hIt) 
    {
        // Get the pointer to the current half edge
        std::shared_ptr<dcel::Halfedge> halfedge = *hIt;

        // Get the visited dynamic attribute of the half edge
        bool visited;
        assert(halfedge->Attributes().Get<bool>(ATTRIB_VISITED, visited));

        // If the current half edge has been visited then continue with the next half edge
        if (visited) 
        {
            continue;
        }

        // 
        size_t v0, v1;
        assert(halfedge->start->Attributes().Get<size_t>(ATTRIB_INDEX, v0));
        assert(halfedge->twin->start->Attributes().Get<size_t>(ATTRIB_INDEX, v1));

        // 
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, v0);
        line->GetPointIds()->SetId(1, v1);

        // Insert the line into the lines object
        lines->InsertNextCell(line);

        // Set the current half edge and its twin as visited
        halfedge->Attributes().Set<bool>(ATTRIB_VISITED, true);
        halfedge->twin->Attributes().Set<bool>(ATTRIB_VISITED, true);
    }

    // Remove the visited dynamic attribute from the half edges of the tessellation
    domain.RemoveHalfedgesAttribute(ATTRIB_VISITED);

	// Initialize the object to store the triangles representing the faces of the geometric 
	// domain. Each face is triangulated with respect of the start vertex from the incident half 
	// edge of each face
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();

	// Traverse through the faces of the tessellation, generate the respective triangles and 
	// insert them in the triangles object
	for (auto it = domain.faces.begin(); it != domain.faces.end(); ++it) 
	{
		// Get the pointer to the current face
		std::shared_ptr<dcel::Face> face = *it;

		// Get the index of the start coordinate of the incident half edge of the face. This is the
		// common vertex for all triangles that represent the face
		size_t v0;
		assert(face->halfedge->start->Attributes().Get<size_t>(ATTRIB_INDEX, v0));

		// Get the pointer to the next half edge to the incident half edge of the face
		std::shared_ptr<dcel::Halfedge> currentHalfedge = face->halfedge->next;

		// Traverse through the half edges of the face. Stop at the previous half edge to the 
		// incident half edge of the face
		do 
		{
			// Get the indices of the end points of the current half edge
			size_t v1, v2;
			assert(currentHalfedge->start->Attributes().Get<size_t>(ATTRIB_INDEX, v1));
			assert(currentHalfedge->twin->start->Attributes().Get<size_t>(ATTRIB_INDEX, v2));

			// Initialize the current triangle and define the indices of its points
			vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
			triangle->GetPointIds()->SetId(0, v0);
			triangle->GetPointIds()->SetId(1, v1);
			triangle->GetPointIds()->SetId(2, v2);

			// Insert the current triangle into the triangles object
			triangles->InsertNextCell(triangle);

			// Move to the next half edge of the face
			currentHalfedge = currentHalfedge->next;

		} while (currentHalfedge != face->halfedge->previous);
	}

	// Remove the index dynamic attribute from the vertices of the tessellation
	domain.RemoveVerticesAttribute(ATTRIB_INDEX);

	// Initialize and define the polydata
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	polydata->SetPolys(triangles);
    polydata->SetLines(lines);

	// Initialize and define the mapper
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(polydata);

	// Initialize and define the tessellation actor
	m_vtkTessellationActor = vtkActor::New();
	m_vtkTessellationActor->SetMapper(mapper);
    m_vtkTessellationActor->GetProperty()->BackfaceCullingOn();
    m_vtkTessellationActor->GetProperty()->SetLineWidth(4);

	// Add the tessellation actor to the renderer (if viewing)
	if (m_viewTessellation) 
	{
		m_vtkRenderer->AddActor(m_vtkTessellationActor);
	}
}*/

void MainTigerWindow::InitTessellationActors(const VF & vf)
{
    // 
    ClearTessellationActors();

    size_t nVertices = vf.countVertices();

    // Initialize the object to store the points (vertices) of the tessellation
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(nVertices);

    // Initialize an array to store the coordinates of the vertices and generate the points
    double P[3] = { 0.0, 0.0, 0.0 };

    // 
    Eigen::Vector3d C = Eigen::Vector3d::Zero();

    // Traverse through the vertices of the tessellation and store their coordinates in the 
    // points object
    for (size_t i = 0; i < nVertices; i += 1)
    {
        C << vf.Vertex(i);
        
        P[0] = C.x();
        P[1] = C.y();
        P[2] = C.z();

        points->InsertPoint(i, P);
    }

    // Initialize an unordered set to track the edges
    std::set<std::tuple<size_t, size_t>> edges;

    // Initialize the object to store the lines representing the wireframe of the geometric
    // domain
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    lines->SetNumberOfCells(vf.countEdges());

    // Initialize the object to store the triangles representing the faces of the geometric 
    // domain. Each face is triangulated with respect of the start vertex from the incident half 
    // edge of each face
    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    triangles->SetNumberOfCells(vf.countTriangles());

    size_t nFaces = vf.countFaces();

    // Traverse through the edges of the tessellation and add their information to the lines 
    // object
    for (size_t i = 0; i < nFaces; i += 1)
    {
        const std::vector<size_t> & face = vf.face(i);

        size_t nFaceIndices = face.size() - 1;

        for (size_t j = 0; j < nFaceIndices; j += 1) 
        {
            std::tuple<size_t, size_t> key(std::min(face[j], face[j + 1]), std::max(face[j], face[j + 1]));

            if (edges.find(key) != edges.end()) 
            {
                continue;
            }

            // 
            vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
            line->GetPointIds()->SetId(0, face[j]);
            line->GetPointIds()->SetId(1, face[j + 1]);

            // Insert the line into the lines object
            lines->InsertNextCell(line);

            edges.insert(key);
        }

        // Check for the last edge
        std::tuple<size_t, size_t> key(std::min(face[0], face[nFaceIndices]), std::max(face[0], face[nFaceIndices]));

        if (edges.find(key) == edges.end())
        {
            // 
            vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
            line->GetPointIds()->SetId(0, face[0]);
            line->GetPointIds()->SetId(1, face[nFaceIndices]);

            lines->InsertNextCell(line);

            edges.insert(key);
        }
        
        for (size_t j = 1; j < nFaceIndices; j += 1) 
        {
            vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
            triangle->GetPointIds()->SetId(0, face[0]);
            triangle->GetPointIds()->SetId(1, face[j]);
            triangle->GetPointIds()->SetId(2, face[j + 1]);

            triangles->InsertNextCell(triangle);
        }
    }

    // Initialize and define the polydata for the edges of the tessellation
    vtkSmartPointer<vtkPolyData> edgesPolydata = vtkSmartPointer<vtkPolyData>::New();
    edgesPolydata->SetPoints(points);
    edgesPolydata->SetLines(lines);

    vtkSmartPointer<vtkPolyDataMapper> edgesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    edgesMapper->SetInputData(edgesPolydata);

    // 
    m_vtkTessellationEdgesActor = vtkActor::New();
    m_vtkTessellationEdgesActor->SetMapper(edgesMapper);
    m_vtkTessellationEdgesActor->GetProperty()->SetLineWidth(4);
    m_vtkTessellationEdgesActor->GetProperty()->SetColor(0.3, 0.3, 0.3);
    m_vtkTessellationEdgesActor->SetVisibility(m_viewTessellationGeometry);

    // Initialize and define the polydata for the faces of the tessellation
    vtkSmartPointer<vtkPolyData> facesPolydata = vtkSmartPointer<vtkPolyData>::New();
    facesPolydata->SetPoints(points);
    facesPolydata->SetPolys(triangles);
    //polydata->SetLines(lines);

    // Initialize and define the mapper
    vtkSmartPointer<vtkPolyDataMapper> facesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    facesMapper->SetInputData(facesPolydata);

    // Initialize and define the tessellation actor
    m_vtkTessellationTilesActor = vtkActor::New();
    m_vtkTessellationTilesActor->SetMapper(facesMapper);
    m_vtkTessellationTilesActor->GetProperty()->BackfaceCullingOn();
    m_vtkTessellationTilesActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
    m_vtkTessellationTilesActor->SetVisibility(m_viewTessellationGeometry);

    // Initialize and define the actor for the bounding box
    m_vtkTessellationBoundBoxActor = InitAxisAlignedBoundingBoxActor(vf);
    m_vtkTessellationBoundBoxActor->SetVisibility(m_viewTessellationBoundingBox);

    // Add the actors to the renderer
    m_vtkRenderer->AddActor(m_vtkTessellationEdgesActor);
    m_vtkRenderer->AddActor(m_vtkTessellationTilesActor);
    m_vtkRenderer->AddActor(m_vtkTessellationBoundBoxActor);
}

void MainTigerWindow::InitGridActor(size_t width, size_t height, size_t widthSegments, size_t heightSegments)
{
	// Calculate the number of points required for drawing the grid
	size_t nPoints = ((widthSegments + 1) * 2) + ((heightSegments - 1) * 2);

	// 
	size_t nHorizontalLines = heightSegments + 1;

	// 
	size_t nVerticalLines = widthSegments + 1;

	// Calculate the number of lines required for drawing the grid
	//size_t nLines = (widthSegments + 1) + (heightSegments + 1);

	// 
	double halfWidth = width / 2.0;
	double halfHeight = height / 2.0;
	double widthStep = width / (double)widthSegments;
	double heightStep = height / (double)heightSegments;

	// Initialize the object to store the points of the grid
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(nPoints);

	// Initialize the object to store the lines of the grid
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

	// 
	for (size_t w = 0; w <= widthSegments; w += 1) 
	{
		// 
		size_t index = w * 2;

		// 
		double P1[3] = { -halfWidth + (widthStep * (double)w), 0.0,  halfHeight };
		double P2[3] = { -halfWidth + (widthStep * (double)w), 0.0, -halfHeight };

		// 
		points->SetPoint(index + 0, P1);
		points->SetPoint(index + 1, P2);

		// Initialize the object representing the current line and store it in the lines object
		vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
		line->GetPointIds()->SetId(0, index + 0);
		line->GetPointIds()->SetId(1, index + 1);
		lines->InsertNextCell(line);
	}

	// Initialize the object representing the top line of the grid and store it in the lines object
	vtkSmartPointer<vtkLine> topLine = vtkSmartPointer<vtkLine>::New();
	topLine->GetPointIds()->SetId(0, 0);
	topLine->GetPointIds()->SetId(1, widthSegments * 2);
	lines->InsertNextCell(topLine);

	// 
	for (size_t h = 1; h < heightSegments; h += 1) 
	{
		// 
		size_t index = (widthSegments * 2) + (h * 2);

		// 
		double P1[3] = { -halfWidth, 0.0, halfHeight - (heightStep * (double)h) };
		double P2[3] = {  halfWidth, 0.0, halfHeight - (heightStep * (double)h) };

		// 
		points->SetPoint(index + 0, P1);
		points->SetPoint(index + 1, P2);

		// Initialize the object representing the current line and store it in the lines object
		vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
		line->GetPointIds()->SetId(0, index + 0);
		line->GetPointIds()->SetId(1, index + 1);
		lines->InsertNextCell(line);
	}

	// Initialize the object representing the bottom line of the grid and store it in the lines 
	// object
	vtkSmartPointer<vtkLine> bottomLine = vtkSmartPointer<vtkLine>::New();
	bottomLine->GetPointIds()->SetId(0, 1);
	bottomLine->GetPointIds()->SetId(1, (widthSegments * 2) + 1);
	lines->InsertNextCell(bottomLine);

	// Initialize and define the polydata
	vtkSmartPointer<vtkPolyData> linesPolydata = vtkSmartPointer<vtkPolyData>::New();
	linesPolydata->SetPoints(points);
	linesPolydata->SetLines(lines);

	// Initialize and define the mapper
	vtkSmartPointer<vtkPolyDataMapper> linesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	linesMapper->SetInputData(linesPolydata);

	// Initialize and define the grid actor
	m_vtkGridActor = vtkActor::New();
	m_vtkGridActor->SetMapper(linesMapper);
    m_vtkGridActor->GetProperty()->SetColor(0.8, 0.8, 0.8);
    m_vtkGridActor->SetVisibility(m_viewGrid);

    // Add the actor to the renderer
    m_vtkRenderer->AddActor(m_vtkGridActor);
}

void MainTigerWindow::InitInterfacePolygonsActor(const std::shared_ptr<InterfacePolygons> intfs)
{
    // Clear the interface polygons actor
    ClearInteracePolygonsActor();

    // Initialize the object to store the points (vertices) of the interface polygons
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(intfs->CountVertices());

    // Initialize the object to store the triangles from the interface polygon. Each interface 
    // polygon is triangulated with respect to one of its vertices
    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    triangles->SetNumberOfCells(intfs->CountTriangles());

    // 
    Eigen::Vector3d C = Eigen::Vector3d::Zero();

    // 
    double P[3] = { 0.0, 0.0, 0.0 };

    // 
    size_t insertedPoints = 0, nVertices = 0, vIdx = 0, nFaces = 0, fIdx = 0, nIndices = 0, i = 0;

    // Traverse through the interface polygons. Add their vertices into the points object, 
    // triangulate their face and insert it into the triangles objecttriangulate their face
    for (auto it = intfs->GetInterfaces().begin(); it != intfs->GetInterfaces().end(); ++it)
    {
        // Get the pointer to the current interface polygon
        std::shared_ptr<VF> intf = *it;

        // If there it no pointer to an interface polygon then continue to the next one
        if (!intf) 
        {
            continue;
        }

        // Get the number of vertices of the geometry
        nVertices = intf->countVertices();

        // Traverse through the vertices of the geometry and insert them into the points object
        for (vIdx = 0; vIdx < nVertices; vIdx += 1) 
        {
            // Get the reference to the coordinates of the current vertex
            C << intf->Vertex(vIdx);

            // Store the coordinate values of the vertex
            P[0] = C.x();
            P[1] = C.y();
            P[2] = C.z();

            // Insert the point into the points object at its respective index
            points->InsertPoint(insertedPoints + vIdx, P);
        }

        // Get the number of faces of the geometry
        nFaces = intf->countFaces();

        // Traverse through the faces of the geometry, triangulate them and insert the triangles 
        // into the triangles object. Each face is triangulated with respect to the start vertex of
        // the incident half edge of each face
        for (fIdx = 0; fIdx < nFaces; fIdx += 1) 
        {
            // Get the pointer to the current face
            const std::vector<size_t> & face = intf->face(fIdx);

            // 
            nIndices = face.size();

            for (i = 1; i < nIndices - 1; i += 1) 
            {
                // Initialize and define the current triangle
                vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                triangle->GetPointIds()->SetId(0, insertedPoints + face[0]);
                triangle->GetPointIds()->SetId(1, insertedPoints + face[i]);
                triangle->GetPointIds()->SetId(2, insertedPoints + face[i + 1]);

                // Insert the triangle into the triangles object
                triangles->InsertNextCell(triangle);
            }
        }

        // Update the number of inserted points
        insertedPoints += nVertices;
    }

    // 
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetPolys(triangles);

    // 
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polyData);

    // Initialize and define the actor for rendering the interface polygons
    m_vtkInterfacePolygonsActor = vtkActor::New();
    m_vtkInterfacePolygonsActor->SetMapper(mapper);
    m_vtkInterfacePolygonsActor->GetProperty()->SetColor(1.0, 1.0, 0.0);
    m_vtkInterfacePolygonsActor->SetVisibility(m_viewInterfacePolygons);

    // Add the actor to the renderer
    m_vtkRenderer->AddActor(m_vtkInterfacePolygonsActor);

    // Update the current interface polygon type
    m_currentInterfacePolygonsType = INTERFACES::PLAIN;

    std::cout << "#Interfaces = " << intfs->CountInterfaces() << std::endl;
}

void MainTigerWindow::InitInterfacePolygonsActor(
    const std::shared_ptr<InterfacePolygons> intfs,
    const std::shared_ptr<EquilibriumAnalysis::Result> results, 
    INTERFACES type)
{
    assert(intfs);
    assert(results);
    assert(type != INTERFACES::PLAIN);

    // Clear the interface polygons actor
    ClearInteracePolygonsActor();

    size_t nInterfaceVertices = intfs->CountVertices();

    // Initialize the object to store the points (vertices) of the interface polygons
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(nInterfaceVertices);

    // Initialize the object to store the color values of the points. Color values go from 0 to 
    // 1. They Are normalized values from the force components from an Equilibrium Analysis result
    vtkSmartPointer<vtkDoubleArray> pointColors = vtkSmartPointer<vtkDoubleArray>::New();
    pointColors->SetNumberOfValues(nInterfaceVertices);

    // Initialize the object to store the triangles from the interface polygon. Each interface 
    // polygon is triangulated with respect to one of its vertices
    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    triangles->SetNumberOfCells(intfs->CountTriangles());

    // 
    Eigen::Vector3d C = Eigen::Vector3d::Zero();

    // 
    double P[3] = { 0.0, 0.0, 0.0 };

    // Keep count of the number of inserted vertices in the points object
    size_t insertedPoints = 0;

    // Traverse through the interface polygons. Add their vertices into the points object, 
    // triangulate their face and insert it into the triangles object triangulate their face
    for (auto it = intfs->GetInterfaces().begin(); it != intfs->GetInterfaces().end(); ++it)
    {
        // Get the pointer to the current interface polygon
        std::shared_ptr<VF> intf = *it;

        // If no pointer to an interface polygon then continue with the next one
        if (!intf) 
        {
            continue;
        }

        // Get the number of vertices of the geometry
        //size_t nVertices = intf->Geometry().vertices.size();

        // 
        bool intfInEquilibrium = false;

        // 
        size_t intfEquilibriumIndex = 0;

        // If the current interface has the ATTRIB_IN_EQUILIBRIUM then get its value
        if (intf->Attributes().Has(ATTRIB_IN_EQUILIBRIUM))
        {
            assert(intf->Attributes().Get<bool>(ATTRIB_IN_EQUILIBRIUM, intfInEquilibrium));

            if (intfInEquilibrium) 
            {
                assert(intf->Attributes().Get<size_t>(ATTRIB_EQUILIBRIUM_INDEX, intfEquilibriumIndex));
            }
        }

        size_t nVertices = intf->countVertices();

        // Get the pointer to the incident half edge of the face of the interface polygon
        //std::shared_ptr<dcel::Halfedge> halfedge = intf->Geometry().faces[0]->halfedge;

        for (size_t vIdx = 0; vIdx < nVertices; vIdx += 1)
        {
            // Get the reference to the coordinates of the current vertex
            //const Eigen::Vector3d & C = halfedge->start->Coords();
            C << intf->Vertex(vIdx);

            // Store the coordinate values of the vertex
            P[0] = C.x();
            P[1] = C.y();
            P[2] = C.z();

            // Insert the point into the points object at its respective index
            points->InsertPoint(insertedPoints + vIdx, P);

            // Initialize the force value of the vertex
            double force = -1e-08;

            // If the current interface is in the equilibrium analysis model then get the required
            // force magnitude from the results object
            if (intfInEquilibrium)
            {
                // Define the key for the respective force 
                std::tuple<size_t, size_t> key = std::make_tuple(intfEquilibriumIndex, vIdx);

                switch (type) 
                {
                    case INTERFACES::COMPRESSION: 
                    {
                        force = results->forces.find(key)->second.compression;
                        break;
                    }

                    case INTERFACES::TENSION:
                    {
                        force = results->forces.find(key)->second.tension;
                        break;
                    }

                    case INTERFACES::UTANGENTIAL:
                    {
                        force = results->forces.find(key)->second.uTangential;
                        break;
                    }

                    case INTERFACES::VTANGENTIAL:
                    {
                        force = results->forces.find(key)->second.vTangential;
                        break;
                    }
                }
            }

            // Get the force value associated to the vertex. Then, normalize its value with respect
            // to the maximum force value from the results and insert it into the points color 
            // object
            pointColors->SetValue(insertedPoints + vIdx, force);
        }

        // Get the number of faces of the geometry
        size_t nFaces = intf->countFaces();

        // Traverse through the faces of the geometry, triangulate them and insert the triangles 
        // into the triangles object. Each face is triangulated with respect to the start vertex of
        // the incident half edge of each face
        for (size_t fIdx = 0; fIdx < nFaces; fIdx += 1)
        {
            // Get the pointer to the current face
            const std::vector<size_t> & face = intf->face(fIdx);

            // 
            size_t nIndices = face.size();

            for (size_t i = 1; i < nIndices - 1; i += 1)
            {
                // Initialize and define the current triangle
                vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                triangle->GetPointIds()->SetId(0, insertedPoints + face[0]);
                triangle->GetPointIds()->SetId(1, insertedPoints + face[i]);
                triangle->GetPointIds()->SetId(2, insertedPoints + face[i + 1]);

                // Insert the triangle into the triangles object
                triangles->InsertNextCell(triangle);
            }
        }

        // Update the number of inserted points
        insertedPoints += nVertices;
    }

    // 
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetPolys(triangles);
    polyData->GetPointData()->SetScalars(pointColors);

    // 
    double min = 0, max = 0;
    results->GetMinMaxForces((EquilibriumAnalysis::Force::TYPE)type, min, max);
    //results->GetCTMinMaxForces(min, max);
    
    double step = (max - min) / 4.0;

    // 
    vtkSmartPointer<vtkColorTransferFunction> colorFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
    //colorFunction->AddRGBPoint(-1.00, 1.0, 1.0, 1.0);
    //colorFunction->AddRGBPoint( 0.00, 0.0, 0.0, 1.0);
    //colorFunction->AddRGBPoint( 0.25, 0.0, 1.0, 1.0);
    //colorFunction->AddRGBPoint( 0.50, 0.0, 1.0, 0.0);
    //colorFunction->AddRGBPoint( 0.75, 1.0, 1.0, 0.0);
    //colorFunction->AddRGBPoint( 1.00, 1.0, 0.0, 0.0);
    colorFunction->AddRGBPoint(-1e-8, 0.5, 0.5, 0.5);
    colorFunction->AddRGBPoint( min, 0.0, 0.0, 1.0);
    colorFunction->AddRGBPoint( min + step, 0.0, 1.0, 1.0);
    colorFunction->AddRGBPoint( (min + max) / 2.0, 0.0, 1.0, 0.0);
    colorFunction->AddRGBPoint( max - step, 1.0, 1.0, 0.0);
    colorFunction->AddRGBPoint( max, 1.0, 0.0, 0.0);



    //colorFunction->AddRGBPoint(-1.0,                             1.0, 1.0, 1.0);
    //colorFunction->AddRGBPoint(result.minTension,                0.0, 0.0, 1.0);
    //colorFunction->AddRGBPoint(result.minTension + step,         1.0, 1.0, 0.0);
    //colorFunction->AddRGBPoint(result.minTension + (2.0 * step), 0.0, 1.0, 0.0);
    //colorFunction->AddRGBPoint(result.minTension + (3.0 * step), 0.0, 1.0, 1.0);
    //colorFunction->AddRGBPoint(result.maxTension,                1.0, 0.0, 0.0);
    colorFunction->Build();

    // 
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polyData);
    mapper->SetLookupTable(colorFunction);

    // Initialize and define the actor for rendering the interface polygons
    m_vtkInterfacePolygonsActor = vtkActor::New();
    m_vtkInterfacePolygonsActor->SetMapper(mapper);
    m_vtkInterfacePolygonsActor->GetProperty()->SetColor(1.0, 1.0, 0.0);
    m_vtkInterfacePolygonsActor->SetVisibility(m_viewInterfacePolygons);

    // Initialize and define the actor for rendering the interface polygon force bar
    m_vtkInterfacePolygonsBarActor = vtkScalarBarActor::New();
    m_vtkInterfacePolygonsBarActor->SetLookupTable(mapper->GetLookupTable());
    m_vtkInterfacePolygonsBarActor->SetTitle(GetInterfaceTypeName(type).c_str());
    m_vtkInterfacePolygonsBarActor->SetNumberOfLabels(6);
    m_vtkInterfacePolygonsBarActor->SetVisibility(m_viewInterfacePolygons);
    m_vtkInterfacePolygonsBarActor->GetTitleTextProperty()->SetColor(0.7, 0.7, 0.7);
    m_vtkInterfacePolygonsBarActor->GetLabelTextProperty()->SetColor(0.7, 0.7, 0.7);
    //m_vtkInterfacePolygonsBarActor->GetTitleTextProperty()->SetShadow(true);

    // Add the actors to the renderer
    m_vtkRenderer->AddActor(m_vtkInterfacePolygonsActor);
    m_vtkRenderer->AddActor2D(m_vtkInterfacePolygonsBarActor);

    // Update the current interface polygon type
    m_currentInterfacePolygonsType = type;
}

void MainTigerWindow::SetInterfacePolygons(bool render)
{
    // Initialize the vector to store the interface polygons between the pieces of the assembly
    std::vector<std::shared_ptr<VF>> interfaces;

    // If there are interface polygons then set them to the workspace and initialize the respective
    // actor. Otherwise, indicate there are no interface polygons
    if (m_workspace->GetAssembly()->CalculateInterfacePolygons(m_workspace->GetTessellation(), interfaces))
    {
        // Set the interface polygons to the workspace
        m_workspace->SetInterfacePolygons(interfaces);

        // Initialize the actor for rendering the interface polygons
        InitInterfacePolygonsActor(m_workspace->GetInterfacePolygons());

        // 
        if (render)
        {
            qvtkWidget->renderWindow()->Render();
            //qvtkWidget->renderWindow()->Render();
        }
    }
    else
    {
        QMessageBox::warning(
            this,
            "No interface polygons",
            "There are no interface polygons between the pieces");
    }
}

void MainTigerWindow::ViewInterfacePolygons(INTERFACES type)
{
    if (m_viewInterfacePolygons && type == m_currentInterfacePolygonsType) 
    {
        return;
    }

    m_viewInterfacePolygons = true;

    switch (m_currentInterfacePolygonsType) 
    {
        case INTERFACES::PLAIN: 
        {
            m_currentInterfacePolygonsType = INTERFACES::PLAIN;

            this->actionViewPlainInterfaces->setChecked(true);
            this->actionViewCompressionForces->setChecked(false);
            this->actionViewTensionForces->setChecked(false);
            this->actionViewUTangentialForces->setChecked(false);
            this->actionViewVTangentialForces->setChecked(false);

            break;
        }

        case INTERFACES::COMPRESSION: 
        {
            m_currentInterfacePolygonsType = INTERFACES::COMPRESSION;

            this->actionViewPlainInterfaces->setChecked(false);
            this->actionViewCompressionForces->setChecked(true);
            this->actionViewTensionForces->setChecked(false);
            this->actionViewUTangentialForces->setChecked(false);
            this->actionViewVTangentialForces->setChecked(false);

            break;
        }

        case INTERFACES::TENSION: 
        {
            m_currentInterfacePolygonsType = INTERFACES::TENSION;

            this->actionViewPlainInterfaces->setChecked(false);
            this->actionViewCompressionForces->setChecked(false);
            this->actionViewTensionForces->setChecked(true);
            this->actionViewUTangentialForces->setChecked(false);
            this->actionViewVTangentialForces->setChecked(false);

            break;
        }

        case INTERFACES::UTANGENTIAL: 
        {
            m_currentInterfacePolygonsType = INTERFACES::UTANGENTIAL;

            this->actionViewPlainInterfaces->setChecked(false);
            this->actionViewCompressionForces->setChecked(false);
            this->actionViewTensionForces->setChecked(false);
            this->actionViewUTangentialForces->setChecked(true);
            this->actionViewVTangentialForces->setChecked(false);

            break;
        }

        case INTERFACES::VTANGENTIAL: 
        {
            m_currentInterfacePolygonsType = INTERFACES::VTANGENTIAL;

            this->actionViewPlainInterfaces->setChecked(false);
            this->actionViewCompressionForces->setChecked(false);
            this->actionViewTensionForces->setChecked(false);
            this->actionViewUTangentialForces->setChecked(false);
            this->actionViewVTangentialForces->setChecked(true);

            break;
        }
    }

    if (m_vtkInterfacePolygonsActor) 
    {
        m_vtkInterfacePolygonsActor->VisibilityOn();
    }

    if (type != INTERFACES::PLAIN && m_vtkInterfacePolygonsBarActor) 
    {
        m_vtkInterfacePolygonsBarActor->VisibilityOn();
    }
}

void MainTigerWindow::ViewTessellationGeometry()
{
    if (m_viewTessellationGeometry) 
    {
        return;
    }

    m_viewTessellationGeometry = true;
    this->actionViewTessellationGeometry->setChecked(true);

    if (m_vtkTessellationEdgesActor) 
    {
        m_vtkTessellationEdgesActor->VisibilityOn();
        m_vtkTessellationTilesActor->VisibilityOn();
    }
}

void MainTigerWindow::on_actionNewArchimedeanSolidGeometry_triggered()
{
}

void MainTigerWindow::on_actionNewBendedSquareGeometry_triggered()
{
    // Open the dialog
    BendedSquareParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Clear the workspace
    ClearWorkspace();

    // Get the parameters from the dialog
    double length = dialog.GetLength();
    double radius = dialog.GetRadius();
    size_t lengthSegments = dialog.GetLengthSegments();
    size_t radialSegments = dialog.GetRadialSegments();
    QString bendDirection = dialog.GetBendDirection();

    // 
    int dir = (bendDirection == "Upwards") ? 1 : -1;

    // 
    VF vf = geometries::BendedSquare(radius, length, radialSegments, lengthSegments, dir);

    // Initialize a new workspace and set its tessellation
    m_workspace = std::make_shared<Workspace>();
    m_workspace->SetTessellation(vf);

    // Initialize the actors that renders the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // If the tessellation geometry is not visible then ask whether or not to show it
    if (!m_viewTessellationGeometry)
    {
        AskViewTessellationGeometry();
    }

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
    //qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionNewConeGeometry_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionNewConeGeometry_triggered");
}

void MainTigerWindow::on_actionNewCyclideGeometry_triggered()
{
    //QMessageBox::information(this, APP_TITLE, "void on_actionNewCyclideGeometry_triggered");

    CyclideParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected) 
    {
        return;
    }

    double a = dialog.getA();
    double b = dialog.getB();
    double c = dialog.getC();
    double d = dialog.getD();
    size_t Rs = dialog.getMajorRadialSegments();
    size_t rs = dialog.getMinorRadialSegments();

    //std::cout << "a = " << a << std::endl;
    //std::cout << "b = " << b << std::endl;
    //std::cout << "c = " << c << std::endl;
    //std::cout << "d = " << d << std::endl;
    //std::cout << "Rs = " << Rs<< std::endl;
    //std::cout << "rs = " << rs << std::endl;

    // Clear the workspace
    ClearWorkspace();

    // 
    VF vf = geometries::Cyclide(a, b, c, d, Rs, rs);

    // Initialize a new workspace and set its tessellation
    m_workspace = std::make_shared<Workspace>();
    m_workspace->SetTessellation(vf);

    // Initialize the actors that renders the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // If the tessellation geometry is not visible then ask whether or not to show it
    if (!m_viewTessellationGeometry)
    {
        AskViewTessellationGeometry();
    }

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionNewCylinderGeometry_triggered()
{
    // Open the dialog
    CylinderParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Clear the workspace
    ClearWorkspace();

    // Get the parameters from the dialog
    double length = dialog.GetLength();
    double radius = dialog.GetRadius();
    size_t lengthSegments = dialog.GetLengthSegments();
    size_t radialSegments = dialog.GetRadialSegments();
    QString axis = dialog.GetAxis();

    // 
    const Eigen::Vector3d C(0.0, 0.0, 0.0);
    
    // 
    const Eigen::Vector3d K = (axis == "X") ? Eigen::Vector3d(1.0, 0.0, 0.0) : 
        ((axis == "Y") ? Eigen::Vector3d(0.0, 1.0, 0.0) : Eigen::Vector3d(0.0, 0.0, 1.0));

    // 
    VF vf = geometries::Cylinder(C, K, radius, length, radialSegments, lengthSegments);

    // Initialize a new workspace and set its tessellation
    m_workspace = std::make_shared<Workspace>();
    m_workspace->SetTessellation(vf);

    // Initialize the actors that renders the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // If the tessellation geometry is not visible then ask whether or not to show it
    if (!m_viewTessellationGeometry) 
    {
        AskViewTessellationGeometry();
    }

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionNewEquilateralTriangleGeometry_triggered()
{
    // Open the dialog
    EquilateralTriangleParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Clear the workspace
    ClearWorkspace();

    // Get the parameters from the dialog
    double length = dialog.GetLength();
    QString plane = dialog.GetPlane();

    Eigen::Vector3d U = Eigen::Vector3d::Zero(), V = Eigen::Vector3d::Zero();
    UVFromPlane(plane, U, V);

    // Get the geometry
    VF vf = geometries::EquilateralTriangle(length, U, V);

    // Initialize a new workspace and set its tessellation
    m_workspace = std::make_shared<Workspace>();
    m_workspace->SetTessellation(vf);

    // Initialize the actors that renders the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // If the tessellation geometry is not visible then ask whether or not to show it
    if (!m_viewTessellationGeometry)
    {
        AskViewTessellationGeometry();
    }

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionNewNsidedRegularPolygonGeometry_triggered()
{
    // Open the dialog
    RegularPolygonParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Clear the workspace
    ClearWorkspace();

    // Get the parameters from the dialog
    size_t nSides = dialog.GetSides();
    double length = dialog.GetLength();
    QString plane = dialog.GetPlane();

    Eigen::Vector3d U = Eigen::Vector3d::Zero(), V = Eigen::Vector3d::Zero();
    UVFromPlane(plane, U, V);

    // Get the geometry
    VF vf = geometries::Polygon(nSides, length, U, V);

    // Initialize a new workspace and set its tessellation
    m_workspace = std::make_shared<Workspace>();
    m_workspace->SetTessellation(vf);

    // Initialize the actors that renders the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // If the tessellation geometry is not visible then ask whether or not to show it
    if (!m_viewTessellationGeometry)
    {
        AskViewTessellationGeometry();
    }

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionNewObjFileGeometry_triggered()
{
    // Show the open file dialog
    QString fileName = QFileDialog::getOpenFileName(
        this, 
        "Open OBJ file", 
        m_lastDirectory, 
        "Wavefront OBJ (*.obj)");

    // Exit the function if the user doesn't select a file
    if (fileName.isNull()) 
    {
        return;
    }

    // Update the location of the last directory
    m_lastDirectory = fileName;

    // Initialize the objects that store the content from the OBJ file
    tinyobj::attrib_t objAttributes;
    std::vector<tinyobj::shape_t> objShapes;
    std::vector<tinyobj::material_t> objMaterials;
    std::string objWarning;
    std::string objError;

    // Load the OBJ file
    bool ret = tinyobj::LoadObj(
        &objAttributes, 
        &objShapes, 
        &objMaterials, 
        &objWarning, 
        &objError, 
        fileName.toStdString().c_str(), 
        nullptr, 
        false);

    // Show any warning message
    if (!objWarning.empty()) 
    {
        QMessageBox::warning(this, "Warning", QString::fromStdString(objWarning));
    }

    // Show any error message and exit the function (if something went wrong then let's be safe and
    // do not generate any tessellation
    if (!objError.empty() || !ret) 
    {
        QMessageBox::critical(this, "Error", QString::fromStdString(objError));
        return;
    }

    // Stop the process if there are more than one shape in the OBJ file. At this moment we only 
    // support single shaped tessellations
    if (objShapes.size() > 1) 
    {
        QMessageBox::critical(this, "Error", "File contains more than one geometry.");
        return;
    }

    // Get the number of vertices and the number of faces from the shape
    size_t nVertices = objAttributes.vertices.size() / 3;
    size_t nFaces = objShapes[0].mesh.num_face_vertices.size();

    // Initialize the object to store the vertex coordinates and vertex indices of the geometry
    VF vf(nVertices, nFaces);
    
    // Traverse through the vertices from the geomnetry, get their coordinate values and store them
    // in the VF object
    for (size_t i = 0; i < nVertices; i += 1) 
    {
        double x = objAttributes.vertices[(i * 3) + 0];
        double y = objAttributes.vertices[(i * 3) + 1];
        double z = objAttributes.vertices[(i * 3) + 2];

        //vf.V[i] << x, y, z;
        vf.addVertex(x, y, z);
    }

    // Initialize an index offset for traversing the face indices properly
    size_t indexOffset = 0;

    // Traverse through the faces and store them in the VF object
    for (size_t i = 0; i < nFaces; i += 1) 
    {
        // Get the number of vertices for the current face
        size_t nFaceVertices = objShapes[0].mesh.num_face_vertices[i];

        // 
        //vf.F[i].resize(nFaceVertices);
        std::vector<size_t> indices(nFaceVertices);

        // 
        for (size_t j = 0; j < nFaceVertices; j += 1) 
        {
            //vf.F[i][j] = objShapes[0].mesh.indices[indexOffset + j].vertex_index;
            indices[j] = objShapes[0].mesh.indices[indexOffset + j].vertex_index;
        }

        vf.addFace(indices);

        indexOffset += nFaceVertices;
    }

    // Clear the workspace
    ClearWorkspace();

    // Initialize a new workspace and set its tessellation
    m_workspace = std::make_shared<Workspace>();
    m_workspace->SetTessellation(vf);

    // Initialize the actors that renders the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // If the tessellation geometry is not visible then ask whether or not to show it
    if (!m_viewTessellationGeometry)
    {
        AskViewTessellationGeometry();
    }

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionNewParaboloidGeometry_triggered()
{
    // Open the dialog
    ParaboloidParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Clear the workspace
    ClearWorkspace();

    // Get the parameters from the dialog
    double width = dialog.GetWidth();
    double height = dialog.GetHeight();;
    size_t ws = dialog.GetWidthSegments();
    size_t hs = dialog.GetHeightSegments();
    QString plane = dialog.GetPlane();
    double A = dialog.GetA();
    double B = dialog.GetB();
    QString type = dialog.GetType();

    Eigen::Vector3d U = Eigen::Vector3d::Zero(), V = Eigen::Vector3d::Zero();
    UVFromPlane(plane, U, V);

    // Get the geometry
    VF vf = (type == "Elliptic") ?
        geometries::EllipticParaboloid(A, B, width, height, ws, hs, U, V) :
        geometries::HyperbolicParaboloid(A, B, width, height, ws, hs, U, V);

    // Initialize a new workspace and set its tessellation
    m_workspace = std::make_shared<Workspace>();
    m_workspace->SetTessellation(vf);

    // Initialize the actors that renders the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // If the tessellation geometry is not visible then ask whether or not to show it
    if (!m_viewTessellationGeometry)
    {
        AskViewTessellationGeometry();
    }

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionNewPlatonicSolidGeometry_triggered()
{
    // Open the dialog
    PlatonicSolidParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Clear the workspace
    ClearWorkspace();

    // Get the parameters from the dialog
    QString solid = dialog.GetSolid();
    double radius = dialog.GetRadius();

    // Get the geometry
    VF vf = (solid == "Tetrahedron") ? geometries::Tetrahedron(radius) : 
        ((solid == "Cube") ? geometries::Hexahedron(radius) : 
        ((solid == "Octahedron") ? geometries::Octahedron(radius) : 
        ((solid == "Dodecahedron") ? geometries::Dodecahedron(radius) : 
        geometries::Icosahedron(radius))));

    // Initialize a new workspace and set its tessellation
    m_workspace = std::make_shared<Workspace>();
    m_workspace->SetTessellation(vf);

    // Initialize the actors that renders the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // If the tessellation geometry is not visible then ask whether or not to show it
    if (!m_viewTessellationGeometry)
    {
        AskViewTessellationGeometry();
    }

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionNewRectangleGeometry_triggered()
{
    // Open the dialog
    RectangleParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Clear the workspace
    ClearWorkspace();

    // Get the parameter values
    double width = dialog.GetWidth();
    double height = dialog.GetHeight();
    QString plane = dialog.GetPlane();

    Eigen::Vector3d U = Eigen::Vector3d::Zero(), V = Eigen::Vector3d::Zero();
    UVFromPlane(plane, U, V);

    // Get the geometry
    VF vf = geometries::Rectangle(width, height, U, V);

    // Initialize a new workspace and set its tessellation
    m_workspace = std::make_shared<Workspace>();
    m_workspace->SetTessellation(vf);

    // Initialize the actors that renders the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // If the tessellation geometry is not visible then ask whether or not to show it
    if (!m_viewTessellationGeometry)
    {
        AskViewTessellationGeometry();
    }

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionCameraBack_triggered() 
{
}

void MainTigerWindow::on_actionCameraBottom_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionCameraBottom_triggered");
}

void MainTigerWindow::on_actionCameraFront_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionCameraFront_triggered");
}

void MainTigerWindow::on_actionCameraLeft_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionCameraLeft_triggered");
}

void MainTigerWindow::on_actionCameraManual_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionCameraManual_triggered");
}

void MainTigerWindow::on_actionCameraReset_triggered()
{
	// Reset the eye transformation matrix by setting it as an identity matrix
	//m_vtkCamera->GetViewTransformMatrix()->Identity();

    // Reset the location and focus of the camera
	m_vtkCamera->SetPosition(3, 3, 3);
	m_vtkCamera->SetFocalPoint(0, 0, 0);

	// Call Render to show changes
	qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionCameraRight_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionCameraRight_triggered");
}

void MainTigerWindow::on_actionCameraTop_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionCameraTop_triggered");
}

void MainTigerWindow::on_actionEditAssemblyClipAdaptiveTileOffset_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements for clipping the pieces. Show a warning message if a requirement is
    // missing
    if (!m_workspace->CheckClippingRequirements(msg))
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // Open the adaptive tile offset clipping parameters dialog
    AdaptiveTileOffsetClippingParametersDialog dialog(this);
    
    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Get the selected function values
    std::string topFunction = dialog.GetTopFunction().toStdString();
    std::string bottomFunction = dialog.GetBottomFunction().toStdString();

    // Clear the actors for rendering the assembly and the interface polygons
    ClearAssemblyActors();
    ClearInteracePolygonsActor();

    // Clip the assembly using adaptive tile offset and set it to the workspace
    m_workspace->SetAssemblyUsingAdaptiveTileOffsetClippedBlocks(topFunction, bottomFunction);

    // Initialize the actors for rendering the assembly
    InitAssemblyActors(m_workspace->GetTessellation(), m_workspace->GetAssembly());

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionEditAssemblyClipTileOffset_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements for clipping the pieces. Show a warning message if a requirement is
    // missing
    if (!m_workspace->CheckClippingRequirements(msg))
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // Open the tile offset clipping parameters dialog
    TileOffsetClippingParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Get the truncation parameter value
    double extrados = dialog.GetExtrados();
    double intrados = dialog.GetIntrados();

    // Validate the values Are correct
    if (intrados == 0.0 || extrados == 0.0)
    {
        QMessageBox::warning(this, "Invalid Value", "Values cannot be zero.");
        return;
    }

    // Clear the actors for rendering the assembly and the interface polygons
    ClearAssemblyActors();
    ClearInteracePolygonsActor();

    // 
    m_workspace->SetAssemblyUsingTileOffsetClippedBlocks(extrados, intrados);

    // Initialize the actors for rendering the assembly
    InitAssemblyActors(m_workspace->GetTessellation(), m_workspace->GetAssembly());

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionEditFlipTessellation_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements for translating the centroid of the tessellation to the origin. 
    // Show a warning message if a requirement is missing
    if (!m_workspace->CheckFlipTessellationRequirements(msg))
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // Erase the face centers from the workspace (edge directions, pieces and interface polygons 
    // are removed as well)
    m_workspace->EraseTileCenters();
    m_workspace->GetTessellation()->ClearDCEL();

    // Delete the actors that contain information about the workspace
    DeleteWorkspaceActors();

    // 
    m_workspace->GetTessellation()->Geometry().flip();

    // Initialize the actors for rendering the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionEditTessellationCentroidToOrigin_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements for translating the centroid of the tessellation to the origin. 
    // Show a warning message if a requirement is missing
    if (!m_workspace->CheckCentroidToOriginRequirements(msg))
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // Get the centroid of the tessellation and invert it
    Eigen::Vector3d C = m_workspace->GetTessellation()->Geometry().centroid(true);
    C *= -1.0;

    // Erase the face centers from the workspace (edge directions, pieces and interface polygons 
    // are removed as well)
    m_workspace->EraseTileCenters();
    m_workspace->GetTessellation()->ClearDCEL();

    // Delete the actors that contain information about the workspace
    DeleteWorkspaceActors();

    // 
    m_workspace->GetTessellation()->Geometry().Translate(C);

    // Initialize the actors for rendering the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionEditTileSubdivision_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements for subdividing the tiles of the tessellation. Show a warning message
    // if a requirement is missing
    if (!m_workspace->CheckTileSubdivisionRequirements(msg))
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // Open the face subdivision parameters dialog
    FaceSubdivisionParametersDialog dialog;

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Get the parameter values
    QString center = dialog.GetCenter();
    QString type = dialog.GetType();

    // Subdivide the faces of the tessellation as indicated
    m_workspace->SubdivideTessellationTiles(Tessellation::GetTileSubdivisionType(type.toStdString()));

    // Get the vertex coordinates and vertex indices of the subdivided tessellation
    VF vf = m_workspace->GetTessellation()->Geometry();

    // Clear the content of the workspace
    ClearWorkspace();

    // Initialize the tessellation using the subdivided geometry
    m_workspace = std::make_shared<Workspace>();
    m_workspace->SetTessellation(vf);

    // Initialize the actors for rendering the faces of the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionEditDualTessellation_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements for calculating the dual of the the tessellation. Show a warning 
    // message if a requirement is missing
    if (!m_workspace->CheckDualTessellationRequirements(msg))
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // 
    dcel::DCEL geom(m_workspace->GetTessellation()->Geometry());
    VF dual = geom.Dual();

    // Initialize the vertex coordinates and vertex indices of the dual of the tessellation
    //VF dual = m_workspace->GetTessellation()->Geometry().Dual();

    // If the tessellation has no dual then show a warning message and exit the function
    //if (!m_workspace->GetTessellation()->Geometry().Dual(dual)) 
    if (dual.countVertices() == 0 || dual.countFaces() == 0)
    {
        QMessageBox::warning(this, "Calculation Incompleted", "tessellation has no dual.");
        return;
    }

    // Clear the content of the workspace (scaling the tessellation makes all elements invalid)
    ClearWorkspace();

    // Initialize the workspace and set the geometry of the tessellation
    m_workspace = std::make_shared<Workspace>();
    m_workspace->SetTessellation(dual);

    // Initialize the actors for rendering the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionEditNormalizeVertices_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements for normalizing the tessellation. Show a warning message if a 
    // requirement is missing
    if (!m_workspace->CheckVertexNormalizationRequirements(msg))
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // Open the scale tessellation parameters dialog
    NormalizeTessellationParametersDialog dialog;

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Get the normalization parameter values
    double radius = dialog.GetRadius();

    // Erase the face centers from the workspace (edge directions, pieces and interface polygons 
    // are removed as well)
    m_workspace->EraseTileCenters();
    m_workspace->GetTessellation()->ClearDCEL();

    // Delete the actors that contain information about the workspace
    DeleteWorkspaceActors();

    // Normalize the vertices of the geometry of the tessellation
    m_workspace->GetTessellation()->Geometry().NormalizeVertices();
    m_workspace->GetTessellation()->Geometry().Scale(radius);

    // Initialize the actors for rendering the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionEditRestoreAssembly_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionEditRestoreAssembly_triggered");
}

void MainTigerWindow::on_actionEditRotateTessellation_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements for scaling the tessellation. Show a warning message if a 
    // requirement is missing
    if (!m_workspace->CheckRotateTessellationRequirements(msg))
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // Open the rotate tessellation parameters dialog
    RotateTessellationParametersDialog dialog;

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Get the rotation parameters, angle must be transformed to radians
    double X = dialog.GetX();
    double Y = dialog.GetY();
    double Z = dialog.GetZ();
    double angle = dialog.GetAngle() * utils::PI / 180.0;

    // Define the axis vector
    Eigen::Vector3d K(X, Y, Z);
    K.normalize();

    // Erase the face centers from the workspace (edge directions, pieces and interface polygons 
    // are removed as well)
    m_workspace->EraseTileCenters();
    m_workspace->GetTessellation()->ClearDCEL();

    // Delete the actors that contain information about the workspace
    DeleteWorkspaceActors();

    // Rotate the geometry of the tessellation
    m_workspace->GetTessellation()->Geometry().Rotate(K, angle);

    // Initialize the actors for rendering the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionEditScaleTessellation_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements for scaling the tessellation. Show a warning message if a 
    // requirement is missing
    if (!m_workspace->CheckScaleTessellationRequirements(msg)) 
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // Open the scale tessellation parameters dialog
    ScaleTessellationParametersDialog dialog;

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Get the scale parameters dialog
    double factor = dialog.GetFactor();

    // Erase the face centers from the workspace (edge directions, pieces and interface polygons 
    // are removed as well)
    m_workspace->EraseTileCenters();
    m_workspace->GetTessellation()->ClearDCEL();

    // Delete the actors that contain information about the workspace
    DeleteWorkspaceActors();

    // Scale the geometry of the tessellation
    m_workspace->GetTessellation()->Geometry().Scale(factor);

    // Initialize the actors for rendering the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

/*void MainTigerWindow::on_actionEditTruncatePieces_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements for truncating the pieces. Show a warning message if a requirement is
    // missing
    if (!m_workspace->CheckTruncatePiecesRequirements(msg))
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // Open the face truncate pieces parameters dialog
    TruncatePiecesParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Get the truncation parameter value
    double intrados = dialog.GetIntrados();
    double extrados = dialog.GetExtrados();

    // Validate the values Are correct
    if (intrados == 0.0 || extrados == 0.0) 
    {
        QMessageBox::warning(this, "Invalid Value", "Values cannot be zero.");
        return;
    }

    // Clear the actors for rendering the assembly and the interface polygons
    ClearAssemblyActor();
    ClearInteracePolygonsActor();

    // Truncate the pieces and set them as the assembly
    m_workspace->SetAssemblyUsingTruncatedPieces(intrados, extrados);

    // Initialize the actor for rendering the assembly
    InitAssemblyActor(m_workspace->GetAssembly());
    
    // Set the interface polygons
    SetInterfacePolygons(false);

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}*/

void MainTigerWindow::on_actionNewSquareGeometry_triggered()
{
	// Open the dialog
	SquareParametersDialog dialog(this);

	// Exit the function if the user canceled the dialog
	if (dialog.exec() == QDialog::Rejected) 
	{
		return;
	}

	// Clear the workspace
	ClearWorkspace();

	// Get the length and plane attributes for a squAre geometry
	double length = dialog.GetLength();
	QString plane = dialog.GetPlane();

    Eigen::Vector3d U = Eigen::Vector3d::Zero(), V = Eigen::Vector3d::Zero();
    UVFromPlane(plane, U, V);

	// Get the geometry
    VF vf = geometries::Square(length, U, V);

	// Initialize a new workspace and set its tessellation
	m_workspace = std::make_shared<Workspace>();
	m_workspace->SetTessellation(vf);

    // Initialize the actors that renders the tessellation
	InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // If the tessellation geometry is not visible then ask whether or not to show it
    if (!m_viewTessellationGeometry)
    {
        AskViewTessellationGeometry();
    }

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionNewSquareGridGeometry_triggered()
{
    // Open the dialog
    SquareGridParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Clear the workspace
    ClearWorkspace();

    // Get the parameters from the dialog
    double width = dialog.GetWidth();
    double height = dialog.GetHeight();
    size_t ws = dialog.GetWidthSegments();
    size_t hs = dialog.GetHeightSegments();
    QString plane = dialog.GetPlane();

    Eigen::Vector3d U = Eigen::Vector3d::Zero(), V = Eigen::Vector3d::Zero();
    UVFromPlane(plane, U, V);

    // Get the geometry
    VF vf = geometries::SquareGrid(width, height, ws, hs, U, V);

    // Initialize a new workspace and set its tessellation
    m_workspace = std::make_shared<Workspace>();
    m_workspace->SetTessellation(vf);

    // Initialize the actors that renders the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // If the tessellation geometry is not visible then ask whether or not to show it
    if (!m_viewTessellationGeometry)
    {
        AskViewTessellationGeometry();
    }

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionNewTorusGeometry_triggered()
{
    // Open the dialog
    TorusParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Clear the workspace
    ClearWorkspace();

    // Get the parameters from the dialog
    //QString axis = dialog.GetAxis();
    size_t majorSegments = dialog.getMajorRadialSegments();
    double majorRadius = dialog.getMajorRadius();
    size_t minorSegments = dialog.getMinorRadialSegments();
    double minorRadius = dialog.getMinorRadius();

    // 
    //const Eigen::Vector3d C(0.0, 0.0, 0.0);

    // 
    //const Eigen::Vector3d K = (axis == "X") ? Eigen::Vector3d(1.0, 0.0, 0.0) :
    //    ((axis == "Y") ? Eigen::Vector3d(0.0, 1.0, 0.0) : Eigen::Vector3d(0.0, 0.0, 1.0));

    // 
    VF vf = geometries::Torus(majorRadius, minorRadius, majorSegments, minorSegments);

    // Initialize a new workspace and set its tessellation
    m_workspace = std::make_shared<Workspace>();
    m_workspace->SetTessellation(vf);

    // Initialize the actors that renders the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // If the tessellation geometry is not visible then ask whether or not to show it
    if (!m_viewTessellationGeometry)
    {
        AskViewTessellationGeometry();
    }

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionNewTruncatedConeGeometry_triggered()
{
    // Open the dialog
    TruncatedConeParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Clear the workspace
    ClearWorkspace();

    // Get the parameters from the dialog
    QString axis = dialog.GetAxis();
    double bottomRadius = dialog.GetBottomRadius();
    double length = dialog.GetLength();
    size_t lengthSegments = dialog.GetLengthSegments();
    size_t radialSegments = dialog.GetRadialSegments();
    double topRadius = dialog.GetTopRadius();

    // 
    const Eigen::Vector3d C(0.0, 0.0, 0.0);

    // 
    const Eigen::Vector3d K = (axis == "X") ? Eigen::Vector3d(1.0, 0.0, 0.0) :
        ((axis == "Y") ? Eigen::Vector3d(0.0, 1.0, 0.0) : Eigen::Vector3d(0.0, 0.0, 1.0));

    // 
    VF vf = geometries::TruncatedCone(C, K, bottomRadius, topRadius, length, radialSegments, lengthSegments);

    // Initialize a new workspace and set its tessellation
    m_workspace = std::make_shared<Workspace>();
    m_workspace->SetTessellation(vf);

    // Initialize the actors that renders the tessellation
    InitTessellationActors(m_workspace->GetTessellation()->Geometry());

    // If the tessellation geometry is not visible then ask whether or not to show it
    if (!m_viewTessellationGeometry)
    {
        AskViewTessellationGeometry();
    }

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionNewWorkspace_triggered()
{
	// 
	if (m_workspace) 
	{
		m_workspace->Clear();
		m_workspace = nullptr;
		std::cout << "on_actionNewWorkspace_triggered Workspace cleAred" << std::endl;
	}

	std::cout << "on_actionNewWorkspace_triggered New Workspace" << std::endl;
}

void MainTigerWindow::on_actionSaveAbaqusFile_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionSaveAbaqusFile_triggered");
}

void MainTigerWindow::on_actionSaveGeogebraJsTessellationFile_triggered()
{
    // Exit the function if there is no tessellation
    if (!m_workspace->GetTessellation()) 
    {
        QMessageBox::critical(this, APP_TITLE, "There is no tessellation in the workspace.");
        return;
    }

    // Get the name of the file
    QString filename = QFileDialog::getSaveFileName(
        this, tr("Save JavaScript File"), "Tessellation.js", tr("JavaScript (*.js)"));

    // Exit the function if there is no file name
    if (filename.trimmed() == "") 
    {
        return;
    }

    // Write the JavaScript file with the commands for generating the tessellation in GeoGebra
    m_workspace->WriteTessellationGeogebraJsFile(filename.toStdString());

    // Indicate the file was generated successfully
    QMessageBox::information(this, APP_TITLE, filename.trimmed() + " saved successfully.");
}

void MainTigerWindow::on_actionSaveObjAllAssemblyFile_triggered()
{
    // Exit the function if there is no assembly
    if (!m_workspace->GetAssembly()) 
    {
        QMessageBox::critical(this, APP_TITLE, "There is no assembly in the workspace.");
        return;
    }

    // Get the selected filename and trim it
    QString filename = QFileDialog::getSaveFileName(
        this, tr("Save OBJ File"), "assembly.obj", tr("WaveFront (*.obj)")).trimmed();

    // Exit the function if there is no file name
    if (filename == "") 
    {
        return;
    }

    // Write the OBJ file with the geometry of the assembly
    m_workspace->GetAssembly()->WriteObjFile(filename.toStdString());

    // Indicate the file was generated successfully
    QMessageBox::information(this, APP_TITLE, filename + " saved successfully.");
}

void MainTigerWindow::on_actionSaveObjAllInterfacePolygonsFile_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionSaveObjAllInterfacePolygonsFile_triggered");
}

void MainTigerWindow::on_actionSaveObjAssemblyFilePerBlock_triggered()
{
    // Exit the function if there is no assembly
    if (!m_workspace->GetAssembly())
    {
        QMessageBox::critical(this, APP_TITLE, "There is no assembly in the workspace.");
        return;
    }

    // Get the path to the current directory (we will let the user to select 
    // the directory in a future update)
    QString filedir = QDir::currentPath();

    // Get the reference to the vector with the blocks
    const std::vector<std::shared_ptr<Block>> & blocks = m_workspace->GetAssembly()->GetBlocks();

    size_t blockIdx = 0;

    // Traverse through the blocks and write their respective .obj file
    for (auto it = blocks.begin(); it != blocks.end(); ++it) 
    {
        // Get the pointer to the current block
        std::shared_ptr<Block> block = *it;

        // Build the path to the file for the current block
        std::string filepath = filedir.toStdString() + "/block_" + std::to_string(blockIdx++) + ".obj";

        block->Geometry().WriteObjFile(filepath);
        // std::cout << filepath << " saved" << std::endl;
    }

    QMessageBox::information(this, APP_TITLE, "Files saved in " + filedir + " successfully.");
}

void MainTigerWindow::on_actionSaveObjTessellationFile_triggered()
{
    // Exit the function if there is no tessellation
    if (!m_workspace->GetTessellation())
    {
        QMessageBox::critical(this, APP_TITLE, "There is no tessellation in the workspace.");
        return;
    }

    // Get the name of the file, then trim it
    QString filename = QFileDialog::getSaveFileName(
        this, tr("Save OBJ File"), "Tessellation.obj", tr("WaveFront (*.obj)")).trimmed();

    // Exit the function if there is no file name
    if (filename == "")
    {
        return;
    }

    // Write the OBJ file with the geometry of the tessellation
    m_workspace->WriteTessellationObjFile(filename.toStdString());

    // Indicate the file was generated successfully
    QMessageBox::information(this, APP_TITLE, filename + " saved successfully.");
}

void MainTigerWindow::on_actionSaveObjWorkspaceFile_triggered()
{
    //QMessageBox::information(this, APP_TITLE, "on_actionSaveObjWorkspaceFile_triggered");

    // Get the selected filename and trim it
    QString filename = QFileDialog::getSaveFileName(
        this, tr("Save OBJ File"), "workspace.obj", tr("WaveFront (*.obj)")).trimmed();

    // Exit the function if there is no file name
    if (filename == "")
    {
        return;
    }

    // Remove the extension from the filename
    //TODO

    // Initialize the OBJ exporter, get the content from the scene, and save it in a OBJ file
    vtkSmartPointer<vtkOBJExporter> exporter = vtkSmartPointer<vtkOBJExporter>::New();
    exporter->SetRenderWindow(qvtkWidget->renderWindow());
    exporter->SetFilePrefix(filename.toStdString().c_str());
    exporter->Write();

    QMessageBox::information(this, APP_TITLE, "OBJ file saved successfully.");
}

void MainTigerWindow::on_actionTicCalculateInterfaces_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements to calculate the interface polygons between the blocks. Show a 
    // warning message if a requirement is missing
    if (!m_workspace->CheckInterfacePolygonsRequirements(msg)) 
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // Initialize the vector to store the interface polygons between the blocks of the assembly
    std::vector<std::shared_ptr<VF>> interfaces;

    // Calculate the interface polygons
    size_t nInterfaces = m_workspace->GetAssembly()->CalculateInterfacePolygons(m_workspace->GetTessellation(), interfaces);

    // If the number of interfaces is not the same as the number of internal edges in the 
    // tessellation then redraw the blocks since there are new disabled blocks
    if (nInterfaces != m_workspace->GetTessellation()->DCEL()->NumberOfInternalEdges()) 
    {
        InitAssemblyActors(m_workspace->GetTessellation(), m_workspace->GetAssembly());
    }

    // Set the interface polygons to the workspace
    m_workspace->SetInterfacePolygons(interfaces);

    // Indicate this type of interface polygons are plain
    m_currentInterfacePolygonsType = INTERFACES::PLAIN;

    // Toggle the view interface polygons options adequately
    this->actionViewPlainInterfaces->setChecked(m_viewInterfacePolygons);
    this->actionViewCompressionForces->setChecked(false);
    this->actionViewTensionForces->setChecked(false);
    this->actionViewUTangentialForces->setChecked(false);
    this->actionViewVTangentialForces->setChecked(false);

    // Initialize the actor that displays the interface polygons
    InitInterfacePolygonsActor(m_workspace->GetInterfacePolygons());
    
    // If the tessellation is not visible then ask whether or not to show it
    if (!m_viewInterfacePolygons)
    {
        AskViewInterfacePolygons();
    }

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionSaveVrmlAssemblyFile_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionSaveVrmlAssemblyFile_triggered");
}

void MainTigerWindow::on_actionSaveVrmlAssemblyFilePerBlock_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionSaveVrmlAssemblyFilePerBlock_triggered");
}

void MainTigerWindow::on_actionSaveVrmlInterfacePolygonsFile_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionSaveVrmlInterfacePolygonsFile_triggered");
}

void MainTigerWindow::on_actionSaveVrmlTessellationFile_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionSaveVrmlTessellationFile_triggered");
}

void MainTigerWindow::on_actionSaveVrmlWorkspaceFile_triggered()
{
    //QMessageBox::information(this, APP_TITLE, "on_actionSaveVrmlWorkspaceFile_triggered");

    // Get the selected filename and trim it
    QString filename = QFileDialog::getSaveFileName(
        this, tr("Save VRML File"), "workspace.wrl", tr("VRML Format (*.wrl)")).trimmed();

    // Exit the function if there is no file name
    if (filename == "")
    {
        return;
    }

    // Initialize the VRML exporter, set the render window, and save its content in the VRML file
    vtkSmartPointer<vtkVRMLExporter> exporter = vtkSmartPointer<vtkVRMLExporter>::New();
    exporter->SetRenderWindow(qvtkWidget->renderWindow());
    exporter->SetFileName(filename.toStdString().c_str());
    exporter->Write();

    QMessageBox::information(this, APP_TITLE, filename + " saved successfully.");
}

void MainTigerWindow::on_actionTicCalculateOverlapping_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionTicCalculateOverlapping_triggered");
}

void MainTigerWindow::on_actionTicEquilibriumAnalysis_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements to run the equilibrium analysis. Show a warning message if a 
    // requirement is missing
    if (!m_workspace->CheckEquilibriumAnalysisRequirements(msg))
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // Open the equilibrium analysis parameters dialog
    EquilibriumAnalysisParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Get the equilibrium analysis parameter values
    double density = dialog.GetDensity();
    double friction = dialog.GetFrictionCoefficient();
    QString lengthUnit = dialog.GetLengthUnit();
    double cWeight = dialog.GetCompressionWeight();
    double tWeight = dialog.GetTensionWeight();
    double uWeight = dialog.GetUTangentialWeight();
    double vWeight = dialog.GetVTangentialWeight();

    bool normalize = dialog.GetNormalize();
    bool verbose = dialog.GetVerbose();
    bool files = dialog.GetFiles();

    // Run the equilibrium analysis
    if (m_workspace->RunEquilibriumAnalysis(
        density,
        friction,
        lengthUnit.toStdString(),
        verbose,
        files,
        cWeight,
        tWeight,
        uWeight,
        vWeight)) 
    {
        std::shared_ptr<EquilibriumAnalysis::Result> results = m_workspace->GetEquilibriumResults();

        // Normalize the force magnitudes if indicated. This normalization divides all force 
        // magnitudes by the density value
        if (normalize) 
        {
            results->Normalize(density);
            std::cout << "Normalized forces" << std::endl;
        }

        if (verbose)
        {
            results->Write();
        }

        // Write the force ranges
        results->WriteForceRanges();

        switch (m_currentInterfacePolygonsType) 
        {
            case INTERFACES::PLAIN: 
            {
                m_currentInterfacePolygonsType = INTERFACES::COMPRESSION;

                this->actionViewPlainInterfaces->setChecked(false);
                this->actionViewCompressionForces->setChecked(m_viewInterfacePolygons);
                this->actionViewTensionForces->setChecked(false);
                this->actionViewUTangentialForces->setChecked(false);
                this->actionViewVTangentialForces->setChecked(false);

                break;
            }

            case INTERFACES::COMPRESSION:
            {
                this->actionViewPlainInterfaces->setChecked(false);
                this->actionViewCompressionForces->setChecked(m_viewInterfacePolygons);
                this->actionViewTensionForces->setChecked(false);
                this->actionViewUTangentialForces->setChecked(false);
                this->actionViewVTangentialForces->setChecked(false);

                break;
            }

            case INTERFACES::TENSION:
            {
                this->actionViewPlainInterfaces->setChecked(false);
                this->actionViewCompressionForces->setChecked(false);
                this->actionViewTensionForces->setChecked(m_viewInterfacePolygons);
                this->actionViewUTangentialForces->setChecked(false);
                this->actionViewVTangentialForces->setChecked(false);

                break;
            }

            case INTERFACES::UTANGENTIAL:
            {
                this->actionViewPlainInterfaces->setChecked(false);
                this->actionViewCompressionForces->setChecked(false);
                this->actionViewTensionForces->setChecked(false);
                this->actionViewUTangentialForces->setChecked(m_viewInterfacePolygons);
                this->actionViewVTangentialForces->setChecked(false);

                break;
            }

            case INTERFACES::VTANGENTIAL:
            {
                this->actionViewPlainInterfaces->setChecked(false);
                this->actionViewCompressionForces->setChecked(false);
                this->actionViewTensionForces->setChecked(false);
                this->actionViewUTangentialForces->setChecked(false);
                this->actionViewVTangentialForces->setChecked(m_viewInterfacePolygons);

                break;
            }
        }

        // Initialize the interface polygons actor usign the result values for coloring the vertices
        InitInterfacePolygonsActor(
            m_workspace->GetInterfacePolygons(), 
            results, 
            m_currentInterfacePolygonsType);

        if (!m_viewInterfacePolygons) 
        {
            AskViewInterfacePolygons();
        }

        // Call Render to show changes
        qvtkWidget->renderWindow()->Render();
    }
    else 
    {
        QMessageBox::critical(this, "Unfeasible System", "System has no solution.");
    }
}

void MainTigerWindow::on_actionTicFixBottomBlocks_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements 
    if (!m_workspace->CheckFixBottomBlocksRequirements(msg))
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    m_workspace->EraseInterfaces();


    Eigen::Vector3d min = Eigen::Vector3d::Zero(), max = Eigen::Vector3d::Zero(), up(0, 1, 0);

    m_workspace->GetTessellation()->Geometry().axisAlignedBoundingBox(min, max);

    toolkit::Plane plane(up, min);

    m_workspace->GetAssembly()->DisableIntersectedBlocks(plane);
    

    // Clear the actors for rendering the assembly and the inteface polygons
    ClearAssemblyActors();
    ClearInteracePolygonsActor();

    // Generate the TIC using the Height-Bisection method. Then, set it as the assembly in the 
    // workspace
    //m_workspace->SetAssemblyUsingHeightBisectionMethod(topHeight, boundary);

    // Initialize the actors for rendering the assembly
    InitAssemblyActors(m_workspace->GetTessellation(), m_workspace->GetAssembly());

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionTicHeightBisectionMethod_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements for generating a TIC using the Height-Bisection method. Show a 
    // warning message if a requirement is missing
    if (!m_workspace->CheckHeightBisectionRequirements(msg))
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // Open the Height-Bisection parameters dialog
    HeightBisectionParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Get the Height-Bisection parameter values
    double topHeight = dialog.GetTopHeight();
    double bottomHeight = dialog.GetBottomHeight();
    bool boundary = dialog.GetBoundary();

    // Clear the actors for rendering the assembly and the inteface polygons
    ClearAssemblyActors();
    ClearInteracePolygonsActor();

    // Generate the TIC using the Height-Bisection method. Then, set it as the assembly in the 
    // workspace
    m_workspace->SetAssemblyUsingHeightBisectionMethod(topHeight, boundary);

    // Initialize the actors for rendering the assembly
    InitAssemblyActors(m_workspace->GetTessellation(), m_workspace->GetAssembly());

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionTicSetCentersAndDirections_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements for setting the centers on the faces of the tessellation. Show a 
    // warning message if a requirement is missing
    if (!m_workspace->CheckSetFaceCentersRequirements(msg))
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // 
    if (!m_workspace->CheckSetEdgeDirectionsRequirements(msg))
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // Open the face center parameters dialog
    CenterDirectionParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // 
    //QString centerValue = dialog.GetCenter();
    int direction = dialog.GetDirection() == "Inside" ? 1 : -1;

    // 
    m_workspace->SetTessellationTileCenters();
    m_workspace->SetTessellationEdgeDirections(direction);

    // 
    InitTileCentersActor(m_workspace->GetTessellation()->DCEL());
    InitEdgeDirectionsActor(m_workspace->GetTessellation()->DCEL());

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionTicTiltingAngleMethod_triggered()
{
    // Initialize the object to store the checking requirements result
    Message msg;

    // Check the requirements for generating a TIC using the Tilting Angle method. Show a warning 
    // message if a requirement is missing
    if (!m_workspace->CheckTiltingAngleRequirements(msg))
    {
        QMessageBox::warning(this, APP_TITLE, msg.getText());
        return;
    }

    // Open the Tilting Angle parameters dialog
    TiltingAngleParametersDialog dialog(this);

    // Exit the function if the user canceled the dialog
    if (dialog.exec() == QDialog::Rejected)
    {
        return;
    }

    // Get the Tilting Angle parameter values
    double angle = dialog.GetAngle();
    QString unit = dialog.GetUnit();

    // 
    if (unit == "Radians")
    {
        angle *= 180.0 / utils::PI;
    }

    // Clear the actors for rendering the assembly and the interface polygons
    ClearAssemblyActors();
    ClearInteracePolygonsActor();

    // Generate the TIC using the Tilting Angle method. Then, set it as the assembly in the 
    // workspace
    m_workspace->SetAssemblyUsingTiltingAngleMethod(angle);

    // Initialize the actors for rendering the assembly
    InitAssemblyActors(m_workspace->GetTessellation(), m_workspace->GetAssembly());

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionViewAssemblyBoundingBoxes_triggered()
{
    QMessageBox::information(this, APP_TITLE, "on_actionViewAssemblyBoundingBoxes_triggered");
}

void MainTigerWindow::on_actionViewAssemblyGeometry_triggered()
{
    m_viewAssemblyGeometry = this->actionViewAssemblyGeometry->isChecked();

    if (m_vtkAssemblyFacesActor) 
    {
        m_vtkAssemblyFacesActor->SetVisibility(m_viewAssemblyGeometry);
        m_vtkAssemblyEdgesActor->SetVisibility(m_viewAssemblyGeometry);
    }
    
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionViewAxes_triggered()
{
	m_viewAxes = this->actionViewAxes->isChecked();

    if (m_vtkAxesActor) 
    {
        m_vtkAxesActor->SetVisibility(m_viewAxes);
    }
    
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionViewCompressionForces_triggered()
{
    m_viewInterfacePolygons = this->actionViewCompressionForces->isChecked();

    this->actionViewPlainInterfaces->setChecked(false);
    this->actionViewTensionForces->setChecked(false);
    this->actionViewUTangentialForces->setChecked(false);
    this->actionViewVTangentialForces->setChecked(false);

    if (m_vtkInterfacePolygonsActor)
    {
        if (m_currentInterfacePolygonsType == INTERFACES::COMPRESSION)
        {
            m_vtkInterfacePolygonsActor->SetVisibility(m_viewInterfacePolygons);
            m_vtkInterfacePolygonsBarActor->SetVisibility(m_viewInterfacePolygons);
        }
        else if(m_workspace->HasEquilibriumResults())
        {
            InitInterfacePolygonsActor(
                m_workspace->GetInterfacePolygons(), 
                m_workspace->GetEquilibriumResults(), 
                INTERFACES::COMPRESSION);

            m_vtkInterfacePolygonsActor->SetVisibility(m_viewInterfacePolygons);
            m_vtkInterfacePolygonsBarActor->SetVisibility(m_viewInterfacePolygons);
        }

        qvtkWidget->renderWindow()->Render();
    }
}

void MainTigerWindow::on_actionViewEdgeDirections_triggered()
{
	m_viewEdgeDirections = this->actionViewEdgeDirections->isChecked();

    if (m_vtkEdgeDirectionsActor) 
    {
        m_vtkEdgeDirectionsActor->SetVisibility(m_viewEdgeDirections);
    }
    
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionViewGrid_triggered()
{
    m_viewGrid = this->actionViewGrid->isChecked();

    if (m_vtkGridActor) 
    {
        m_vtkGridActor->SetVisibility(m_viewGrid);
    }
    
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionViewTessellationBoundingBox_triggered()
{
    m_viewTessellationBoundingBox = this->actionViewTessellationBoundingBox->isChecked();

    if (m_vtkTessellationBoundBoxActor) 
    {
        m_vtkTessellationBoundBoxActor->SetVisibility(m_viewTessellationBoundingBox);
    }

    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionViewTessellationGeometry_triggered()
{
	m_viewTessellationGeometry = this->actionViewTessellationGeometry->isChecked();

    if (m_vtkTessellationEdgesActor) 
    {
        m_vtkTessellationEdgesActor->SetVisibility(m_viewTessellationGeometry);
        m_vtkTessellationTilesActor->SetVisibility(m_viewTessellationGeometry);
    }
    
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionViewTileCenters_triggered()
{
    m_viewTileCenters = this->actionViewTileCenters->isChecked();

    if (m_vtkTileCentersActor) 
    {
        m_vtkTileCentersActor->SetVisibility(m_viewTileCenters);
    }

    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionViewUTangentialForces_triggered()
{
    m_viewInterfacePolygons = this->actionViewUTangentialForces->isChecked();

    this->actionViewPlainInterfaces->setChecked(false);
    this->actionViewCompressionForces->setChecked(false);
    this->actionViewTensionForces->setChecked(false);
    this->actionViewVTangentialForces->setChecked(false);

    if (m_vtkInterfacePolygonsActor)
    {
        if (m_currentInterfacePolygonsType == INTERFACES::UTANGENTIAL)
        {
            m_vtkInterfacePolygonsActor->SetVisibility(m_viewInterfacePolygons);
            m_vtkInterfacePolygonsBarActor->SetVisibility(m_viewInterfacePolygons);
        }
        else if (m_workspace->HasEquilibriumResults())
        {
            InitInterfacePolygonsActor(
                m_workspace->GetInterfacePolygons(),
                m_workspace->GetEquilibriumResults(),
                INTERFACES::UTANGENTIAL);

            m_vtkInterfacePolygonsActor->SetVisibility(m_viewInterfacePolygons);
            m_vtkInterfacePolygonsBarActor->SetVisibility(m_viewInterfacePolygons);
        }

        qvtkWidget->renderWindow()->Render();
    }
}

void MainTigerWindow::on_actionViewVTangentialForces_triggered()
{
    m_viewInterfacePolygons = this->actionViewVTangentialForces->isChecked();

    this->actionViewPlainInterfaces->setChecked(false);
    this->actionViewCompressionForces->setChecked(false);
    this->actionViewTensionForces->setChecked(false);
    this->actionViewUTangentialForces->setChecked(false);

    if (m_vtkInterfacePolygonsActor)
    {
        if (m_currentInterfacePolygonsType == INTERFACES::VTANGENTIAL)
        {
            m_vtkInterfacePolygonsActor->SetVisibility(m_viewInterfacePolygons);
            m_vtkInterfacePolygonsBarActor->SetVisibility(m_viewInterfacePolygons);
        }
        else if (m_workspace->HasEquilibriumResults())
        {
            InitInterfacePolygonsActor(
                m_workspace->GetInterfacePolygons(),
                m_workspace->GetEquilibriumResults(),
                INTERFACES::VTANGENTIAL);

            m_vtkInterfacePolygonsActor->SetVisibility(m_viewInterfacePolygons);
            m_vtkInterfacePolygonsBarActor->SetVisibility(m_viewInterfacePolygons);
        }

        qvtkWidget->renderWindow()->Render();
    }
}

void MainTigerWindow::on_actionViewHideAll_triggered()
{
    if (m_viewAssemblyGeometry)
    {
        if (m_vtkAssemblyEdgesActor) 
        {
            m_vtkAssemblyEdgesActor->VisibilityOff();
            m_vtkAssemblyFacesActor->VisibilityOff();
        }
        
        m_viewAssemblyGeometry = false;
        this->actionViewAssemblyGeometry->setChecked(false);
    }

    if (m_viewAssemblyBlocksBoundingBoxes) 
    {
        m_viewAssemblyBlocksBoundingBoxes = false;
        this->actionViewAssemblyBoundingBoxes->setChecked(false);
    }

    if (m_viewAxes) 
    {
        if (m_vtkAxesActor) 
        {
            m_vtkAxesActor->VisibilityOff();
        }
        
        m_viewAxes = false;
        this->actionViewAxes->setChecked(false);
    }

    if (m_viewEdgeDirections) 
    {
        if (m_vtkEdgeDirectionsActor) 
        {
            m_vtkEdgeDirectionsActor->VisibilityOff();
        }
        
        m_viewEdgeDirections = false;
        this->actionViewEdgeDirections->setChecked(false);
    }

    if (m_viewEdgeRotatedVectors) 
    {
        m_viewEdgeRotatedVectors = false;
        //this->actionViewEdgeRotatedVectors->setChecked(false);
    }

    if (m_viewTileCenters) 
    {
        if (m_vtkTileCentersActor) 
        {
            m_vtkTileCentersActor->VisibilityOff();
        }
        
        m_viewTileCenters = false;
        this->actionViewTileCenters->setChecked(false);
    }

    if (m_viewTessellationGeometry) 
    {
        if (m_vtkTessellationEdgesActor) 
        {
            m_vtkTessellationEdgesActor->VisibilityOff();
            m_vtkTessellationTilesActor->VisibilityOff();
        }
        
        m_viewTessellationGeometry = false;
        this->actionViewTessellationGeometry->setChecked(false);
    }

    if (m_viewTessellationBoundingBox) 
    {
        if (m_vtkTessellationBoundBoxActor) 
        {
            m_vtkTessellationBoundBoxActor->VisibilityOff();
        }

        m_viewTessellationBoundingBox = false;
        this->actionViewTessellationBoundingBox->setChecked(false);
    }

    if (m_viewGrid) 
    {
        if (m_vtkGridActor) 
        {
            m_vtkGridActor->VisibilityOff();
        }
        
        m_viewGrid = false;
        this->actionViewGrid->setChecked(false);
    }

    if (m_viewInterfacePolygons) 
    {
        if (m_vtkInterfacePolygonsActor) 
        {
            m_vtkInterfacePolygonsActor->VisibilityOff();
        }

        if (m_vtkInterfacePolygonsBarActor)
        {
            m_vtkInterfacePolygonsBarActor->VisibilityOff();
        }
        
        m_viewInterfacePolygons = false;
        
        this->actionViewPlainInterfaces->setChecked(false);
        this->actionViewCompressionForces->setChecked(false);
        this->actionViewTensionForces->setChecked(false);
        this->actionViewUTangentialForces->setChecked(false);
        this->actionViewVTangentialForces->setChecked(false);
    }

    if (m_viewSectionHeights) 
    {
        // Remove actor
        m_viewSectionHeights = false;
        //this->actionViewSectionHeights->setChecked(false);
    }

    // Call Render to show changes
    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::on_actionViewPlainInterfaces_triggered()
{
    m_viewInterfacePolygons = this->actionViewPlainInterfaces->isChecked();
    
    this->actionViewCompressionForces->setChecked(false);
    this->actionViewTensionForces->setChecked(false);
    this->actionViewUTangentialForces->setChecked(false);
    this->actionViewVTangentialForces->setChecked(false);

    if (m_vtkInterfacePolygonsBarActor) 
    {
        m_vtkInterfacePolygonsBarActor->VisibilityOff();
    }

    if (m_vtkInterfacePolygonsActor) 
    {
        if (m_currentInterfacePolygonsType == INTERFACES::PLAIN) 
        {
            m_vtkInterfacePolygonsActor->SetVisibility(m_viewInterfacePolygons);
        }
        else 
        {
            InitInterfacePolygonsActor(m_workspace->GetInterfacePolygons());

            m_vtkInterfacePolygonsActor->SetVisibility(m_viewInterfacePolygons);
        }

        qvtkWidget->renderWindow()->Render();
    }
}

void MainTigerWindow::on_actionViewTensionForces_triggered()
{
    m_viewInterfacePolygons = this->actionViewTensionForces->isChecked();

    this->actionViewPlainInterfaces->setChecked(false);
    this->actionViewCompressionForces->setChecked(false);
    this->actionViewUTangentialForces->setChecked(false);
    this->actionViewVTangentialForces->setChecked(false);

    if (m_vtkInterfacePolygonsActor)
    {
        if (m_currentInterfacePolygonsType == INTERFACES::TENSION)
        {
            m_vtkInterfacePolygonsActor->SetVisibility(m_viewInterfacePolygons);
            m_vtkInterfacePolygonsBarActor->SetVisibility(m_viewInterfacePolygons);
        }
        else if (m_workspace->HasEquilibriumResults())
        {
            InitInterfacePolygonsActor(
                m_workspace->GetInterfacePolygons(),
                m_workspace->GetEquilibriumResults(),
                INTERFACES::TENSION);

            m_vtkInterfacePolygonsActor->SetVisibility(m_viewInterfacePolygons);
            m_vtkInterfacePolygonsBarActor->SetVisibility(m_viewInterfacePolygons);
        }

        qvtkWidget->renderWindow()->Render();
    }
}

/*void MainTigerWindow::on_actionViewInterfacePolygons_triggered()
{
    m_viewInterfacePolygons = this->actionViewInterfacePolygons->isChecked();

    if (m_vtkInterfacePolygonsActor) 
    {
        m_vtkInterfacePolygonsActor->SetVisibility(m_viewInterfacePolygons);
    }
    
    if (m_vtkInterfacePolygonsBarActor) 
    {
        m_vtkInterfacePolygonsBarActor->SetVisibility(m_viewInterfacePolygons);
    }
    
    qvtkWidget->renderWindow()->Render();
}*/

void MainTigerWindow::on_actionViewShowAll_triggered()
{
    if (!m_viewAssemblyGeometry)
    {
        if (m_vtkAssemblyFacesActor) 
        {
            m_vtkAssemblyFacesActor->VisibilityOn();
            m_vtkAssemblyEdgesActor->VisibilityOn();
        }
        
        m_viewAssemblyGeometry = true;
        this->actionViewAssemblyGeometry->setChecked(true);
    }

    if (!m_viewAssemblyBlocksBoundingBoxes) 
    {
        m_viewAssemblyBlocksBoundingBoxes = true;
        this->actionViewAssemblyBoundingBoxes->setChecked(true);
    }

    if (!m_viewAxes)
    {
        if (m_vtkAxesActor) 
        {
            m_vtkAxesActor->VisibilityOn();
        }
        
        m_viewAxes = true;
        this->actionViewAxes->setChecked(true);
    }

    if (!m_viewEdgeDirections)
    {
        if (m_vtkEdgeDirectionsActor) 
        {
            m_vtkEdgeDirectionsActor->VisibilityOn();
        }
        
        m_viewEdgeDirections = true;
        this->actionViewEdgeDirections->setChecked(true);
    }

    if (!m_viewEdgeRotatedVectors)
    {
        m_viewEdgeRotatedVectors = true;
        //this->actionViewEdgeRotatedVectors->setChecked(true);
    }

    if (!m_viewTileCenters)
    {
        if (m_vtkTileCentersActor) 
        {
            m_vtkTileCentersActor->VisibilityOn();
        }
        
        m_viewTileCenters = true;
        this->actionViewTileCenters->setChecked(true);
    }

    ViewTessellationGeometry();

    if (!m_viewTessellationBoundingBox) 
    {
        if (m_vtkTessellationBoundBoxActor) 
        {
            m_vtkTessellationBoundBoxActor->VisibilityOn();
        }

        m_viewTessellationBoundingBox = true;
        this->actionViewTessellationBoundingBox->setChecked(true);
    }

    if (!m_viewGrid)
    {
        if (m_vtkGridActor) 
        {
            m_vtkGridActor->VisibilityOn();
        }
        
        m_viewGrid = true;
        this->actionViewGrid->setChecked(true);
    }

    if (!m_viewInterfacePolygons)
    {
        if (m_vtkInterfacePolygonsActor) 
        {
            m_vtkInterfacePolygonsActor->VisibilityOn();
        }

        if (m_vtkInterfacePolygonsBarActor)
        {
            m_vtkInterfacePolygonsBarActor->VisibilityOn();
        }
        
        m_viewInterfacePolygons = true;

        switch (m_currentInterfacePolygonsType) 
        {
            case INTERFACES::PLAIN: 
            {
                this->actionViewPlainInterfaces->setChecked(true);
                break;
            }

            case INTERFACES::COMPRESSION:
            {
                this->actionViewCompressionForces->setChecked(true);
                break;
            }

            case INTERFACES::TENSION:
            {
                this->actionViewTensionForces->setChecked(true);
                break;
            }

            case INTERFACES::UTANGENTIAL:
            {
                this->actionViewUTangentialForces->setChecked(true);
                break;
            }

            case INTERFACES::VTANGENTIAL:
            {
                this->actionViewVTangentialForces->setChecked(true);
                break;
            }
        }
    }

    if (!m_viewSectionHeights)
    {
        m_viewSectionHeights = true;
        //this->actionViewSectionHeights->setChecked(true);
    }

    qvtkWidget->renderWindow()->Render();
}

void MainTigerWindow::slotExit()
{
	qApp->exit();
}

bool MainTigerWindow::ValidateWorkspacesDirectory(bool mkdir) const
{
    QDir dir(m_workspacesDir);

    return dir.exists() ? true : (mkdir ? dir.mkdir(m_workspacesDir) : false);
}

TilePickerStyle::TilePickerStyle() : m_tigerWindow(nullptr)
{
}

TilePickerStyle::~TilePickerStyle()
{
    m_tigerWindow = nullptr;
}

void TilePickerStyle::OnLeftButtonDown()
{
    assert(m_tigerWindow);

    // Get the location of the click in window coordinates
    int* position = this->GetInteractor()->GetEventPosition();
    std::cout << "Left click at (" << position[0] << ", " << position[1] << ")" << std::endl;

    //vtkIdType id = m_tigerWindow->GetPickedTileCellId(position[0], position[1]);
    //std::cout << "ID = " << id << std::endl;

    // Initialize the cell picker and define it
    //vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
    //picker->SetTolerance(0.0005);

    // 
    //picker->Pick(position[0], position[1], 0, this->GetDefaultRenderer());

    //double* scenePosition = picker->GetPickPosition();
    //std::cout << "Picked cell " << picker->GetCellId() << std::endl;

    //if (picker->GetCellId() != -1) 
    //{
    //    std::cout << "Pick position is " << scenePosition[0] << ", " << scenePosition[1] << ", " << scenePosition[2] << std::endl;
    //}

    // Forward the event
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}

void TilePickerStyle::OnRightButtonDown()
{
    // Get the location of the click in window coordinates
    int* position = this->GetInteractor()->GetEventPosition();
    std::cout << "Right click at (" << position[0] << ", " << position[1] << ")" << std::endl;

    // Initialize the cell picker and define it
    //vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
    //picker->SetTolerance(0.0005);

    // 
    //picker->Pick(position[0], position[1], 0, this->GetDefaultRenderer());

    //double* scenePosition = picker->GetPickPosition();
    //std::cout << "Picked cell " << picker->GetCellId() << std::endl;

    //if (picker->GetCellId() != -1) 
    //{
    //    std::cout << "Pick position is " << scenePosition[0] << ", " << scenePosition[1] << ", " << scenePosition[2] << std::endl;
    //}

    // Forward the event
    vtkInteractorStyleTrackballCamera::OnRightButtonDown();
}

void TilePickerStyle::SetTigerWindow(MainTigerWindow * tigerWindow)
{
    m_tigerWindow = tigerWindow;
}
