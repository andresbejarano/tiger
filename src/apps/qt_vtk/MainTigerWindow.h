#ifndef _MAIN_TIGER_WINDOW_H_
#define _MAIN_TIGER_WINDOW_H_

#pragma once

#include <QMainWindow>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkFollower.h>
#include <vtkWeakPointer.h>
#include <vtkScalarBarActor.h>
#include <tiger/workspace/workspace.h>
#include "ui_MainTigerWindow.h"

//
// The class representing the main TIGER Qt VTK app window. This class handles 
// both Qt UI and VTK elements.
//
class MainTigerWindow : public QMainWindow, private Ui::MainTigerWindow
{
  Q_OBJECT

private:

    // The different content types for the interface polygons actor
    enum class INTERFACES 
    { 
        PLAIN, 
        COMPRESSION, 
        TENSION, 
        UTANGENTIAL, 
        VTANGENTIAL 
    };

    // The type of current interface polygon
    INTERFACES m_currentInterfacePolygonsType;

    // The path to the last open directory by the user
    QString m_lastDirectory;

    // Indicates if the blocks bounding boxes of the assembly are visible
    bool m_viewAssemblyBlocksBoundingBoxes;

    // Indicates if the indices of the blocks are visible
    bool m_viewAssemblyBlockIndices;

	// Indicates if the geometry of the assembly is visible
	bool m_viewAssemblyGeometry;

	// Indicates if the basis vectors representing the (X, Y, Z) axes are 
    // visible
	bool m_viewAxes;

	// Indicates if the edge directions are visible
	bool m_viewEdgeDirections;

	// Indicates if the edge rotated vectors are visible
	bool m_viewEdgeRotatedVectors;

	// Indicates if the tile centes are visible
	bool m_viewTileCenters;

    // Indicates if the bounding box of the tessellation is visible
    bool m_viewTessellationBoundingBox;

	// Indicates if the geometry of the tessellation is visible
	bool m_viewTessellationGeometry;

	// Indicates if the XZ reference grid is visible
	bool m_viewGrid;

	// Indicates if the interface polygons are visible
	bool m_viewInterfacePolygons;

	// Indicates if the section height parameter values are visible
	bool m_viewSectionHeights;

    // The actor to render the bounding boxes of the blocks of the assembly
    vtkWeakPointer<vtkActor> m_vtkAssemblyBlocksBoundingBoxesActor;

    // The vector with the actors to render the block indices of the assembly
    std::vector<vtkWeakPointer<vtkFollower>> m_vtkAssemblyBlockIndicesActor;

    // The actor to render the edges of the assembly (blocks)
    vtkWeakPointer<vtkActor> m_vtkAssemblyEdgesActor;

    // The actor to render the faces of the assembly (blocks)
    vtkWeakPointer<vtkActor> m_vtkAssemblyFacesActor;

    // The actor to render the reference axes
    vtkWeakPointer<vtkActor> m_vtkAxesActor;

    // The camera of the VTK scene
    vtkWeakPointer<vtkCamera> m_vtkCamera;

    // The actor to render the edge directions in the tessellation
    vtkWeakPointer<vtkActor> m_vtkEdgeDirectionsActor;

    // The actor to render the reference grid
    vtkWeakPointer<vtkActor> m_vtkGridActor;

    // The actor to render the interface polygons
    vtkWeakPointer<vtkActor> m_vtkInterfacePolygonsActor;

    // The actor to render the scalar bar for the force magnitudes on the 
    // interface polygons
    vtkWeakPointer<vtkScalarBarActor> m_vtkInterfacePolygonsBarActor;

    // The actor to render the bounding box of the tessellation
    vtkWeakPointer<vtkActor> m_vtkTessellationBoundBoxActor;

    // The actor to render the edges of the tessellation
    vtkWeakPointer<vtkActor> m_vtkTessellationEdgesActor;

    // The actor to render the tiles of the tessellation
    vtkWeakPointer<vtkActor> m_vtkTessellationTilesActor;

    // The actor to render the center of the tiles of the tessellation
    vtkWeakPointer<vtkActor> m_vtkTileCentersActor;

    // The renderer of the VTK scene
    vtkWeakPointer<vtkRenderer> m_vtkRenderer;

	// The pointer to the workspace
	std::shared_ptr<Workspace> m_workspace;

    // The path to the directory to store the workspaces information
    QString m_workspacesDir;

public:
	
	//
	// Constructor of the class.
	//
	MainTigerWindow();
	
	//
	// Destructor of the class.
	//
	~MainTigerWindow();

	//
	// Indicates if the basis vectors representing the (X, Y, Z) axes are 
    // visible.
	// @return bool The view state of the axes.
	//
	bool AreAxesVisible() const;

	//
	// Indicates if the edge directions are visible.
	// @return bool The view state of the edge directions.
	//
	bool AreEdgeDirectionsVisible() const;

	//
	// Indicates if the edge rotated vectors are visible.
	// @return bool The view state of the edge rotated vectors.
	//
	bool AreEdgeRotatedVectorsVisible() const;

	//
	// Indicates if the face centers are visible.
	// @return bool The view state of the face centers.
	//
	bool AreFaceCentersVisible() const;

	//
	// Indicates if the piece interface polygons are visible.
	// @return bool The view state of the interface polygons.
	//
	bool AreInterfacePolygonsVisible() const;

	//
	// Indicates if the section height parameter values are visible.
	// @return bool The view state of the section height parameter values.
	//
	bool AreSectionHeightsVisible() const;

private:

    //
    //
    void AskViewInterfacePolygons();

    //
    // Asks whether to show the tessellation geometryor not. If yes then it 
    // visualizes the tessellation geometry actor (if any). This function does 
    // not call the render function.
    //
    void AskViewTessellationGeometry();

    //
    // Clears the actors to render the assembly. There are three actors: 
    // an actor for the faces
    // an actor for the edges
    // an actor for the block indices
    //
    void ClearAssemblyActors();

    //
    // Clears the actors to render the tessellation. There are two actors: one 
    // for the tiles and one for the edges.
    //
    void ClearTessellationActors();

    //
    // Clears the actor to render the interface polygons.
    //
    void ClearInteracePolygonsActor();

public:

	//
	// Clears the content of the workspace, including the actors associated to 
    // its elements.
	//
	void ClearWorkspace();

private:

    //
    // Returns the number of workspaces stored in the workspaces directory.
    // @return size_t The number of workspaces.
    //
    size_t CountWorkspaces() const;

    //
    // Deletes the actors associated to the elements of a workspace.
    //
    void DeleteWorkspaceActors();

private:

    //
    // @return std::string
    //
    std::string GetInterfaceTypeName(INTERFACES type) const;

    //
    // Initializes the actors to render the information of the assembly. 
    // NOTE: Call this function only after adding the renderer to the VTK 
    // render window of the Qt VTK widget.
    // @param const std::shared_ptr<Assembly> assembly The pointer to the 
    // assembly.
    //
    void InitAssemblyActors(
        const std::shared_ptr<Tessellation> tessellation, 
        const std::shared_ptr<Assembly> assembly);

    // 
    // 
    // 
    void InitAssemblyBlockIndicesActor();

	//
	// Initializes the actor to render the axes vectors. Reference axes point 
    // along X, Y and Z directions. NOTE: Call this function only after adding
    // the renderer to the VTK render window of the Qt VTK widget.
	//
	void InitAxesActor();

    //
    // Initialized a vtkActor with the edges of the axis aligned bounding box 
    // of a geometry.
    // @param VF & vf The reference to the geometry.
    // @return vtkWeakPointer<vtkActor> The pointer to the actor.
    //
    vtkWeakPointer<vtkActor> InitAxisAlignedBoundingBoxActor(const VF & vf);

    //
    // Initializes the actor to render the centers of the tiles of the 
    // tessellation. Each tile center is represented using an octahedron. 
    // NOTE: Call this function only after adding the renderer to the VTK 
    // render window of the Qt VTK widget.
    // @param const dcel::DCEL & domain The reference to the geometry of the 
    // tessellation. The tiles must have the ATTRIB_CENTER dynamic attribute.
    // @param double radius The radius of the octahedra.
    //
    void InitTileCentersActor(
        const std::shared_ptr<dcel::DCEL> dcel, 
        double radius = 0.1);

    //
    // Initializes the actor to render the edge directions in the 
    // tessellation. 
    // NOTE: Call this function only after adding the renderer to the VTK 
    // render window of the Qt VTK widget.
	// @param const dcel::DCEL & domain The reference to the geometry of the 
    // tessellation.
    // @param double length The length of the lines representing the edge 
    // directions.
    //
    void InitEdgeDirectionsActor(
        const std::shared_ptr<dcel::DCEL> dcel, 
        double length = 0.5);

	//
	// Initializes the actor to render the tessellation. NOTE: Call this 
    // function only after adding the renderer to the VTK render window of the
    // Qt VTK widget.
	// @param const dcel::DCEL & domain The reference to the geometry of the 
    // tessellation.
	//
	//void InitTessellationActor(const dcel::DCEL & domain);

    //
    // Initializes the actors to render the content of the tessellation. 
    // NOTE: Call this function only after adding the renderer to the VTK 
    // render window of the Qt VTK widget.
    // @param const VF & vf The reference to the geometry of the tessellation.
    //
    void InitTessellationActors(const VF & vf);

	//
	// Initializes the actor to render the reference grid. Reference grid lies
    // on the XY plane. 
	// NOTE: Call this function only after adding the renderer to the VTK 
    // render window of the Qt VTK widget.
	// @param size_t width
	// @param size_t height
	// @param size_t widthSegments
	// @param size_t heightSegments
	//
	void InitGridActor(
		size_t width = 2, 
		size_t height = 2, 
		size_t widthSegments = 20, 
		size_t heightSegments = 20);

    //
    // Initializes the actor to render the interface polygons. NOTE: Call this
    // function only after adding the renderer to the VTK render window of the
    // Qt VTK widget.
    // @param const std::shared_ptr<InterfacePolygons> intfs The pointer to 
    // the interface polygons.
    //
    void InitInterfacePolygonsActor(
        const std::shared_ptr<InterfacePolygons> intfs);

    //
    // @param const std::shared_ptr<InterfacePolygons> intfs
    // const std::shared_ptr<EquilibriumAnalysis::Result> results
    // @param INTERFACES type
    //
    void InitInterfacePolygonsActor(
        const std::shared_ptr<InterfacePolygons> intfs,
        const std::shared_ptr<EquilibriumAnalysis::Result> results, 
        INTERFACES type,
        bool useForceCaps = false,
        double minForceCap = 1.0,
        double maxForceCap = 1.0);

    //
    //
    void SetInterfacePolygons(bool render = false);

    //
    // Adds the interface polygons actor (if any) to the renderer and 
    // indicates the interfaces are being shown. This function does not call 
    // the render function.
    // @param INTERFACES type
    //
    void ViewInterfacePolygons(INTERFACES type);

    //
    // Shows the tessellation geometry actor (if any). This function does not 
    // call the render function.
    //
    void ViewTessellationGeometry();

public slots:

    //
    //
    void on_actionCameraBack_triggered();

    //
    //
    void on_actionCameraBottom_triggered();

    //
    //
    void on_actionCameraFront_triggered();

    //
    //
    void on_actionCameraLeft_triggered();

    //
    //
    void on_actionCameraManual_triggered();

    //
	// Slot function called when the user clicks the camera reset menu item.
	//
	void on_actionCameraReset_triggered();

    //
    //
    void on_actionCameraRight_triggered();

    //
    //
    void on_actionCameraTop_triggered();

    //
    // Slot function called when the user clicks the tessellation dual menu 
    // item.
    //
    void on_actionEditDualTessellation_triggered();

    //
    // Slot function executed when the user clicks the edit assembly adaptive 
    // tile offset menu item.
    //
    void on_actionEditAssemblyClipAdaptiveTileOffset_triggered();

    //
    // Slot function called when the user clicks the edit assembly clip tile 
    // offset menu item.
    //
    void on_actionEditAssemblyClipTileOffset_triggered();

    //
    //
    void on_actionEditFlipTessellation_triggered();

    //
    // Slot function called when the user clicks the center to origin menu 
    // item.
    //
    void on_actionEditTessellationCentroidToOrigin_triggered();

    //
    // Slot function called when the user clicks the face subdivision menu 
    // item.
    //
    void on_actionEditTileSubdivision_triggered();

    //
    // Slot function called when the user clicks the normalize tessellation 
    // menu item.
    //
    void on_actionEditNormalizeVertices_triggered();

    //
    //
    void on_actionEditRestoreAssembly_triggered();

    //
    // Slot function called when the user clicks the rotate tessellation menu 
    // item.
    //
    void on_actionEditRotateTessellation_triggered();

    //
    // Slot function called when the user clicks the scale tessellation menu 
    // item.
    //
    void on_actionEditScaleTessellation_triggered();

    //
    //
    void on_actionNewArchimedeanSolidGeometry_triggered();

    // 
    // Slot function called when the user clicks the new barrel vault geometry
    // menu item.
    // 
    void on_actionNewBarrelVaultGeometry_triggered();

    //
    // Slot function called when the user clicks the new bended square 
    // geometry menu item.
    //
    void on_actionNewBendedSquareGeometry_triggered();

    //
    //
    void on_actionNewConeGeometry_triggered();

    //
    //
    //
    void on_actionNewCyclideGeometry_triggered();

    //
    // Slot function called when the user clicks the new cylinder geometry 
    // menu item.
    //
    void on_actionNewCylinderGeometry_triggered();

    //
    // Slot function called when the user clicks the new equiliateral triangle
    // geometry menu item.
    //
    void on_actionNewEquilateralTriangleGeometry_triggered();

    //
    // Slot function executed when the user clicks the new N-sided regular 
    // polygon geometry menu item.
    //
    void on_actionNewNsidedRegularPolygonGeometry_triggered();

    //
    // Slot function executed when the user clicks the new tessellation from 
    // an OBJ file menu item.
    //
    void on_actionNewObjFileGeometry_triggered();

    //
    // Slot function called when the user clicks the new elliptic paraboloid 
    // menu item.
    //
    void on_actionNewParaboloidGeometry_triggered();

    //
    // Slot function called when the user clicks the new Platonic solid 
    // geometry menu item.
    //
    void on_actionNewPlatonicSolidGeometry_triggered();

    //
    // Slot function called when the user clicks the new Rectangle geometry 
    // menu.
    //
    void on_actionNewRectangleGeometry_triggered();

	//
    // Slot function called when the user clicks the new square geometry menu 
    // item.
	//
	void on_actionNewSquareGeometry_triggered();

    //
    // Slot function called when the user clicks the new square grid geometry 
    // menu item.
    //
    void on_actionNewSquareGridGeometry_triggered();

    //
    // Slot function called when the user clicks the new torus geometry menu 
    // item.
    //
    void on_actionNewTorusGeometry_triggered();

    //
    // Slot function called when the user clicks the new truncated cone 
    // geometry menu item.
    //
    void on_actionNewTruncatedConeGeometry_triggered();

	//
	//
	void on_actionNewWorkspace_triggered();

    //
    //
    void on_actionSaveAbaqusFile_triggered();

    //
    //
    void on_actionSaveGeogebraJsTessellationFile_triggered();

    //
    //
    void on_actionSaveObjAllAssemblyFile_triggered();

    //
    //
    void on_actionSaveObjAllInterfacePolygonsFile_triggered();

    //
    //
    void on_actionSaveObjAssemblyFilePerBlock_triggered();

    //
    //
    void on_actionSaveObjTessellationFile_triggered();

    //
    //
    void on_actionSaveObjWorkspaceFile_triggered();

    //
    //
    void on_actionSaveVrmlAssemblyFile_triggered();

    //
    //
    void on_actionSaveVrmlAssemblyFilePerBlock_triggered();

    //
    //
    void on_actionSaveVrmlInterfacePolygonsFile_triggered();

    //
    //
    void on_actionSaveVrmlTessellationFile_triggered();

    //
    //
    void on_actionSaveVrmlWorkspaceFile_triggered();

    // 
    // 
    // 
    void on_actionTicAddBlockForceTorqueLoads_triggered();

    // 
    // 
    // 
    void on_actionTicAddBlockWeightLoads_triggered();

    // 
    // 
    // 
    void on_actionTicCalculateInterfaces_triggered();

    // 
    // 
    // 
    void on_actionTicCalculateOverlapping_triggered();

    // 
    // 
    // 
    //void on_actionTicEnableDisableAllBlocks_triggered();

    // 
    // Slot function called when the user clicks the enable/disable specific 
    // blocks menu item.
    // 
    void on_actionTicEnableDisableSpecificBlocks_triggered();

    // 
    // 
    // 
    void on_actionTicFixBottomBlocks_triggered();

    //
    // Slot function called when the user clicks the TIC Height-Bisection 
    // method menu item.
    //
    void on_actionTicHeightBisectionMethod_triggered();

    // 
    // Slot function called when the user clicks the TIC reset block loads 
    // menu item.
    // 
    void on_actionTicResetBlockLoads_triggered();

    //
    // Slot function called when the user clicks the TIC set center and 
    // direction menu item.
    //
    void on_actionTicSetCentersAndDirections_triggered();

    //
    // Slot function executed when the user clicks the Equilibrium Analysis 
    // menu item.
    //
    void on_actionTicStaticEquilibriumAnalysis_triggered();

    //
    // Slot function called when the user clicks the TIC Tilting Angle method 
    // menu item.
    //
    void on_actionTicTiltingAngleMethod_triggered();

    // 
    // Slot function called when the user toggles the visualization of the 
    // block indices.
    // 
    void on_actionViewAssemblyBlockIndices_triggered();

    //
    //
    void on_actionViewAssemblyBoundingBoxes_triggered();

	//
	// Slot function called when the user toggles the visualization of the 
    // assembly.
	//
	void on_actionViewAssemblyGeometry_triggered();

	//
	// Slot function called when the user toggles the visualization of the 
    // (X, Y, Z) axes.
	//
	void on_actionViewAxes_triggered();

    //
    // Slot function called when the user clicks the cap forces menu item.
    //
    void on_actionViewCapInterfaceForces_triggered();

    //
    //
    void on_actionViewCompressionForces_triggered();

	//
	// Slot function called when the user toggles the visualization of the 
    // edge directions.
	//
	void on_actionViewEdgeDirections_triggered();

    //
    // Slot function called when the user toggles the visualization of the 
    // reference grid.
    //
    void on_actionViewGrid_triggered();

    //
    // Slot function called when the user hides all workspace elements.
    //
    void on_actionViewHideAll_triggered();

    //
    //
    void on_actionViewPlainInterfaces_triggered();

    //
    //
    void on_actionViewTensionForces_triggered();

    //
    //
    void on_actionViewTessellationBoundingBox_triggered();

	//
	// Slot function called when the user toggles the visualization of the 
    // tessellation.
	//
	void on_actionViewTessellationGeometry_triggered();

    //
    // Slot function called when the user toggles the visualization of the 
    // face centers.
    //
    void on_actionViewTileCenters_triggered();

    //
    //
    void on_actionViewUTangentialForces_triggered();

    //
    //
    void on_actionViewVTangentialForces_triggered();
	
	//
	// Slot function called when the user shows all workspace elements.
	//
	void on_actionViewShowAll_triggered();
	
	//
	//
	virtual void slotExit();

 private:

     //
     // Checks the directory for storing the information of the workspaces 
     // exists.
     // @param bool mkdir Indicates whether to make the directory if it 
     // doesn't exist.
     // @return bool Indicates whether the workspaces directory exists or not.
     //
     bool ValidateWorkspacesDirectory(bool mkdir = false) const;

};

#endif
