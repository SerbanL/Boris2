#pragma once

#include "CompileFlags.h"
#if GRAPHICS == 1

#include "BorisLib.h"
#include "TextObject.h"
#include "BorisGraphics.h"
#include "WinSpaces.h"
#include "PhysQRep.h"
#include "Display_Defs.h"



////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// MESH WINDOW (graphical display)

class BorisMeshWindow :
	public WinSpace,
	public ProgramState<BorisMeshWindow,
	std::tuple<PhysQRep, DBL3, std::string, std::string, std::string, std::string, bool>, std::tuple<>>
{
	friend class BorisDisplay;

private:

	//the computed graphical representation of a physical quantity (e.g. representation of magnetization, temperature etc.)
	PhysQRep physQRep;

	//Text Format for drawing text in mesh window
	IDWriteTextFormat *pMeshWinTextFormat = nullptr;

	//when mouse hovers over mesh values are picked and set here to be displayed
	DBL3 meshPosition;
	std::string meshName, typeName, mouse_mesh_info_position, mouse_mesh_info_value;
	bool displayMeshInfo = false;

private:

	void ZoomNotchUp(void);
	void ZoomNotchDn(void);
	//rather than changing detail level using mouse wheel, set it directly using this
	void SetDetailLevel(double detail_level);
	double GetDetailLevel(void);

	//Set mesh display render threshold values for faster rendering when we have many cells
	void SetRenderThresholds(INT3 renderthresholds);
	INT3 GetRenderThresholds(void);

protected:

	void AdjustWindowDimensions(D2D1_RECT_F delta) { ChangeWindowSides(delta); }

	//compute the physical quantity representation for the given physical quantity
	void UpdatePhysQRep(std::vector<PhysQ> physQ);

	//adjust mesh display camera fov and averaging stencil dimension for default physical quantity representation
	void AutoSetMeshDisplaySettings(std::vector<PhysQ> physQ);

	double Get_PhysQRep_DetailLevel(void) { return physQRep.get_detail_level(); }

	void DrawWindow(void);
	void DrawWindow_Quick(void) { DrawWindow(); }

	ActionOutcome NewMessage(AC_ aCode, INT2 mouse, std::string data = "");

	void SetDefaultCursor(void) { SetCursor(LoadCursor(nullptr, IDC_ARROW)); }

	//image_cropping specify normalized cropping within the mesh image, as left, bottom, right, top : 0, 0 point is left, bottom of screen as far as user is concerned.
	bool SaveMeshImage(std::string fileName, DBL4 image_cropping);

	bool SaveImage(std::string fileName, std::vector<PhysQ> physQ);

public:

	//Setup new window space with rectangle specified by ratios of the entire screen
	BorisMeshWindow(D2D1_RECT_F ratios, INT2 winId);
	~BorisMeshWindow();

	void RepairObjectState(void) {}
};

#endif
