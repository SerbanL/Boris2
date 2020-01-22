#pragma once

#include "CompileFlags.h"
#if GRAPHICS == 1

#include "BorisLib.h"
#include "TextObject.h"
#include "BorisGraphics.h"
#include "WinSpaces.h"
#include "PhysQRep.h"
#include "Display_Defs.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// MESH WINDOW (graphical display)

class BorisMeshWindow :
	public WinSpace,
	public ProgramState<BorisMeshWindow,
	tuple<PhysQRep, DBL3, string, string, string, string, bool>, tuple<>>
{
	friend class BorisDisplay;

private:

	//the computed graphical representation of a physical quantity (e.g. representation of magnetization, temperature etc.)
	PhysQRep physQRep;

	//Text Format for drawing text in mesh window
	IDWriteTextFormat *pMeshWinTextFormat = nullptr;

	//when mouse hovers over mesh values are picked and set here to be displayed
	DBL3 meshPosition;
	string meshName, typeName, mouse_mesh_info_position, mouse_mesh_info_value;
	bool displayMeshInfo = false;

private:

	void ZoomNotchUp(void);
	void ZoomNotchDn(void);

protected:

	void AdjustWindowDimensions(D2D1_RECT_F delta) { ChangeWindowSides(delta); }

	//compute the physical quantity representation for the given physical quantity
	void UpdatePhysQRep(vector<PhysQ> physQ);

	//adjust mesh display camera fov and averaging stencil dimension for default physical quantity representation
	void AutoSetMeshDisplaySettings(vector<PhysQ> physQ);

	double Get_PhysQRep_DetailLevel(void) { return physQRep.get_detail_level(); }

	void DrawWindow(void);
	void DrawWindow_Quick(void) { DrawWindow(); }

	ActionOutcome NewMessage(AC_ aCode, INT2 mouse, string data = "");

	void SetDefaultCursor(void) { SetCursor(LoadCursor(nullptr, IDC_ARROW)); }

	//image_cropping specify normalized cropping within the mesh image, as left, bottom, right, top : 0, 0 point is left, bottom of screen as far as user is concerned.
	bool SaveMeshImage(string fileName, DBL4 image_cropping);

public:

	//Setup new window space with rectangle specified by ratios of the entire screen
	BorisMeshWindow(D2D1_RECT_F ratios, INT2 winId);
	~BorisMeshWindow();

	void RepairObjectState(void) {}
};

#endif
