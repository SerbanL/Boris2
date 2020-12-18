#pragma once

//Contains BorisDisplay when running in non-graphics mode (basic text console only)

#include "CompileFlags.h"
#if GRAPHICS == 0

#include <string>
#include <vector>
#include <iostream>

#include "TextFormatting.h"
#include "BorisLib.h"

#include <SFML/System.hpp>
#include <SFML/Graphics/Texture.hpp>

#pragma comment(lib, "sfml-system.lib")
#pragma comment(lib, "sfml-graphics.lib")



////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Basic functionality in console-only mode. Portable. Still need to shape meshes using mask files in console-only mode.

class BorisGraphics :
	public ProgramState<BorisGraphics,
	std::tuple<float, float, float, float, float, float, float, float, float, float, float, int, int>, std::tuple<>>
{
private:

	float camX, camY, camZ, camDistance, fovDeg;
	float camUx, camUy, camUz;
	float view_shiftX, view_shiftY, view_shiftZ;
	int original_wndWidth, original_wndHeight;

public:

	BorisGraphics(void) :
		ProgramStateNames(this,
			{
				VINFO(camX), VINFO(camY), VINFO(camZ), VINFO(camDistance), VINFO(fovDeg),
				VINFO(camUx), VINFO(camUy), VINFO(camUz), VINFO(view_shiftX), VINFO(view_shiftY), VINFO(view_shiftZ),
				VINFO(original_wndWidth), VINFO(original_wndHeight)
			}, {})
	{}

	//from image file extract a raw bitmap (BYTE array with 4 BYTE-sized entries as B-G-R-A for each pixel) to specified pixel size (so image file is rescaled to specified size)
	void GetBitmapFromImage(std::string fileName, std::vector<unsigned char>& bitmap, INT2 n_plane);

	bool MakeVideoFromFileSequence(std::string directory, std::vector<std::string>& fileNames, unsigned int fps, double scaling, int quality)
	{
		std::cout << "Not available in console-only mode" << std::endl;

		return false;
	}

	void RepairObjectState(void) {}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 

struct PhysQRepComponent :
	public ProgramState<PhysQRepComponent,
	std::tuple<bool, double, int, bool, std::string, std::string, std::string>, std::tuple<>>
{
private:

	bool scale_to_magnitude = true;
	double exponent = 0.2;
	int obSelector = 0;
	bool drawFrame = true;
	std::string unit;
	std::string typeName;
	std::string meshName;

public:

	PhysQRepComponent(void) :
		ProgramStateNames(this,
			{
				VINFO(scale_to_magnitude), VINFO(exponent), VINFO(obSelector), VINFO(drawFrame), VINFO(unit), VINFO(typeName), VINFO(meshName)
			}, {})
	{}

	void RepairObjectState(void) {}
};

class PhysQRep :
	public ProgramState<PhysQRep,
	std::tuple<double, Rect, double, std::string, std::string, DBL2, std::string, std::vector<PhysQRepComponent>>, std::tuple<>>
{
private:

	double m_to_l = 1e9;
	Rect focusRect;
	double detail_level = 5e-9;
	std::string meshName_focused;
	std::string typeName_focused;
	DBL2 minmax_focusedmeshValues = DBL2();
	std::string unit;
	std::vector<PhysQRepComponent> physQRep;

public:

	PhysQRep(void) :
		ProgramStateNames(this, { VINFO(m_to_l), VINFO(focusRect), VINFO(detail_level), VINFO(meshName_focused), VINFO(typeName_focused), VINFO(minmax_focusedmeshValues), VINFO(unit), VINFO(physQRep) }, {})
	{}

	void RepairObjectState(void) {}
};

class BorisMeshWindow :
	public ProgramState<BorisMeshWindow,
	std::tuple<PhysQRep, DBL3, std::string, std::string, std::string, std::string, bool>, std::tuple<>>

{
private:

	PhysQRep physQRep;

	//when mouse hovers over mesh values are picked and set here to be displayed
	DBL3 meshPosition;
	std::string meshName, typeName, mouse_mesh_info_position, mouse_mesh_info_value;
	bool displayMeshInfo = false;

public:

	BorisMeshWindow(void) :
		ProgramStateNames(this, { VINFO(physQRep), VINFO(meshPosition), VINFO(meshName), VINFO(typeName), VINFO(mouse_mesh_info_position), VINFO(mouse_mesh_info_value), VINFO(displayMeshInfo) }, {})
	{}

	void RepairObjectState(void) {}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Boris Display container and director object. All access to display elements and drawing done through here.

class BorisDisplay :
	public ProgramState<BorisDisplay,
	std::tuple<BorisGraphics*, BorisMeshWindow*>, std::tuple<>>
{
private: //Private data

	BorisGraphics *pBG = nullptr;
	BorisMeshWindow *pbMeshWin = nullptr;

public:  //Public data

private: //Private methods	

public:

	BorisDisplay() :
		ProgramStateNames(this, { VINFO_KEEPPTR(pBG), VINFO_KEEPPTR(pbMeshWin) }, {})
	{
		pBG = new BorisGraphics();
		pbMeshWin = new BorisMeshWindow();
	}

	~BorisDisplay() 
	{
		delete pBG;
		delete pbMeshWin;
	}

	void RepairObjectState(void) {}
	
	//-----------------------------------------------Thread-safe access points to BorisDisplay

	void ClearDataBox(void) {}

	//display various message types in the console
	void DisplayConsoleMessage(std::string text) { std::cout << text << std::endl; }
	void DisplayConsoleError(std::string text) { std::cout << text << std::endl; }
	void DisplayConsoleWarning(std::string text) { std::cout << text << std::endl; }
	void DisplayConsoleListing(std::string text) { std::cout << text << std::endl; }

	//display message after stripping out the formatting
	void DisplayFormattedConsoleMessage(std::string text);

	//add a new entry in the data box as an interative object (use formatted text)
	void NewDataBoxField(std::string formattedText) {}
	//update the data box at given lineIdx to the given value, passed as a std::string
	void UpdateDataBoxField(int lineIdx, std::string value_string) {}

	//set the console line entry text
	void SetConsoleEntryLineText(std::string text) { std::cout << text; }

	bool SaveMeshImage(std::string fileName, DBL4 image_cropping = DBL4(0, 0, 1, 1))
	{ 
		std::cout << "Not available in console-only mode" << std::endl; 
		return false;
	}

	//Set mesh display detail level directly (N/A in non-graphics mode)
	void SetDetailLevel(double detail_level) {}
	double GetDetailLevel(void) { return 0.0; }

	//Set mesh display render threshold values for faster rendering when we have many cells
	void SetRenderThresholds(INT3 renderthresholds) {}
	INT3 GetRenderThresholds(void) { return INT3(); }

	//allows access to Boris Graphics public methods.
	BorisGraphics* BGMethods(void) { return pBG; }

	void WaitForDisplayEnd(void) {}
};

#endif