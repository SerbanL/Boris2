//Defines graphics methods built on D3D

#pragma once

#include "CompileFlags.h"
#if OPERATING_SYSTEM == OS_WIN

#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <d2d1.h>
#include <wincodec.h>
#include <dwrite.h>
#include <D3d9types.h>

#include "BorisLib.h"

#include "D3D.h"
#include "Boris_Enums_Defs.h"

#pragma comment(lib, "WindowsCodecs.lib")
#pragma comment (lib, "d2d1.lib")
#pragma comment (lib, "dwrite.lib")

//Constants using in drawing the coordinate system
#define COORDSYSTLENGTHPIXELS	80
#define COORDSYSTCAMDIST	38

//3D view options
#define ROTATIONSENS	10.0F
#define CAMERAZSENS	1000.0F
#define CAMERAROTSENS	5.0F
#define MINFOV	5.0F
#define MAXFOV	170.0F

//scaling factor for fov change using the mouse wheel
#define ZOOMNOTCH 1.15F



//Text object outline enum: types of outline to be used for text objects
enum TOO_ {TOO_NONE = 0, TOO_ROUND, TOO_SQUARE};

//cell display object - objects in objCol (ObjectBufferCollection). These must be loaded in objCol in the same order as they appear in the enum.
enum CDO_ {CDO_NONE = -1, CDO_CONE = 0, CDO_ARROW, CDO_HALFARROW, CDO_CUBE, CDO_HALFCUBE, CDO_CUBEFRAME, CDO_DODECA, CDO_AXISX, CDO_AXISY, CDO_AXISZ};

////////////////////////////////////////////////////////////////////////////////////////
//Text formatting specifier

struct FormatSpecifier {

	D2D1_COLOR_F textColor, bgrndColor;
	
	bool italic, bold;
	
	TOO_ textOutline;

	DWRITE_TEXT_ALIGNMENT textAlignment;
	DWRITE_PARAGRAPH_ALIGNMENT paragraphAlignment;

	FormatSpecifier(void) {

		textColor = D2D1::ColorF(D2D1::ColorF::White);
		bgrndColor = D2D1::ColorF(0,0);

		italic = false;
		bold = false;

		textOutline = TOO_NONE;

		textAlignment = DWRITE_TEXT_ALIGNMENT_LEADING;
		paragraphAlignment = DWRITE_PARAGRAPH_ALIGNMENT_NEAR;
	}

	FormatSpecifier(D2D1_COLOR_F textColor) {

		this->textColor = textColor;
		bgrndColor = D2D1::ColorF(0,0);

		italic = false;
		bold = false;

		textOutline = TOO_NONE;

		textAlignment = DWRITE_TEXT_ALIGNMENT_LEADING;
		paragraphAlignment = DWRITE_PARAGRAPH_ALIGNMENT_NEAR;
	}

	FormatSpecifier(D2D1_COLOR_F textColor, D2D1_COLOR_F bgrndColor, bool italic, bool bold, TOO_ textOutline) {

		this->textColor = textColor;
		this->bgrndColor = bgrndColor;
		
		this->italic = italic;
		this->bold = bold;
		
		this->textOutline = textOutline;

		textAlignment = DWRITE_TEXT_ALIGNMENT_LEADING;
		paragraphAlignment = DWRITE_PARAGRAPH_ALIGNMENT_NEAR;
	}

	bool operator==(const FormatSpecifier &rhs) {

		 return ( SameColor(textColor, rhs.textColor) && SameColor(bgrndColor, rhs.bgrndColor) &&
				  italic == rhs.italic && bold == rhs.bold && textOutline == rhs.textOutline &&
				  textAlignment == rhs.textAlignment && paragraphAlignment == rhs.paragraphAlignment);
	}
};

////////////////////////////////////////////////////////////////////////////////////////

class BorisGraphics : 
	public D3D,
	public ProgramState<BorisGraphics,
	std::tuple<float, float, float, float, float, float, float, float, float, float, float, int, int>, std::tuple<>>
{
private:

	ID2D1SolidColorBrush *pVariableColorBrush;

	//Text settings. textMetrics.height must always contain the height of a single line of text.
	IDWriteTextFormat *pNormalTextFormat, *pBoldTextFormat, *pItalicTextFormat, *pBoldItalicTextFormat;

	//default font and font sizes. set at initialization, cannot be changed later
	const std::string defFont;
	const float defFontSize;

	float fontPixHeight;

	//for monospaced fonts the width is fixed - calculated once when font is set
	float monospacedfontPixWidth;

	//collection of objects used to draw geometrical shapes
	ObjectBufferCollection objCol;

	//"original" main window width and height : these values are saved in simulation files.
	//When loading a simulation file (in particular on a screen with different resolution), then you can compare these values with the current wndWidth and wndHeight values.
	//Any quantities which have been saved as pixels can then be adjusted to be displayed correctly at the loaded resolution
	int original_wndWidth, original_wndHeight;

private:

	//setting up methods
	void Setup2DGraphicsResources(void);
	HRESULT Setup3DGraphicsResources(void);

	void Initialize(void);
	void Cleanup(void);

	//Load into objCol the geometrical object specified in the given VIN file, and using the given pixel and vertex compiled shaders (give filenames)
	HRESULT LoadVINFile(std::string filename, std::string VSFile, std::string PSFile, ObjectBufferCollection &objCol);

public:

	BorisGraphics(HWND hWnd, std::string defFont, float defFontSize);
	~BorisGraphics();

	void RepairObjectState(void);

	//---------------------------------------------------

	//Brush and text format creation methods
	HRESULT CreateSolidBrush(D2D1_COLOR_F color, ID2D1SolidColorBrush **ppBrush);
	HRESULT CreateTextFormat(std::string fontName, FLOAT fontSize, DWRITE_FONT_WEIGHT fontWeight, DWRITE_FONT_STYLE fontStyle, IDWriteTextFormat **ppTextFormat);

	//Create geometry from given path points
	HRESULT CreatePath(D2D1_POINT_2F *ppathpoints, UINT32 pointsCount, ID2D1PathGeometry **ppPathGeometry);

	void CalculateFontPixelsHeight(void);
	void CalculateMonospacedFontPixelsWidth(void);

	//---------------------------------------------------

	//Special camera conntrol routines
	void RotateCameraAboutOrigin(float dAzim, float dPolar);
	void AdjustCameraDistanceFromOrigin(float dZ);
	void RotateCameraView(float dAngle);
	
	//Calculate FOV required to fit the given number of logical units along the viewing width. viewDims is the view port size in pixels.
	float SetCameraFOVToFitUnits(float units, D2D1_SIZE_F viewDims);
	
	//reset the rendering transforms
	void ResetTransformation(void);

	void FillRectangle(D2D1_RECT_F rect, D2D1_COLOR_F color);
	void FillPromptRectangle(D2D1_RECT_F rect);

	//Draw text methods using full formatting specifier
	void DrawTextLine(std::string text, D2D1_RECT_F textRect, FormatSpecifier fs);
	//Simpler text draw method with externally provided IDWriteTextFormat, but fewer options : just the text, top-left position and color
	void DrawSimpleTextLine(std::string text, FLT2 position, D2D1_COLOR_F textColor, IDWriteTextFormat** ppTextFormat);

	//Draw stored 3D object
	void DrawCBObject(CDO_ objColSelector, XMMATRIX Rotation, XMMATRIX Scaling, XMMATRIX Translation, XMFLOAT4 Color) { D3D::DrawCBObject(Rotation, Scaling, Translation, Color, objCol, objColSelector); }
	void DrawCBObject(CDO_ objColSelector, XMFLOAT4 Color) { D3D::DrawCBObject(XMMatrixIdentity(), XMMatrixIdentity(), XMMatrixIdentity(), Color, objCol, objColSelector); }

	//draw a batch of objects with precomputed transforms (CBObjectTransform). Draw using selected object from ObjectBufferCollection objCol
	void DrawCBObjectBatch(std::vector<CBObjectTransform>& transformBatch, CDO_ objColSelector);

	void DrawFrameObject(CDO_ objColSelector, XMMATRIX Rotation, XMMATRIX Scaling, XMMATRIX Translation, XMFLOAT4 Color) { D3D::DrawFrameObject(Rotation, Scaling, Translation, Color, objCol, objColSelector); }
	void DrawFrameObject(CDO_ objColSelector, XMFLOAT4 Color) { D3D::DrawFrameObject(XMMatrixIdentity(), XMMatrixIdentity(), XMMatrixIdentity(), Color, objCol, objColSelector); }

	//Draw the coordinate system in the bottom-left corner of the given spaceRect
	void DrawCoordinateSystem(D2D1_RECT_F spaceRect);

	//Get dimensions in physical pixels of formatted text std::string
	float GetFontStringPixelsWidth(const std::string& str, const FormatSpecifier& fs);
	float GetMonospacedFontStringPixelsWidth(const std::string& str);

	float GetFontPixelsHeight(void) { return fontPixHeight; }
	float GetMonospacedFontPixelsWidth(void) { return monospacedfontPixWidth; }

	//save currently displayed image (in specified rectangle) to file
	bool SaveScreenToFile(std::string fileName, D2D1_RECT_F capture_rect);

	//Get camera properties
	FLT3 GetCameraPosition(void) { return FLT3(camX, camY, camZ); }
	FLT3 GetCameraUp(void) { return FLT3(camUx, camUy, camUz); }
	FLT3 GetViewShift(void) { return FLT3(view_shiftX, view_shiftY, view_shiftZ); }
	float GetFOV(void) { return fovDeg; }
	float GetCameraDistance(void) { return camDistance; }

};

//Wrapper for BorisGraphics static pointer (only one needed, hence static). If an object needs access to graphics, derive its class from GraphicalObject. Only the top container object needs to call the full constructor, the rest can just call the empty one.
class GraphicalObject {

private:

	bool thisObjectHasAllocatedMemory;

protected:

	static BorisGraphics *pBG;

public:

	//default "empty" constructor
	GraphicalObject() { thisObjectHasAllocatedMemory = false; }

	//full constructor - only the top level graphical object (e.g. a display coordinator object) should call this, the rest can just call the empty constructor.
	GraphicalObject(HWND hWnd, std::string defFont, float defFontSize) {

		if(!pBG) {

			pBG = new BorisGraphics(hWnd, defFont, defFontSize);
			thisObjectHasAllocatedMemory = true;
		}
		else thisObjectHasAllocatedMemory = false;
	}

	virtual ~GraphicalObject() { 
		
		if(thisObjectHasAllocatedMemory) {

			if(pBG) delete pBG;
			pBG = nullptr;
		}
	}
};

#endif