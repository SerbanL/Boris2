//Defines abstract Boris window spaces.
//Implementations and top display coordinator object should go in a separate file

#pragma once

#include "TextObject.h"
#include "BorisGraphics.h"
#include "BorisTypes.h"

//default text color options
//See: https://msdn.microsoft.com/en-us/library/windows/desktop/dd370907(v=vs.85).aspx

#define MESSAGECOLOR D2D1::ColorF(D2D1::ColorF::Green)  //rgba = (0, 0.5, 0, 1)
#define ERRORCOLOR D2D1::ColorF(D2D1::ColorF::Red)		//rgba = (1, 0, 0, 1)
#define LISTINGCOLOR D2D1::ColorF(D2D1::ColorF::White)  //rgba = (1, 1, 1, 1)
#define USERCOLOR D2D1::ColorF(D2D1::ColorF::Yellow)	//rgba = (1, 1, 0, 1)
#define BGRNDTEXTCOLOR D2D1::ColorF(0,0)				//rgba = (0, 0, 0, 0)
#define PROMPTCOLOR D2D1::ColorF(13882323, 0.7)			//Light Gray: rgba = (211/255, 211/255, 211/255, 1). Thus single rgba value as UINT32 : 13882323 = 256^2 * 211 + 256 * 211 + 211
#define POOPUPCOLOR D2D1::ColorF(11119017, 0.7)			//Gray: rgba = (255/255, 255/255, 255/255, 1). Thus single rgba value as UINT32 : 13882323 = 256^2 * 255 + 256 * 255 + 255
#define RESIZEFRAMECOLOR D2D1::ColorF(D2D1::ColorF::Orange)

#define DOUBLECLICKMS	160								//maximum time in ms between clicks to be considered a double click

#define RESIZEBORDER	10								//Number of pixels of resizing border
#define RESIZEFRAMEBORDER	3							//Number of pixels of displayed resizing frame

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Main window action codes (mouse actions, keyboard actions etc). Many of these are directly translated from WM_ codes in WndProc, some are not, so it's easier to work with these.
enum AC_ {AC_NONE = 0,
	      AC_MOUSELEFTDOWN, AC_MOUSERIGHTDOWN, AC_MOUSERIGHTUP, AC_MOUSELEFTUP, AC_MOUSEMIDDLEDOWN, AC_MOUSEMIDDLEUP, AC_CLICK, AC_DOUBLECLICK, AC_MOUSEWHEEL, 
		  AC_MOUSEMOVE, 
		  AC_DROPFILES, 
		  AC_KEYBOARD, 
		  AC_KBDENTER, AC_KBDBACK, AC_KBDINS, AC_KBDDEL, AC_KBDLEFT, AC_KBDRIGHT, AC_KBDDN, AC_KBDUP, AC_KBDEND, AC_KBDHOME
		 };

////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Action outcome codes, usually generated after processing AC_ messages, but can also be used to pass messages around
enum AO_ {AO_NONE = 0, AO_NOTHANDLED, AO_TEXTRETURNED, AO_SETTOPMOST, AO_REFRESH, AO_MAKENEWWINDOW, AO_DESTROYWINDOW, AO_INTERACTSPACES, AO_WINRESIZED};

//struct for passing action outcome messages after being handled in NewMessage method
struct ActionOutcome {

	//Multiple action outcome codes could be generated, but this should always have size >= 1
	vector<AO_> aoCodes;

	//text generated as a result of AO_TEXTRETURNED
	string text;

	//mouse coordinates where action occured
	INT2 mouse;

	//the winId of the window space which generated the action outcome
	INT2 winId;

	ActionOutcome(void) {

		text = "";
		mouse = INT2(0,0);
		aoCodes.push_back(AO_NONE);
	}

	ActionOutcome(AO_ aoCode) {

		text = "";
		mouse = INT2(0,0);
		this->aoCodes.push_back(aoCode);
	}

	bool IsCodeSet(AO_ aoCode) { 
		
		if(aoCodes.size() == 1) 
			return (aoCode == aoCodes[0]);
			
		for(int i = 0; i < aoCodes.size(); i++) {
		
			if(aoCode == aoCodes[i]) return true;
		}

		return false;
	}

	void SetCode(INT2 winId, INT2 mouse, AO_ aoCode) { aoCodes.clear(); aoCodes.push_back(aoCode); this->winId = winId; this->mouse = mouse; }

	void SetCodes(INT2 winId, INT2 mouse, int numCodes, ...) {

		va_list codes;

		va_start(codes, numCodes);

		SetCode(winId, mouse, va_arg(codes, AO_));

		for(int i = 1; i < numCodes; i++)
			AddCode(winId, mouse, va_arg(codes, AO_));

		va_end(codes);

		this->winId = winId;
		this->mouse = mouse;
	}

	void AddCode(INT2 winId, INT2 mouse, AO_ aoCode) { 
		
		if(aoCodes[0] != AO_NONE && aoCodes[0] != AO_NOTHANDLED)
			aoCodes.push_back(aoCode);
		else aoCodes[0] = aoCode;
		
		this->winId = winId; 
		this->mouse = mouse;
	}
	
	AO_ GetCode(int idx) { return aoCodes[idx]; }

	size_t NumCodes(void) { return aoCodes.size(); }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

//Resizing types - these values are purposely set as they are (see WindowResizing)
enum RSZ_ { RSZ_NONE = 0, RSZ_N = 1, RSZ_E = 2, RSZ_NE = 3, RSZ_S = 4, RSZ_SE = 6, RSZ_W = 7, RSZ_NW = 8, RSZ_SW = 11 };

class WindowResizing {

private:

	int resizeType;

	//resize enable/disable flags for window sides
	bool allowLeft, allowRight, allowTop, allowBottom;

public:

	WindowResizing(void) { resizeType = RSZ_NONE; EnableResizing(false); }

	void EnableResizing(bool flag = true) { allowLeft = flag; allowRight = flag; allowTop = flag; allowBottom = flag; }
	void EnableResizing(bool left, bool right, bool top, bool bottom) { allowLeft = left; allowRight = right; allowTop = top; allowBottom = bottom; }

	bool IsResizable(void) { return (allowLeft | allowRight | allowTop | allowBottom); }

	bool IsResizing(void) { return (resizeType != RSZ_NONE); }

	void ResetResize(void) { resizeType = RSZ_NONE; }

	//Set resizing cursor and value for resizeType indicator if mouse entered resizing area. Return true if it is no longer in resizing area (but was before).
	bool IsMouseInResizingArea(INT2 mouse, D2D1_RECT_F spaceRect) {

		UniformResizeRect(spaceRect, RESIZEBORDER);

		if(!IsInside(spaceRect, mouse)) {

			ResetResize();

			//mouse is on resize border (assuming here mouse is inside the original spaceRect). Now find the type of resize.
			if(allowRight && mouse.i >= spaceRect.right) resizeType += RSZ_E;
			if(allowLeft && mouse.i <= spaceRect.left) resizeType += RSZ_W;
			if(allowTop && mouse.j <= spaceRect.top) resizeType += RSZ_N;
			if(allowBottom && mouse.j >= spaceRect.bottom) resizeType += RSZ_S;

			//Now set resizing cursor
			switch(resizeType) {

			case RSZ_S:
			case RSZ_N:
				SetCursor(LoadCursor(NULL, IDC_SIZENS));
				break;

			case RSZ_W:
			case RSZ_E:
				SetCursor(LoadCursor(NULL, IDC_SIZEWE));
				break;

			case RSZ_NE:
			case RSZ_SW:
				SetCursor(LoadCursor(NULL, IDC_SIZENESW));
				break;

			case RSZ_NW:
			case RSZ_SE:
				SetCursor(LoadCursor(NULL, IDC_SIZENWSE));
				break;

			default:
				break;
			}
		}
		else if(IsResizing()) { ResetResize(); return false; }

		return true;
	}

	D2D1_RECT_F GetResizingDeltas(INT2 mouseDelta) {
		
		D2D1_RECT_F delta = D2D1::RectF(0,0,0,0);

		if(resizeType == RSZ_NE || resizeType == RSZ_E || resizeType == RSZ_SE) delta.right += (float)mouseDelta.i;
		if(resizeType == RSZ_NW || resizeType == RSZ_W || resizeType == RSZ_SW) delta.left += (float)mouseDelta.i;
		if(resizeType == RSZ_NW || resizeType == RSZ_N || resizeType == RSZ_NE) delta.top += (float)mouseDelta.j;
		if(resizeType == RSZ_SW || resizeType == RSZ_S || resizeType == RSZ_SE) delta.bottom += (float)mouseDelta.j;

		return delta;
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct WinSpaceSettings {

	//the window space rectangle in pixels, together with dimensions (width and height - redundant info but useful)
	D2D1_RECT_F spaceRect;
	D2D1_SIZE_F spaceRectDims;

	//The window spaces are normally drawn in reverse order as they appear in bWin (index 0 is on top). Some special window spaces can be set to be drawn always on the bottom or on top.
	bool stickyDisplayBottom, stickyDisplayTop;

public:

	WinSpaceSettings(void) {};

	WinSpaceSettings(D2D1_RECT_F spaceRect, D2D1_SIZE_F spaceRectDims, bool stickyDisplayBottom, bool stickyDisplayTop) {

		this->spaceRect = spaceRect;
		this->spaceRectDims = spaceRectDims;

		this->stickyDisplayBottom = stickyDisplayBottom;
		this->stickyDisplayTop = stickyDisplayTop;
	}

	void RecallSettings(D2D1_RECT_F &spaceRect, D2D1_SIZE_F &spaceRectDims, bool &stickyDisplayBottom, bool &stickyDisplayTop) {

		spaceRect = this->spaceRect;
		spaceRectDims = this->spaceRectDims;

		stickyDisplayBottom = this->stickyDisplayBottom;
		stickyDisplayTop = this->stickyDisplayTop;
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Abstract class for Boris windows spaces. All the basic common functionality and properties of actual window spaces implementations must go here.

class BorisWinSpace : public GraphicalObject {

private: //Private data

protected: //Protected data

	//windows id number: major and minor. major will typically be the type of window (use an enum), minor will be incremented when multiple windows of same type are added.
	INT2 winId;

	//the window space rectangle in pixels, together with dimensions (width and height - redundant info but useful)
	D2D1_RECT_F spaceRect;
	D2D1_SIZE_F spaceRectDims;

	//background color
	D2D1_COLOR_F bgrndColor, resizingFrameColor;

	int dragStartX, dragStartY;
	bool mouseLeftDown, mouseMiddleDown, mouseRightDown;

	//flag to indicate if this window, as a top level window, has been entered (might need to change mouse cursor etc)
	bool windowEntered;

	//enable maximising or not?
	bool allowMaximising;

	WindowResizing resizing;

	//used to toggle between maximized and normal size, if enabled.
	bool windowMaximised;

	//The window spaces are normally drawn in reverse order as they appear in bWin (index 0 is on top). Some special window spaces can be set to be drawn always on the bottom or on top.
	bool stickyDisplayBottom, stickyDisplayTop;

	//a value storing system ticks in ms, used to time some actions (e.g. custom left double-click implementation - simply call IsDoubleClick() on AC_MOUSELEFTDOWN message where this is needed)
	DWORD dblClickTimeCounter;

	//used to store settings before a change, if we need to come back to them (e.g. maximize / restore type toggling).
	WinSpaceSettings saveWinSpaceSettings;

private: //Private methods

protected: //Protected methods

	//Windows event occured in this space: process it. This is shared by all derived classes and defines common responses, in addition to particular responses in NewMessage
	ActionOutcome NewMessage_CommonResponses(AC_ aCode, INT2 mouse, char param);

	bool IsDoubleClick(AC_ &aCode);

	//set window size using ratios of screen dimensions / change window size by adjusting sides using delta values
	void ResizeWindow(D2D1_RECT_F ratios);
	void ChangeWindowSides(D2D1_RECT_F delta);

	//maximise / recall window size toggling
	void MaximiseWindowSize(void);
	void RestoreWindowSize(void);

	//toggle between the MaximiseWindowSize and RestoreWindowSize depending on the windowMaximised flag setting
	void Toggle_Maximise_Normal(INT2 mouse, ActionOutcome &actionResult);

	void DrawResizingFrame(void);

public: //Public methods

	//Setup new window space with rectangle specified by ratios of the entire screen
	BorisWinSpace(D2D1_RECT_F ratios, INT2 winId);
	virtual ~BorisWinSpace() = 0;

	//Draw the window content. Do not call it directly, should only be called by Refresh method in BorisDisplay.
	virtual void DrawWindow(void) = 0;

	//Windows event occured in this space: process it.
	virtual ActionOutcome NewMessage(AC_ aCode, INT2 mouse, char param) = 0;

	virtual void SetDefaultCursor(void) = 0;

	//This should be implemented to call ChangeWindowSides method. ChangeWindowSides also appears in objects derived from this class and adds extra functionality on top. By implementing this, the right version of ChangeWindowSides will be called.
	virtual void AdjustWindowDimensions(D2D1_RECT_F delta) = 0;

	void WindowLeft(void);
	void WindowEntered(void);

	//set background color for this space
	void SetBackground(D2D1_COLOR_F bgrndColor) { this->bgrndColor = bgrndColor; }

	bool IsMouseInWindow(INT2 mouse) { return IsInside(spaceRect, mouse); }

	//handle various flags
	void SetStickyDisplayBottom(bool state = true) { stickyDisplayBottom = state; if(stickyDisplayBottom) stickyDisplayTop = false; }
	void SetStickyDisplayTop(bool state = true) { stickyDisplayTop = state; if(stickyDisplayTop) stickyDisplayBottom = false; }
	bool IsStickyDisplayBottom(void) { return stickyDisplayBottom; }
	bool IsStickyDisplayTop(void) { return stickyDisplayTop; }

	void IsMaximisable(bool flag) { allowMaximising = flag; }
	void EnableResizing(bool flag) { resizing.EnableResizing(flag); }
	void EnableResizing(bool left, bool right, bool top, bool bottom) { resizing.EnableResizing(left, right, top, bottom); }

	void GetSpaceTopleftCoordinates(float &X, float &Y) { X = spaceRect.left; Y = spaceRect.top; }
	D2D1_RECT_F GetSpaceRect(void) { return spaceRect; }

	INT2 GetwinId(void) { return winId; }
};

inline BorisWinSpace::~BorisWinSpace() {}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Derived from the abstract BorisWinSpace, but does not implement it. Defines basic functionality for displaying text.
//Derive window spaces implementations that need this functionality from this.

class BorisTextDisplay : public BorisWinSpace {

private: //private data

	vector<string> formattingVector;

protected: //protected data

	//Text is displayed using TextObject objects. Each line can have a number of these (thus a vector of vectors). Each principal entry corresponds to a text line exactly, and is guaranteed to fit in a window line when displayed.
	//this vector should not be empty
	TextLines textLines;

	//the top line index to display (starting at the top of the window space)
	//the coordinates of stored text objects are not directly screen coordinates, but these can be calculated using the space rect and topLine index.
	//Text object coordinates are relative to the space rect, and each new line is added with Y coordinate incrementing the last stored line.
	int topLine;

	//Action handler functionoid: pass this to interactive text objects. When the text objects are interacted with, they will use this functionoid to call the handler routine, using their id and message parameter
	SimTOFunct *pActionHandler;

	//allow text lines to wrap when exceeding width of window?
	bool allowTextLineWrapping;

private: //private methods

	//Indicate if given line index will overflow the space rect for the current topLine index.
	//Also calculate the required Y screen offset for drawing text objects for the current topLine index (add to the text object coordinate to obtain actual screen coordinate)
	bool LineOverflowsSpace(int lineIdx, float *pscreenYOffset = NULL);

protected: //protected methods

	ActionOutcome NewMessage_CommonResponses(AC_ aCode, INT2 mouse, char param);

	//Set topLine so that the last line will not overflow the space rect. If it does, set topLine so that the last line is as far down in the console as possible
	void ResetTopLineIndex(void);

	//draw the textLines content. return true if lines had to be recalculated to fit the window width
	bool DrawTextLines(bool doDraw = true);

	//set window size using ratios of screen dimensions / change window size by adjusting sides using delta values
	void ResizeWindow(D2D1_RECT_F ratios);
	void ChangeWindowSides(D2D1_RECT_F delta);

	//maximise / recall window size toggling : overloads the base method, calls it then adds extra settings for text
	void MaximiseWindowSize(void);
	void RestoreWindowSize(void);

	//toggle between the MaximiseWindowSize and RestoreWindowSize depending on the windowMaximised flag setting
	void Toggle_Maximise_Normal(INT2 mouse, ActionOutcome &actionResult);

public: //protected methods

	BorisTextDisplay(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pActionHandler);
	virtual ~BorisTextDisplay() = 0;

	void SetTextLineWrapping(bool flag) { allowTextLineWrapping = flag; if(allowTextLineWrapping) textLines.SetWidthLimit(spaceRectDims.width); else textLines.SetWidthLimit(0); }

	int LastLine(void) { return textLines.LastLine(); }

	//get text line from lineIdx with coordinates as it appears on screen
	TextLine GetRawTextLine_ScreenCoords(int lineIdx);

	//get text line index from given mouse coordinates
	int GetMouseClickLineIdx(INT2 mouse);

	//Set text in window space, erasing any existing text
	void SetText(TextLine textLine) { textLines.set(textLine); }

	//Set text in window space, erasing any existing text
	void SetText(string text, FormatSpecifier fs);

	//Add text using text formatting (either on a new line or at the given line index as an insertion, pushing down all existing lines from that index onwards)
	//<b> </b> : bold
	//<i> </i> : italic
	//[on], [os], [or] : text outline (TOO_NONE, TOO_SQUARE, TOO_ROUND)
	//[tcr,g,b,a/tc] : text color specified using r,g,b,a
	//[bcr,g,b,a/bc] : background color
	//[io0,0,textId/io] : interactive object settings: first two numbers identify the major and minor TextObject type (used to identify the type of action generated when interacting with this TextObject). Last field contains a text identifier.
	//</io> use this to specify interactive object end
	void NewFormattedTextLine(int lineIdx, string text);
};