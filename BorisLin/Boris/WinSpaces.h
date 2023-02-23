//Defines abstract Boris window spaces.
//Implementations and top display coordinator object should go in a separate file

#pragma once

#include "CompileFlags.h"
#if GRAPHICS == 1

#include "TextObject.h"
#include "BorisGraphics.h"
#include "BorisLib.h"

#include "ColorDefs.h"

#define DOUBLECLICKMS	250								//maximum time in ms between clicks to be considered a double click
#define HOVERTIME		300								//minimum time in ms to wait before showing something when mouse is hovering
#define HOVERTIME_LONG	1000							//minimum time in ms to wait before showing something when mouse is hovering

#define RESIZEBORDER	10								//Number of pixels of resizing border
#define RESIZEFRAMEBORDER	3							//Number of pixels of displayed resizing frame



////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Main window action codes (mouse actions, keyboard actions etc). Many of these are directly translated from WM_ codes in WndProc, some are not, so it's easier to work with these.
enum AC_ 
{
	AC_NONE = 0, 
	AC_ALLWINDOWSLEFT,
	AC_MOUSELEFTDOWN, AC_SHIFT_MOUSELEFTDOWN, AC_MOUSERIGHTDOWN, AC_MOUSERIGHTUP, AC_MOUSELEFTUP, AC_MOUSEMIDDLEDOWN, AC_MOUSEMIDDLEUP, AC_CLICK, AC_DOUBLECLICK, AC_MOUSEWHEEL,
	AC_MOUSEMOVE, AC_HOVERCHECK,
	AC_DROPFILES, 
	AC_KEYBOARD, 
	AC_KBDENTER, AC_KBDBACK, AC_KBDINS, AC_KBDDEL, AC_KBDLEFT, AC_KBDSHIFTLEFT, AC_KBDRIGHT, AC_KBDSHIFTRIGHT, AC_KBDDN, AC_KBDUP, AC_KBDEND, AC_KBDSHIFTEND, AC_KBDHOME, AC_KBDSHIFTHOME, AC_KBDPGDN, AC_KBDPGUP, AC_KBDESC,
	AC_CTRL_C, AC_CTRL_V,
	AC_INTERACTOBJECTS, AC_INTERACTOBJECTWITHWINDOW, AC_DROPINTERACTOBJECTS, AC_DROPINTERACTOBJECTWITHWINDOW,
	AC_CONSOLECOMMAND, AC_CONSOLEENTRY, AC_CONSOLECOMMAND_ENTRY, AC_CONSOLECOMMAND_NOPARAMS_ENTRY,
	AC_POPUPEDITTEXTBOXRETURNEDTEXT
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Action outcome codes, usually generated after processing AC_ messages, but can also be used to pass messages around
enum AO_ {
	AO_NOTHING = 0, AO_NOTHANDLED, 
	AO_TEXTRETURNED, AO_MESSAGERETURNED, AO_MESHFOCUS2, AO_FILEDROPPEDINCONSOLE, AO_FILEDROPPEDINMESH,
	AO_SETTOPMOST, 
	AO_REFRESH, AO_REFRESHWINDOW, AO_DELAYEDREFRESH, AO_DRAW, AO_DRAWWINDOW,
	AO_STARTINTERACTION, AO_ENDINTERACTION, AO_HOVERCHECK, AO_SHOWHOVERINFO, AO_DESTROYWINDOW, AO_CHECKIOINTERACTION, AO_CHECKIODROPINTERACTION, AO_WINRESIZED,
	AO_STARTPOPUPEDITBOX, AO_POPUPEDITBOXRETURNEDTEXT,
	AO_RECALCULATEMESHDISPLAY, 
	AO_ADDCONSOLECOORDINATE
};

//ActionOutcome object for passing action outcome messages after being handled in NewMessage method
class ActionOutcome {

	//Multiple action outcome codes could be generated, but this should always have size >= 1
	std::vector<AO_> aoCodes;

public:

	//text generated as a result of AO_TEXTRETURNED or AO_MESSAGERETURNED
	std::string text;

	//mouse coordinates where action occured
	INT2 mouse;

	//the winId of the window space which generated the action outcome
	INT2 winId;

private:

	//helpers for SetCodes with parameter pack
	void AddCodes(AO_ code) { aoCodes.push_back(code); }
	template <typename... Type> void AddCodes(AO_ aoCode, Type... aoCodes) {

		AddCodes(aoCode);
		AddCodes((AO_)adCodes...);
	}

public:

	ActionOutcome(void) {

		text = "";
		mouse = INT2(0,0);
		aoCodes.push_back(AO_NOTHING);
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

	//set a single code
	void SetCode(INT2 winId, INT2 mouse, AO_ aoCode) { aoCodes.clear(); aoCodes.push_back(aoCode); this->winId = winId; this->mouse = mouse; }

	//set multiple codes using a parameter pack
	template <typename... Type> void SetCodes(INT2 winId, INT2 mouse, AO_ aoCode, Type... aoCodes) {

		SetCode(winId, mouse, aoCode);
		AddCodes((AO_)aoCodes...);
	}

	//Add a code: aoCodes always has size >= 1. If no code set before aoCodes[0] is AO_NOTHING.
	void AddCode(INT2 winId, INT2 mouse, AO_ aoCode) { 
		
		if(aoCodes[0] != AO_NOTHING && aoCodes[0] != AO_NOTHANDLED)
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
enum RSZ_ { RSZ_NONE = 0, RSZ_N = 1, RSZ_E = 2, RSZ_NE = 3, RSZ_S = 4, RSZ_SE = 6, RSZ_W = 7, RSZ_NW = 8, RSZ_SW = 11, RSZ_JUSTENTERED, RSZ_JUSTLEFT };

class WinSpace;
typedef void (WinSpace::*CursorMethod)(void);

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

	bool IsResizeCursorSet(void) { return (resizeType != RSZ_NONE); }

	void ResetResize(void) { resizeType = RSZ_NONE; }

	//Set resizing cursor and value for resizeType indicator if mouse entered resizing area. Return RSZ_JUSTENTERED if it was not a resize cursor before but now is, RSZ_JUSTLEFT for the opposite case.
	int MouseNewEntryInResizingArea(INT2 mouse, D2D1_RECT_F spaceRect) {

		if(!IsResizable()) return RSZ_NONE;

		UniformResizeRect(spaceRect, RESIZEBORDER);

		if(!IsInside(spaceRect, mouse)) {

			bool isnewResize = !IsResizeCursorSet();
			int newresizeType = RSZ_NONE;

			//mouse is on resize border (assuming here mouse is inside the original spaceRect). Now find the type of resize.
			if(allowRight && mouse.i >= spaceRect.right) newresizeType += RSZ_E;
			if(allowLeft && mouse.i <= spaceRect.left) newresizeType += RSZ_W;
			if(allowTop && mouse.j <= spaceRect.top) newresizeType += RSZ_N;
			if(allowBottom && mouse.j >= spaceRect.bottom) newresizeType += RSZ_S;

			if(resizeType != newresizeType) {

				//change resize cursor type
				resizeType = newresizeType;

				//Now set resizing cursor
				switch(resizeType) {

				case RSZ_S:
				case RSZ_N:
					SetCursor(LoadCursor(nullptr, IDC_SIZENS));
					break;

				case RSZ_W:
				case RSZ_E:
					SetCursor(LoadCursor(nullptr, IDC_SIZEWE));
					break;

				case RSZ_NE:
				case RSZ_SW:
					SetCursor(LoadCursor(nullptr, IDC_SIZENESW));
					break;

				case RSZ_NW:
				case RSZ_SE:
					SetCursor(LoadCursor(nullptr, IDC_SIZENWSE));
					break;

				default:
					break;
				}
			}

			if(isnewResize) return RSZ_JUSTENTERED;
			else return resizeType;
		}
		else if(IsResizeCursorSet()) { 
			
			ResetResize(); 
			return RSZ_JUSTLEFT;
		}

		return RSZ_NONE;
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

class WinSpace : public GraphicalObject {

private: //Private data

protected: //Protected data

	//windows id number: major and minor. major will typically be the type of window (use an enum), minor will be incremented when multiple windows of same type are added.
	INT2 winId;

	//the window space rectangle in pixels, together with dimensions (width and height - redundant info but useful)
	D2D1_RECT_F spaceRect;
	D2D1_SIZE_F spaceRectDims;

	//background color
	D2D1_COLOR_F bgrndColor = D2D1::ColorF(D2D1::ColorF::White, 1.0f);
	D2D1_COLOR_F resizingFrameColor = RESIZEFRAMECOLOR;

	INT2 dragStart = INT2(0);
	bool mouseLeftDown = false, mouseMiddleDown = false, mouseRightDown = false;

	//flag to indicate if this window, as a top level window, has been entered (might need to change mouse cursor etc)
	bool windowEntered = false;

	//enable maximising or not?
	bool allowMaximising = false;

	WindowResizing resizing;

	//used to toggle between maximized and normal size, if enabled.
	bool windowMaximised = false;

	//The window spaces are normally drawn in reverse order as they appear in bWin (index 0 is on top). Some special window spaces can be set to be drawn always on the bottom or on top.
	bool stickyDisplayBottom = false, stickyDisplayTop = false;

	//a value storing system ticks in ms, used to time some actions (e.g. custom left double-click implementation - simply call IsDoubleClick() on AC_MOUSELEFTDOWN message where this is needed)
	DWORD dblClickTimeCounter, resizeHoverTimeCounter, hoverTimeCounter;

	//used to store settings before a change, if we need to come back to them (e.g. maximize / restore type toggling).
	WinSpaceSettings saveWinSpaceSettings;

private: //Private methods

protected: //Protected methods

	//Windows event occured in this space: process it. This is shared by all derived classes and defines common responses, in addition to particular responses in NewMessage
	ActionOutcome NewMessage_CommonResponses(AC_ aCode, INT2 mouse, std::string data);

	bool IsDoubleClick(AC_ &aCode);
	bool IsHoveringOverResizeArea(void);

	//set window size using ratios of screen dimensions / change window size by adjusting sides using delta values
	void ResizeWindow(D2D1_RECT_F ratios);
	void ChangeWindowSides(D2D1_RECT_F delta);

	//maximise / recall window size toggling
	void MaximiseWindowSize(void);
	void RestoreWindowSize(void);

	//toggle between the MaximiseWindowSize and RestoreWindowSize depending on the windowMaximised flag setting
	void Toggle_Maximise_Normal(INT2 mouse, ActionOutcome &actionResult);

	void DrawFrame(void);

	void DrawResizingFrame(void);

	//set background color for this space
	void SetBackground(D2D1_COLOR_F bgrndColor) { this->bgrndColor = bgrndColor; }

	void IsMaximisable(bool flag) { allowMaximising = flag; }
	void EnableResizing(bool flag) { resizing.EnableResizing(flag); }
	void EnableResizing(bool left, bool right, bool top, bool bottom) { resizing.EnableResizing(left, right, top, bottom); }

	void GetSpaceTopleftCoordinates(float &X, float &Y) { X = spaceRect.left; Y = spaceRect.top; }

public: //Public methods

	//Setup new window space with rectangle specified by ratios of the entire screen
	WinSpace(D2D1_RECT_F ratios, INT2 winId);
	virtual ~WinSpace() {}

	//Draw the window content. Do not call it directly, should only be called by Refresh method in BorisDisplay.
	virtual void DrawWindow(void) = 0;

	//Similar to DrawWindow but do not refresh interacitve objects
	virtual void DrawWindow_Quick(void) = 0;

	//Windows event occured in this space: process it.
	virtual ActionOutcome NewMessage(AC_ aCode, INT2 mouse, std::string data = "") = 0;

	virtual void SetDefaultCursor(void) = 0;

	//This should be implemented to call ChangeWindowSides method. ChangeWindowSides also appears in objects derived from this class and adds extra functionality on top. By implementing this, the right version of ChangeWindowSides will be called.
	virtual void AdjustWindowDimensions(D2D1_RECT_F delta) = 0;

	bool IsMouseInWindow(INT2 mouse) { return IsInside(spaceRect, mouse); }

	void WindowLeft(void);
	void WindowEntered(void);

	D2D1_RECT_F GetSpaceRect(void) { return spaceRect; }

	INT2 GetwinId(void) { return winId; }

	//handle various flags
	void SetStickyDisplayBottom(bool state = true) { stickyDisplayBottom = state; if(stickyDisplayBottom) stickyDisplayTop = false; }
	void SetStickyDisplayTop(bool state = true) { stickyDisplayTop = state; if(stickyDisplayTop) stickyDisplayBottom = false; }
	bool IsStickyDisplayBottom(void) { return stickyDisplayBottom; }
	bool IsStickyDisplayTop(void) { return stickyDisplayTop; }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Derived from the abstract WinSpace, but does not implement it. Defines basic functionality for displaying text.
//Derive window spaces implementations that need this functionality from this.

class TextDisplay : public WinSpace {

private: //private data

protected: //protected data

	//Text is displayed using TextObject objects. Each line can have a number of these (thus a vector of vectors). Each principal entry corresponds to a text line exactly, and is guaranteed to fit in a window line when displayed.
	//this vector should not be empty
	TextLines textLines;

	//the top line index to display (starting at the top of the window space)
	//the coordinates of stored text objects are not directly screen coordinates, but these can be calculated using the space rect and topLine index.
	//Text object coordinates are relative to the space rect, and each new line is added with Y coordinate incrementing the last stored line.
	int topLine;

	//allow text lines to wrap when exceeding width of window?
	bool allowTextLineWrapping;

private: //private methods

	//Indicate if given line index will overflow the space rect for the current topLine index.
	//Also calculate the required Y screen offset for drawing text objects for the current topLine index (add to the text object coordinate to obtain actual screen coordinate)
	bool LineOverflowsSpace(int lineIdx, float *pscreenYOffset = nullptr);

protected: //protected methods

	ActionOutcome NewMessage_CommonResponses(AC_ aCode, INT2 mouse, std::string data);

	//Set topLine so that the last line will not overflow the space rect. If it does, set topLine so that the last line is as far down in the console as possible
	void ResetTopLineIndex(void);

	//draw the textLines content. return true if lines had to be recalculated to fit the window width
	bool DrawTextLines(bool doDraw = true);

	//Faster version of Draw where we don't refresh any interactive objects but draw them in their current state. Also doesn't re-check for re-alignment of text objects, assumes everything is correct.
	//Thus this method purely draws the screen from current state, which is assumed to be in the correct.
	//use this whenever you know there cannot be any changes to interactive objects or any other settings (e.g. window dimensions)
	//when there are alot of interactive objects on screen the refresh rate can drop significantly making the interface sluggish, so use the full Draw method sparingly)
	void DrawTextLines_Quick(void);

	//set window size using ratios of screen dimensions / change window size by adjusting sides using delta values
	void ResizeWindow(D2D1_RECT_F ratios);
	void ChangeWindowSides(D2D1_RECT_F delta);

	//maximise / recall window size toggling : overloads the base method, calls it then adds extra settings for text
	void MaximiseWindowSize(void);
	void RestoreWindowSize(void);

	//toggle between the MaximiseWindowSize and RestoreWindowSize depending on the windowMaximised flag setting
	void Toggle_Maximise_Normal(INT2 mouse, ActionOutcome &actionResult);

	//if text line wrapping is set then text lines will not exceed the space rectangle, but will be split into multiple lines, each fitting the space rectangle (this is now called a split paragraph).
	void SetTextLineWrapping(bool flag) { allowTextLineWrapping = flag; if(allowTextLineWrapping) textLines.SetWidthLimit(spaceRectDims.width); else textLines.SetWidthLimit(0); }

	bool InLastPara(int lineIdx) { return (lineIdx >= LastPara() && lineIdx <= LastLine()); }

	//for given prompter character index in last paragraph, convert to line index and character index in that line.
	INT2 EntryLineIndex(int promptPos) { return textLines.GetParaLine_CharIdx(INT2(LastLine(), promptPos)); }

	//Set text in window space, erasing any existing text
	void SetText(const TextLine& textLine) { textLines.set(textLine); }

	//Set text in window space, erasing any existing text
	void SetText(TextLines& __textLines) 
	{ 
		ClearWindow();

		for (int idx = 0; idx < __textLines.size(); idx++) {

			textLines.push(__textLines[idx]);
		}
	}

	//Set text in window space, erasing any existing text
	void SetText(std::string text, FormatSpecifier fs);

	//Add text using text formatting (either on a new line or at the given line index as an insertion, pushing down all existing lines from that index onwards)
	void NewFormattedTextLine(int lineIdx, std::string text);

	void ClearWindow(void) { textLines.clear(); topLine = 0; }

public: //public methods

	TextDisplay(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pActionHandler);
	virtual ~TextDisplay() {}

	//index of last line
	int LastLine(void) { return textLines.LastLine(); }
	//index of start line of last paragraph
	int LastPara(void) { return textLines.LastPara(); }

	//get text line from lineIdx with coordinates as it appears on screen (same with just a single text object, but returned as a text line)
	TextLine GetRawTextLine_ScreenCoords(int lineIdx);
	TextLine GetRawTextLine_ScreenCoords(INT2 objIdx);

	TextObject& GetTextObjectRef(INT2 line_obj_idx) { return textLines[line_obj_idx.major][line_obj_idx.minor]; }

	//get text line index from given mouse coordinates
	int GetMouseClickLineIdx(INT2 mouse);

	//get text line and object index from given mouse coordinates
	INT2 GetMouseClickObjectIdx(INT2 mouse);

	//get text line and object index from given mouse coordinates
	INT3 GetMouseClickFullIdx(INT2 mouse);
};

#endif