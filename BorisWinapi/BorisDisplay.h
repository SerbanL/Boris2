//Defines display organisation. Uses BorisGraphics (built on D3D) to actually draw things.
//
//Classes: 
//BorisDisplay is a container object.
//WinSpace is an abstract class setting out the general structure and minimum components of a window space (e.g. console, mesh display, data box etc.). 
//Specific spaces must implement this class and they are contained in BorisDisplay.

#pragma once

#include "CompileFlags.h"
#if GRAPHICS == 1

#include "BorisLib.h"
#include "TextObject.h"
#include "BorisGraphics.h"
#include "WinSpaces.h"
#include "PhysQRep.h"

//default text font size
#define FONT	"Consolas"		//nice monospaced font (default for Visual Studio)
#define TEXTSIZE 15
#define MESHTEXTSIZE	30

//number of pixels to shift the hover info text box away from mouse position (so the mouse is inside the window when created)
#define HOVERINFOTEXTBOX_MOUSESHIFT 1

//enum for different types of window spaces : these are the winId major identifier
enum WIN_ { WIN_ALL = -1, WIN_CONSOLE = 0, WIN_DATABOX, WIN_MESHWINDOW, WIN_IOIPOPUPTEXTBOX, WIN_HOVERINFOTEXTBOX };

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// CONSOLE (display text and allows user input)

class BorisConsole : public TextDisplay {

	friend class BorisDisplay;

private: //Private data

	//save user commands
	TextLines savedCommands;

	//index of command recall (index for savedCommands).
	int recallIdx;

	//prompter position in text on entry line (index of character)
	int promptPos;

	//text selection start and end indexes of characters on entry line (inclusive). Note, it is possible selectStart > selectEnd, these are just the start and end indexes of a text selection action.
	int selectStart, selectEnd;

	//prompter rectangle user to draw it
	D2D1_RECT_F promptRect;

	//rectangle of text selection box
	D2D1_RECT_F selectRect;

	bool mouseOverInteractiveObject;

public:  //Public data

private: //Private methods

		 //set prompter on a new line
	void SetNewLinePrompter(void);

	//calculate prompter rectangle at the current position
	void SetPrompterRect(bool recalculateTopLine = true);

	//decrement/increment prompter position on entry line. Return true if prompter position has changed.
	bool DecrementPrompter(void);
	bool IncrementPrompter(void);
	void SetPrompterEnd(void);
	void SetPrompterHome(void);

	//delete text on entry line between selectStart and selectEnd characters, inclusive, or if text selection not available just delete at prompter
	void DeleteTextFromEntryLine(void);

	void ResetTextSelector(void) { selectEnd = 0; selectStart = selectEnd; }

	//Enter pressed, return command line text in command (if any). Return true if there was text on line, otherwise false. If true then set a new line prompter.
	bool NewCommandEntered(string &command);

	//Text on the entry line must always be of the form: first word bold, further words italic non-bold. Must have USERCOLOR textcolor and BGRNDTEXTCOLOR with no outline.
	string FormatEntryLineText(string text);

	//make sure formatting of entry line is correct 
	void SetFormattingonEntryLine(void);

protected: //Protected methods

	void DrawWindow(void);

	//Implementation of message handling routine
	ActionOutcome NewMessage(AC_ aCode, INT2 mouse, string data = "");

	void SetDefaultCursor(void) { SetCursor(LoadCursor(nullptr, IDC_IBEAM)); }

	void AdjustWindowDimensions(D2D1_RECT_F delta) { ChangeWindowSides(delta); }

	//New user entry at prompter position
	void NewUserEntry(string text);

	//add new simple text line (fixed formatting)
	void NewTextLine(string text, FormatSpecifier fs);

	//overloads the TextDisplay base method to add extra functionality : new text line with variable formatting
	void NewFormattedTextLine(string text);

	//insert text at prompter by adding to respective TextObject.
	void InsertTextatPrompter(string text);

	//set text as the entry line text (formatted type)
	void SetEntryLineText(string text);

public:

	//Setup new window space with rectangle specified by ratios of the entire screen
	BorisConsole(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pConsoleActionHandler);
	~BorisConsole() {}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// DATABOX WINDOW

class BorisTextBox : public TextDisplay {

	friend class BorisDisplay;

private: //Private data

public:  //Public data

private: //Private methods

protected: //Protected methods

	void DrawWindow(void);

	//Implementation of message handling routine
	ActionOutcome NewMessage(AC_ aCode, INT2 mouse, string data = "");

	void SetDefaultCursor(void) { SetCursor(LoadCursor(nullptr, IDC_HAND)); }

	void AdjustWindowDimensions(D2D1_RECT_F delta) { ChangeWindowSides(delta); }

public:  //Public methods

	//Setup new window space with rectangle specified by ratios of the entire screen
	BorisTextBox(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pActionHandler);
	~BorisTextBox() {}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// IOI Popup TextBox (carries information from one interactive object to another : inter-object interaction)

class BorisIOIPopupTextBox : public TextDisplay {

	friend class BorisDisplay;

private: //Private data

	//the extracted, or popped-up paragraph
	TextLine poppedParagaph;

	//the space which spawned this Popup
	INT2 parent_winId;

	//the Id of the interactive object in the parent window which spawned this popup (as a result of an interaction with it) - by dragging this popup to another interactive object, an inter-object interaction can be implemented (in the action handler)
	INT2 IO_Id;

	//the index (line, object) of the IO_Id object in the parent window
	INT2 spawningObject_index;

public:  //Public data

private: //Private methods

protected: //Protected methods

	void DrawWindow(void);

	//Implementation of message handling routine
	ActionOutcome NewMessage(AC_ aCode, INT2 mouse, string data = "");

	void SetDefaultCursor(void) { SetCursor(LoadCursor(nullptr, IDC_HAND)); }

	void AdjustWindowDimensions(D2D1_RECT_F delta) { ChangeWindowSides(delta); }

public:  //Public methods

	//Setup new window space with rectangle specified by ratios of the entire screen
	BorisIOIPopupTextBox(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pActionHandler);
	~BorisIOIPopupTextBox() {}

	void ResetActionHandler(void);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Hover Info TextBox (displays fixed text and is deleted as soon as the mouse leaves it)

class BorisHoverInfoTextBox : public TextDisplay {

	friend class BorisDisplay;

protected:

	void DrawWindow(void);

	//Implementation of message handling routine
	ActionOutcome NewMessage(AC_ aCode, INT2 mouse, string data = "");

	void SetDefaultCursor(void) { SetCursor(LoadCursor(nullptr, IDC_HAND)); }

	void AdjustWindowDimensions(D2D1_RECT_F delta) { ChangeWindowSides(delta); }

public: 

	//Setup new window space with rectangle specified by ratios of the entire screen
	BorisHoverInfoTextBox(D2D1_RECT_F ratios, INT2 winId);
	~BorisHoverInfoTextBox() {}
};

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Boris Display container and director object. All access to display elements and drawing done through here.

class BorisDisplay :
	public GraphicalObject,
	public Threads<BorisDisplay>,
	public ProgramState<BorisDisplay,
	tuple<BorisGraphics*, BorisMeshWindow*>, tuple<>>
{
private: //Private data

	//entire screen rectangle
	D2D1_RECT_F screenRect;

	//background color of screen
	D2D1_COLOR_F bgrndColor;

	//list of boris window spaces. The order of these can change. Index 0 specifies top-most window (drawn last, i.e. on top of everything else), etc.
	vector_lut<WinSpace*> bWin;

	//also need pointers to special window spaces which must always exist: Mesh Display, Console, Data Box (so specific methods can be accessed easily in these classes derived from the abstract WinSpace)
	BorisMeshWindow *pbMeshWin;
	BorisConsole *pbConsole;
	BorisTextBox *pbDataBox;

	//pairs of windows with linked sides : these sides have coordinates which are maintained during resizing. The order of winIds in the pair is important.
	vector<PAIR> winLinkages_RightLeft, winLinkages_BottomTop;

	//need thread-safe access to BorisDisplay - guard all external points of access (public methods)
	mutex displayMutex;

public:  //Public data

private: //Private methods

	//Bring the window space currently at currPos in bWin vector to the top.
	void BringWindowToFront(INT2 winId);
	void BringWindowToFront(int currPos);

	//Create new window of type specified by WIN_ identifier, dimensions given by ratios of screen dimensions. The window space is added to bWin vector, but also returrned by this method.
	WinSpace* AddWinSpace(WIN_ winIdMajor, D2D1_RECT_F ratios, SimTOFunct *pActionHandler = nullptr);

	//Delete window space with given winId unique identifier
	void DelWinSpace(INT2 winId);

	//make a popup text box which carries a paragraph of text (this will be a representation of an entry in a data structure) from the parent window. 
	//Also carry an interactive object Id (this is the object which called for this popup to be made). By dragging this popup to another interactive object, an inter-object interaction is achieved.
	void MakeIOIPopupTextBox(INT2 parent_winId, INT2 mouse);

	//For the given popup try to interact it with an interactive object
	void PopupTextBoxInteraction(INT2 popupId, INT2 mouse);

	//For the given popup try to interact it with an interactive object, then destroy the popup
	void PopupTextBoxDropInteraction(INT2 popupId, INT2 mouse);

	//Make a text box which displays some fixed info. The window is deleted as soon as the mouse leaves it.
	void MakeHoverInfoTextBox(INT2 mouse, string formatted_text_info_string);

	//link 2 windows on their left - right sides, or bottom - top sides The order is important.
	void LinkWindowSides_RightLeft(INT2 winId_RightSide, INT2 winId_LeftSide) { winLinkages_RightLeft.push_back(PAIR(winId_RightSide, winId_LeftSide)); }
	void LinkWindowSides_BottomTop(INT2 winId_BottomSide, INT2 winId_TopSide) { winLinkages_BottomTop.push_back(PAIR(winId_BottomSide, winId_TopSide)); }
	void ClearWindowLinkages(void) { winLinkages_RightLeft.resize(0); winLinkages_BottomTop.resize(0); }

	//based on the above linkages, set all sides which are linked to winId
	void RecalculateWindowLinks(INT2 winId);

	void SendHoverCheck_ThreadSafe(INT2 mouse) { displayMutex.lock(); NewMessage(AC_HOVERCHECK, mouse); displayMutex.unlock(); }

	//-----------------------------------------------Advanced display methods

	//run this on a thread to animate the transition between mesh view of *pphysQ with start and end settings. *parameter varies between 0 and 1 and increases based on start time and duration
	void AnimateMeshViewChange(vector<PhysQ>* pphysQ, PhysQRepSettings start, PhysQRepSettings end, double* parameter, DWORD start_time_ms, DWORD duration_ms);

	//-----------------------------------------------Mirror public methods - accessed through their thread-safe callers from outside BorisDisplay; available here for internal use.

	void Refresh(int winIdMajor = WIN_ALL);

	//Process message received from WndProc
	ActionOutcome NewMessage(AC_ aCode, INT2 mouse, string data = "");

public:

	BorisDisplay(HWND hWnd, SimTOFunct *pConsoleActionHandler);

	~BorisDisplay();

	void RepairObjectState(void) {}

	double Get_MeshDisplay_DetailLevel(void) { return pbMeshWin->Get_PhysQRep_DetailLevel(); }

	//-----------------------------------------------Thread-safe access points to BorisDisplay

	void Refresh_ThreadSafe(void) { displayMutex.lock(); Refresh(); displayMutex.unlock(); }

	void ClearScreen(void) { displayMutex.lock(); pbConsole->ClearWindow(); pbDataBox->ClearWindow(); Refresh(); displayMutex.unlock(); }
	void ClearDataBox(void) { displayMutex.lock(); pbDataBox->ClearWindow(); Refresh(); displayMutex.unlock(); }

	//Process message received from WndProc
	ActionOutcome NewMessage_ThreadSafe(AC_ aCode, INT2 mouse, string data = "") { displayMutex.lock(); ActionOutcome ao = NewMessage(aCode, mouse, data); displayMutex.unlock(); return ao; }

	//display various message types in the console
	void DisplayConsoleMessage(string text);
	void DisplayConsoleError(string text);
	void DisplayConsoleListing(string text);
	void DisplayFormattedConsoleMessage(string text);

	//add a new entry in the data box as an interative object (use formatted text)
	void NewDataBoxField(string formattedText);
	//update the data box at given lineIdx to the given value, passed as a string
	void UpdateDataBoxField(int lineIdx, string value_string);

	//set the console line entry text
	void SetConsoleEntryLineText(string text) { displayMutex.lock(); pbConsole->SetEntryLineText(text); Refresh(WIN_CONSOLE); displayMutex.unlock(); }

	//adjust display for default view of physical quantity representation - animate view change
	void AutoSetMeshDisplaySettings(vector<PhysQ> physQ);

	void AutoSetMeshDisplaySettings_KeepOrientation(vector<PhysQ> physQ);

	//adjust display for default view of physical quantity representation - no animation used, just set.
	void AutoSetMeshDisplaySettings_Sudden(vector<PhysQ> physQ) { displayMutex.lock(); pbMeshWin->AutoSetMeshDisplaySettings(physQ); displayMutex.unlock(); }

	//Calculate graphical representation of given physical quantity and refresh screen
	//only update if mutex can be locked right away, as this is used in performance-critical parts
	void UpdateMeshDisplay(vector<PhysQ> physQ) { if (displayMutex.try_lock()) { pbMeshWin->UpdatePhysQRep(physQ); displayMutex.unlock(); } }

	//image_cropping specify normalized cropping within the mesh image, as left, bottom, right, top : 0, 0 point is left, bottom of screen as far as user is concerned.
	bool SaveMeshImage(string fileName, DBL4 image_cropping = DBL4(0, 0, 1, 1)) { displayMutex.lock(); bool success = pbMeshWin->SaveMeshImage(fileName, image_cropping); displayMutex.unlock(); return success; }

	//allows access to Boris Graphics public methods.
	BorisGraphics* BGMethods(void) { return pBG; }
};

#endif