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

#include "Display_BorisConsole.h"
#include "Display_BorisTextBox.h"
#include "Display_BorisIOIPopupTextBox.h"
#include "Display_BorisIOIPopupEditBox.h"
#include "Display_BorisHoverInfoTextBox.h"
#include "Display_BorisMeshWindow.h"

#include "Display_Defs.h"

//enum for different types of window spaces : these are the winId major identifier
enum WIN_ { WIN_ALL = -1, WIN_NONE, WIN_CONSOLE, WIN_DATABOX, WIN_MESHWINDOW, WIN_IOIPOPUPTEXTBOX, WIN_IOIPOPUPEDITBOX, WIN_HOVERINFOTEXTBOX };



////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Boris Display container and director object. All access to display elements and drawing done through here.

class BorisDisplay :
	public GraphicalObject,
	public Threads<BorisDisplay>,
	public ProgramState<BorisDisplay,
	std::tuple<BorisGraphics*, BorisMeshWindow*>, std::tuple<>>
{
private:

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
	std::vector<PAIR> winLinkages_RightLeft, winLinkages_BottomTop;

	//need thread-safe access to BorisDisplay - guard all external points of access (public methods)
	std::mutex displayMutex;

private:

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

	//make a popup edit box. After editing the text send it back to the original interactive object to use.
	void MakeIOIPopupEditBox(INT2 parent_winId, INT2 mouse);

	//For the given popup try to interact it with an interactive object
	void PopupTextBoxInteraction(INT2 popupId, INT2 mouse);

	//For the given popup try to interact it with an interactive object, then destroy the popup
	void PopupTextBoxDropInteraction(INT2 popupId, INT2 mouse);

	//A popup edit box has returned text (e.g. after pressing enter)
	void PopupEditBoxReturnedText(INT2 popupId, INT2 mouse);

	//Make a text box which displays some fixed info. The window is deleted as soon as the mouse leaves it.
	void MakeHoverInfoTextBox(INT2 mouse, std::string formatted_text_info_string);

	//link 2 windows on their left - right sides, or bottom - top sides The order is important.
	void LinkWindowSides_RightLeft(INT2 winId_RightSide, INT2 winId_LeftSide) { winLinkages_RightLeft.push_back(PAIR(winId_RightSide, winId_LeftSide)); }
	void LinkWindowSides_BottomTop(INT2 winId_BottomSide, INT2 winId_TopSide) { winLinkages_BottomTop.push_back(PAIR(winId_BottomSide, winId_TopSide)); }
	void ClearWindowLinkages(void) { winLinkages_RightLeft.resize(0); winLinkages_BottomTop.resize(0); }

	//based on the above linkages, set all sides which are linked to winId
	void RecalculateWindowLinks(INT2 winId);

	void SendHoverCheck_ThreadSafe(INT2 mouse) { displayMutex.lock(); NewMessage(AC_HOVERCHECK, mouse); displayMutex.unlock(); }

	//-----------------------------------------------Advanced display methods

	//run this on a thread to animate the transition between mesh view of *pphysQ with start and end settings. *parameter varies between 0 and 1 and increases based on start time and duration
	void AnimateMeshViewChange(std::vector<PhysQ>* pphysQ, PhysQRepSettings start, PhysQRepSettings end, double* parameter, DWORD start_time_ms, DWORD duration_ms);

	//-----------------------------------------------Mirror public methods - accessed through their thread-safe callers from outside BorisDisplay; available here for internal use.

	//refresh display, including recalculating interactive objects
	void Refresh(int winIdMajor = WIN_ALL);

	//draw display but do not recalculate interactive objects
	void Draw(int winIdMajor = WIN_ALL);

	//Process message received from WndProc
	ActionOutcome NewMessage(AC_ aCode, INT2 mouse, std::string data = "");

public:

	BorisDisplay(HWND hWnd, SimTOFunct *pConsoleActionHandler);

	~BorisDisplay();

	void RepairObjectState(void) {}

	double Get_MeshDisplay_DetailLevel(void) { return pbMeshWin->Get_PhysQRep_DetailLevel(); }

	//-----------------------------------------------Thread-safe access points to BorisDisplay

	void Refresh_ThreadSafe(void) { displayMutex.lock(); Refresh(); displayMutex.unlock(); }
	void Draw_ThreadSafe(void) { displayMutex.lock(); Draw(); displayMutex.unlock(); }

	void ClearScreen(void) { displayMutex.lock(); pbConsole->ClearWindow(); pbDataBox->ClearWindow(); Refresh(); displayMutex.unlock(); }
	void ClearDataBox(void) { displayMutex.lock(); pbDataBox->ClearWindow(); Refresh(); displayMutex.unlock(); }

	//Process message received from WndProc
	ActionOutcome NewMessage_ThreadSafe(AC_ aCode, INT2 mouse, std::string data = "") { displayMutex.lock(); ActionOutcome ao = NewMessage(aCode, mouse, data); displayMutex.unlock(); return ao; }

	//display various message types in the console
	void DisplayConsoleMessage(std::string text);
	void DisplayConsoleError(std::string text);
	void DisplayConsoleWarning(std::string text);
	void DisplayConsoleListing(std::string text);
	void DisplayFormattedConsoleMessage(std::string text);

	//add a new entry in the data box as an interative object (use formatted text)
	void NewDataBoxField(std::string formattedText);
	//update the data box at given lineIdx to the given value, passed as a std::string
	void UpdateDataBoxField(int lineIdx, std::string value_string);

	//set the console line entry text
	void SetConsoleEntryLineText(std::string text) { displayMutex.lock(); pbConsole->SetEntryLineText(text); Draw(WIN_CONSOLE); displayMutex.unlock(); }

	//adjust display for default view of physical quantity representation - animate view change
	void AutoSetMeshDisplaySettings(std::vector<PhysQ> physQ);

	void AutoSetMeshDisplaySettings_KeepOrientation(std::vector<PhysQ> physQ);

	//adjust display for default view of physical quantity representation - no animation used, just set.
	void AutoSetMeshDisplaySettings_Sudden(std::vector<PhysQ> physQ) { displayMutex.lock(); pbMeshWin->AutoSetMeshDisplaySettings(physQ); displayMutex.unlock(); }

	//Calculate graphical representation of given physical quantity and refresh screen
	//only update if std::mutex can be locked right away, as this is used in performance-critical parts
	void UpdateMeshDisplay_Auxiliary(std::vector<PhysQ> physQ) { if (displayMutex.try_lock()) { pbMeshWin->UpdatePhysQRep(physQ); Refresh(); displayMutex.unlock(); } }
	void UpdateMeshDisplay(std::vector<PhysQ> physQ, bool asynchronous)
	{ 
		if (asynchronous) single_call_launch<std::vector<PhysQ>>(&BorisDisplay::UpdateMeshDisplay_Auxiliary, physQ, THREAD_DISPLAY);
		else UpdateMeshDisplay_Auxiliary(physQ);
	}

	//Calculate graphical representation of given physical quantity and draw screen (so doesn't update interactive objects, just draws)
	//only update if std::mutex can be locked right away, as this is used in performance-critical parts
	void UpdateMeshDisplay_Quick_Auxiliary(std::vector<PhysQ> physQ) { if (displayMutex.try_lock()) { pbMeshWin->UpdatePhysQRep(physQ); Draw(); displayMutex.unlock(); } }
	void UpdateMeshDisplay_Quick(std::vector<PhysQ> physQ, bool asynchronous) 
	{ 
		if (asynchronous) single_call_launch<std::vector<PhysQ>>(&BorisDisplay::UpdateMeshDisplay_Quick_Auxiliary, physQ, THREAD_DISPLAY);
		else UpdateMeshDisplay_Quick_Auxiliary(physQ);
	}

	//image_cropping specify normalized cropping within the mesh image, as left, bottom, right, top : 0, 0 point is left, bottom of screen as far as user is concerned.
	bool SaveMeshImage(std::string fileName, DBL4 image_cropping = DBL4(0, 0, 1, 1)) { displayMutex.lock(); bool success = pbMeshWin->SaveMeshImage(fileName, image_cropping); displayMutex.unlock(); return success; }

	bool SaveImage(std::string fileName, std::vector<PhysQ> physQ) { displayMutex.lock(); bool success = pbMeshWin->SaveImage(fileName, physQ); displayMutex.unlock(); return success; }

	//Set mesh display detail level directly
	void SetDetailLevel(double detail_level) { displayMutex.lock(); pbMeshWin->SetDetailLevel(detail_level); displayMutex.unlock(); }
	double GetDetailLevel(void) { return pbMeshWin->GetDetailLevel();	}

	//camera control using methods rather than mouse messages
	void RotateCameraAboutOrigin(float dAzim, float dPolar) { displayMutex.lock(); pBG->RotateCameraAboutOrigin(dAzim, dPolar); displayMutex.unlock(); }
	void AdjustCameraDistanceFromOrigin(float dZ) { displayMutex.lock(); pBG->AdjustCameraDistanceFromOrigin(dZ); displayMutex.unlock(); }
	void RotateCameraView(float dAngle) { displayMutex.lock(); pBG->RotateCameraView(dAngle); displayMutex.unlock(); }
	void Shift3DOriginPixelPosition(float dX, float dY) { displayMutex.lock(); pBG->Shift3DOriginPixelPosition(dX, dY); displayMutex.unlock(); }

	//Set mesh display render threshold values for faster rendering when we have many cells
	void SetRenderThresholds(INT3 renderthresholds) { displayMutex.lock(); pbMeshWin->SetRenderThresholds(renderthresholds); displayMutex.unlock(); }
	INT3 GetRenderThresholds(void) { return pbMeshWin->GetRenderThresholds(); }

	//allows access to Boris Graphics public methods.
	BorisGraphics* BGMethods(void) { return pBG; }

	//Lock and unlock the display std::mutex : calling this to ensure that any code that follows is not executed during a display refresh or before any other public BorisDisplay method has finished.
	void WaitForDisplayEnd(void) { displayMutex.lock(); displayMutex.unlock(); }
};

#endif
