#pragma once

#include "CompileFlags.h"
#if GRAPHICS == 1

#include "BorisLib.h"
#include "TextObject.h"
#include "BorisGraphics.h"
#include "WinSpaces.h"



////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// IOI Popup TextBox (carries information from one interactive object to another : inter-object interaction)

class BorisIOIPopupTextBox : 
	public TextDisplay 
{
	friend class BorisDisplay;

private:

	//the extracted, or popped-up paragraph
	TextLine poppedParagaph;

	//the space which spawned this Popup
	INT2 parent_winId;

	//the Id of the interactive object in the parent window which spawned this popup (as a result of an interaction with it) - by dragging this popup to another interactive object, an inter-object interaction can be implemented (in the action handler)
	INT2 IO_Id;

	//the index (line, object) of the IO_Id object in the parent window
	INT2 spawningObject_index;

protected:

	void DrawWindow(void);
	void DrawWindow_Quick(void) { DrawWindow(); }

	//Implementation of message handling routine
	ActionOutcome NewMessage(AC_ aCode, INT2 mouse, std::string data = "");

	void SetDefaultCursor(void) { SetCursor(LoadCursor(nullptr, IDC_HAND)); }

	void AdjustWindowDimensions(D2D1_RECT_F delta) { ChangeWindowSides(delta); }

public:

	//Setup new window space with rectangle specified by ratios of the entire screen
	BorisIOIPopupTextBox(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pActionHandler);
	~BorisIOIPopupTextBox() {}

	void ResetActionHandler(void);
};

#endif
