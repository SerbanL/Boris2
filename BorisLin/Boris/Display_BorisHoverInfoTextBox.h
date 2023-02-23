#pragma once

#include "CompileFlags.h"
#if GRAPHICS == 1

#include "BorisLib.h"
#include "TextObject.h"
#include "BorisGraphics.h"
#include "WinSpaces.h"



////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Hover Info TextBox (displays fixed text and is deleted as soon as the mouse leaves it)

class BorisHoverInfoTextBox : 
	public TextDisplay 
{
	friend class BorisDisplay;

protected:

	void DrawWindow(void);
	void DrawWindow_Quick(void) { DrawWindow(); }

	//Implementation of message handling routine
	ActionOutcome NewMessage(AC_ aCode, INT2 mouse, std::string data = "");

	void SetDefaultCursor(void) { SetCursor(LoadCursor(nullptr, IDC_HAND)); }

	void AdjustWindowDimensions(D2D1_RECT_F delta) { ChangeWindowSides(delta); }

public:

	//Setup new window space with rectangle specified by ratios of the entire screen
	BorisHoverInfoTextBox(D2D1_RECT_F ratios, INT2 winId);
	~BorisHoverInfoTextBox() {}
};

#endif
