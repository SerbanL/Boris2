#include "stdafx.h"
#include "Display_BorisHoverInfoTextBox.h"

#if GRAPHICS == 1

////////////////////////////////////////////////////////////////////////////////////////////////////////////HOVER INFO TEXT BOX

BorisHoverInfoTextBox::BorisHoverInfoTextBox(D2D1_RECT_F ratios, INT2 winId) : 
	TextDisplay(ratios, winId, nullptr) 
{
}

void BorisHoverInfoTextBox::DrawWindow(void) 
{
	//clear the drawing area
	pBG->FillRectangle(spaceRect, bgrndColor);

	//Draw text lines
	DrawTextLines();
}

ActionOutcome BorisHoverInfoTextBox::NewMessage(AC_ aCode, INT2 mouse, std::string data) 
{
	ActionOutcome actionResult;

	//First implement responses shared by all derived classes
	actionResult = TextDisplay::NewMessage_CommonResponses(aCode, mouse, data);

	//Next implement particular responses of this derived class
	switch (aCode) {

	case AC_MOUSEMOVE:
	{
		if (!IsMouseInWindow(mouse))
			actionResult.SetCode(winId, mouse, AO_DESTROYWINDOW);
	}
	break;

	default:
		break;
	}

	return actionResult;
}

#endif