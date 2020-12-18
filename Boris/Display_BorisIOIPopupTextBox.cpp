#include "stdafx.h"
#include "Display_BorisIOIPopupTextBox.h"

#if GRAPHICS == 1

////////////////////////////////////////////////////////////////////////////////////////////////////////////IOI POPUP TEXT BOX

BorisIOIPopupTextBox::BorisIOIPopupTextBox(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pActionHandler) : 
	TextDisplay(ratios, winId, pActionHandler) 
{
}

void BorisIOIPopupTextBox::DrawWindow(void) 
{
	//clear the drawing area
	pBG->FillRectangle(spaceRect, bgrndColor);

	//Draw text lines
	DrawTextLines();
}

ActionOutcome BorisIOIPopupTextBox::NewMessage(AC_ aCode, INT2 mouse, std::string data) 
{
	ActionOutcome actionResult;

	//First implement responses shared by all derived classes
	actionResult = TextDisplay::NewMessage_CommonResponses(aCode, mouse, data);

	//Next implement particular responses of this derived class
	switch (aCode) {

	case AC_MOUSEMOVE:
	{
		if (mouseLeftDown) {

			ShiftRect(spaceRect, FLT2(mouse - dragStart));
			dragStart = mouse;

			//this space has been moved : check for interactions with other objects
			actionResult.AddCode(winId, mouse, AO_CHECKIOINTERACTION);
		}

		if (!IsMouseInWindow(mouse))
			actionResult.SetCode(winId, mouse, AO_DESTROYWINDOW);
	}
	break;

	case AC_MOUSELEFTUP:
		//before destroying window, check if the popup was dropped on an object it's supposed to interact with
		actionResult.SetCode(winId, mouse, AO_CHECKIODROPINTERACTION);
		break;

	default:
		break;
	}

	return actionResult;
}

void BorisIOIPopupTextBox::ResetActionHandler(void) 
{
	for (int lineIdx = 0; lineIdx <= textLines.LastLine(); lineIdx++) {
		for (int objIdx = 0; objIdx <= textLines[lineIdx].LastElem(); objIdx++) {

			textLines[lineIdx][objIdx].ResetActionHandler();
		}
	}
}

#endif