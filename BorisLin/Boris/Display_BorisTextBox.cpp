#include "stdafx.h"
#include "Display_BorisTextBox.h"

#if GRAPHICS == 1

////////////////////////////////////////////////////////////////////////////////////////////////////////////DATABOX WINDOW

BorisTextBox::BorisTextBox(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pActionHandler) : 
	TextDisplay(ratios, winId, pActionHandler) 
{
}

void BorisTextBox::DrawWindow(void) 
{
	//clear the drawing area
	pBG->FillRectangle(spaceRect, bgrndColor);

	//Draw text lines
	DrawTextLines();

	DrawResizingFrame();
}

ActionOutcome BorisTextBox::NewMessage(AC_ aCode, INT2 mouse, std::string data) 
{
	ActionOutcome actionResult;

	//First implement responses shared by all derived classes
	actionResult = TextDisplay::NewMessage_CommonResponses(aCode, mouse, data);

	//Next implement particular responses of this derived class
	switch (aCode) {

	case AC_MOUSEWHEEL:
	{
		if (windowEntered) {

			actionResult.SetCode(winId, mouse, AO_REFRESH);

			int wheelDirection = ToNum(data);
			topLine -= wheelDirection;
			if (topLine < 0) topLine = 0;
			if (topLine > LastLine()) topLine = LastLine();
		}
		else actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
	}
	break;

	case AC_MOUSELEFTDOWN:
	{
		if (windowEntered) {

			if (IsDoubleClick(aCode)) {

				//Double click occured
				Toggle_Maximise_Normal(mouse, actionResult);
			}
			else {

				INT2 click = GetMouseClickObjectIdx(mouse);

				//clicked on a text object
				if (click.i >= 0 && click.j >= 0) {

					if (textLines[click.i][click.j].IsActionSet()) {

						//it's an interactive object. Send object interaction code. If the object then asks for an inter-object interaction to start, then send this signal further (AO_STARTINTERACTION)
						if (textLines[click.i][click.j].ObjectInteraction(aCode) == AO_STARTINTERACTION)
							actionResult.AddCode(winId, mouse, AO_STARTINTERACTION);
					}
				}
			}
		}
	}
	break;

	case AC_MOUSERIGHTDOWN:
	{
		if (windowEntered) {

			INT2 click = GetMouseClickObjectIdx(mouse);
			if (click.i >= 0 && click.j >= 0) {

				if (textLines[click.i][click.j].IsActionSet()) {

					//it's an interactive object. Send object interaction code : this will ask for this entry line to be deleted, by deleting entry in Simulation::dataBoxList
					textLines[click.i][click.j].ObjectInteraction(aCode);
				}
			}
		}
	}
	break;

	default:
		break;
	}

	return actionResult;
}

#endif