#include "stdafx.h"
#include "WinSpaces.h"

#if GRAPHICS == 1

////////////////////////////////////////////////////////////////////////////////////////////////////////////BORIS WIN SPACE BASE CLASS

WinSpace::WinSpace(D2D1_RECT_F ratios, INT2 winId) : GraphicalObject() {

	this->winId = winId;

	//Space rectangle in pixels
	ResizeWindow(ratios);
}

void WinSpace::WindowLeft(void) {

	windowEntered = false;
	mouseLeftDown = false;
	mouseMiddleDown = false;
	mouseRightDown = false;

	resizing.ResetResize();
}

void WinSpace::WindowEntered(void) {

	windowEntered = true;
	SetDefaultCursor();
}

void WinSpace::ResizeWindow(D2D1_RECT_F ratios) {

	//Space rectangle in pixels
	spaceRect.bottom = (pBG->wndHeight) * ratios.bottom;
	spaceRect.left = (pBG->wndWidth) * ratios.left;
	spaceRect.top = (pBG->wndHeight) * ratios.top;
	spaceRect.right = (pBG->wndWidth) * ratios.right;

	spaceRectDims.width = spaceRect.right - spaceRect.left;
	spaceRectDims.height = spaceRect.bottom - spaceRect.top;
}

void WinSpace::ChangeWindowSides(D2D1_RECT_F delta) {

	if (spaceRect.left + delta.left < spaceRect.right)
		spaceRect.left += delta.left;

	if (spaceRect.right + delta.right > spaceRect.left)
		spaceRect.right += delta.right;

	if (spaceRect.bottom + delta.bottom > spaceRect.top)
		spaceRect.bottom += delta.bottom;

	if (spaceRect.top + delta.top < spaceRect.bottom)
		spaceRect.top += delta.top;

	spaceRectDims.width = spaceRect.right - spaceRect.left;
	spaceRectDims.height = spaceRect.bottom - spaceRect.top;
}

void WinSpace::MaximiseWindowSize(void) {

	saveWinSpaceSettings = WinSpaceSettings(spaceRect, spaceRectDims, stickyDisplayBottom, stickyDisplayTop);

	ResizeWindow(D2D1::RectF(0, 0, 1, 1));

	SetStickyDisplayBottom(false);

	windowMaximised = true;
}

void WinSpace::RestoreWindowSize(void) {

	saveWinSpaceSettings.RecallSettings(spaceRect, spaceRectDims, stickyDisplayBottom, stickyDisplayTop);

	windowMaximised = false;
}

void WinSpace::DrawFrame(void) 
{
	pBG->FillRectangle(D2D1::RectF(spaceRect.left - RESIZEFRAMEBORDER, spaceRect.top - RESIZEFRAMEBORDER, spaceRect.left, spaceRect.bottom + RESIZEFRAMEBORDER), FRAMECOLOR);
	pBG->FillRectangle(D2D1::RectF(spaceRect.right, spaceRect.top - RESIZEFRAMEBORDER, spaceRect.right + RESIZEFRAMEBORDER, spaceRect.bottom + RESIZEFRAMEBORDER), FRAMECOLOR);
	pBG->FillRectangle(D2D1::RectF(spaceRect.left, spaceRect.top - RESIZEFRAMEBORDER, spaceRect.right, spaceRect.top), FRAMECOLOR);
	pBG->FillRectangle(D2D1::RectF(spaceRect.left, spaceRect.bottom, spaceRect.right, spaceRect.bottom + RESIZEFRAMEBORDER), FRAMECOLOR);
}

void WinSpace::DrawResizingFrame(void) {

	if (resizing.IsResizeCursorSet()) {

		if (mouseLeftDown || IsHoveringOverResizeArea()) {

			pBG->FillRectangle(D2D1::RectF(spaceRect.left, spaceRect.top, spaceRect.left + RESIZEFRAMEBORDER, spaceRect.bottom), resizingFrameColor);
			pBG->FillRectangle(D2D1::RectF(spaceRect.right - RESIZEFRAMEBORDER, spaceRect.top, spaceRect.right, spaceRect.bottom), resizingFrameColor);
			pBG->FillRectangle(D2D1::RectF(spaceRect.left, spaceRect.top, spaceRect.right, spaceRect.top + RESIZEFRAMEBORDER), resizingFrameColor);
			pBG->FillRectangle(D2D1::RectF(spaceRect.left, spaceRect.bottom - RESIZEFRAMEBORDER, spaceRect.right, spaceRect.bottom), resizingFrameColor);
		}
	}
}

void WinSpace::Toggle_Maximise_Normal(INT2 mouse, ActionOutcome &actionResult) {

	if (allowMaximising) {

		if (!windowMaximised) {

			MaximiseWindowSize();
			actionResult.AddCode(winId, mouse, AO_SETTOPMOST);
		}
		else RestoreWindowSize();

		actionResult.AddCode(winId, mouse, AO_REFRESH);
	}
}

ActionOutcome WinSpace::NewMessage_CommonResponses(AC_ aCode, INT2 mouse, std::string data) {

	//Windows event occured in this space: process it. This is called by all derived classes and defines common responses, so only add behavior here that all classes should implement. Further responses to same windows message may follow after this.
	//
	//If you want to force other windows to try to handle a message then set AO_NOTHANDLED code.
	//
	//Messages which set AO_NOTHANDLED code and do nothing else, are meant to be handled only by a particular window. In this case do not handle those messages anywhere else, except in that window. The program will get there
	//eventually, irrespective of where the mouse is or window position is.
	//
	//It is possible to be here without the mouse being inside this window, so check for mouse position if you only want to handle the message if the mouse is here. 
	//Normally you would check windowEntered flag for this, except if you want to set this window as the top window. It could be another window is the top window, in which case windowEntered = false here.
	//
	//In NewMessage routines you should always check for windowEntered flag before handling a message, unless you want to be able to handle a message even when another window is on top.
	//In this case you need to use the AO_NOTHANDLED code here for that message (as described above), irrespective of which window this is. The intended particular NewMessage routine will respond to it, as no other routine will handle the message.

	ActionOutcome actionResult;

	switch (aCode) {

	case AC_MOUSELEFTDOWN:
	{
		if (IsMouseInWindow(mouse)) {

			if (!windowEntered) {
				//window just entered ... ask for this window to become the top window
				actionResult.SetCode(winId, mouse, AO_SETTOPMOST);
				//set default prompter and flag
				WindowEntered();
			}

			dragStart = mouse;
			mouseLeftDown = true;
		}
		else actionResult.SetCode(winId, mouse, AO_NOTHANDLED);	//not in this window, so give other windows a chance to handle this
	}
	break;

	case AC_MOUSELEFTUP:
	{
		if (windowEntered) {

			mouseLeftDown = false;

			if (resizing.IsResizeCursorSet()) {

				resizing.ResetResize();
				SetDefaultCursor();
			}

			actionResult.SetCode(winId, mouse, AO_REFRESH);
		}
		else actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
	}
	break;

	case AC_MOUSERIGHTDOWN:
	{
		if (IsMouseInWindow(mouse)) {

			if (!windowEntered) {

				actionResult.SetCode(winId, mouse, AO_SETTOPMOST);
				WindowEntered();
			}

			dragStart = mouse;
			mouseRightDown = true;
		}
		else actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
	}
	break;

	case AC_MOUSERIGHTUP:
	{
		if (windowEntered) {

			mouseRightDown = false;
		}
		else actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
	}
	break;

	case AC_MOUSEMIDDLEDOWN:
	{
		if (IsMouseInWindow(mouse)) {

			if (!windowEntered) {

				actionResult.SetCode(winId, mouse, AO_SETTOPMOST);
				WindowEntered();
			}

			dragStart = mouse;
			mouseMiddleDown = true;
		}
		else actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
	}
	break;

	case AC_MOUSEMIDDLEUP:
	{
		if (windowEntered) {

			mouseMiddleDown = false;
		}
		else actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
	}
	break;

	case AC_MOUSEMOVE:
	{
		if (mouseLeftDown && resizing.IsResizeCursorSet()) {

			ChangeWindowSides(resizing.GetResizingDeltas(mouse - dragStart));

			dragStart = mouse;

			actionResult.SetCodes(winId, mouse, AO_REFRESH, AO_WINRESIZED);
		}
		else {

			if (IsMouseInWindow(mouse)) {

				//set window entered if not set already
				if (!windowEntered) {

					actionResult.SetCode(winId, mouse, AO_SETTOPMOST);
					WindowEntered();
				}

				//check for resizing area
				if (!mouseLeftDown) {

					//set resizing cursor as required
					switch (resizing.MouseNewEntryInResizingArea(mouse, spaceRect)) {

					case RSZ_JUSTENTERED:
						//call for a timed delayed refresh - if mouse is still in resizing area then resizing frame will show up
						resizeHoverTimeCounter = GetSystemTickCount();
						actionResult.AddCode(winId, mouse, AO_DELAYEDREFRESH);
						break;

					case RSZ_JUSTLEFT:
						SetDefaultCursor();
						actionResult.AddCode(winId, mouse, AO_REFRESH);
						break;

					default:
						break;

					}
				}
			}
			else {

				actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
			}
		}
	}
	break;

	case AC_KEYBOARD:
	case AC_KBDENTER:
	case AC_KBDEND:
	case AC_KBDSHIFTEND:
	case AC_KBDHOME:
	case AC_KBDSHIFTHOME:
	case AC_KBDDEL:
	case AC_KBDBACK:
	case AC_KBDUP:
	case AC_KBDDN:
	case AC_KBDLEFT:
	case AC_KBDSHIFTLEFT:
	case AC_KBDRIGHT:
	case AC_KBDSHIFTRIGHT:
	case AC_KBDPGUP:
	case AC_KBDPGDN:
	case AC_KBDESC:
	case AC_CTRL_C:
	case AC_CTRL_V:
	case AC_DROPFILES:
	case AC_SHIFT_MOUSELEFTDOWN:
		//These codes can only be handled by particular windows, so respond to them there. Notice windowEntered flag not checked.
		actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		break;

	default:
		break;
	}

	return actionResult;
}

bool WinSpace::IsDoubleClick(AC_ &aCode) {

	if (GetSystemTickCount() - dblClickTimeCounter <= DOUBLECLICKMS) {

		aCode = AC_DOUBLECLICK;
		return true;
	}
	else {

		dblClickTimeCounter = GetSystemTickCount();
		return false;
	}
}

bool WinSpace::IsHoveringOverResizeArea(void) {

	return (GetSystemTickCount() - resizeHoverTimeCounter >= HOVERTIME);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////TEXT DISPLAY DERIVED ABSTRACT CLASS

TextDisplay::TextDisplay(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pActionHandler) : 
	WinSpace(ratios, winId), textLines(pActionHandler) 
{

	topLine = 0;

	allowTextLineWrapping = true;

	textLines.set(TextLine("", FormatSpecifier(USERCOLOR, BGRNDTEXTCOLOR, false, true, TOO_NONE)));

	ResizeWindow(ratios);
}

ActionOutcome TextDisplay::NewMessage_CommonResponses(AC_ aCode, INT2 mouse, std::string data) {

	ActionOutcome actionResult = WinSpace::NewMessage_CommonResponses(aCode, mouse, data);

	switch (aCode) {

	case AC_MOUSEMOVE:
	{
		if (mouseLeftDown && resizing.IsResizeCursorSet()) {

			ChangeWindowSides(resizing.GetResizingDeltas(mouse - dragStart));

			dragStart = mouse;

			actionResult.AddCode(winId, mouse, AO_REFRESH);
		}
	}
	break;
	}

	return actionResult;
}

bool TextDisplay::DrawTextLines(bool doDraw) 
{
	//Draw text lines in view. If lines had to be recalculated so they fit in window width, then reset top-line index
	return textLines.Draw(topLine, spaceRectDims.height, spaceRect.left, spaceRect.top, doDraw);
}

//Faster version of Draw where we don't refresh any interactive objects but draw them in their current state. Also doesn't re-check for re-alignment of text objects, assumes everything is correct.
//Thus this method purely draws the screen from current state, which is assumed to be in the correct state.
//use this whenever you know there cannot be any changes to interactive objects or any other settings (e.g. window dimensions)
//when there are alot of interactive objects on screen the refresh rate can drop significantly making the interface sluggish, so use the full Draw method sparingly)
void TextDisplay::DrawTextLines_Quick(void)
{
	//Draw text lines in view.
	textLines.Draw_Quick(topLine, spaceRectDims.height, spaceRect.left, spaceRect.top);
}

void TextDisplay::ResizeWindow(D2D1_RECT_F ratios) {

	WinSpace::ResizeWindow(ratios);

	if (allowTextLineWrapping) textLines.SetWidthLimit(spaceRectDims.width);
	DrawTextLines(false);
}

void TextDisplay::ChangeWindowSides(D2D1_RECT_F delta) {

	WinSpace::ChangeWindowSides(delta);

	if (allowTextLineWrapping) textLines.SetWidthLimit(spaceRectDims.width);
	DrawTextLines(false);
}

void TextDisplay::MaximiseWindowSize(void) {

	WinSpace::MaximiseWindowSize();

	if (allowTextLineWrapping) textLines.SetWidthLimit(spaceRectDims.width);
	DrawTextLines(false);
}

void TextDisplay::RestoreWindowSize(void) {

	WinSpace::RestoreWindowSize();

	if (allowTextLineWrapping) textLines.SetWidthLimit(spaceRectDims.width);
	DrawTextLines(false);
}

void TextDisplay::Toggle_Maximise_Normal(INT2 mouse, ActionOutcome &actionResult) {

	if (allowMaximising) {

		if (!windowMaximised) {

			MaximiseWindowSize();
			actionResult.AddCode(winId, mouse, AO_SETTOPMOST);
		}
		else RestoreWindowSize();

		actionResult.AddCode(winId, mouse, AO_REFRESH);
	}
}

bool TextDisplay::LineOverflowsSpace(int lineIdx, float *pscreenYOffset) {

	float screenYOffset = spaceRect.top;
	if (pscreenYOffset)
		*pscreenYOffset = screenYOffset;

	if (lineIdx > textLines.LastLine()) return true;

	//find Y coordinate offset which must be used when drawing text objects: TextObject coordinates are relative to the space rect, and the topLine line must always display at the top of the space rect
	screenYOffset -= textLines[topLine].top();
	if (pscreenYOffset)
		*pscreenYOffset = screenYOffset;

	//check if this line's actual screen coordinates will overflow the space rect
	if (textLines[lineIdx].bottom() + screenYOffset > spaceRect.bottom) return true;

	return false;
}

void TextDisplay::ResetTopLineIndex(void) {

	//Set topLine so that the new line will not overflow the space rect

	int lastLine = LastLine();

	if (LineOverflowsSpace(lastLine)) {

		for (topLine = lastLine; topLine >= 0; topLine--) {

			if (LineOverflowsSpace(lastLine)) {

				topLine++;
				break;
			}
		}

		//NOTE: the possibility of topLine = -1 should never occur normally. It could happen if the space window is so small that a single line cannot be displayed (no text can be displayed then)
		if (topLine < 0) topLine = lastLine;
	}
}

int TextDisplay::GetMouseClickLineIdx(INT2 mouse) {

	INT2 click = textLines.GetTextObjectIdx(mouse - INT2((int)spaceRect.left, (int)spaceRect.top), topLine);

	return click.i;
}

INT2 TextDisplay::GetMouseClickObjectIdx(INT2 mouse) {

	INT2 click = textLines.GetTextObjectIdx(mouse - INT2((int)spaceRect.left, (int)spaceRect.top), topLine);

	return click;
}

INT3 TextDisplay::GetMouseClickFullIdx(INT2 mouse) {

	INT3 click = textLines.GetFullIdx(mouse - INT2((int)spaceRect.left, (int)spaceRect.top), topLine);

	return click;
}

TextLine TextDisplay::GetRawTextLine_ScreenCoords(int lineIdx) {

	//get text line from lineIdx with coordinates as it appears on screen

	TextLine textLine;
	textLine = textLines[lineIdx];

	textLine.set_top_left(spaceRect.left, spaceRect.top + textLines[lineIdx].top() - textLines[topLine].top());

	return textLine;
}

TextLine TextDisplay::GetRawTextLine_ScreenCoords(INT2 objIdx) 
{
	TextLine textLine_singleObject;
	TextObject textObject = textLines[objIdx.i][objIdx.j];
	textLine_singleObject.push(textObject);

	textLine_singleObject.set_top_left(spaceRect.left + textLines[objIdx.i][objIdx.j].left(), spaceRect.top + textLines[objIdx.i].top() - textLines[topLine].top());

	return textLine_singleObject;
}

void TextDisplay::SetText(std::string text, FormatSpecifier fs) {

	topLine = 0;

	textLines.set(TextLine(text, fs));
}

void TextDisplay::NewFormattedTextLine(int lineIdx, std::string text) 
{
	textLines.insert_paragraph(lineIdx, textLines.BuildFormattedTextLine(text, lineIdx));
}

#endif