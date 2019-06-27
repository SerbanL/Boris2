#include "stdafx.h"
#include "BorisWinSpaces.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////BORIS WIN SPACE BASE CLASS

BorisWinSpace::BorisWinSpace(D2D1_RECT_F ratios, INT2 winId) : GraphicalObject() {

	this->winId = winId;

	//Space rectangle in pixels
	ResizeWindow(ratios);

	//default background color
	bgrndColor = D2D1::ColorF(D2D1::ColorF::White, 1.0f);
	resizingFrameColor = RESIZEFRAMECOLOR;

	dragStartX = 0; dragStartY = 0;
	mouseLeftDown = false;
	mouseMiddleDown = false;
	mouseRightDown = false;

	//various flags
	stickyDisplayBottom = false; 
	stickyDisplayTop = false;
	windowEntered = false;
	allowMaximising = false;

	windowMaximised = false;
}

void BorisWinSpace::WindowLeft(void) { 
	
	windowEntered = false;
	mouseLeftDown = false;
	mouseMiddleDown = false;
	mouseRightDown = false;

	resizing.ResetResize();
}

void BorisWinSpace::WindowEntered(void) { 
	
	windowEntered = true; 
	SetDefaultCursor();
}

void BorisWinSpace::ResizeWindow(D2D1_RECT_F ratios) {

	//Space rectangle in pixels
	spaceRect.bottom = (pBG->wndHeight) * ratios.bottom;
	spaceRect.left = (pBG->wndWidth) * ratios.left;
	spaceRect.top = (pBG->wndHeight) * ratios.top;
	spaceRect.right = (pBG->wndWidth) * ratios.right;

	spaceRectDims.width = spaceRect.right - spaceRect.left;
	spaceRectDims.height = spaceRect.bottom - spaceRect.top;
}

void BorisWinSpace::ChangeWindowSides(D2D1_RECT_F delta) {
	
	if(spaceRect.left + delta.left < spaceRect.right)
		spaceRect.left += delta.left;

	if(spaceRect.right + delta.right > spaceRect.left)
		spaceRect.right += delta.right;

	if(spaceRect.bottom + delta.bottom > spaceRect.top)
		spaceRect.bottom += delta.bottom;
	
	if(spaceRect.top + delta.top < spaceRect.bottom)
		spaceRect.top += delta.top;

	spaceRectDims.width = spaceRect.right - spaceRect.left;
	spaceRectDims.height = spaceRect.bottom - spaceRect.top;
}

void BorisWinSpace::MaximiseWindowSize(void) {

	saveWinSpaceSettings = WinSpaceSettings(spaceRect, spaceRectDims, stickyDisplayBottom, stickyDisplayTop);

	ResizeWindow(D2D1::RectF(0, 0, 1, 1));

	SetStickyDisplayTop();
	SetStickyDisplayBottom(false);

	windowMaximised = true;
}

void BorisWinSpace::RestoreWindowSize(void) {

	saveWinSpaceSettings.RecallSettings(spaceRect, spaceRectDims, stickyDisplayBottom, stickyDisplayTop);

	windowMaximised = false;
}

void BorisWinSpace::DrawResizingFrame(void) {

	if(resizing.IsResizing()) {

		if(mouseLeftDown) {

			pBG->FillRectangle(D2D1::RectF(spaceRect.left, spaceRect.top, spaceRect.left + RESIZEFRAMEBORDER, spaceRect.bottom), resizingFrameColor);
			pBG->FillRectangle(D2D1::RectF(spaceRect.right - RESIZEFRAMEBORDER, spaceRect.top, spaceRect.right, spaceRect.bottom), resizingFrameColor);
			pBG->FillRectangle(D2D1::RectF(spaceRect.left, spaceRect.top, spaceRect.right, spaceRect.top + RESIZEFRAMEBORDER), resizingFrameColor);
			pBG->FillRectangle(D2D1::RectF(spaceRect.left, spaceRect.bottom - RESIZEFRAMEBORDER, spaceRect.right, spaceRect.bottom), resizingFrameColor);
		}
	}
}

void BorisWinSpace::Toggle_Maximise_Normal(INT2 mouse, ActionOutcome &actionResult) {

	if(allowMaximising) {

		if(!windowMaximised) { 
			
			MaximiseWindowSize(); 
			actionResult.AddCode(winId, mouse, AO_SETTOPMOST); 
		}
		else RestoreWindowSize();

		actionResult.AddCode(winId, mouse, AO_REFRESH);
	}
}

ActionOutcome BorisWinSpace::NewMessage_CommonResponses(AC_ aCode, INT2 mouse, char param) {

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

	switch(aCode) {

	case AC_MOUSELEFTDOWN:
		{
			if(IsMouseInWindow(mouse)) {
				
				if(!windowEntered) {
					//window just entered ... ask for this window to become the top window
					actionResult.SetCode(winId, mouse, AO_SETTOPMOST);
					//set default prompter and flag
					WindowEntered();
				}

				dragStartX = mouse.i; dragStartY = mouse.j;
				mouseLeftDown = true;
			}
			else actionResult.SetCode(winId, mouse, AO_NOTHANDLED);	//not in this window, so give other windows a chance to handle this
		}
		break;

	case AC_MOUSELEFTUP:
		{
			if(windowEntered) { 
				
				mouseLeftDown = false;
				resizing.ResetResize();
				actionResult.SetCode(winId, mouse, AO_REFRESH);
			}
			else actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		}
		break;

	case AC_MOUSERIGHTDOWN:
		{
			if(IsMouseInWindow(mouse)) { 
				
				if(!windowEntered) {
					
					actionResult.SetCode(winId, mouse, AO_SETTOPMOST);
					WindowEntered();
				}

				dragStartX = mouse.i; dragStartY = mouse.j;
				mouseRightDown = true;
			}
			else actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		}
		break;

	case AC_MOUSERIGHTUP:
		{
			if(windowEntered) { 
				
				mouseRightDown = false;
			}
			else actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		}
		break;

	case AC_MOUSEMIDDLEDOWN:
		{
			if(IsMouseInWindow(mouse)) { 
				
				if(!windowEntered) {
					
					actionResult.SetCode(winId, mouse, AO_SETTOPMOST);
					WindowEntered();
				}

				dragStartX = mouse.i; dragStartY = mouse.j;
				mouseMiddleDown = true;
			}
			else actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		}
		break;

	case AC_MOUSEMIDDLEUP:
		{
			if(windowEntered) { 
				
				mouseMiddleDown = false;
			}
			else actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		}
		break;

	case AC_MOUSEMOVE:
		{
			if(mouseLeftDown && resizing.IsResizing()) {

				ChangeWindowSides( resizing.GetResizingDeltas( mouse - INT2(dragStartX, dragStartY) ) );

				dragStartX = mouse.i;
				dragStartY = mouse.j;

				actionResult.AddCode(winId, mouse, AO_REFRESH);
				actionResult.AddCode(winId, mouse, AO_WINRESIZED);
			}
			else {

				if(IsMouseInWindow(mouse)) {
				
					//set window entered if not set already
					if(!windowEntered) {
					
						actionResult.SetCode(winId, mouse, AO_SETTOPMOST);
						WindowEntered();
					}

					//check for resizing area
					if(!mouseLeftDown) {

						//set resizing cursor as required
						if(!resizing.IsMouseInResizingArea(mouse, spaceRect)) {

							actionResult.AddCode(winId, mouse, AO_REFRESH);
							SetDefaultCursor();
						}
					}
				}
				else actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
			}
		}
		break;

	case AC_KEYBOARD:
		//This can only be handled by the console, so respond to it there. Notice windowEntered flag not checked.
		actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		break;

	case AC_KBDENTER:
		//This can only be handled by the console, so respond to it there
		actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		break;

	case AC_KBDEND:
		//This can only be handled by the console, so respond to it there
		actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		break;

	case AC_KBDHOME:
		//This can only be handled by the console, so respond to it there
		actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		break;

	case AC_KBDDEL:
		//This can only be handled by the console, so respond to it there
		actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		break;

	case AC_KBDBACK:
		//This can only be handled by the console, so respond to it there
		actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		break;

	case AC_KBDUP:
		//This can only be handled by the console, so respond to it there
		actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		break;

	case AC_KBDDN:
		//This can only be handled by the console, so respond to it there
		actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		break;

	case AC_KBDLEFT:
		//This can only be handled by the console, so respond to it there
		actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		break;

	case AC_KBDRIGHT:
		//This can only be handled by the console, so respond to it there
		actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
		break;

	default:
		break;
	}

	return actionResult;
}

bool BorisWinSpace::IsDoubleClick(AC_ &aCode) { 
	
	if(GetTickCount() - dblClickTimeCounter <= DOUBLECLICKMS) { 
		
		aCode = AC_DOUBLECLICK; 
		return true; 
	} 
	else { 
		
		dblClickTimeCounter = GetTickCount(); 
		return false; 
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////TEXT DISPLAY DERIVED ABSTRACT CLASS

BorisTextDisplay::BorisTextDisplay(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pActionHandler) : BorisWinSpace(ratios,winId) {

	topLine = 0;

	allowTextLineWrapping = true;

	textLines.set( TextLine("", FormatSpecifier(USERCOLOR, BGRNDTEXTCOLOR, false, true, TOO_NONE)) );

	//construct text formatting vector from the global const char* array defined in BorisGraphics.h
	for(int i = 0; i < TF_COUNT; i++) {

		formattingVector.push_back(formattingStrings[i]);
	}

	this->pActionHandler = pActionHandler;

	ResizeWindow(ratios);
}

BorisTextDisplay::~BorisTextDisplay() { 
	
	if(pActionHandler) delete pActionHandler;
}

ActionOutcome BorisTextDisplay::NewMessage_CommonResponses(AC_ aCode, INT2 mouse, char param) {

	ActionOutcome actionResult = BorisWinSpace::NewMessage_CommonResponses(aCode, mouse, param);

	switch(aCode) {

	case AC_MOUSEMOVE:
		{
			if(mouseLeftDown && resizing.IsResizing()) {

				ChangeWindowSides( resizing.GetResizingDeltas( mouse - INT2(dragStartX, dragStartY) ) );

				dragStartX = mouse.i;
				dragStartY = mouse.j;

				actionResult.AddCode(winId, mouse, AO_REFRESH);
			}
		}
		break;
	}

	return actionResult;
}

bool BorisTextDisplay::DrawTextLines(bool doDraw) {

	//Draw text lines in view. If liens had to be recalculated so they fit in window width, then reset top-line index
	return ( textLines.Draw(topLine, spaceRectDims.height, spaceRect.left, spaceRect.top, doDraw) );
}

void BorisTextDisplay::ResizeWindow(D2D1_RECT_F ratios) {
	
	BorisWinSpace::ResizeWindow(ratios);

	if(allowTextLineWrapping) textLines.SetWidthLimit(spaceRectDims.width);
	DrawTextLines(false);
}

void BorisTextDisplay::ChangeWindowSides(D2D1_RECT_F delta) {
	
	BorisWinSpace::ChangeWindowSides(delta);

	if(allowTextLineWrapping) textLines.SetWidthLimit(spaceRectDims.width);
	DrawTextLines(false);
}

void BorisTextDisplay::MaximiseWindowSize(void) {

	BorisWinSpace::MaximiseWindowSize();

	if(allowTextLineWrapping) textLines.SetWidthLimit(spaceRectDims.width);
	DrawTextLines(false);
}

void BorisTextDisplay::RestoreWindowSize(void) {

	BorisWinSpace::RestoreWindowSize();
	
	if(allowTextLineWrapping) textLines.SetWidthLimit(spaceRectDims.width);
	DrawTextLines(false);
}

void BorisTextDisplay::Toggle_Maximise_Normal(INT2 mouse, ActionOutcome &actionResult) {

	if(allowMaximising) {

		if(!windowMaximised) { 
			
			MaximiseWindowSize(); 
			actionResult.AddCode(winId, mouse, AO_SETTOPMOST); 
		}
		else RestoreWindowSize();

		actionResult.AddCode(winId, mouse, AO_REFRESH);
	}
}

bool BorisTextDisplay::LineOverflowsSpace(int lineIdx, float *pscreenYOffset) {

	float screenYOffset = spaceRect.top;
	if(pscreenYOffset)
		*pscreenYOffset = screenYOffset;

	if(lineIdx > textLines.LastLine()) return true;

	//find Y coordinate offset which must be used when drawing text objects: TextObject coordinates are relative to the space rect, and the topLine line must always display at the top of the space rect
	screenYOffset -= textLines[topLine].top();
	if(pscreenYOffset)
		*pscreenYOffset = screenYOffset;

	//check if this line's actual screen coordinates will overflow the space rect
	if(textLines[lineIdx].bottom() + screenYOffset > spaceRect.bottom) return true;

	return false;
}

void BorisTextDisplay::ResetTopLineIndex(void) {

	//Set topLine so that the new line will not overflow the space rect

	int lastLine = LastLine();

	if(LineOverflowsSpace(lastLine)) {
		
		for(topLine = lastLine; topLine >= 0; topLine--) {

			if(LineOverflowsSpace(lastLine)) {

				topLine++;
				break;
			}
		}

		//NOTE: the possibility of topLine = -1 should never occur normally. It could happen if the space window is so small that a single line cannot be displayed (no text can be displayed then)
		if(topLine < 0) topLine = lastLine;
	}
}

int BorisTextDisplay::GetMouseClickLineIdx(INT2 mouse) {

	INT2 click = textLines.GetTextObjectIdx(mouse - INT2((int)spaceRect.left, (int)spaceRect.top), topLine);

	if(click.j >= 0) return click.i;
	else return -1;
}

TextLine BorisTextDisplay::GetRawTextLine_ScreenCoords(int lineIdx) { 

	//get text line from lineIdx with coordinates as it appears on screen

	TextLine textLine; 
	textLine = textLines[lineIdx]; 

	textLine.set_top_left(spaceRect.left, spaceRect.top + textLines[lineIdx].top() - textLines[topLine].top());

	return textLine; 
}

void BorisTextDisplay::SetText(string text, FormatSpecifier fs) {

	topLine = 0;

	textLines.set( TextLine(text, fs) );
}

//Display new line using text formatting: 
	//<b> </b> : bold
	//<i> </i> : italic
	//[on], [os], [or] : text outline (TOO_NONE, TOO_SQUARE, TOO_ROUND)
	//[tcr,g,b,a/tc] : text color specified using r,g,b,a
	//[bcr,g,b,a/bc] : background color
	//[io0,0,textId/io] : interactive object settings: first two numbers identify the major and minor TextObject type (used to identify the type of action generated when interacting with this TextObject). Last field contains a text identifier.
	//</io> use this to specify interactive object end
void BorisTextDisplay::NewFormattedTextLine(int lineIdx, string text) {

	TextLine newLine;

	//Default settings : modified by formatting
	FormatSpecifier fs;

	SimTOFunct *pActionHandler_ = NULL;

	string remainingText;

	while(true) {

		//search for the first occurence of a format specifier
		int formCode = SplitOnFormatSpecifier(text, remainingText, formattingVector);

		if(text.length()) newLine.push( TextObject(text, fs, pActionHandler_) );

		//If no remaining text, then there's nothing left to do
		if(!remainingText.length()) break;

		//change formatting depending on the format specifier found. There may be multiple specifiers one after another.
		while(true) {

			//make changes to formatting if any required
			switch(formCode) {

			case TF_BOLDSTART:
				fs.bold = true;
				break;

			case TF_BOLDEND:
				fs.bold = false;
				break;

			case TF_ITALICSTART:
				fs.italic = true;
				break;

			case TF_ITALICEND:
				fs.italic = false;
				break;

			case TF_ONONE:
				fs.textOutline = TOO_NONE;
				break;

			case TF_OSQUARE:
				fs.textOutline = TOO_SQUARE;
				break;

			case TF_OROUND:
				fs.textOutline = TOO_ROUND;
				break;

			case TF_TCSTART:
				{
					string save_text = remainingText;	//save remaining text in case the passed text does not have the correct formatting sequence
					text = remainingText;
					formCode = SplitOnFormatSpecifier(text, remainingText, formattingVector);
					if(formCode == TF_TCEND) {

						vector<string> rgba = split(text, ',');
						if(rgba.size() == 4) {
							D3DCOLORVALUE newColor = { ToNum<FLOAT>(rgba[0]),ToNum<FLOAT>(rgba[1]),ToNum<FLOAT>(rgba[2]),ToNum<FLOAT>(rgba[3]) };
							fs.textColor = newColor;
						}
						else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
					}
					else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
				}
				break;

			case TF_BCSTART:
				{
					string save_text = remainingText;	//save remaining text in case the passed text does not have the correct formatting sequence
					text = remainingText;
					formCode = SplitOnFormatSpecifier(text, remainingText, formattingVector);
					if(formCode == TF_BCEND) {

						vector<string> rgba = split(text, ',');
						if(rgba.size() == 4) {
							D3DCOLORVALUE newColor = { ToNum<FLOAT>(rgba[0]),ToNum<FLOAT>(rgba[1]),ToNum<FLOAT>(rgba[2]),ToNum<FLOAT>(rgba[3]) };
							fs.bgrndColor = newColor;
						}
						else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
					}
					else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
				}
				break;

			case TF_SETACTIONSPECSTART:
				{
					string save_text = remainingText;	//save remaining text in case the passed text does not have the correct formatting sequence
					text = remainingText;
					formCode = SplitOnFormatSpecifier(text, remainingText, formattingVector);
					if(formCode == TF_SETACTIONSPECEND) {

						vector<string> ids = split(text, ',');
						if(ids.size() >= 3) {
						
							pActionHandler_ = pActionHandler;
							if(pActionHandler_) {

								pActionHandler_->set_iop( InteractiveObjectProperties( INT2(ToNum<int>(ids[0]), ToNum<int>(ids[1])), ids[2]) );
							}
						}
						else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
					}
					else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
				}
				break;

			case TF_SETACTIONEND:
				pActionHandler_ = NULL;
				break;

			default:
				break;
			}

			//check to see if other formatting specifiers follow immediately.
			text = remainingText;
			string remainingText_;
			formCode = SplitOnFormatSpecifier(text, remainingText_, formattingVector);

			//if text contains something then no immediate specifiers found, otherwise take another iteration to process formCode and check for further specifiers.
			if(text.length() > 0) break;
			else remainingText = remainingText_;
		}

		//more text to process
		text = remainingText;
	}

	textLines.insert(lineIdx, newLine);
}