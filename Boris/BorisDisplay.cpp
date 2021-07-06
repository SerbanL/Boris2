#include "stdafx.h"
#include "BorisDisplay.h"

#if GRAPHICS == 1

////////////////////////////////////////////////////////////////////////////////////////////////////////////BORIS DISPLAY

BorisDisplay::BorisDisplay(HWND hWnd, SimTOFunct *pConsoleActionHandler) : 
	GraphicalObject(hWnd, FONT, TEXTSIZE),
	ProgramStateNames(this, { VINFO_KEEPPTR(pBG), VINFO_KEEPPTR(pbMeshWin) }, {})
{
	RECT rc;
	GetClientRect(hWnd, &rc);
	RECTtoD2D1RECT(rc, screenRect);

	//default background color
	bgrndColor = D2D1::ColorF(D2D1::ColorF::Gray, 1.0f);

	//Console
	pbConsole = dynamic_cast<BorisConsole*>(AddWinSpace(WIN_CONSOLE, D2D1::RectF(0, 0, 0.85, 0.2), pConsoleActionHandler));
	pbConsole->SetBackground(D2D1::ColorF(D2D1::ColorF::Black, 1.0f));
	pbConsole->IsMaximisable(true);
	pbConsole->EnableResizing(false, true, false, true);

	//Databox
	pbDataBox = dynamic_cast<BorisTextBox*>(AddWinSpace(WIN_DATABOX, D2D1::RectF(0.85, 0, 1, 0.2), pConsoleActionHandler));
	pbDataBox->SetBackground(D2D1::ColorF(D2D1::ColorF::LightGray, 1.0f));
	pbDataBox->IsMaximisable(true);
	pbDataBox->EnableResizing(true, false, false, true);
	pbDataBox->SetTextLineWrapping(false);

	//Mesh Display
	pbMeshWin = dynamic_cast<BorisMeshWindow*>(AddWinSpace(WIN_MESHWINDOW, D2D1::RectF(0.0, 0.2, 1, 1)));
	pbMeshWin->SetBackground(D2D1::ColorF(D2D1::ColorF::White, 1.0f));
	pbMeshWin->SetStickyDisplayBottom();
	pbMeshWin->IsMaximisable(true);
	pbMeshWin->EnableResizing(false);

	LinkWindowSides_RightLeft(pbConsole->GetwinId(), pbDataBox->GetwinId());
	LinkWindowSides_BottomTop(pbDataBox->GetwinId(), pbMeshWin->GetwinId());
	LinkWindowSides_BottomTop(pbConsole->GetwinId(), pbMeshWin->GetwinId());
}

BorisDisplay::~BorisDisplay() 
{
	for (int i = 0; i < bWin.size(); i++) {

		if (bWin[i]) delete bWin[i];
	}
}

WinSpace* BorisDisplay::AddWinSpace(WIN_ winIdMajor, D2D1_RECT_F ratios, SimTOFunct *pActionHandler) 
{
	WinSpace* newSpace = nullptr;

	int winIdMinor = bWin.push_back(newSpace, winIdMajor);

	switch (winIdMajor) {

	case WIN_CONSOLE:
		newSpace = new BorisConsole(ratios, INT2(winIdMajor, winIdMinor), pActionHandler);
		break;

	case WIN_DATABOX:
		newSpace = new BorisTextBox(ratios, INT2(winIdMajor, winIdMinor), pActionHandler);
		break;

	case WIN_MESHWINDOW:
		newSpace = new BorisMeshWindow(ratios, INT2(winIdMajor, winIdMinor));
		break;

	case WIN_IOIPOPUPTEXTBOX:
		newSpace = new BorisIOIPopupTextBox(ratios, INT2(winIdMajor, winIdMinor), pActionHandler);
		break;

	case WIN_IOIPOPUPEDITBOX:
		newSpace = new BorisIOIPopupEditBox(ratios, INT2(winIdMajor, winIdMinor), pActionHandler);
		break;

	case WIN_HOVERINFOTEXTBOX:
		newSpace = new BorisHoverInfoTextBox(ratios, INT2(winIdMajor, winIdMinor));
		break;

	default:
		break;
	}

	bWin[INT2(winIdMajor, winIdMinor)] = newSpace;

	return newSpace;
}

void BorisDisplay::DelWinSpace(INT2 winId) 
{
	//delete the actual entry : delete the WinSpace object pointed to in bWin
	delete bWin[winId];
	//erase the entry in bWIn for this winId
	bWin.erase(winId);
}

void BorisDisplay::MakeIOIPopupTextBox(INT2 parent_winId, INT2 mouse) 
{
	//the parent must be derived from TextDisplay
	TextDisplay *pParent = dynamic_cast<TextDisplay*>(bWin[parent_winId]);

	INT2 click = pParent->GetMouseClickObjectIdx(mouse);

	//the last text line is the console entry line: don't use it
	if (click.i >= pParent->LastLine() || click.i < 0 || click.j < 0) return;

	//get the text object from parent window (as a text line)
	TextLine popup_textLine = pParent->GetRawTextLine_ScreenCoords(click);

	//the interactive object which called for this popup to be spawned
	TextObject& toRef = pParent->GetTextObjectRef(click);

	//just make sure it's an interactive object
	if (!toRef.IsActionSet()) return;

	D2D1_RECT_F spaceRect = popup_textLine.get_rect();
	spaceRect.bottom += 1;		//need this otherwise the line might not display

	//make the popup window
	BorisIOIPopupTextBox *pnewPopUpBox = dynamic_cast<BorisIOIPopupTextBox*>(AddWinSpace(WIN_IOIPOPUPTEXTBOX, GetNormalisedRect(spaceRect, screenRect.right, screenRect.bottom)));

	//options and settings for it
	pnewPopUpBox->SetTextLineWrapping(false);
	pnewPopUpBox->SetBackground(POPUPCOLOR);
	pnewPopUpBox->SetStickyDisplayTop();
	pnewPopUpBox->SetText(popup_textLine);
	pnewPopUpBox->spawningObject_index = click;

	//don't want interactive objects in this space to actually call handlers, they just need to look like the original popped line
	pnewPopUpBox->ResetActionHandler();
	pnewPopUpBox->parent_winId = parent_winId;
	pnewPopUpBox->IO_Id = toRef.get_iop_Id();

	BringWindowToFront(pnewPopUpBox->GetwinId());

	//the popup data box was created with a mouse left click, so transmit that to it
	pnewPopUpBox->NewMessage(AC_MOUSELEFTDOWN, mouse);

	//CHECK : RESPONSIVITY IMPROVEMENT
	//Draw();
}

//make a popup edit box. After editing the text send it back to the original interactive object to use.
void BorisDisplay::MakeIOIPopupEditBox(INT2 parent_winId, INT2 mouse)
{
	//the parent must be derived from TextDisplay
	TextDisplay *pParent = dynamic_cast<TextDisplay*>(bWin[parent_winId]);

	INT2 click = pParent->GetMouseClickObjectIdx(mouse);

	//the last text line is the console entry line: don't use it
	if (click.i >= pParent->LastLine() || click.i < 0 || click.j < 0) return;

	//get the text object from parent window (as a text line)
	TextLine popup_textLine = pParent->GetRawTextLine_ScreenCoords(click);
	std::string rawText;
	rawText << popup_textLine;

	//the interactive object which called for this popup to be spawned
	TextObject& toRef = pParent->GetTextObjectRef(click);

	//just make sure it's an interactive object
	if (!toRef.IsActionSet()) return;

	D2D1_RECT_F spaceRect = popup_textLine.get_rect();
	spaceRect.bottom += 1;		//need this otherwise the line might not display

	//make the popup window
	BorisIOIPopupEditBox *pnewPopUpBox = dynamic_cast<BorisIOIPopupEditBox*>(AddWinSpace(WIN_IOIPOPUPEDITBOX, GetNormalisedRect(spaceRect, screenRect.right, screenRect.bottom)));

	//options and settings for it
	pnewPopUpBox->SetTextLineWrapping(false);
	pnewPopUpBox->SetBackground(ERRORCOLOR);
	pnewPopUpBox->SetStickyDisplayTop();
	pnewPopUpBox->SetText(rawText, FormatSpecifier());
	pnewPopUpBox->spawningObject_index = click;

	//don't want interactive objects in this space to actually call handlers
	pnewPopUpBox->ResetActionHandler();
	pnewPopUpBox->parent_winId = parent_winId;
	pnewPopUpBox->IO_Id = toRef.get_iop_Id();

	BringWindowToFront(pnewPopUpBox->GetwinId());

	//the popup data box was created with a mouse left click, so transmit that to it
	pnewPopUpBox->NewMessage(AC_MOUSELEFTDOWN, mouse);
}

//Make a text box which displays some fixed info. The window is deleted as soon as the mouse leaves it.
void BorisDisplay::MakeHoverInfoTextBox(INT2 mouse, std::string formatted_text_info_string)
{
	//there may be multiple lines in the foramtted text
	std::vector<std::string> messageLines = split(formatted_text_info_string, "\r", "\n");

	//text info to display - create it from formatted text std::string
	TextLines tlines;
	
	for (int idx = 0; idx < messageLines.size(); idx++) {

		tlines.push(tlines.BuildFormattedTextLine(messageLines[idx], idx));
	}

	//the rectangle to make the window with : this is the rectangle that encompasses the created text lines
	D2D1_RECT_F spaceRect = tlines.get_rect();
	
	//need this otherwise the line might not display	
	spaceRect.bottom += 1;
	
	//shift rectangle to mouse position with a small displacement so the mouse will be inside the window
	spaceRect.top += mouse.j - HOVERINFOTEXTBOX_MOUSESHIFT;
	spaceRect.bottom += mouse.j - HOVERINFOTEXTBOX_MOUSESHIFT;
	spaceRect.left += mouse.i - HOVERINFOTEXTBOX_MOUSESHIFT;
	spaceRect.right += mouse.i - HOVERINFOTEXTBOX_MOUSESHIFT;

	//make hover info text box
	BorisHoverInfoTextBox *pnewHoverInfoTextBox = dynamic_cast<BorisHoverInfoTextBox*>(AddWinSpace(WIN_HOVERINFOTEXTBOX, GetNormalisedRect(spaceRect, screenRect.right, screenRect.bottom)));

	//options and settings for it
	pnewHoverInfoTextBox->SetTextLineWrapping(false);
	pnewHoverInfoTextBox->SetBackground(HOVERINFOCOLOR);
	pnewHoverInfoTextBox->SetStickyDisplayTop();
	
	//set text lines in created text box
	pnewHoverInfoTextBox->SetText(tlines);

	BringWindowToFront(pnewHoverInfoTextBox->GetwinId());

	//CHECK : RESPONSIVITY IMPROVEMENT
	//Draw();
}

void BorisDisplay::PopupTextBoxInteraction(INT2 popupId, INT2 mouse) 
{
	//the popup space
	BorisIOIPopupTextBox *bpupSpace = dynamic_cast<BorisIOIPopupTextBox*>(bWin[popupId]);
	//the window containing the mouse
	TextDisplay *pWindow = nullptr;

	//find first window which contains mouse (except for the popup evidently)
	int bWin_idx;
	for (bWin_idx = 0; bWin_idx < bWin.size(); bWin_idx++) {

		//skip the popup window
		if (bWin[bWin_idx]->GetwinId() == bpupSpace->GetwinId()) continue;

		if (bWin[bWin_idx]->IsMouseInWindow(mouse)) {

			//if bWin[bWin_idx] cannot be cast to TextDisplay* (because *bWin[bWin_idx] is not derived from TextDisplay) then nullptr is set in pWindow (don't use reinterpret_cast here !!)
			pWindow = dynamic_cast<TextDisplay*>(bWin[bWin_idx]);
			break;
		}
	}

	if (!pWindow) return;

	//check interaction with the window containing the mouse (if different from parent window of this popup)
	if (bWin[bWin_idx]->GetwinId() != bpupSpace->parent_winId) {

		//call action handler on the spawning interactive object to check for interactions with this space
		TextDisplay* pParent = dynamic_cast<TextDisplay*>(bWin[bpupSpace->parent_winId]);
		TextObject& toRef = pParent->GetTextObjectRef(bpupSpace->spawningObject_index);
		if (toRef.IsActionSet()) {

			//the interacting "object" now is actually another window space
			toRef.set_iop_interactingObjectId(bWin[bWin_idx]->GetwinId());

			//try to interact objects (implementation done in handler). If handler says interaction has finished then delete this popup.
			if (toRef.ObjectInteraction(AC_INTERACTOBJECTWITHWINDOW) == AO_ENDINTERACTION)
				DelWinSpace(bpupSpace->GetwinId());

			return;
		}
	}

	//textline and object clicked in parent window	
	INT2 click = pWindow->GetMouseClickObjectIdx(mouse);

	//check interaction with another interactive object
	if (click.i >= 0 && click.j >= 0) {

		TextObject& toRef = pWindow->GetTextObjectRef(click);

		//must have action set, but do not check it against taos
		if (toRef.IsActionSet() && !toRef.is_text_alignment_object()) {

			//mouse is over an interactive object (toRef) - send message to action handler to check if these objects should be interacted

			//set interactive object identifier which spawned this popup after being iteracted with
			toRef.set_iop_interactingObjectId(bpupSpace->IO_Id);

			TextObject& toInteracting = pWindow->GetTextObjectRef(bpupSpace->spawningObject_index);

			//now call action handler with code - implementation of any interaction done there. If handler says interaction has finished then delete this popup.
			if (toRef.ObjectInteraction(AC_INTERACTOBJECTS, &toInteracting) == AO_ENDINTERACTION)
				DelWinSpace(bpupSpace->GetwinId());
		}
	}
}

void BorisDisplay::PopupTextBoxDropInteraction(INT2 popupId, INT2 mouse) 
{
	//the popup space
	BorisIOIPopupTextBox *bpupSpace = dynamic_cast<BorisIOIPopupTextBox*>(bWin[popupId]);
	//the window containing the mouse
	TextDisplay *pWindow = nullptr;

	//find first window which contains mouse (except for the popup evidently)
	int bWin_idx;
	for (bWin_idx = 0; bWin_idx < bWin.size(); bWin_idx++) {

		//skip the popup window
		if (bWin[bWin_idx]->GetwinId() == bpupSpace->GetwinId()) continue;

		if (bWin[bWin_idx]->IsMouseInWindow(mouse)) {

			//if bWin[bWin_idx] cannot be cast to TextDisplay* (because *bWin[bWin_idx] is not derived from TextDisplay) then nullptr is set in pWindow (don't use reinterpret_cast here !!)
			pWindow = dynamic_cast<TextDisplay*>(bWin[bWin_idx]);
			break;
		}
	}

	if (!pWindow) {

		//delete this popup before returning
		DelWinSpace(bpupSpace->GetwinId());

		return;
	}

	//check interaction with the window containing the mouse (if different from parent window of this popup)
	if (bWin[bWin_idx]->GetwinId() != bpupSpace->parent_winId) {

		//call action handler on the spawning interactive object to check for interactions with this space
		TextDisplay* pParent = dynamic_cast<TextDisplay*>(bWin[bpupSpace->parent_winId]);
		TextObject& toRef = pParent->GetTextObjectRef(bpupSpace->spawningObject_index);
		if (toRef.IsActionSet()) {

			//the interacting "object" now is actually another window space
			toRef.set_iop_interactingObjectId(bWin[bWin_idx]->GetwinId());

			//try to interact objects (implementation done in handler), then delete this popup.
			toRef.ObjectInteraction(AC_DROPINTERACTOBJECTWITHWINDOW);
			DelWinSpace(bpupSpace->GetwinId());

			return;
		}
	}

	//textline and object clicked in parent window	
	INT2 click = pWindow->GetMouseClickObjectIdx(mouse);

	//check interaction with another interactive object
	if (click.i >= 0 && click.j >= 0) {

		TextObject& toRef = pWindow->GetTextObjectRef(click);

		//must have action set, but do not check it against taos
		if (toRef.IsActionSet() && !toRef.is_text_alignment_object()) {

			//mouse is over an interactive object (toRef) - send message to action handler to check if these objects should be interacted

			//set interactive object identifier which spawned this popup after being iteracted with
			toRef.set_iop_interactingObjectId(bpupSpace->IO_Id);
			//now call action handler with code - implementation of any interaction done there, then delete this popup.
			toRef.ObjectInteraction(AC_DROPINTERACTOBJECTS);
			DelWinSpace(bpupSpace->GetwinId());

			return;
		}
	}

	//delete it
	DelWinSpace(bpupSpace->GetwinId());
}

void BorisDisplay::PopupEditBoxReturnedText(INT2 popupId, INT2 mouse)
{
	//the popup space
	BorisIOIPopupEditBox *bpupSpace = dynamic_cast<BorisIOIPopupEditBox*>(bWin[popupId]);
	
	//the interactive object which generated this edit box
	TextDisplay *pParent = dynamic_cast<TextDisplay*>(bWin[bpupSpace->parent_winId]);
	TextObject& toRef = pParent->GetTextObjectRef(bpupSpace->spawningObject_index);

	//before ending the popup edit box signal to the parent interactive object a text has been returned - action can be taken in the handler, e.g. updating a value from the text in toRef.
	//Pass this popup edit text box as the object to interact with as it carries potential text to transfer to the parent window
	toRef.ObjectInteraction(AC_POPUPEDITTEXTBOXRETURNEDTEXT, &dynamic_cast<TextDisplay*>(bpupSpace)->GetTextObjectRef(INT2(0)));

	//delete it
	DelWinSpace(bpupSpace->GetwinId());
}

void BorisDisplay::RecalculateWindowLinks(INT2 winId) 
{
	D2D1_RECT_F rect = bWin[winId]->GetSpaceRect();

	for (int idx = 0; idx < (int)winLinkages_RightLeft.size(); idx++) {

		//right to left linkage
		if (winLinkages_RightLeft[idx].first == winId) {

			INT2 winId_linked = winLinkages_RightLeft[idx].second;
			D2D1_RECT_F rect_Linked = bWin[winId_linked]->GetSpaceRect();

			float delta = rect.right - rect_Linked.left;

			if (IsNZ(delta)) {

				bWin[winId_linked]->AdjustWindowDimensions(D2D1::RectF(delta, 0, 0, 0));

				RecalculateWindowLinks(winId_linked);
			}
		}

		//left to right linkage
		if (winLinkages_RightLeft[idx].second == winId) {

			INT2 winId_linked = winLinkages_RightLeft[idx].first;
			D2D1_RECT_F rect_Linked = bWin[winId_linked]->GetSpaceRect();

			float delta = rect.left - rect_Linked.right;

			if (IsNZ(delta)) {

				bWin[winId_linked]->AdjustWindowDimensions(D2D1::RectF(0, 0, delta, 0));

				RecalculateWindowLinks(winId_linked);
			}
		}
	}

	for (int idx = 0; idx < (int)winLinkages_BottomTop.size(); idx++) {

		//bottom to top linkage
		if (winLinkages_BottomTop[idx].first == winId) {

			INT2 winId_linked = winLinkages_BottomTop[idx].second;
			D2D1_RECT_F rect_Linked = bWin[winId_linked]->GetSpaceRect();

			float delta = rect.bottom - rect_Linked.top;

			if (IsNZ(delta)) {

				bWin[winId_linked]->AdjustWindowDimensions(D2D1::RectF(0, delta, 0, 0));

				RecalculateWindowLinks(winId_linked);
			}
		}

		//top to bottom linkage
		if (winLinkages_BottomTop[idx].second == winId) {

			INT2 winId_linked = winLinkages_BottomTop[idx].first;
			D2D1_RECT_F rect_Linked = bWin[winId_linked]->GetSpaceRect();

			float delta = rect.top - rect_Linked.bottom;

			if (IsNZ(delta)) {

				bWin[winId_linked]->AdjustWindowDimensions(D2D1::RectF(0, 0, 0, delta));

				RecalculateWindowLinks(winId_linked);
			}
		}
	}

	Refresh();
}

void BorisDisplay::BringWindowToFront(INT2 winId) 
{
	BringWindowToFront(bWin.get_index_from_id(winId));
}

void BorisDisplay::BringWindowToFront(int currPos) 
{
	INT2 winId = bWin[currPos]->GetwinId();

	//set the calling window as the top window, unless it is a sticky bottom window
	if (!bWin[currPos]->IsStickyDisplayBottom() && !bWin[0]->IsStickyDisplayTop()) {

		bWin.move(currPos);
	}

	//must all tell all the other windows the mouse is somewhere else. 
	//This is because of sticky bottom windows, which will be prevented from calling WindowLeft() themselves, as other higher up windows will handle the message before.
	for (int i = 0; i < bWin.size(); i++) {

		if (bWin[i]->GetwinId() != winId) bWin[i]->WindowLeft();
	}

	Refresh();
}

//refresh display, including recalculating interactive objects
void BorisDisplay::Refresh(int winIdMajor)
{
	if (!bWin.size()) return;

	//Need the { pBG->BeginD3DDraw(); ... pBG->EndD3DDraw(); } pair. Graphical updates should not be called from anywhere else.

	pBG->BeginD3DDraw();
	
	//reset tranformation
	pBG->ResetTransformation();

	if (winIdMajor == WIN_ALL) {

		//Draw all windows : call takes the form Refresh();

		//clear the drawing area
		pBG->FillRectangle(screenRect, bgrndColor);

		//Draw all windows
		for (int i = (int)bWin.size() - 1; i >= 0; i--) {

			bWin[i]->DrawWindow();
		}
	}
	else {

		//Draw all windows of this type
		for (int i = 0; i < bWin.size(); i++) {

			INT2 winId = bWin[i]->GetwinId();
			if (winId.i == winIdMajor) bWin[i]->DrawWindow();
		}
	}

	pBG->EndD3DDraw();
}

//draw display but do not recalculate interactive objects
void BorisDisplay::Draw(int winIdMajor)
{
	if (!bWin.size()) return;

	//Need the { pBG->BeginD3DDraw(); ... pBG->EndD3DDraw(); } pair. Graphical updates should not be called from anywhere else.

	pBG->BeginD3DDraw();

	//reset tranformation
	pBG->ResetTransformation();

	if (winIdMajor == WIN_ALL) {

		//Draw all windows : call takes the form Refresh();

		//clear the drawing area
		pBG->FillRectangle(screenRect, bgrndColor);

		//Draw all windows
		for (int i = (int)bWin.size() - 1; i >= 0; i--) {

			bWin[i]->DrawWindow_Quick();
		}
	}
	else {

		//Draw all windows of this type
		for (int i = 0; i < bWin.size(); i++) {

			INT2 winId = bWin[i]->GetwinId();
			if (winId.i == winIdMajor) bWin[i]->DrawWindow_Quick();
		}
	}

	pBG->EndD3DDraw();
}

ActionOutcome BorisDisplay::NewMessage(AC_ aCode, INT2 mouse, std::string data) 
{
	ActionOutcome actionResult;

	//Dispatch message to top-most window first (first in bWin vector).
	//The top-most window will decide if it will relinquish top-most position (e.g. if a mouse click occured, check if we need to make another window the top one).

	if (aCode == AC_ALLWINDOWSLEFT) {

		//The display area was left : signal this to all windows and set cursor back to arrow default.
		for (int i = 0; i < bWin.size(); i++)
			bWin[i]->WindowLeft();

		SetCursor(LoadCursor(nullptr, IDC_ARROW));
	}
	else if (bWin.size()) {

		//Starting from the top, check for the first window to handle this message
		for (int i = 0; i < bWin.size(); i++) {

			actionResult = bWin[i]->NewMessage(aCode, mouse, data);
			if (!actionResult.IsCodeSet(AO_NOTHANDLED)) break;
		}

		//nothing handled this code: mouse is not in any window. Make sure every AC_ code is handled by some window
		if (actionResult.IsCodeSet(AO_NOTHANDLED)) {

			for (int i = 0; i < bWin.size(); i++)
				bWin[i]->WindowLeft();

			SetCursor(LoadCursor(nullptr, IDC_ARROW));
			Draw();
		}

		//check for returned action code
		for (int i = 0; i < actionResult.NumCodes(); i++) {

			switch (actionResult.GetCode(i)) {

			case AO_SETTOPMOST:
			{
				//set the calling window as the top window (unless it is a sticky bottom window)
				BringWindowToFront(actionResult.winId);
			}
			break;

			case AO_STARTINTERACTION:
			{
				//start inter-object interaction: call routine to make a popup textbox for inter-object interaction
				MakeIOIPopupTextBox(actionResult.winId, mouse);
			}
			break;

			case AO_STARTPOPUPEDITBOX:
			{
				//start edit box generated by an interactive object: call routine to make a popup edit box for it
				MakeIOIPopupEditBox(actionResult.winId, mouse);
			}
			break;

			case AO_POPUPEDITBOXRETURNEDTEXT:
			{
				PopupEditBoxReturnedText(actionResult.winId, mouse);
				Refresh();
			}
			break;

			case AO_HOVERCHECK:
			{
				delayed_call_launch<INT2>(&BorisDisplay::SendHoverCheck_ThreadSafe, mouse, int((float)HOVERTIME_LONG*1.25));
			}
			break;

			case AO_SHOWHOVERINFO:
			{
				//show a pop-up window with text from actionResult.text - this window will be active as long as the mouse is within it
				MakeHoverInfoTextBox(mouse, actionResult.text);
			}
			break;

			case AO_DESTROYWINDOW:
			{
				DelWinSpace(actionResult.winId);
				Refresh();
			}
			break;

			case AO_CHECKIOINTERACTION:
				PopupTextBoxInteraction(actionResult.winId, mouse);
				Refresh();
				break;

			case AO_CHECKIODROPINTERACTION:
				//popup text box was dropped (mouse left up) - check for interaction then destroy the popup
				PopupTextBoxDropInteraction(actionResult.winId, mouse);
				Refresh();
				break;

			case AO_REFRESH:
				Refresh();
				break;

			case AO_REFRESHWINDOW:
				Refresh(actionResult.winId.major);
				break;

			case AO_DRAW:
				Draw();
				break;

			case AO_DRAWWINDOW:
				Draw(actionResult.winId.major);
				break;

			case AO_DELAYEDREFRESH:
				//Wait for 25% longer than the hover time before refreshing to check for hovering : error margin
				delayed_call_launch(&BorisDisplay::Refresh_ThreadSafe, int((float)HOVERTIME*1.25));
				break;

			case AO_WINRESIZED:
				RecalculateWindowLinks(actionResult.winId);
				break;

			case AO_ADDCONSOLECOORDINATE:
				pbConsole->NewSpacedUserEntry(actionResult.text);
				Draw();
				break;

			default:
				break;
			}
		}
	}

	return actionResult;
}

void BorisDisplay::DisplayConsoleMessage(std::string text) 
{
	//for external calls only
	displayMutex.lock();

	//Check for newline and carriage return separators (\n, \r)
	std::vector<std::string> messageLines = split(text, "\r", "\n");

	//add them line by line
	for (int i = 0; i < messageLines.size(); i++)
		pbConsole->NewTextLine(messageLines[i], FormatSpecifier(MESSAGECOLOR));

	Refresh(WIN_CONSOLE);

	displayMutex.unlock();
}

void BorisDisplay::DisplayConsoleError(std::string text)
{
	//for external calls only
	displayMutex.lock();

	//Check for newline and carriage return separators (\n, \r)
	std::vector<std::string> messageLines = split(text, "\r", "\n");

	//add them line by line
	for (int i = 0; i < messageLines.size(); i++)
		pbConsole->NewTextLine(messageLines[i], FormatSpecifier(ERRORCOLOR));

	Refresh(WIN_CONSOLE);

	displayMutex.unlock();
}

void BorisDisplay::DisplayConsoleWarning(std::string text)
{
	//for external calls only
	displayMutex.lock();

	//Check for newline and carriage return separators (\n, \r)
	std::vector<std::string> messageLines = split(text, "\r", "\n");

	//add them line by line
	for (int i = 0; i < messageLines.size(); i++)
		pbConsole->NewTextLine(messageLines[i], FormatSpecifier(WARNINGCOLOR));

	Refresh(WIN_CONSOLE);

	displayMutex.unlock();
}

void BorisDisplay::DisplayConsoleListing(std::string text) 
{
	//for external calls only
	displayMutex.lock();

	//Check for newline and carriage return separators (\n, \r)
	std::vector<std::string> messageLines = split(text, "\r", "\n");

	//add them line by line
	for (int i = 0; i < messageLines.size(); i++)
		pbConsole->NewTextLine(messageLines[i], FormatSpecifier(LISTINGCOLOR));

	Refresh(WIN_CONSOLE);

	displayMutex.unlock();
}

void BorisDisplay::DisplayFormattedConsoleMessage(std::string text) 
{
	//for external calls only
	displayMutex.lock();

	//Check for newline and carriage return separators (\n, \r)
	std::vector<std::string> messageLines = split(text, "\r", "\n");

	//add them line by line
	for (int i = 0; i < messageLines.size(); i++)
		pbConsole->NewFormattedTextLine(messageLines[i]);

	Refresh(WIN_CONSOLE);

	displayMutex.unlock();
}

void BorisDisplay::NewDataBoxField(std::string formattedText) 
{
	//for external calls only
	displayMutex.lock();

	pbDataBox->NewFormattedTextLine(pbDataBox->LastLine(), formattedText);

	Refresh(WIN_DATABOX);

	displayMutex.unlock();
}

void BorisDisplay::UpdateDataBoxField(int lineIdx, std::string value_string) 
{
	//for external calls only
	displayMutex.lock();

	//the value is in the second object (it's not an interacting object). The first one is IOI_DATABOXFIELDLABEL, an interactive object.
	if(lineIdx < pbDataBox->textLines.size() && pbDataBox->textLines[lineIdx].size() >= 2)
		pbDataBox->textLines[lineIdx][1].set(value_string);

	//do not refresh here - this is typically called in a full data box update loop, so only need one refresh at the end : done by the caller

	displayMutex.unlock();
}

//adjust display for default view of physical quantity representation
void BorisDisplay::AutoSetMeshDisplaySettings(std::vector<PhysQ> physQ) 
{ 
	displayMutex.lock(); 

	PhysQRepSettings start_settings = pbMeshWin->physQRep.get_current_settings();
	PhysQRepSettings end_settings(physQ, pBG->wndWidth, pbMeshWin->spaceRect);

	//interpolation parameter to vary from 0 to 1
	double parameter = 0;
	//the refresh rate
	double refresh_ms = ANIMATIONREFRESH_MS;
	//animation duration at the refresh rate
	double duration_ms = FOCUSCHANGEDURATION_MS;
	
	//animate start to end view - use blocking thread as animation is fairly short. If non-blocking would have to solve the complication with the std::mutex, not worth it in this case.
	set_blocking_thread(THREAD_TIMEDREFRESH);
	timed_call_launch<std::vector<PhysQ>*, PhysQRepSettings, PhysQRepSettings, double*, DWORD, DWORD>
		(&BorisDisplay::AnimateMeshViewChange, &physQ, start_settings, end_settings, &parameter, GetSystemTickCount(), duration_ms, refresh_ms, duration_ms, THREAD_TIMEDREFRESH);
		
	//make sure the end view is actually set
	pbMeshWin->physQRep.CalculateRepresentation_NewSettings(physQ, end_settings);

	displayMutex.unlock(); 
}

void BorisDisplay::AutoSetMeshDisplaySettings_KeepOrientation(std::vector<PhysQ> physQ)
{
	displayMutex.lock();

	PhysQRepSettings start_settings = pbMeshWin->physQRep.get_current_settings();
	PhysQRepSettings end_settings(physQ, pBG->wndWidth, pbMeshWin->spaceRect);

	//keep all settings apart from focusRect (set to new focus) and viewShift (set viewshift to default)
	Rect newfocusRect = end_settings.focusRect;
	float viewShiftX = end_settings.view_shiftX;
	float viewShiftY = end_settings.view_shiftY;
	float viewShiftZ = end_settings.view_shiftZ;
	
	end_settings = start_settings;
	end_settings.focusRect = newfocusRect;
	end_settings.view_shiftX = viewShiftX;
	end_settings.view_shiftY = viewShiftY;
	end_settings.view_shiftZ = viewShiftZ;

	//interpolation parameter to vary from 0 to 1
	double parameter = 0;
	//the refresh rate
	double refresh_ms = ANIMATIONREFRESH_MS;
	//animation duration at the refresh rate
	double duration_ms = FOCUSCHANGEDURATION_MS/2;

	//animate start to end view - use blocking thread as animation is fairly short. If non-blocking would have to solve the complication with the std::mutex, not worth it in this case.
	set_blocking_thread(THREAD_TIMEDREFRESH);
	timed_call_launch<std::vector<PhysQ>*, PhysQRepSettings, PhysQRepSettings, double*, DWORD, DWORD>
		(&BorisDisplay::AnimateMeshViewChange, &physQ, start_settings, end_settings, &parameter, GetSystemTickCount(), duration_ms, refresh_ms, duration_ms, THREAD_TIMEDREFRESH);

	//make sure the end view is actually set
	pbMeshWin->physQRep.CalculateRepresentation_NewSettings(physQ, end_settings);

	displayMutex.unlock();
}

void BorisDisplay::AnimateMeshViewChange(std::vector<PhysQ>* pphysQ, PhysQRepSettings start, PhysQRepSettings end, double* parameter, DWORD start_time_ms, DWORD duration_ms)
{
	PhysQRepSettings newSettings = start * (1 - *parameter) + end * (*parameter);

	pbMeshWin->physQRep.CalculateRepresentation_NewSettings(*pphysQ, newSettings);
	Refresh();

	//parameter must vary between zero and 1.

	//use tanh method : start slow, middle fast, end slow - looks better than the linear method
	//alpha varies between zero and 1
	double alpha = double(GetSystemTickCount() - start_time_ms) / duration_ms;
	//x varies between 0 and 2*TANHCUTOFF
	double x = TANHCUTOFF * 2.0 * alpha;
	//(tanh(x) + 1) / 2 varies in interval (0, 1)
	*parameter = (tanh(x) + 1.0) / 2.0;

	//Linear method : vary linearly with time (don't use this as it doesn't look so good)
	//*parameter = double(GetSystemTickCount() - start_time_ms) / duration_ms;
}

#endif