#include "stdafx.h"
#include "BorisDisplay.h"

#if GRAPHICS == 1

////////////////////////////////////////////////////////////////////////////////////////////////////////////CONSOLE WINDOW

BorisConsole::BorisConsole(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pConsoleActionHandler) : TextDisplay(ratios, winId, pConsoleActionHandler) {

	recallIdx = 0;

	promptPos = 0;
	ResetTextSelector();

	topLine = 0;

	SetPrompterRect();

	mouseOverInteractiveObject = false;
}

void BorisConsole::SetNewLinePrompter(void) {

	//new line with blank text and default settings for user input
	textLines.push(TextLine("", FormatSpecifier(USERCOLOR, BGRNDTEXTCOLOR, false, true, TOO_NONE)));

	//starting prompter position on input line
	promptPos = 0;

	//Recalculate prompter rectangle
	SetPrompterRect();
}

void BorisConsole::SetPrompterRect(bool recalculateTopLine) {

	//if prompter position went out of range set it to the end
	if (promptPos >= textLines.para_length(LastPara())) promptPos = textLines.para_length(LastPara());

	//Need to set topLine so that the new line will not overflow the space rect
	if (recalculateTopLine) ResetTopLineIndex();

	INT2 elIdx = EntryLineIndex(promptPos);

	promptRect = textLines[elIdx.i].get_char_rect(elIdx.j);

	ShiftRect(promptRect, spaceRect.left, spaceRect.top - textLines[topLine].top());

	//finally, calculate text selection box if enabled
	if (selectStart == selectEnd) return;

	int charStart, charEnd;
	charStart = (selectStart <= selectEnd ? selectStart : selectEnd);
	charEnd = (selectStart <= selectEnd ? selectEnd : selectStart);

	INT2 select_Sidx = EntryLineIndex(charStart);
	INT2 select_Eidx = EntryLineIndex(charEnd);

	selectRect = textLines[select_Sidx.i].get_char_rect(select_Sidx.j);
	D2D1_RECT_F selectRectE = textLines[select_Eidx.i].get_char_rect(select_Eidx.j);
	selectRect.right = selectRectE.right;
	selectRect.bottom = selectRectE.bottom;

	ShiftRect(selectRect, spaceRect.left, spaceRect.top - textLines[topLine].top());
}

bool BorisConsole::DecrementPrompter(void) {

	bool changed = true;

	if (promptPos > 0) {

		promptPos--;
		SetPrompterRect();
	}
	else changed = false;

	return changed;
}

bool BorisConsole::IncrementPrompter(void) {

	bool changed = true;

	if (promptPos < textLines.para_length(LastPara())) {

		promptPos++;
		SetPrompterRect();
	}
	else changed = false;

	return changed;
}

void BorisConsole::SetPrompterEnd(void) {

	promptPos = textLines.para_length(LastPara());

	SetPrompterRect();
}

void BorisConsole::SetPrompterHome(void) {

	promptPos = 0;
	SetPrompterRect();
}

void BorisConsole::DeleteTextFromEntryLine(void) {

	if (selectStart == selectEnd) {

		//no selection, just delete character at prompter
		INT2 elIdx = EntryLineIndex(promptPos);
		textLines[elIdx.i].delchar(elIdx.j);
	}
	else {

		int charStart, charEnd;
		charStart = (selectStart <= selectEnd ? selectStart : selectEnd);
		charEnd = (selectStart <= selectEnd ? selectEnd : selectStart);

		//delete text selection
		INT2 select_SIdx = EntryLineIndex(charStart);
		INT2 select_EIdx = EntryLineIndex(charEnd);

		textLines[select_SIdx.i].deltext(select_SIdx.j, select_EIdx.j + 1 - select_SIdx.j);

		promptPos = (selectStart < selectEnd ? selectStart : selectEnd);

		//cancel text selection now
		ResetTextSelector();
	}
}

bool BorisConsole::NewCommandEntered(string &command) {

	//cancel text selection now (if any)
	ResetTextSelector();

	command = textLines.get_paragraph(LastLine());

	if (command.length()) {

		//if line not blank then save it before creating a new line
		savedCommands.push(textLines[LastLine()]);

		//recallIdx follows last element, otherwise stays put
		if (recallIdx == savedCommands.LastLine() - 1) recallIdx++;

		SetNewLinePrompter();

		return true;
	}
	else return false;
}

string BorisConsole::FormatEntryLineText(string text) {

	//USERCOLOR (Yellow) with BGRNDCOLOR (no color)
	string formattedText = "<b>[tc1,1,0,1/tc][bc0,0,0,0/bc]";

	vector<string> words = split(text, " ");

	//if last character is a separator then make an additional entry in words (blank)
	if (text[text.length() - 1] == ' ')
		words.push_back("");

	if (words.size() > 1) {

		formattedText += words[0];

		formattedText += "</b> <i>[tc1,1,0,1/tc][bc0,0,0,0/bc]";

		words.erase(words.begin());
		formattedText += combine(words, " ");
	}
	else formattedText += text;

	return formattedText;
}

void BorisConsole::SetFormattingonEntryLine(void) 
{
	int lineIndex = textLines.LastPara();

	textLines.replace_paragraph(
		lineIndex,
		textLines.BuildFormattedTextLine(FormatEntryLineText(textLines.get_paragraph(LastPara())), lineIndex));
}

void BorisConsole::DrawWindow(void) {

	//clear the drawing area
	pBG->FillRectangle(spaceRect, bgrndColor);

	DrawTextLines();
	SetPrompterRect(false);

	//Draw prompter if it's within the space window and text selection not enabled
	if (selectStart == selectEnd && promptRect.top >= spaceRect.top && promptRect.bottom <= spaceRect.bottom && promptRect.right <= spaceRect.right)
		pBG->FillRectangle(promptRect, PROMPTCOLOR);

	//Draw text selection box if it's within the space window, and is enabled
	if (selectStart != selectEnd && selectRect.top >= spaceRect.top && selectRect.bottom <= spaceRect.bottom && selectRect.right <= spaceRect.right)
		pBG->FillRectangle(selectRect, PROMPTCOLOR);

	DrawResizingFrame();
}

ActionOutcome BorisConsole::NewMessage(AC_ aCode, INT2 mouse, string data) 
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

			SetPrompterRect(false);
		}
		else actionResult.SetCode(winId, mouse, AO_NOTHANDLED);
	}
	break;

	case AC_KEYBOARD:
	{
		//This message, and others like it are only meant to be handled here. This is why windowEntered flag is not checked. No other windows will handle this message, and AO_NOTHANDLED flag was set in NewMessage_CommonResponses
		actionResult.SetCodes(winId, mouse, AO_TEXTRETURNED, AO_REFRESHWINDOW);
		NewUserEntry(data);

		//get text entered so far to check for autocompletion
		actionResult.text = textLines.get_paragraph(LastLine());
	}
	break;

	case AC_KBDLEFT:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);

		ResetTextSelector();
		DecrementPrompter();
	}
	break;

	case AC_KBDSHIFTLEFT:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);

		if (selectStart == selectEnd) {

			selectStart = promptPos;
			selectEnd = promptPos;
		}

		if (DecrementPrompter()) {

			selectEnd--;
			SetPrompterRect();
		}
	}
	break;

	case AC_KBDRIGHT:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);

		ResetTextSelector();
		IncrementPrompter();
	}
	break;

	case AC_KBDSHIFTRIGHT:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);

		if (selectStart == selectEnd) {

			selectStart = promptPos;
			selectEnd = promptPos;
		}

		if (IncrementPrompter()) {

			selectEnd++;
			SetPrompterRect();
		}
	}
	break;

	case AC_KBDDN:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);

		if (savedCommands.IsNotEmpty()) {

			if (recallIdx < savedCommands.LastLine()) recallIdx++;

			textLines.replace_paragraph(LastLine(), savedCommands[recallIdx]);

			SetPrompterEnd();
		}
	}
	break;

	case AC_KBDUP:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);

		if (savedCommands.IsNotEmpty()) {

			if (recallIdx > 0) recallIdx--;

			textLines.replace_paragraph(LastLine(), savedCommands[recallIdx]);

			SetPrompterEnd();
		}
	}
	break;

	case AC_KBDBACK:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);

		if (DecrementPrompter()) DeleteTextFromEntryLine();

		SetPrompterRect();
	}
	break;

	case AC_KBDDEL:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);

		DeleteTextFromEntryLine();

		SetPrompterRect();
	}
	break;

	case AC_KBDENTER:
	{
		if (NewCommandEntered(actionResult.text)) {
			actionResult.SetCodes(winId, mouse, AO_MESSAGERETURNED, AO_REFRESH);
		}
	}
	break;

	case AC_KBDEND:

		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);
		ResetTextSelector();
		SetPrompterEnd();
		break;

	case AC_KBDSHIFTEND:

		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);
		selectStart = promptPos;
		selectEnd = selectStart;
		SetPrompterEnd();
		selectEnd = promptPos;
		SetPrompterRect();
		break;

	case AC_KBDHOME:

		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);
		ResetTextSelector();
		SetPrompterHome();
		break;

	case AC_KBDSHIFTHOME:

		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);
		selectStart = promptPos;
		selectEnd = selectStart;
		SetPrompterHome();
		selectEnd = promptPos;
		SetPrompterRect();
		break;

	case AC_KBDPGUP:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);

		int linesUp = int(spaceRectDims.height / textLines[topLine].height());
		topLine -= linesUp;
		if (topLine < 0) topLine = 0;
	}
	break;

	case AC_KBDPGDN:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);

		int linesUp = int(spaceRectDims.height / textLines[topLine].height());
		topLine += linesUp;
		if (topLine > LastLine()) topLine = LastLine();
	}
	break;

	case AC_KBDESC:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);
		ResetTextSelector();
		SetEntryLineText("");
	}
	break;

	case AC_CTRL_C:
	{
		string entryLineText = textLines.get_paragraph(LastLine());
		if (selectStart != selectEnd) SetClipboardText(entryLineText.substr(selectStart, selectEnd + 1 - selectStart));
		else SetClipboardText(entryLineText.substr(promptPos, 1));
	}
	break;

	case AC_CTRL_V:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);
		InsertTextatPrompter(data);
	}
	break;

	case AC_MOUSEMOVE:
	{
		if (windowEntered) {

			//INT2 click = GetMouseClickObjectIdx(mouse);
			INT3 click = GetMouseClickFullIdx(mouse);

			if (click.i >= 0 && click.j >= 0) {

				//if dragging the mouse with left button down, then select text in entry line
				if (mouseLeftDown && InLastPara(click.i)) {

					selectEnd = textLines.GetParaCharIdx(INT2(click.i, click.k));

					SetPrompterRect();
					actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);
				}

				//If mouse is hovering over text object with action set then set hand cursor
				if (textLines[click.i][click.j].IsActionSet() && !mouseOverInteractiveObject) {

					SetCursor(LoadCursor(nullptr, IDC_HAND));
					mouseOverInteractiveObject = true;

					//start hover timer : if hovering long enough over interactive object a text info box may popup if required
					//AO_HOVERCHECK starts a delayed call launch, which sends AC_HOVERCHECK message after a set hover time
					//currently disabled hover info
					//hoverTimeCounter = GetTickCount();
					//actionResult.AddCode(winId, mouse, AO_HOVERCHECK);
				}

				if (!textLines[click.i][click.j].IsActionSet() && mouseOverInteractiveObject) {

					SetDefaultCursor();
					mouseOverInteractiveObject = false;
				}
			}
			else if (mouseOverInteractiveObject) { SetDefaultCursor(); mouseOverInteractiveObject = false; }
		}
	}
	break;
	
	//note: AC_HOVERCHECK is issued from AC_MOUSEMOVE (by setting result AO_HOVERCHECK, which after a delays launches AC_HOVERCHECK) but is currently disabled - info is obtained using shift-click only.
	case AC_HOVERCHECK:
	{
		if (windowEntered) {

			INT3 click = GetMouseClickFullIdx(mouse);

			if (click.i >= 0 && click.j >= 0) {

				if (textLines[click.i][click.j].IsActionSet() && mouseOverInteractiveObject) {

					//mouse is over an interactive object. If a hover is detected then attempt to show hover info.
					if (GetTickCount() - hoverTimeCounter >= HOVERTIME) {

						InteractiveObjectActionOutcome interactionOutcome = textLines[click.i][click.j].ObjectInteraction(aCode);
						if (interactionOutcome == AO_SHOWHOVERINFO) {

							actionResult.SetCode(winId, mouse, AO_SHOWHOVERINFO);
							actionResult.text = interactionOutcome.text;
						}
					}
				}
			}
		}
	}
	break;
	
	case AC_SHIFT_MOUSELEFTDOWN:
	{
		if (windowEntered) {

			INT3 click = GetMouseClickFullIdx(mouse);

			if (click.i >= 0 && click.j >= 0) {

				if (textLines[click.i][click.j].IsActionSet()) {

					//mouse is over an interactive object - attempt to show info
					InteractiveObjectActionOutcome interactionOutcome = textLines[click.i][click.j].ObjectInteraction(AC_HOVERCHECK);
					if (interactionOutcome == AO_SHOWHOVERINFO) {

						actionResult.SetCode(winId, mouse, AO_SHOWHOVERINFO);
						actionResult.text = interactionOutcome.text;
					}
				}
			}
		}
	}
	break;

	case AC_MOUSELEFTDOWN:
	{
		if (windowEntered) {

			//AddCode not SetCode as we need to keep the AO_SETTOPMOST code set by NewMessage_CommonResponses
			actionResult.AddCode(winId, mouse, AO_REFRESH);

			bool doubleClick = IsDoubleClick(aCode);
			if (doubleClick) aCode = AC_DOUBLECLICK;

			INT3 click = GetMouseClickFullIdx(mouse);

			//clicked on a text object
			if (click.i >= 0 && click.j >= 0) {

				if (textLines[click.i][click.j].IsActionSet()) {

					//it's an interactive object. Send object interaction code. If the object then asks for an inter-object interaction to start, then send this signal further (AO_STARTINTERACTION)
					if (textLines[click.i][click.j].ObjectInteraction(aCode) == AO_STARTINTERACTION)
						actionResult.AddCode(winId, mouse, AO_STARTINTERACTION);
				}
				else if (InLastPara(click.i) && click.k >= 0) {

					if (!doubleClick) {

						//clicked on entry line (the last paragraph)
						promptPos = textLines.GetParaCharIdx(INT2(click.i, click.k));

						//could be the start of a text selection
						selectStart = promptPos;
						selectEnd = selectStart;
					}
					else {

						//double clicked on text on entry line : select clicked word
						string entryLineText = textLines.get_paragraph(LastLine());

						pair<int, int> selectionIndexes = get_word_indexes(entryLineText, promptPos);
						selectStart = selectionIndexes.first;
						selectEnd = selectionIndexes.second;
					}

					SetPrompterRect();
				}
				else if (doubleClick) {

					//Double click occured
					Toggle_Maximise_Normal(mouse, actionResult);
					SetPrompterRect();
				}
			}
			else if (doubleClick) {

				//Double click occured
				Toggle_Maximise_Normal(mouse, actionResult);
				SetPrompterRect();
			}
		}
	}
	break;

	case AC_MOUSERIGHTDOWN:
	{
		if (windowEntered) {

			//AddCode not SetCode as we need to keep the AO_SETTOPMOST code set by NewMessage_CommonResponses
			actionResult.AddCode(winId, mouse, AO_REFRESH);

			//cancel text selection now (if any)
			ResetTextSelector();

			INT3 click = GetMouseClickFullIdx(mouse);

			if (click.i >= 0 && click.j >= 0) {

				if (textLines[click.i][click.j].IsActionSet()) {

					textLines[click.i][click.j].ObjectInteraction(aCode);
				}
			}
		}
	}
	break;

	case AC_DROPFILES:
	{
		if (IsMouseInWindow(mouse)) {

			actionResult.SetCode(winId, mouse, AO_FILEDROPPEDINCONSOLE);
			//name of the file in text
			actionResult.text = data;
		}
	}
	break;

	default:
		break;
	}

	return actionResult;
}

void BorisConsole::NewUserEntry(string text) {

	//if selection enabled then replace with enetered text : first delete selection, then insert text
	if (selectStart != selectEnd) DeleteTextFromEntryLine();

	INT2 elIdx = EntryLineIndex(promptPos);

	textLines[elIdx.i].insert(elIdx.j, text);
	promptPos += (int)text.length();

	//Now make sure the entry line text is in the required formatting
	SetFormattingonEntryLine();

	//might need to split entry lines into multiple lines to fit window width - DrawTextLines(false) will do this without graphically updating the screen
	DrawTextLines(false);

	SetPrompterRect();
}

void BorisConsole::NewTextLine(string text, FormatSpecifier fs) {

	//If some text has already been entered on the entry line, then insert the new message before, pushing the entry line down one
	textLines.insert_paragraph(LastLine(), TextLine(text, fs));

	//might need to split entry lines into multiple lines to fit window width - DrawTextLines(false) will do this without graphically updating the screen
	DrawTextLines(false);

	//entry line was pushed down one, so need to recalculate prompter rectangle
	SetPrompterRect();
}

void BorisConsole::InsertTextatPrompter(string text) {

	if (!text.length()) return;

	//cancel text selection now (if any)
	ResetTextSelector();

	//Stop insertion with newline characters : messes up the entry line and not needed.
	vector<string> messageLines = split(text, "\r", "\n");

	INT2 elIdx = EntryLineIndex(promptPos);

	//Insert text at prompter position using current formatting
	textLines[elIdx.i].insert(elIdx.j, messageLines[0]);

	//insert any further text lines using default formatting
	for (int i = 1; i < messageLines.size(); i++) {
		textLines.push(TextLine(messageLines[i], FormatSpecifier(USERCOLOR, BGRNDTEXTCOLOR, true, false, TOO_NONE)));
	}

	//Now make sure the entry line text is in the required formatting
	SetFormattingonEntryLine();

	//might need to split entry lines into multiple lines to fit window width - DrawTextLines(false) will do this without graphically updating the screen
	DrawTextLines(false);

	SetPrompterRect();
}

//This sets the default formatting for the entry line text - first word bold, further words italic, USERCOLOR, BGRNDTEXTCOLOR, no outline.
void BorisConsole::SetEntryLineText(string text) {

	//Stop insertion with newline characters : messes up the entry line and not needed.
	vector<string> messageLines = split(text, "\r", "\n");

	int lineIndex = textLines.LastPara();

	if (messageLines.size())
		textLines.replace_paragraph(lineIndex, textLines.BuildFormattedTextLine(FormatEntryLineText(messageLines[0]), lineIndex));
	else
		textLines.replace_paragraph(lineIndex, textLines.BuildFormattedTextLine(FormatEntryLineText(""), lineIndex));

	//might need to split entry lines into multiple lines to fit window width - DrawTextLines(false) will do this without graphically updating the screen
	DrawTextLines(false);

	SetPrompterEnd();
}

void BorisConsole::NewFormattedTextLine(string text) 
{
	//Display text on last line: existing console entry text will be pushed down
	TextDisplay::NewFormattedTextLine(LastLine(), text);

	//might need to split entry lines into multiple lines to fit window width - DrawTextLines(false) will do this without graphically updating the screen
	DrawTextLines(false);

	//entry line was pushed down one, so need to recalculate prompter rectangle
	SetPrompterRect();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////DATABOX WINDOW

BorisTextBox::BorisTextBox(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pActionHandler) : TextDisplay(ratios, winId, pActionHandler) {

}

void BorisTextBox::DrawWindow(void) {

	//clear the drawing area
	pBG->FillRectangle(spaceRect, bgrndColor);

	//Draw text lines
	DrawTextLines();

	DrawResizingFrame();
}

ActionOutcome BorisTextBox::NewMessage(AC_ aCode, INT2 mouse, string data) {

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////IOI POPUP TEXT BOX

BorisIOIPopupTextBox::BorisIOIPopupTextBox(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pActionHandler) : TextDisplay(ratios, winId, pActionHandler) {

}

void BorisIOIPopupTextBox::DrawWindow(void) {

	//clear the drawing area
	pBG->FillRectangle(spaceRect, bgrndColor);

	//Draw text lines
	DrawTextLines();
}

ActionOutcome BorisIOIPopupTextBox::NewMessage(AC_ aCode, INT2 mouse, string data) {

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

void BorisIOIPopupTextBox::ResetActionHandler(void) {

	for (int lineIdx = 0; lineIdx <= textLines.LastLine(); lineIdx++) {
		for (int objIdx = 0; objIdx <= textLines[lineIdx].LastElem(); objIdx++) {

			textLines[lineIdx][objIdx].ResetActionHandler();
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////HOVER INFO TEXT BOX

BorisHoverInfoTextBox::BorisHoverInfoTextBox(D2D1_RECT_F ratios, INT2 winId) : TextDisplay(ratios, winId, nullptr) {

}

void BorisHoverInfoTextBox::DrawWindow(void) {

	//clear the drawing area
	pBG->FillRectangle(spaceRect, bgrndColor);

	//Draw text lines
	DrawTextLines();
}

ActionOutcome BorisHoverInfoTextBox::NewMessage(AC_ aCode, INT2 mouse, string data) {

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////MESH WINDOW

BorisMeshWindow::BorisMeshWindow(D2D1_RECT_F ratios, INT2 winId) : 
	WinSpace(ratios, winId),
	ProgramStateNames(this, { VINFO(physQRep), VINFO(meshPosition), VINFO(meshName), VINFO(typeName), VINFO(mouse_mesh_info_position), VINFO(mouse_mesh_info_value), VINFO(displayMeshInfo) }, {})
{
	//font format : bold
	pBG->CreateTextFormat(FONT, MESHTEXTSIZE, DWRITE_FONT_WEIGHT_BOLD, DWRITE_FONT_STYLE_NORMAL, &pMeshWinTextFormat);
	pMeshWinTextFormat->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_LEADING);
	pMeshWinTextFormat->SetParagraphAlignment(DWRITE_PARAGRAPH_ALIGNMENT_NEAR);
	pMeshWinTextFormat->SetWordWrapping(DWRITE_WORD_WRAPPING_NO_WRAP);						//wrapping done manually, don't need in-built one
}

BorisMeshWindow::~BorisMeshWindow() 
{ 
	SafeRelease(&pMeshWinTextFormat);
}

void BorisMeshWindow::DrawWindow(void) {

	//clear the drawing area
	pBG->FillRectangle(spaceRect, bgrndColor);

	//draw coordinate system
	pBG->DrawCoordinateSystem(spaceRect);

	//draw the physical quantity set to be displayed in the mesh
	physQRep.DrawPhysQRep();

	string text = " Focused Mesh : " + physQRep.get_focused_meshName() + " - " + physQRep.get_focused_typeName();
	pBG->DrawSimpleTextLine(text, FLT2(spaceRect.left, spaceRect.top + MESHTEXTSIZE), D2D1::ColorF(D2D1::ColorF::Gray), &pMeshWinTextFormat);

	text = " Minimum  : " + physQRep.get_min_value_string();
	pBG->DrawSimpleTextLine(text, FLT2(spaceRect.left, spaceRect.top + 2 * MESHTEXTSIZE), D2D1::ColorF(D2D1::ColorF::Blue), &pMeshWinTextFormat);

	text = " Maximum  : " + physQRep.get_max_value_string();
	pBG->DrawSimpleTextLine(text, FLT2(spaceRect.left, spaceRect.top + 3 * MESHTEXTSIZE), D2D1::ColorF(D2D1::ColorF::Red), &pMeshWinTextFormat);

	if (displayMeshInfo) {

		text = " Position : " + mouse_mesh_info_position;
		pBG->DrawSimpleTextLine(text, FLT2(spaceRect.left, spaceRect.top + 4 * MESHTEXTSIZE), D2D1::ColorF(D2D1::ColorF::Green), &pMeshWinTextFormat);

		text = " Value    : " + mouse_mesh_info_value;
		pBG->DrawSimpleTextLine(text, FLT2(spaceRect.left, spaceRect.top + 5 * MESHTEXTSIZE), D2D1::ColorF(D2D1::ColorF::Green), &pMeshWinTextFormat);
	}

	//draw resizing frame if needed
	DrawResizingFrame();
}

ActionOutcome BorisMeshWindow::NewMessage(AC_ aCode, INT2 mouse, string data) {

	ActionOutcome actionResult;

	//First implement responses shared by all derived classes
	actionResult = NewMessage_CommonResponses(aCode, mouse, data);

	//Next implement particular responses of this derived class
	switch (aCode) {

	case AC_MOUSEMIDDLEDOWN:
		if (windowEntered) SetCursor(LoadCursor(nullptr, IDC_CROSS));		//cross cursor: indicates rotation/zoom mode
		break;

	case AC_MOUSEMIDDLEUP:
		if (windowEntered) SetDefaultCursor();
		break;

	case AC_MOUSELEFTDOWN:
		if (windowEntered) {

			if (IsDoubleClick(aCode)) {

				//Double click occured

				if (physQRep.GetMouseInfo(INT2(mouse.x, mouse.y), &meshName, &typeName)) {

					actionResult.text = meshName;
					actionResult.SetCode(winId, mouse, AO_MESHFOCUS2);
				}
			}
			else {

				SetCursor(LoadCursor(nullptr, IDC_HAND));		//hand cursor: indicates rotation/zoom mode
			}
		}
		break;

	case AC_MOUSELEFTUP:
		if (windowEntered) SetDefaultCursor();
		break;

	case AC_MOUSERIGHTDOWN:
		if (windowEntered) SetCursor(LoadCursor(nullptr, IDC_CROSS));		//cross cursor: indicates rotation/zoom mode
		break;

	case AC_MOUSERIGHTUP:
		if (windowEntered) SetDefaultCursor();
		break;

	case AC_MOUSEMOVE:
	{
		if (windowEntered) {

			//rotate camera view about its own axis (x mouse movement), zoom (y mouse movement)
			if (mouseRightDown) {

				pBG->AdjustCameraDistanceFromOrigin((mouse.j - dragStart.y)*CAMERAZSENS / pBG->wndHeight);
				pBG->RotateCameraView((float)(mouse.i - dragStart.x)*CAMERAROTSENS / pBG->wndWidth);
				dragStart = mouse;
				actionResult.AddCode(winId, mouse, AO_REFRESH);

				physQRep.SetAlphaBlend();
			}

			//rotate camera about origin
			else if (mouseMiddleDown) {

				pBG->RotateCameraAboutOrigin((float)(dragStart.x - mouse.i)*ROTATIONSENS / pBG->wndWidth, (float)(dragStart.y - mouse.j)*ROTATIONSENS / pBG->wndHeight);
				dragStart = mouse;
				actionResult.AddCode(winId, mouse, AO_REFRESH);

				physQRep.SetAlphaBlend();
			}

			//shift camera for mesh display view
			else if (mouseLeftDown) {

				INT2 cam_shift = (mouse - dragStart);
				pBG->Shift3DOriginPixelPosition(cam_shift.x, cam_shift.y);
				dragStart = mouse;
				actionResult.AddCode(winId, mouse, AO_REFRESH);

				physQRep.SetAlphaBlend();
			}
			else {

				//see if mouse is hovering over a mesh
				DBL3 newMeshPosition;
				string unit;
				double meshValue;
				
				if (physQRep.GetMouseInfo(INT2(mouse.x, mouse.y), &meshName, &typeName, &newMeshPosition, &meshValue, &unit)) {

					displayMeshInfo = true;
					actionResult.AddCode(winId, mouse, AO_REFRESH);

					//yes it is - is there new information?
					if (meshPosition != newMeshPosition) {

						//get position within mesh and corresponding value to be displayed
						meshPosition = newMeshPosition;
						mouse_mesh_info_position = ToString(meshPosition, "m");
						mouse_mesh_info_value = ToString(meshValue, unit) + string(" (") + meshName + " - " + typeName + string(")");
					}
				}
				else {

					if(displayMeshInfo) actionResult.AddCode(winId, mouse, AO_REFRESH);
					displayMeshInfo = false;
				}
			}
		}
		else displayMeshInfo = false;
	}
	break;

	case AC_MOUSEWHEEL:
	{
		if (windowEntered) {

			actionResult.SetCode(winId, mouse, AO_RECALCULATEMESHDISPLAY);

			int wheelDirection = ToNum(data);

			if (wheelDirection > 0) ZoomNotchUp();
			else ZoomNotchDn();
		}
	}
	break;

	case AC_DROPFILES:
	{
		if (IsMouseInWindow(mouse)) {

			actionResult.SetCode(winId, mouse, AO_FILEDROPPEDINMESH);
			//name of the file in text
			actionResult.text = data;
		}
	}
	break;

	default:
		break;
	}

	return actionResult;
}

void BorisMeshWindow::AutoSetMeshDisplaySettings(vector<PhysQ> physQ) 
{
	physQRep.CalculateRepresentation_AutoSettings(physQ, spaceRect);
}

void BorisMeshWindow::UpdatePhysQRep(vector<PhysQ> physQ) 
{
	physQRep.CalculateRepresentation(physQ);
}

bool BorisMeshWindow::SaveMeshImage(string fileName)
{
	return pBG->SaveScreenToFile(fileName, spaceRect);
}

void BorisMeshWindow::ZoomNotchUp(void) 
{
	physQRep.adjust_detail_level(DETAILNOTCHUP);
}

void BorisMeshWindow::ZoomNotchDn(void) 
{
	physQRep.adjust_detail_level(DETAILNOTCHDN);
}

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

BorisDisplay::~BorisDisplay() {

	for (int i = 0; i < bWin.size(); i++) {

		if (bWin[i]) delete bWin[i];
	}
}

WinSpace* BorisDisplay::AddWinSpace(WIN_ winIdMajor, D2D1_RECT_F ratios, SimTOFunct *pActionHandler) {

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

	case WIN_HOVERINFOTEXTBOX:
		newSpace = new BorisHoverInfoTextBox(ratios, INT2(winIdMajor, winIdMinor));
		break;

	default:
		break;
	}

	bWin[INT2(winIdMajor, winIdMinor)] = newSpace;

	return newSpace;
}

void BorisDisplay::DelWinSpace(INT2 winId) {

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

	Refresh();
}

//Make a text box which displays some fixed info. The window is deleted as soon as the mouse leaves it.
void BorisDisplay::MakeHoverInfoTextBox(INT2 mouse, string formatted_text_info_string)
{
	//there may be multiple lines in the foramtted text
	vector<string> messageLines = split(formatted_text_info_string, "\r", "\n");

	//text info to display - create it from formatted text string
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

	Refresh();
}

void BorisDisplay::PopupTextBoxInteraction(INT2 popupId, INT2 mouse) {

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
			//now call action handler with code - implementation of any interaction done there. If handler says interaction has finished then delete this popup.
			if (toRef.ObjectInteraction(AC_INTERACTOBJECTS) == AO_ENDINTERACTION)
				DelWinSpace(bpupSpace->GetwinId());
		}
	}
}

void BorisDisplay::PopupTextBoxDropInteraction(INT2 popupId, INT2 mouse) {

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

void BorisDisplay::RecalculateWindowLinks(INT2 winId) {

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

void BorisDisplay::BringWindowToFront(INT2 winId) {

	BringWindowToFront(bWin.get_index_from_id(winId));
}

void BorisDisplay::BringWindowToFront(int currPos) {

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

void BorisDisplay::Refresh(int winIdMajor) {

	if (!bWin.size()) return;

	//Need the { pBG->BeginD3DDraw(); ... pBG->EndD3DDraw(); } pair. Graphical updates should not be called from anywhere else.

	pBG->BeginD3DDraw();

	//reset tranformation
	pBG->ResetTransformation();

	if (winIdMajor < 0) {

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

ActionOutcome BorisDisplay::NewMessage(AC_ aCode, INT2 mouse, string data) {

	ActionOutcome actionResult;

	//Dispatch message to top-most window first (first in bWin vector).
	//The top-most window will decide if it will relinquish top-most position (e.g. if a mouse click occured, check if we need to make another window the top one).

	if (bWin.size()) {

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
			Refresh();
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

			case AO_DELAYEDREFRESH:
				//Wait for 25% longer than the hover time before refreshing to check for hovering : error margin
				delayed_call_launch(&BorisDisplay::Refresh_ThreadSafe, int((float)HOVERTIME*1.25));
				break;

			case AO_WINRESIZED:
				RecalculateWindowLinks(actionResult.winId);
				break;

			default:
				break;
			}
		}
	}

	return actionResult;
}

void BorisDisplay::DisplayConsoleMessage(string text) 
{
	//for external calls only
	displayMutex.lock();

	//Check for newline and carriage return separators (\n, \r)
	vector<string> messageLines = split(text, "\r", "\n");

	//add them line by line
	for (int i = 0; i < messageLines.size(); i++)
		pbConsole->NewTextLine(messageLines[i], FormatSpecifier(MESSAGECOLOR));

	Refresh(WIN_CONSOLE);

	displayMutex.unlock();
}

void BorisDisplay::DisplayConsoleError(string text)
{
	//for external calls only
	displayMutex.lock();

	//Check for newline and carriage return separators (\n, \r)
	vector<string> messageLines = split(text, "\r", "\n");

	//add them line by line
	for (int i = 0; i < messageLines.size(); i++)
		pbConsole->NewTextLine(messageLines[i], FormatSpecifier(ERRORCOLOR));

	Refresh(WIN_CONSOLE);

	displayMutex.unlock();
}

void BorisDisplay::DisplayConsoleListing(string text) 
{
	//for external calls only
	displayMutex.lock();

	//Check for newline and carriage return separators (\n, \r)
	vector<string> messageLines = split(text, "\r", "\n");

	//add them line by line
	for (int i = 0; i < messageLines.size(); i++)
		pbConsole->NewTextLine(messageLines[i], FormatSpecifier(LISTINGCOLOR));

	Refresh(WIN_CONSOLE);

	displayMutex.unlock();
}

void BorisDisplay::DisplayFormattedConsoleMessage(string text) 
{
	//for external calls only
	displayMutex.lock();

	//Check for newline and carriage return separators (\n, \r)
	vector<string> messageLines = split(text, "\r", "\n");

	//add them line by line
	for (int i = 0; i < messageLines.size(); i++)
		pbConsole->NewFormattedTextLine(messageLines[i]);

	Refresh(WIN_CONSOLE);

	displayMutex.unlock();
}

void BorisDisplay::NewDataBoxField(string formattedText) 
{
	//for external calls only
	displayMutex.lock();

	pbDataBox->NewFormattedTextLine(pbDataBox->LastLine(), formattedText);

	Refresh(WIN_DATABOX);

	displayMutex.unlock();
}

void BorisDisplay::UpdateDataBoxField(int lineIdx, string value_string) {

	//for external calls only
	displayMutex.lock();

	//the value is in the second object (it's not an interacting object). The first one is IOI_DATABOXFIELDLABEL, an interactive object.
	if(lineIdx < pbDataBox->textLines.size() && pbDataBox->textLines[lineIdx].size() >= 2)
		pbDataBox->textLines[lineIdx][1].set(value_string);

	//do not refresh here - this is typically called in a full data box update loop, so only need one refresh at the end : done by the caller

	displayMutex.unlock();
}

//adjust display for default view of physical quantity representation
void BorisDisplay::AutoSetMeshDisplaySettings(vector<PhysQ> physQ) 
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
	
	//animate start to end view - use blocking thread as animation is fairly short. If non-blocking would have to solve the complication with the mutex, not worth it in this case.
	set_blocking_thread(THREAD_TIMEDREFRESH);
	timed_call_launch<vector<PhysQ>*, PhysQRepSettings, PhysQRepSettings, double*, DWORD, DWORD>
		(&BorisDisplay::AnimateMeshViewChange, &physQ, start_settings, end_settings, &parameter, GetTickCount(), duration_ms, refresh_ms, duration_ms, THREAD_TIMEDREFRESH);
		
	//make sure the end view is actually set
	pbMeshWin->physQRep.CalculateRepresentation_NewSettings(physQ, end_settings);

	displayMutex.unlock(); 
}

void BorisDisplay::AutoSetMeshDisplaySettings_KeepOrientation(vector<PhysQ> physQ)
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

	//animate start to end view - use blocking thread as animation is fairly short. If non-blocking would have to solve the complication with the mutex, not worth it in this case.
	set_blocking_thread(THREAD_TIMEDREFRESH);
	timed_call_launch<vector<PhysQ>*, PhysQRepSettings, PhysQRepSettings, double*, DWORD, DWORD>
		(&BorisDisplay::AnimateMeshViewChange, &physQ, start_settings, end_settings, &parameter, GetTickCount(), duration_ms, refresh_ms, duration_ms, THREAD_TIMEDREFRESH);

	//make sure the end view is actually set
	pbMeshWin->physQRep.CalculateRepresentation_NewSettings(physQ, end_settings);

	displayMutex.unlock();
}

void BorisDisplay::AnimateMeshViewChange(vector<PhysQ>* pphysQ, PhysQRepSettings start, PhysQRepSettings end, double* parameter, DWORD start_time_ms, DWORD duration_ms)
{
	PhysQRepSettings newSettings = start * (1 - *parameter) + end * (*parameter);

	pbMeshWin->physQRep.CalculateRepresentation_NewSettings(*pphysQ, newSettings);
	Refresh();

	//parameter must vary between zero and 1.

	//use tanh method : start slow, middle fast, end slow - looks better than the linear method
	//alpha varies between zero and 1
	double alpha = double(GetTickCount() - start_time_ms) / duration_ms;
	//x varies between 0 and 2*TANHCUTOFF
	double x = TANHCUTOFF * 2.0 * alpha;
	//(tanh(x) + 1) / 2 varies in interval (0, 1)
	*parameter = (tanh(x) + 1.0) / 2.0;

	//Linear method : vary linearly with time (don't use this as it doesn't look so good)
	//*parameter = double(GetTickCount() - start_time_ms) / duration_ms;
}

#endif