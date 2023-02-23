#include "stdafx.h"
#include "Display_BorisConsole.h"

#if GRAPHICS == 1

////////////////////////////////////////////////////////////////////////////////////////////////////////////CONSOLE WINDOW

BorisConsole::BorisConsole(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pConsoleActionHandler) :
	TextDisplay(ratios, winId, pConsoleActionHandler)
{
	recallIdx = 0;

	promptPos = 0;
	ResetTextSelector();

	topLine = 0;

	SetPrompterRect();

	mouseOverInteractiveObject = false;
}

void BorisConsole::SetNewLinePrompter(void)
{
	//new line with blank text and default settings for user input
	textLines.push(TextLine("", FormatSpecifier(USERCOLOR, BGRNDTEXTCOLOR, false, true, TOO_NONE)));

	//starting prompter position on input line
	promptPos = 0;

	//Recalculate prompter rectangle
	SetPrompterRect();
}

void BorisConsole::SetPrompterRect(bool recalculateTopLine)
{
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

bool BorisConsole::DecrementPrompter(void)
{
	bool changed = true;

	if (promptPos > 0) {

		promptPos--;
		SetPrompterRect();
	}
	else changed = false;

	return changed;
}

bool BorisConsole::IncrementPrompter(void)
{
	bool changed = true;

	if (promptPos < textLines.para_length(LastPara())) {

		promptPos++;
		SetPrompterRect();
	}
	else changed = false;

	return changed;
}

void BorisConsole::SetPrompterEnd(void)
{
	promptPos = textLines.para_length(LastPara());

	SetPrompterRect();
}

void BorisConsole::SetPrompterHome(void)
{
	promptPos = 0;
	SetPrompterRect();
}

void BorisConsole::DeleteTextFromEntryLine(void)
{
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

bool BorisConsole::NewCommandEntered(std::string &command)
{
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

std::string BorisConsole::FormatEntryLineText(std::string text)
{
	//USERCOLOR (Yellow) with BGRNDCOLOR (no color)
	std::string formattedText = "<b>[tc1,1,0,1/tc][bc0,0,0,0/bc]";

	std::vector<std::string> words = split(text, " ");

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

void BorisConsole::DrawWindow(void)
{
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

void BorisConsole::DrawWindow_Quick(void)
{
	//clear the drawing area
	pBG->FillRectangle(spaceRect, bgrndColor);

	DrawTextLines_Quick();
	SetPrompterRect(false);

	//Draw prompter if it's within the space window and text selection not enabled
	if (selectStart == selectEnd && promptRect.top >= spaceRect.top && promptRect.bottom <= spaceRect.bottom && promptRect.right <= spaceRect.right)
		pBG->FillRectangle(promptRect, PROMPTCOLOR);

	//Draw text selection box if it's within the space window, and is enabled
	if (selectStart != selectEnd && selectRect.top >= spaceRect.top && selectRect.bottom <= spaceRect.bottom && selectRect.right <= spaceRect.right)
		pBG->FillRectangle(selectRect, PROMPTCOLOR);

	DrawResizingFrame();
}

ActionOutcome BorisConsole::NewMessage(AC_ aCode, INT2 mouse, std::string data)
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
		actionResult.SetCodes(winId, mouse, AO_TEXTRETURNED, AO_DRAWWINDOW);
		NewUserEntry(data);

		//get text entered so far to check for autocompletion
		actionResult.text = textLines.get_paragraph(LastLine());
	}
	break;

	case AC_KBDLEFT:
	{
		actionResult.SetCode(winId, mouse, AO_DRAWWINDOW);

		ResetTextSelector();
		DecrementPrompter();
	}
	break;

	case AC_KBDSHIFTLEFT:
	{
		actionResult.SetCode(winId, mouse, AO_DRAWWINDOW);

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
		actionResult.SetCode(winId, mouse, AO_DRAWWINDOW);

		ResetTextSelector();
		IncrementPrompter();
	}
	break;

	case AC_KBDSHIFTRIGHT:
	{
		actionResult.SetCode(winId, mouse, AO_DRAWWINDOW);

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
		actionResult.SetCode(winId, mouse, AO_DRAWWINDOW);

		if (DecrementPrompter()) DeleteTextFromEntryLine();

		SetPrompterRect();
	}
	break;

	case AC_KBDDEL:
	{
		actionResult.SetCode(winId, mouse, AO_DRAWWINDOW);

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

		actionResult.SetCode(winId, mouse, AO_DRAWWINDOW);
		ResetTextSelector();
		SetPrompterEnd();
		break;

	case AC_KBDSHIFTEND:

		actionResult.SetCode(winId, mouse, AO_DRAWWINDOW);
		selectStart = promptPos;
		selectEnd = selectStart;
		SetPrompterEnd();
		selectEnd = promptPos;
		SetPrompterRect();
		break;

	case AC_KBDHOME:

		actionResult.SetCode(winId, mouse, AO_DRAWWINDOW);
		ResetTextSelector();
		SetPrompterHome();
		break;

	case AC_KBDSHIFTHOME:

		actionResult.SetCode(winId, mouse, AO_DRAWWINDOW);
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
		actionResult.SetCode(winId, mouse, AO_DRAWWINDOW);
		ResetTextSelector();
		SetEntryLineText("");
	}
	break;

	case AC_CTRL_C:
	{
		std::string entryLineText = textLines.get_paragraph(LastLine());
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

			INT3 click = GetMouseClickFullIdx(mouse);

			if (click.i >= 0 && click.j >= 0) {

				//if dragging the mouse with left button down, then select text in entry line
				if (mouseLeftDown && InLastPara(click.i)) {

					selectEnd = textLines.GetParaCharIdx(INT2(click.i, click.k));

					SetPrompterRect();
					actionResult.SetCode(winId, mouse, AO_DRAWWINDOW);
				}

				//If mouse is hovering over text object with action set then set hand cursor
				if (textLines[click.i][click.j].IsActionSet() && !mouseOverInteractiveObject) {

					SetCursor(LoadCursor(nullptr, IDC_HAND));
					mouseOverInteractiveObject = true;

					//start hover timer : if hovering long enough over interactive object a text info box may popup if required
					//AO_HOVERCHECK starts a delayed call launch, which sends AC_HOVERCHECK message after a set hover time
					//currently disabled hover info
					//hoverTimeCounter = GetSystemTickCount();
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
					if (GetSystemTickCount() - hoverTimeCounter >= HOVERTIME) {

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

					//it's an interactive object. Send object interaction code. If the object then asks for an inter-object interaction to start, then send this signal further (AO_STARTINTERACTION or AO_STARTPOPUPEDITBOX)
					InteractiveObjectActionOutcome ioOutcome = textLines[click.i][click.j].ObjectInteraction(aCode);
					if (ioOutcome == AO_STARTINTERACTION || ioOutcome == AO_STARTPOPUPEDITBOX) {

						actionResult.AddCode(winId, mouse, (AO_)ioOutcome.actionOutcome);
					}
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
						std::string entryLineText = textLines.get_paragraph(LastLine());

						std::pair<int, int> selectionIndexes = get_word_indexes(entryLineText, promptPos);
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

void BorisConsole::NewUserEntry(std::string text)
{
	//if selection enabled then replace with enetered text : first delete selection, then insert text
	if (selectStart != selectEnd) DeleteTextFromEntryLine();

	INT2 elIdx = EntryLineIndex(promptPos);

	textLines[elIdx.i].insert(elIdx.j, text);
	promptPos += (int)text.length();

	//Now make sure the entry line text is in the required formatting
	SetFormattingonEntryLine();

	SetPrompterRect();
}

//very similar to NewUserEntry, but add a space first if any text exists before prompt position
void BorisConsole::NewSpacedUserEntry(std::string text)
{
	//if selection enabled then replace with enetered text : first delete selection, then insert text
	if (selectStart != selectEnd) DeleteTextFromEntryLine();

	INT2 elIdx = EntryLineIndex(promptPos);

	std::string consoleText = textLines.get_paragraph(LastLine());
	if (promptPos > 0 && promptPos - 1 < consoleText.length() && consoleText[promptPos - 1] != ' ') {

		textLines[elIdx.i].insert(elIdx.j, " " + text);
		promptPos += (int)text.length() + 1;
	}
	else {

		textLines[elIdx.i].insert(elIdx.j, text);
		promptPos += (int)text.length();
	}

	//Now make sure the entry line text is in the required formatting
	SetFormattingonEntryLine();

	SetPrompterRect();
}

void BorisConsole::NewTextLine(std::string text, FormatSpecifier fs)
{
	//If some text has already been entered on the entry line, then insert the new message before, pushing the entry line down one
	textLines.insert_paragraph(LastLine(), TextLine(text, fs));

	//might need to split entry lines into multiple lines to fit window width - DrawTextLines(false) will do this without graphically updating the screen
	DrawTextLines(false);

	//entry line was pushed down one, so need to recalculate prompter rectangle
	SetPrompterRect();
}

void BorisConsole::InsertTextatPrompter(std::string text)
{
	if (!text.length()) return;

	//cancel text selection now (if any)
	ResetTextSelector();

	//Stop insertion with newline characters : messes up the entry line and not needed.
	std::vector<std::string> messageLines = split(text, "\r", "\n");

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
void BorisConsole::SetEntryLineText(std::string text)
{
	//Stop insertion with newline characters : messes up the entry line and not needed.
	std::vector<std::string> messageLines = split(text, "\r", "\n");

	int lineIndex = textLines.LastPara();
	
	if (messageLines.size())
		textLines.replace_paragraph(lineIndex, textLines.BuildFormattedTextLine(FormatEntryLineText(messageLines[0]), lineIndex));
	else
		textLines.replace_paragraph(lineIndex, textLines.BuildFormattedTextLine(FormatEntryLineText(""), lineIndex));
		
	//might need to split entry lines into multiple lines to fit window width - DrawTextLines(false) will do this without graphically updating the screen
	DrawTextLines(false);

	SetPrompterEnd();
}

void BorisConsole::NewFormattedTextLine(std::string text)
{
	//Display text on last line: existing console entry text will be pushed down
	TextDisplay::NewFormattedTextLine(LastLine(), text);

	//might need to split entry lines into multiple lines to fit window width - DrawTextLines(false) will do this without graphically updating the screen
	DrawTextLines(false);

	//entry line was pushed down one, so need to recalculate prompter rectangle
	SetPrompterRect();
}

#endif