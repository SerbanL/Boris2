#include "stdafx.h"
#include "Display_BorisIOIPopupEditBox.h"

#if GRAPHICS == 1

////////////////////////////////////////////////////////////////////////////////////////////////////////////IOI POPUP EDIT BOX

BorisIOIPopupEditBox::BorisIOIPopupEditBox(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pActionHandler) :
	TextDisplay(ratios, winId, pActionHandler)
{
	SetDefault();
}

void BorisIOIPopupEditBox::SetPrompterRect(void)
{
	//if prompter position went out of range set it to the end
	if (promptPos >= textLines.para_length(LastPara())) promptPos = textLines.para_length(LastPara());

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

//set edit box size depending on text length
void BorisIOIPopupEditBox::SetBoxSize(void)
{
	INT2 elIdx = EntryLineIndex(promptPos);
	promptRect = textLines[elIdx.i].get_char_rect(elIdx.j);

	spaceRect.right = spaceRect.left + textLines[0][0].width() + promptRect.right - promptRect.left;
}

//set edit box to default settings
void BorisIOIPopupEditBox::SetDefault(void)
{
	promptPos = textLines.para_length(LastPara());
	SetBoxSize();
	ResetTextSelector();
	SetPrompterRect();
}

bool BorisIOIPopupEditBox::DecrementPrompter(void)
{
	bool changed = true;

	if (promptPos > 0) {

		promptPos--;
		SetPrompterRect();
	}
	else changed = false;

	return changed;
}

bool BorisIOIPopupEditBox::IncrementPrompter(void)
{
	bool changed = true;

	if (promptPos < textLines.para_length(LastPara())) {

		promptPos++;
		SetPrompterRect();
	}
	else changed = false;

	return changed;
}

void BorisIOIPopupEditBox::SetPrompterEnd(void)
{
	promptPos = textLines.para_length(LastPara());

	SetPrompterRect();
}

void BorisIOIPopupEditBox::SetPrompterHome(void)
{
	promptPos = 0;
	SetPrompterRect();
}

void BorisIOIPopupEditBox::DeleteTextFromEntryLine(void)
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

	SetBoxSize();
	SetPrompterRect();
}

void BorisIOIPopupEditBox::NewUserEntry(std::string text)
{
	//if selection enabled then replace with enetered text : first delete selection, then insert text
	if (selectStart != selectEnd) DeleteTextFromEntryLine();

	INT2 elIdx = EntryLineIndex(promptPos);

	textLines[elIdx.i].insert(elIdx.j, text);
	promptPos += (int)text.length();

	SetBoxSize();
	SetPrompterRect();
}

void BorisIOIPopupEditBox::InsertTextatPrompter(std::string text)
{
	if (!text.length()) return;

	//cancel text selection now (if any)
	ResetTextSelector();

	//Stop insertion with newline characters : messes up the entry line and not needed.
	std::vector<std::string> messageLines = split(text, "\r", "\n");

	INT2 elIdx = EntryLineIndex(promptPos);

	//Insert text at prompter position using current formatting
	textLines[elIdx.i].insert(elIdx.j, messageLines[0]);

	SetBoxSize();
	SetPrompterRect();
}

void BorisIOIPopupEditBox::DrawWindow(void)
{
	//clear the drawing area
	pBG->FillRectangle(spaceRect, bgrndColor);

	//Draw text lines
	DrawTextLines();

	//Draw prompter if it's within the space window and text selection not enabled
	if (selectStart == selectEnd && promptRect.top >= spaceRect.top && promptRect.bottom <= spaceRect.bottom && promptRect.right <= spaceRect.right)
		pBG->FillRectangle(promptRect, PROMPTCOLOR);

	//Draw text selection box if it's within the space window, and is enabled
	if (selectStart != selectEnd && selectRect.top >= spaceRect.top && selectRect.bottom <= spaceRect.bottom && selectRect.right <= spaceRect.right)
		pBG->FillRectangle(selectRect, PROMPTCOLOR);

	DrawFrame();
}

ActionOutcome BorisIOIPopupEditBox::NewMessage(AC_ aCode, INT2 mouse, std::string data)
{
	ActionOutcome actionResult;

	//First implement responses shared by all derived classes
	actionResult = TextDisplay::NewMessage_CommonResponses(aCode, mouse, data);

	//Next implement particular responses of this derived class
	switch (aCode) {

	case AC_KBDESC:
	{
		actionResult.SetCode(winId, mouse, AO_DESTROYWINDOW);
	}
	break;

	case AC_KBDENTER:
	{
		actionResult.SetCode(winId, mouse, AO_POPUPEDITBOXRETURNEDTEXT);
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

	case AC_KBDBACK:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESH);

		if (DecrementPrompter()) DeleteTextFromEntryLine();

		SetPrompterRect();
	}
	break;

	case AC_KBDDEL:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESH);

		DeleteTextFromEntryLine();

		SetPrompterRect();
	}
	break;

	case AC_KBDEND:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);
		ResetTextSelector();
		SetPrompterEnd();
	}
	break;

	case AC_KBDSHIFTEND:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);
		selectStart = promptPos;
		selectEnd = selectStart;
		SetPrompterEnd();
		selectEnd = promptPos;
		SetPrompterRect();
	}
	break;

	case AC_KBDHOME:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);
		ResetTextSelector();
		SetPrompterHome();
	}
	break;

	case AC_KBDSHIFTHOME:
	{
		actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);
		selectStart = promptPos;
		selectEnd = selectStart;
		SetPrompterHome();
		selectEnd = promptPos;
		SetPrompterRect();
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

	case AC_MOUSELEFTDOWN:
	{
		if (!IsMouseInWindow(mouse))
			actionResult.SetCode(winId, mouse, AO_DESTROYWINDOW);

		if (windowEntered) {

			//AddCode not SetCode as we need to keep the AO_SETTOPMOST code set by NewMessage_CommonResponses
			//actionResult.AddCode(winId, mouse, AO_REFRESH);

			bool doubleClick = IsDoubleClick(aCode);
			if (doubleClick) aCode = AC_DOUBLECLICK;

			INT3 click = GetMouseClickFullIdx(mouse);

			//clicked on a text object
			if (click.i >= 0 && click.j >= 0) {

				if (InLastPara(click.i) && click.k >= 0) {

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
			}
		}
	}
	break;

	case AC_MOUSERIGHTDOWN:
	{
		if (!IsMouseInWindow(mouse))
			actionResult.SetCode(winId, mouse, AO_DESTROYWINDOW);
		else ResetTextSelector();
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
					actionResult.SetCode(winId, mouse, AO_REFRESHWINDOW);
				}
			}
		}
	}
	break;

	case AC_MOUSELEFTUP:
		//handle some special codes explicitly to stop window from being destroyed by default case
		break;

	default:
		//if an action next explicitly handled here is called then destroy window
		actionResult.SetCode(winId, mouse, AO_DESTROYWINDOW);
		break;
	}

	return actionResult;
}

void BorisIOIPopupEditBox::ResetActionHandler(void)
{
	for (int lineIdx = 0; lineIdx <= textLines.LastLine(); lineIdx++) {
		for (int objIdx = 0; objIdx <= textLines[lineIdx].LastElem(); objIdx++) {

			textLines[lineIdx][objIdx].ResetActionHandler();
		}
	}
}

#endif