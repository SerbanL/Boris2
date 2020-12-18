#pragma once

#include "CompileFlags.h"
#if GRAPHICS == 1

#include "BorisLib.h"
#include "TextObject.h"
#include "BorisGraphics.h"
#include "WinSpaces.h"



////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// CONSOLE (display text and allows user input)

class BorisConsole :
	public TextDisplay
{
	friend class BorisDisplay;

private:

	//save user commands
	TextLines savedCommands;

	//index of command recall (index for savedCommands).
	int recallIdx;

	//prompter position in text on entry line (index of character)
	int promptPos;

	//text selection start and end indexes of characters on entry line (inclusive). Note, it is possible selectStart > selectEnd, these are just the start and end indexes of a text selection action.
	int selectStart, selectEnd;

	//prompter rectangle user to draw it
	D2D1_RECT_F promptRect;

	//rectangle of text selection box
	D2D1_RECT_F selectRect;

	bool mouseOverInteractiveObject;

private:

	//set prompter on a new line
	void SetNewLinePrompter(void);

	//calculate prompter rectangle at the current position
	void SetPrompterRect(bool recalculateTopLine = true);

	//decrement/increment prompter position on entry line. Return true if prompter position has changed.
	bool DecrementPrompter(void);
	bool IncrementPrompter(void);
	void SetPrompterEnd(void);
	void SetPrompterHome(void);

	//delete text on entry line between selectStart and selectEnd characters, inclusive, or if text selection not available just delete at prompter
	void DeleteTextFromEntryLine(void);

	void ResetTextSelector(void) { selectEnd = 0; selectStart = selectEnd; }

	//Enter pressed, return command line text in command (if any). Return true if there was text on line, otherwise false. If true then set a new line prompter.
	bool NewCommandEntered(std::string &command);

	//Text on the entry line must always be of the form: first word bold, further words italic non-bold. Must have USERCOLOR textcolor and BGRNDTEXTCOLOR with no outline.
	std::string FormatEntryLineText(std::string text);

	//make sure formatting of entry line is correct 
	void SetFormattingonEntryLine(void);

protected:

	void DrawWindow(void);
	void DrawWindow_Quick(void);

	//Implementation of message handling routine
	ActionOutcome NewMessage(AC_ aCode, INT2 mouse, std::string data = "");

	void SetDefaultCursor(void) { SetCursor(LoadCursor(nullptr, IDC_IBEAM)); }

	void AdjustWindowDimensions(D2D1_RECT_F delta) { ChangeWindowSides(delta); }

	//New user entry at prompter position
	void NewUserEntry(std::string text);

	//very similar to NewUserEntry, but add a space first if any text exists before prompt position
	void NewSpacedUserEntry(std::string text);

	//add new simple text line (fixed formatting)
	void NewTextLine(std::string text, FormatSpecifier fs);

	//overloads the TextDisplay base method to add extra functionality : new text line with variable formatting
	void NewFormattedTextLine(std::string text);

	//insert text at prompter by adding to respective TextObject.
	void InsertTextatPrompter(std::string text);

	//set text as the entry line text (formatted type)
	void SetEntryLineText(std::string text);

public:

	//Setup new window space with rectangle specified by ratios of the entire screen
	BorisConsole(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pConsoleActionHandler);
	~BorisConsole() {}
};

#endif