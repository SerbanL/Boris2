#pragma once

#include "CompileFlags.h"
#if GRAPHICS == 1

#include "BorisLib.h"
#include "TextObject.h"
#include "BorisGraphics.h"
#include "WinSpaces.h"



////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// IOI Popup EditBox (allows editing the text of an interactive object)

class BorisIOIPopupEditBox :
	public TextDisplay
{
	friend class BorisDisplay;

private:

	//the text to edit
	TextLine poppedParagaph;

	//the space which spawned this Popup
	INT2 parent_winId;

	//the Id of the interactive object in the parent window which spawned this popup (this is the interactive object we edit here)
	INT2 IO_Id;

	//the index (line, object) of the IO_Id object in the parent window
	INT2 spawningObject_index;

	/// data related to editing text

	//prompter position in text (index of character)
	int promptPos;

	//text selection start and end indexes of characters on entry line (inclusive). Note, it is possible selectStart > selectEnd, these are just the start and end indexes of a text selection action.
	int selectStart, selectEnd;

	//prompter rectangle user to draw it
	D2D1_RECT_F promptRect;

	//rectangle of text selection box
	D2D1_RECT_F selectRect;

public:

private:

	//calculate prompter rectangle at the current position
	void SetPrompterRect(void);
	//set edit box size depending on text length
	void SetBoxSize(void);
	//set edit box to default settings
	void SetDefault(void);

	void ResetTextSelector(void) { selectEnd = 0; selectStart = selectEnd; }

	//prompter control
	bool DecrementPrompter(void);
	bool IncrementPrompter(void);
	void SetPrompterEnd(void);
	void SetPrompterHome(void);
	//delete text depending on prompter position
	void DeleteTextFromEntryLine(void);
	//insert text at prompter position
	void NewUserEntry(std::string text);
	void InsertTextatPrompter(std::string text);

protected:

	void DrawWindow(void);
	void DrawWindow_Quick(void) { DrawWindow(); }

	//Implementation of message handling routine
	ActionOutcome NewMessage(AC_ aCode, INT2 mouse, std::string data = "");

	void SetDefaultCursor(void) { SetCursor(LoadCursor(nullptr, IDC_HAND)); }

	void AdjustWindowDimensions(D2D1_RECT_F delta) { ChangeWindowSides(delta); }

	virtual void SetText(std::string text, FormatSpecifier fs)
	{
		TextDisplay::SetText(text, fs);
		SetDefault();
	}

public:  //Public methods

	//Setup new window space with rectangle specified by ratios of the entire screen
	BorisIOIPopupEditBox(D2D1_RECT_F ratios, INT2 winId, SimTOFunct *pActionHandler);
	~BorisIOIPopupEditBox() {}

	void ResetActionHandler(void);
};

#endif
