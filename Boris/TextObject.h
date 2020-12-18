#pragma once

#include "CompileFlags.h"
#if GRAPHICS == 1

#include <string>
#include <vector>
#include <memory>
#include <utility>

#include "BorisGraphics.h"
#include "BorisLib.h"
#include "InteractiveObject.h"

#include "TextFormatting.h"

//maximum characters a text object can display (text it actually holds is not limited, but we do want to limit the text displayed on screen).
#define MAX_TEXTOBJECT_LENGTH	50



////////////////////////////////////////////////////////////////////////////////////////////////////////////TEXT OBJECT

class TextObject : 
	public GraphicalObject 
{
	//TextObject should not be used directly from outside, but rather through TextLine and TextLines
	friend class TextLine;
	friend class TextLines;

	//tao pointer pair : the TextObject* itself and TextLine* for the TextObject TextLine
	typedef std::pair<TextLine*, TextObject*> tao_ptr;
	//used for synchronised taos shared vector
	typedef std::shared_ptr< std::vector<tao_ptr> > tao_vec_ptr;

private:

	//This is the full text stored in this object, which may be longer than the displayed text (only interactive objects which are not text alignment objects are limited in display length)
	std::string fulltext;

	//number of characters to display from fulltext : limit to MAX_TEXTOBJECT_LENGTH
	int displayTextLength;

	//text box for drawing
	D2D1_RECT_F toRect;

	//formatting options
	FormatSpecifier fs;

	//Action handler functionoid: called when user interacts with this object (if set), e.g. through a mouse click.
	//A shared_ptr is used to allow for easy garbage collection. When an interactive object is made a new pActionHandler is made. Since it's a shared_ptr we don't have to worry about deleting it in the destructor. This is useful because in the meantime this object could
	//be copied over to another text object, in which case the reference count of pActionHandler increases. With a normal pointer simply deleting pActionHandler in the destructor would leave the other text object pActionHandler dangling.
	std::shared_ptr<SimTOFunct> pActionHandler = nullptr;

	//shared pointer to vector of pointers to text alignment objects (tao), including this one. If empty this is not a text alignment object. This shared_ptr is shared between the synchronised taos. 
	//If one gets deleted then it removes its entry in the vector, thus this information is immediately available to all other synchronised taos. The first tao makes the shared_ptr, the others which synchronise on creation just add themselves to the vector.
	tao_vec_ptr aligned_objects = nullptr;

private:

	//-------------------Draw

	//draw the text object, adding coordinates shifts if specified. Also limit the text rectangle in the DrawTextLine method to given widthLimit if specified.
	void Draw(float dX = 0, float dY = 0, float widthLimit = 0);

	//-------------------Auxiliary methods

	//does the formatted text overflow the given right coordinate?
	bool TextOverflows(float right) { return (toRect.right > right); }

	//split the text of the current object using the separator if its right rectangle coordinate overflows the given value.
	//the nonOverflowingTO object will have same formatting but its text is guaranteed not to overflow. Return true if an overflow was detected.
	bool SplitOnOverflow(TextObject &nonOverflowingTO, float right, char separator);

	//from this object extract the first word (identified using a space separator) to the new object (keep same formatting), including the space separator. If no space found then shift the whole text.
	void ShiftFirstWord(TextObject &oneWordTextObject);

	//from the set full text set the text std::string to display
	void set_display_text_length(void)
	{
		//only limit interactive objects which are not text alignment objects
		if (pActionHandler && !aligned_objects) {

			if (fulltext.length() >= MAX_TEXTOBJECT_LENGTH) displayTextLength = MAX_TEXTOBJECT_LENGTH;
			else displayTextLength = fulltext.length();
		}
		else displayTextLength = fulltext.length();
	}

	//-------------------Coordinates and rectangles set / get methods

	//calculate rectangle for current top-left coordinates (for display purposes, so use text)
	void CalculateRectangle(void) 
	{ 
		float width = (float)pBG->GetMonospacedFontPixelsWidth() * displayTextLength;
		float height = (float)pBG->GetFontPixelsHeight(); 
		toRect = D2D1::RectF(toRect.left, toRect.top, toRect.left + width, toRect.top + height); 
	}

	//calculate rectangle for given top-left coordinates
	void CalculateRectangle(float X, float Y) 
	{ 
		float width = (float)pBG->GetMonospacedFontPixelsWidth() * displayTextLength;
		float height = (float)pBG->GetFontPixelsHeight(); 
		toRect = D2D1::RectF(X, Y, X + width, Y + height); 
	}

	//set top-left coordinates but do not calculate rectangle
	void set_top_left(float X, float Y) { toRect = D2D1::RectF(X, Y, X + width(), Y + height()); }

	void bottom_left(float &X, float &Y) { X = toRect.left; Y = toRect.bottom; }
	void bottom_right(float &X, float &Y) { X = toRect.right; Y = toRect.bottom; }
	void top_right(float &X, float &Y) { X = toRect.right; Y = toRect.top; }
	void top_left(float &X, float &Y) { X = toRect.left; Y = toRect.top; }

	//get the rectangle of substring with given starting index and length
	D2D1_RECT_F substr_rect(int charIdx, int length = 1);
	//as above but get just the right coordinate
	float substr_right(int charIdx, int length = 1);

	//get character index from mouse coordinates (-1 if none). Note, mouse coordinates must be relative the to the top-left coordinates of window in which this TextObject appears
	int get_char_index(INT2 relMouse);

	//-------------------Various properties getters and setters

	int displayedtext_length(void) { return displayTextLength; }
	bool IsBlank(void) { return (fulltext.compare("") == 0); }

	void SetFormatSpecifier(FormatSpecifier fs) { this->fs = fs; set_top_left(toRect.left, toRect.top); }

	//get tao number
	int get_text_alignment_object_number(void) { return pActionHandler->iop.majorId; }
	//add new synchronising tao in shared vector
	void synch_text_alignment_object(tao_ptr addThis) { aligned_objects->push_back(addThis); }

public:  //public methods

	//Default constructor
	TextObject();

	//constructor : no action handler set
	TextObject(std::string text, FormatSpecifier fs);

	//full constructor
	TextObject(std::string text, FormatSpecifier fs, SimTOFunct *pActionHandler_, float X = 0, float Y = 0);

	//full constructor for text alignment object starter. *ptao_line is the TextLine on which this TextObject being created will be located
	TextObject(std::string text, FormatSpecifier fs, SimTOFunct *pActionHandler_, TextLine* ptao_line, float X = 0, float Y = 0);

	//full constructor for text alignment object with synchronisation (synchronise to existing tao: synch_tao). *ptao_line is the TextLine on which this TextObject being created will be located
	TextObject(std::string text, FormatSpecifier fs, SimTOFunct *pActionHandler_, TextLine* ptao_line, TextObject* psynch_tao, float X = 0, float Y = 0);

	//destructor - must manage aligned_objects
	~TextObject();

	//copy constructor for lvalues
	TextObject(const TextObject& copyThis) { *this = copyThis; }

	TextObject& operator=(const TextObject& copyThis);

	//-------------------Text modification methods

	//Add to TextObject text
	template <typename VType> void operator+=(const VType &rhs) 
	{ 
		fulltext += ToString(rhs); 
		set_display_text_length();
		CalculateRectangle(); 
	}

	//use this to insert text in TextObject at the index position
	template <typename VType> void insert(int charIdx, VType &rhs) 
	{ 
		if (!GoodIdx((int)fulltext.length(), charIdx)) return;
		fulltext = fulltext.substr(0, charIdx) + ToString(rhs) + fulltext.substr(charIdx);
		set_display_text_length();
		CalculateRectangle(); 
	}

	//Set TextObject text
	template <typename VType> void set(VType &rhs) 
	{ 
		fulltext = ToString(rhs);
		set_display_text_length();
		CalculateRectangle(); 
	}

	//delete character from end of std::string or from the given index (if index is valid). Return true if deletion results in empty std::string.
	bool delchar(int delIdx);

	//-------------------Text formatting modification methods

	void SetBackgroundColor(D2D1_COLOR_F bgrndColor) { fs.bgrndColor = bgrndColor; }

	//-------------------Text get methods

	//Get textObject text into a std::string
	friend std::string& operator<<(std::string &lhs, const TextObject &rhs) { lhs = rhs.fulltext; return lhs; }

	//as above but using a method
	std::string GetText(void) { return fulltext; }

	//Get textObject text into a std::string by adding to it
	friend std::string& operator+=(std::string &lhs, const TextObject &rhs) { lhs += rhs.fulltext; return lhs; }

	FormatSpecifier GetFormatSpecifier(void) { return fs; }

	//-------------------Comparisons

	//complete comparison : text and formatting for non-interactive objects. Interactive objects are never considered identical.
	bool operator==(const TextObject &rhs) { return (fulltext.compare(rhs.fulltext) == 0 && fs == rhs.fs && pActionHandler == nullptr && rhs.pActionHandler == nullptr); }

	//same type of object ? i.e. same formatting for non-interactive objects. Interactive objects are always considered to be different.
	bool SameType(const TextObject &rhs) { return (fs == rhs.fs && pActionHandler == nullptr && rhs.pActionHandler == nullptr); }

	//-------------------Coordinates and rectangles public get methods

	float left(void) { return toRect.left; }
	float right(void) { return toRect.right; }
	float bottom(void) { return toRect.bottom; }
	float top(void) { return toRect.top; }
	float height(void) { return (toRect.bottom - toRect.top); }
	float width(void) { return (toRect.right - toRect.left); }

	//-------------------Interactive object methods

	//this object has been interacted with as specified in the action code (e.g. a moust click. Note here we don't need to know what actionCode is, just pass it to the handler which will know what to do)
	InteractiveObjectActionOutcome ObjectInteraction(int actionCode) { if (pActionHandler && !aligned_objects) return pActionHandler->Action(actionCode, this); else return (AO_)0; }
	
	//same as above, but when the Action handler is called, instead of passing this TextObject, pass pInteractingObject
	InteractiveObjectActionOutcome ObjectInteraction(int actionCode, TextObject* pInteractingObject) { if (pActionHandler && !aligned_objects) return pActionHandler->Action(actionCode, pInteractingObject); else return (AO_)0; }

	//check state of this object using the state handler (if set, i.e. an interactive object).
	InteractiveObjectStateChange CheckStateChange(void);

	//Check if this object is an interactive object type
	bool IsActionSet(void) { return (pActionHandler != nullptr); }

	//check for tao
	bool is_text_alignment_object(void) { return (aligned_objects != nullptr); }

	void ResetActionHandler(void) { pActionHandler = nullptr; }

	//get interactive object identifier from iop (InteractiveObjectProperties object) stored in pActionHandler
	INT2 get_iop_Id(void) { if (pActionHandler) return INT2(pActionHandler->iop.majorId, pActionHandler->iop.minorId); else return INT2(-1); }
	void set_iop_interactingObjectId(INT2 interactingObjectId) { if (pActionHandler) pActionHandler->iop.interactingObjectId = interactingObjectId; }

	//get/set iop state
	IOS_ get_iop_state(void) { if (pActionHandler) return pActionHandler->iop.state; else return IOS_ON; }
	void set_iop_state(IOS_ state) { if (pActionHandler) pActionHandler->iop.state = state; }

	//get pointer InteractiveObjectProperties
	InteractiveObjectProperties* get_iop_pointer(void) { if (pActionHandler) return &(pActionHandler->iop); else return nullptr; }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////TEXT LINE
//Collection of contiguous TextObjects forming a single text line

//mark split lines
enum SPLIT_ { SPLIT_NO = 0, SPLIT_START, SPLIT_MIDDLE, SPLIT_END };

class TextLine {

	friend class TextLines;

private:

	std::vector<TextObject> textLineVEC;

	SPLIT_ splitLine;

public:

private:

	//-------------------Draw

	//draw the text line, adding coordinates shifts if specified. Limit drawing up to given width limit (if specified)
	void Draw(float dX = 0, float dY = 0, float widthLimit = 0);

	//-------------------Auxiliary methods

	//convert index in textLine overall text into an INT2 index : TextObject index and further index inside text of respective TextObject
	INT2 IndextoINT2Index(int idx);

	//Recalculate rectangles of TextObjects on textLine from given index up (all if invalid index)
	void RecalculateRectangles(int idx = -1);

	//check for tao_line pointers in tao shared vector and replace them with pointer to this. Explanation:
	//Mostly we just use assignment operator line1 = line2, where the intention is line2 will go out of scope before line1. In this case the pointers to line2 are replaced by pointers to line1 by the assignment operator and everything works.
	//In some special cases however we need to build line1 out of line2 by pushing text objects one by one. In this case the line2 pointer never gets replaced and results in a dangling pointer. After line1 has been built use this method on line1 with tao_line being line2.
	//FixTaoLine : replace all tao_line occurences with this
	void FixTaoLine(TextLine* tao_line);

	//-------------------Text modification methods

	//clear back to a blank textline
	void clear() { textLineVEC.resize(0); textLineVEC.push_back(TextObject("", FormatSpecifier())); }

	//-------------------Coordinates and rectangles set / get methods

	//get index of clicked character from mouse coordinates (-1 if none). Note, mouse coordinates must be relative the to the top-left coordinates of window in which this TextObject appears
	int get_char_index(INT2 relMouse);
	//get index of clicked text object on this line
	int get_textobject_index(INT2 relMouse);
	//get both - combine the above methods in one
	INT2 get_char_and_textobject_index(INT2 relMouse);

	//fix the rectangle of the last TextObject on this line (i.e. make it snap to the text) - used on lines which have been split, since their last rectangle is elongated up to the window width.
	TextLine& FixSplitLineRectangle(void);

	//elongate rectangle of last TextObject on this line up to the widthLimit value - used for split lines
	TextLine& SetSplitLineRectangle(float widthLimit);

public:

	TextLine(void);

	TextLine(std::string text, FormatSpecifier fs);

	~TextLine(void);

	//copy constructor for lvalues
	TextLine(const TextLine& copyThis) { *this = copyThis; }

	//assignment operator for lvalues
	TextLine& operator=(const TextLine& copyThis);

	//-------------------Indexing

	//return TextObject reference at given index
	TextObject& operator[](const int &objIdx) { if (GoodIdx(LastElem(), objIdx)) return textLineVEC[objIdx]; else return textLineVEC[LastElem()]; }

	//-------------------Text modification methods

	//Add text to last object on textLine
	template <typename VType> void operator+=(const VType &rhs) { textLineVEC[LastElem()] += rhs; }

	//add two text lines into one (e.g. textLine += addThisLine;)
	friend TextLine& operator+=(TextLine &lhs, TextLine &rhs);

	//Add a TextObject at the end of this line, replacing the last TextObject if blank
	void push(const TextObject &newTextObject);

	//use this to insert text in TextLine at the character index position. Will use formatting settings for the Text Object at that position.
	template <typename VType> void insert(int charIdx, VType &insertText);

	//set text of a given text object on this line
	template <typename VType> void set_textobject(int toIdx, VType &setText) { if (GoodIdx(LastElem(), toIdx)) { textLineVEC[toIdx].set(setText); RecalculateRectangles(toIdx + 1); } }

	//delete character from end of std::string or from the given index (if index is valid). Return true if deletion results in empty textLine std::string.
	bool delchar(int delIdx);

	//delete text from start character to given length in number of characters
	bool deltext(int startIdx, int numChars = -1);

	//delete textobject from end of textLine, or from the given index
	void delobject(int delIdx = -1);

	//-------------------Text get methods

	//Get text of TextLine object into a std::string
	friend std::string& operator<<(std::string &lhs, const TextLine &rhs);

	//Get textObject text into a std::string by adding to it
	friend std::string& operator+=(std::string &lhs, const TextLine &rhs) { std::string text; text << rhs; lhs += text; return lhs; }

	//-------------------Coordinates and rectangles set / get methods

	//set starting top-left coordinates of this line
	void set_top_left(float X, float Y);

	//get various coordinates of this line
	float bottom(void) { return textLineVEC[0].bottom(); }
	float top(void) { return textLineVEC[0].top(); }
	float right(void) { return textLineVEC[LastElem()].right(); }
	float left(void) { return textLineVEC[0].left(); }
	float height(void) { return (textLineVEC[0].bottom() - textLineVEC[0].top()); }
	void bottom_left(float &X, float &Y) { textLineVEC[0].bottom_left(X, Y); }
	void top_left(float &X, float &Y) { textLineVEC[0].top_left(X, Y); }

	//get rectangle of entire text line
	D2D1_RECT_F get_rect(void) { float left, right, top, bottom; textLineVEC[0].top_left(left, top); textLineVEC[LastElem()].bottom_right(right, bottom); return D2D1::RectF(left, top, right, bottom); }

	//return rectangle of character at charIdx
	D2D1_RECT_F get_char_rect(int charIdx);

	//-------------------Various properties getters and setters

	//search for text object on this line and return its starting character index
	int CharIndex_of_TextObject(TextObject* pTO);

	//get text length
	int length(void);

	//get size of textLineVEC
	int size(void) { return (int)textLineVEC.size(); }

	bool IsBlank(void) { return (LastElem() == 0 && textLineVEC[0].IsBlank()); }

	int LastElem(void) { return ((int)textLineVEC.size() - 1); }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////TEXT LINES
//Collection of TextLine objects arranged contiguously along Y direction

class TextLines {

private:

	//Action handler functionoid: pass this to interactive text objects. When the text objects are interacted with, they will use this functionoid to call the handler routine, using their id and message parameter
	SimTOFunct * pActionHandler;

	//main text storage - line by line
	std::vector<TextLine> textLinesVEC;

	//if set then text is wrapped to this limit
	float widthLimit;

public:

private:

	//Recalculate rectangles of TextLine objects from given index up (all if invalid index)
	void RecalculateRectangles(int idx = -1);

	//Check the text line at given index, object by object to see if the state of any interactive object has changed (each in turn calls its State handler). Return true if something changed.
	//Also implement functionality for special state settings - see InteractiveObject.h
	bool CheckLine_InteractiveObjectState(int lineIdx);

	//split given line into multiple lines fitting the widthLimit value above. Do not split interactive objects or whole words.
	std::vector<TextLine> SplitLinetoFitWidth(TextLine& textLine);

	//spplit line at lineIdx into multiple lines which fit to widthLimit
	void FitLinetoWidth(int &lineIdx);

	bool LineNeedsRecalculation(int lineIdx);

	//from given line index, return an index set to the start of a split line (SPLIT_START) if applicable (i.e. if lineIdx line is part of a split line sequence)
	int GetSplitStartIndex(int lineIdx);
	int GetSplitEndIndex(int lineIdx);

	//-------------------Text modification methods

	//delete line at given index, allowing deletion of a line in a split line sequence
	void delline(int lineIdx);

public:

	TextLines(void);
	TextLines(SimTOFunct *pActionHandler);
	~TextLines() { if (pActionHandler) delete pActionHandler; }

	//should be called just after TextLines object is declared - this is held in a TextDisplay window and should use its window width value
	void SetWidthLimit(float widthLimit) { this->widthLimit = widthLimit; }

	//-------------------Draw

	//draw textLines starting from topLine index, drawing all lines up to the given height. The lines are drawn relative to the topLine position, i.e. subtract top Y coordinate of topLine. Add dX and dY shifts.
	//return true if any lines in view had to be recalculated to fit the window width. doDraw flag decides if the lines are actually drawn (doDraw = false is used to test line width overflow, without actually updating the graphical display)
	bool Draw(int &topLine, float height, float dX = 0, float dY = 0, bool doDraw = true);

	//Faster version of Draw where we don't refresh any interactive objects but draw them in their current state. Also doesn't re-check for re-alignment of text objects, assumes everything is correct.
	//Thus this method purely draws the screen from current state, which is assumed to be in the correct.
	//use this whenever you know there cannot be any changes to interactive objects or any other settings (e.g. window dimensions)
	//when there are alot of interactive objects on screen the refresh rate can drop significantly making the interface sluggish, so use the full Draw method sparingly)
	void Draw_Quick(int &topLine, float height, float dX = 0, float dY = 0);

	//-------------------Indexing

	//return TextLine reference at given index
	TextLine& operator[](const int &lineIdx) { if (GoodIdx(LastLine(), lineIdx)) return textLinesVEC[lineIdx]; else return textLinesVEC[LastLine()]; }

	//index with text line index (which will be snapped to start of paragraph) and character index within the paragraph. Returns the text line containing the indexed character
	TextLine& operator[](const INT2 &i2Idx) { return textLinesVEC[GetParaLine_CharIdx(i2Idx).i]; }

	//from lineIdx and paragraph character index given, find paragraph which contains lineIdx and return index of line and character index in that line containing the indexed character
	INT2 GetParaLine_CharIdx(INT2 i2Idx);

	//this goes back the other way : from line index and character index in paragraph containing that line, return character index in paragraph
	int GetParaCharIdx(INT2 i2Idx);

	//-------------------Text modification methods

	//Set a single text line
	void set(const TextLine &setThisTextLine);

	void clear(void) { textLinesVEC.clear(); textLinesVEC.push_back(TextLine()); }

	//append text to last TextObject on last TextLine
	void operator+=(const std::string &rhs);

	//add a text line at the end, replacing the last line if it is blank
	void push(const TextLine &newTextLine);

	//insert text line at given index, snapped to the start of a line sequence
	void insert_paragraph(int lineIdx, const TextLine &insertTextLine);

	//replace paragraph at given index (i.e. replace entire split line sequence which includes the given lineIdx)
	void replace_paragraph(int lineIdx, const TextLine &replaceTextLine);

	//delete paragraph at given index (i.e. delete entire split line sequence which includes the given lineIdx)
	void del_paragraph(int lineIdx);

	//-------------------Text get methods

	//get text of paragraph containing lineIdx
	std::string get_paragraph(int lineIdx);

	//-------------------Indexes retrieval methods

	//from relative coordinates (relative to the window space containing these text lines) and given topLine index, return the text line index
	int GetLineIdx(INT2 relCoord, int topLine);
	//from relative coordinates (relative to the window space containing these text lines) and given topLine index, return the text line and text object index
	INT2 GetTextObjectIdx(INT2 relCoord, int topLine);
	//from relative coordinates (relative to the window space containing these text lines) and given topLine index, return the text line and character index on that line
	INT2 GetCharacterIdx(INT2 relCoord, int topLine);
	//from relative coordinates (relative to the window space containing these text lines) and given topLine index, return the text line, text object and character index on that line : this order, .e. INT3(lineIdx, textobjectIdx, characterIdx)
	INT3 GetFullIdx(INT2 relCoord, int topLine);

	//-------------------Various properties getters and setters

	//last line index
	int LastLine(void) { return ((int)textLinesVEC.size() - 1); }

	//return start line index of last paragraph
	int LastPara(void) { return GetSplitStartIndex(LastLine()); }

	//get length of paragraph containing lineIdx line
	int para_length(int lineIdx);

	int size(void) { return (int)textLinesVEC.size(); }

	bool IsEmpty(void) { return (textLinesVEC.size() == 1 && textLinesVEC[0].length() == 0); }
	bool IsNotEmpty(void) { return !(textLinesVEC.size() == 1 && textLinesVEC[0].length() == 0); }

	//get smallest rectangle encompassing all text lines
	D2D1_RECT_F get_rect(void);

	//-------------------Formatted Text methods

	//Defined in TextFormatting.h

	//Make text using text formatting (either on a new line or at the given line index as an insertion, pushing down all existing lines from that index onwards)
	//<b> </b> : bold

	//<i> </i> : italic

	//[on], [os], [or] : text outline (TOO_NONE, TOO_SQUARE, TOO_ROUND)

	//[tcr,g,b,a/tc] : text color specified using r,g,b,a

	//[bcr,g,b,a/bc] : background color

	//[iomajorId,minorId,auxId,textId/io] : interactive object settings: first two numbers identify the major and minor TextObject type (used to identify the type of action generated when interacting with this TextObject). Next is an auxiliary id number. Last field contains a text identifier.
	//</io> use this to specify interactive object end

	//[anumber/a]	: text alignment object. This is a special type of interactive text object, whose state is checked internally (in CheckLine_InteractiveObjectsState) - not by the state handler method, and has no action assigned to it. 
	//It is used to align text across multiple neighboring lines. This is the first text aligner, see next. The number must be unique on this text line. Number is the majorId in the interactive object.
	//[sanumber/sa] : text alignment object as above. This object will attempt to synchronise with another text aligner object of same number on the line above.
	//Text alignment objects hold a vector of smart pointers to all other text alignment objects to which it is synchronised. 
	//They consist simply of spaces and the number of spaces is modified in the internal state handler (at display update time) so that they all finish at the same character number on each of their respective lines, whilst keeping the number of spaces to a minimum. 
	//This means the text that comes after a text alignment object always starts at the same character index across the text lines holding these synchronised objects (so looks aligned in the console if using a monospaced font). 
	//Can use multiple text alignment objects on the same line, just make sure you number them differently.

	//</c> clear formatting back to default

	//in addition to text also pass the intended lineIndex where this must be placed. This is used by synchronising text alignment objects.
	TextLine BuildFormattedTextLine(std::string text, int lineIndex);
};

#endif
