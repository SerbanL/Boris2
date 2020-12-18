#include "stdafx.h"
#include "TextObject.h"

#if GRAPHICS == 1

////////////////////////////////////////////////////////////////////////////////////////////////////////////TEXT OBJECT

//Default constructor
TextObject::TextObject() : GraphicalObject()
{
	//Default settings
	fulltext = "";
	set_display_text_length();
	CalculateRectangle(0, 0);
}


//almost full constructor : no action handler set
TextObject::TextObject(std::string text_, FormatSpecifier fs) : GraphicalObject()
{
	fulltext = text_;
	set_display_text_length();
	this->fs = fs;
	CalculateRectangle(0, 0);
}

//full constructor
TextObject::TextObject(std::string text_, FormatSpecifier fs, SimTOFunct *pActionHandler_, float X, float Y)
{
	fulltext = text_;
	set_display_text_length();
	this->fs = fs;
	CalculateRectangle(X, Y);

	if (pActionHandler_) {

		pActionHandler = std::shared_ptr<SimTOFunct>(new SimTOFunct(pActionHandler_));
	}
}

//full constructor for text alignment object starter. *ptao_line is the TextLine on which this TextObject being created will be located
TextObject::TextObject(std::string text_, FormatSpecifier fs, SimTOFunct *pActionHandler_, TextLine* ptao_line, float X, float Y)
{
	fulltext = text_;
	set_display_text_length();
	this->fs = fs;
	CalculateRectangle(X, Y);

	if (pActionHandler_) {

		pActionHandler = std::shared_ptr<SimTOFunct>(new SimTOFunct(pActionHandler_));

		//starting tao
		aligned_objects = tao_vec_ptr(new std::vector<tao_ptr>());
		//for now this object is the only tao in the vector
		aligned_objects->push_back(tao_ptr(ptao_line, this));
	}
}

//full constructor for text alignment object with synchronisation (synchronise to existing tao: synch_tao). *ptao_line is the TextLine on which this TextObject being created will be located
TextObject::TextObject(std::string text_, FormatSpecifier fs, SimTOFunct *pActionHandler_, TextLine* ptao_line, TextObject* psynch_tao, float X, float Y)
{
	fulltext = text_;
	set_display_text_length();
	this->fs = fs;
	CalculateRectangle(X, Y);

	if (pActionHandler_) {

		pActionHandler = std::shared_ptr<SimTOFunct>(new SimTOFunct(pActionHandler_));

		//check if we can actually synchronise to the passed (supposedly) tao - if not there's a bug in the code! This textobject will not become a tao then.
		if (psynch_tao->is_text_alignment_object() && pActionHandler_->iop.majorId == psynch_tao->get_text_alignment_object_number()) {

			//add this tao to the shared vector...
			psynch_tao->synch_text_alignment_object(tao_ptr(ptao_line, this));
			//...and get the list as well - ref count increases
			aligned_objects = psynch_tao->aligned_objects;
		}
	}
}

TextObject::~TextObject()
{
	if (aligned_objects) {

		//this is a tao being deleted. Search for its entry in the vector shared with the other synchronised taos (shared through shared_ptr) and erase it.
		for (int idx = 0; idx < aligned_objects->size(); idx++) {

			if ((*aligned_objects)[idx].second == this) {

				aligned_objects->erase(aligned_objects->begin() + idx);
				break;
			}
		}
	}
}

TextObject& TextObject::operator=(const TextObject& copyThis)
{
	if (aligned_objects) {

		//before making the assignment make sure this tao is deleted from the shared vector, otherwise we can end up with dangling pointers
		for (int idx = 0; idx < aligned_objects->size(); idx++) {

			TextObject* to = (*aligned_objects)[idx].second;

			if (to == this) {

				aligned_objects->erase(aligned_objects->begin() + idx);
				break;
			}
		}
	}

	//now can safely make the assignment
	toRect = copyThis.toRect;
	fs = copyThis.fs;
	pActionHandler = copyThis.pActionHandler;
	aligned_objects = copyThis.aligned_objects;
	//only copy full text after the above properties have been copied - important
	fulltext = copyThis.fulltext;
	set_display_text_length();
	//might need to adjust rectangle
	CalculateRectangle();

	//If copyThis is a tao: when making an assignment to_1 = to_2, the intention will be for to_2 to go out of scope before to_1. 
	//In this case we want the address of to_1 to replace that of to_2 in the shared tao vector (note a similar consideration applied for TextLine but that is done in the TextLine class.
	if (aligned_objects) {

		for (int idx = 0; idx < aligned_objects->size(); idx++) {

			if ((*aligned_objects)[idx].second == &copyThis) {

				(*aligned_objects)[idx].second = this;
			}
		}
	}

	return *this;
}

//draw the text object, adding coordinates shifts if specified
void TextObject::Draw(float dX, float dY, float widthLimit)
{
	if (IsNZ(widthLimit) && toRect.right > widthLimit) {

		//width limit will be exceeded - adjust rectangle to this limit and extract the largest substring which fits in it : this is what wil be drawn instead.

		D2D1_RECT_F drawRect = toRect;
		drawRect.right = widthLimit;

		std::string cut_text = fulltext;

		for (int i = 1; i <= fulltext.length(); i++) {

			if (substr_right(0, i) > widthLimit) {

				cut_text = fulltext.substr(0, i - 1);
				break;
			}
		}

		if (cut_text.length() > displayTextLength) {

			//limit displayed text and show ... termination
			pBG->DrawTextLine(cut_text.substr(0, displayTextLength - 3) + std::string("..."), GetShiftedRect(drawRect, dX, dY), fs);
		}
		else pBG->DrawTextLine(cut_text, GetShiftedRect(drawRect, dX, dY), fs);
	}
	else {

		if (fulltext.length() > displayTextLength) {

			//limit displayed text and show ... termination
			pBG->DrawTextLine(fulltext.substr(0, displayTextLength - 3) + std::string("..."), GetShiftedRect(toRect, dX, dY), fs);
		}
		else pBG->DrawTextLine(fulltext, GetShiftedRect(toRect, dX, dY), fs);
	}
}

bool TextObject::SplitOnOverflow(TextObject &nonOverflowingTO, float right, char separator)
{
	//make it same type of object
	nonOverflowingTO = *this;

	//is the current text overflowing?
	if (TextOverflows(right)) {

		//overflow detected: need to extract the largest whole-word std::string which fits in one line
		std::string leftString;
		nonOverflowingTO.set(leftString);

		while (true) {

			//Shift one word at a time from current text to leftString, continue until the largest non-overflowing std::string is extracted
			//It is possible the current text is a single word which overflows. In this case the output TextObject will have a blank std::string.
			if (ShiftSubstring_R2L(leftString, fulltext, separator)) {

				nonOverflowingTO.set(leftString);

				if (nonOverflowingTO.TextOverflows(right)) {

					ShiftSubstring_L2R(leftString, fulltext, separator);
					nonOverflowingTO.set(leftString);
					break;
				}

				continue;
			}

			break;
		}

		//update this text object with the split text.
		set(fulltext);
		set_top_left(nonOverflowingTO.left(), nonOverflowingTO.top());

		//overflow was detected, nonOverflowingTO object will not be overflowing
		return true;
	}

	//no overflow was detected
	return false;
}

void TextObject::ShiftFirstWord(TextObject &oneWordTextObject)
{
	oneWordTextObject = *this;
	std::string oneWordText;
	oneWordTextObject.set(oneWordText);

	if (!ShiftSubstring_R2L(oneWordText, fulltext, ' ')) {
		oneWordTextObject.set(fulltext);
		fulltext = "";
	}
	else oneWordTextObject.set(oneWordText);

	set(fulltext);
	set_top_left(oneWordTextObject.left(), oneWordTextObject.top());
}

bool TextObject::delchar(int delIdx) 
{
	//Delete character from index position (do not delete if index outside range). Return true if deletion results in empty std::string.

	if (GoodIdx((int)fulltext.length() - 1, delIdx)) {

		set(fulltext.substr(0, delIdx) + fulltext.substr(delIdx + 1));
	}

	if (!fulltext.length()) return true;

	return false;
}

D2D1_RECT_F TextObject::substr_rect(int charIdx, int length) 
{
	//get the rectangle of substring starting at index idxS with given length
	//if idxS is not a valid index then just get the rectangle for a space (" ") at the end of the text

	D2D1_RECT_F substringRect;

	if (GoodIdx(displayTextLength - 1, charIdx)) {

		float left = pBG->GetMonospacedFontPixelsWidth() * charIdx + toRect.left;
		float width = pBG->GetMonospacedFontPixelsWidth() * length;

		substringRect = D2D1::RectF(left, toRect.top, left + width, toRect.bottom);
	}
	else {

		float width = (float)pBG->GetMonospacedFontPixelsWidth();

		substringRect = D2D1::RectF(toRect.right, toRect.top, toRect.right + width, toRect.bottom);
	}

	return substringRect;
}

float TextObject::substr_right(int charIdx, int length) 
{
	float right;

	if (GoodIdx(displayTextLength - 1, charIdx)) {

		float left = pBG->GetMonospacedFontPixelsWidth() * charIdx + toRect.left;
		float width = pBG->GetMonospacedFontPixelsWidth() * length;

		right = left + width;
	}
	else {

		float width = (float)pBG->GetMonospacedFontPixelsWidth();

		right = toRect.right + width;
	}

	return right;
}

int TextObject::get_char_index(INT2 relMouse) 
{
	for (int l = 1; l <= displayTextLength; l++) {

		if (relMouse.i <= (int)substr_right(0, l)) {

			return l - 1;
		}
	}

	return -1;
}

InteractiveObjectStateChange TextObject::CheckStateChange(void) 
{
	InteractiveObjectStateChange stateChanged;

	if (pActionHandler) {

		//This is a normal interactive text object
		if (!aligned_objects) {

			stateChanged = pActionHandler->State(this);
		}
		//this is a text alignment object (tao) so do not use the external state handler : internal state handler
		else {

			TextLine* tl = (*aligned_objects)[0].first;
			TextObject* to = (*aligned_objects)[0].second;

			int min_length = to->displayedtext_length();
			int char_idx = tl->CharIndex_of_TextObject(to) + to->displayedtext_length();
			
			//check if there is a mismatch in ending character indexes.
			for (int idx = 1; idx < aligned_objects->size(); idx++) {

				tl = (*aligned_objects)[idx].first;
				to = (*aligned_objects)[idx].second;

				if (char_idx != (tl->CharIndex_of_TextObject(to) + to->displayedtext_length())) {

					pActionHandler->iop.state = IOS_MISALIGNEDTAO;
					return stateChanged(true);
				}

				min_length = (min_length < to->displayedtext_length() ? min_length : to->displayedtext_length());
			}
			
			//also return IOS_MISALIGNEDTAO if the minimum length across all synchronised taos is not met
			
			if (min_length > 0) {

				pActionHandler->iop.state = IOS_MISALIGNEDTAO;
				return stateChanged(true);
			}
		}
	}

	return stateChanged;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////TEXT LINE

TextLine::TextLine(void) 
{
	splitLine = SPLIT_NO;

	//Make a default blank entry
	textLineVEC.push_back(TextObject("", FormatSpecifier()));
}

TextLine::TextLine(std::string text, FormatSpecifier fs) 
{
	splitLine = SPLIT_NO;

	textLineVEC.push_back(TextObject(text, fs));
}

TextLine::~TextLine()
{
	//no need to do anything here : TextObject destructors will be called for all TextObjects in textLineVEC, so any cleaning delegated to them
}

TextLine& TextLine::operator=(const TextLine& copyThis) 
{
	//before assigning copyThis to *this we must make sure taos in this object have been removed first. If we simply do the assignment they'll become dangling taos (i.e. their stored pointers will become dangling).
	//This can happen in some surprising situations, e.g. this TextLine is in a vector and we delete its entry - it is very much possible its destructor is never called!!!
	//Instead what happens std::vector assigns another TextLine at this address by calling this assignment operator, then the assigned TextLine destructor is called instead.
	//We must then make sure any taos here are removed first, and the new taos are assigned properly.

	//search for taos and delete them if they hold a pointer to this textline
	for (int to_idx = 0; to_idx < textLineVEC.size(); to_idx++) {

		//found a tao - search the shared vector for this TextLine pointer
		if (textLineVEC[to_idx].is_text_alignment_object()) {

			int tao_idx = 0;
			while (tao_idx < textLineVEC[to_idx].aligned_objects->size()) {

				TextLine* tl = (*(textLineVEC[to_idx].aligned_objects))[tao_idx].first;

				//found a tao on this TextLine holding a valid pointer : erase it
				if (tl == this) {

					textLineVEC[to_idx].aligned_objects->erase(textLineVEC[to_idx].aligned_objects->begin() + tao_idx);
				}
				else tao_idx++;
			}
		}
	}

	//can now safely make the assignment
	textLineVEC = copyThis.textLineVEC;
	splitLine = copyThis.splitLine;
	
	//check all text objects on this line to see if they are tao (text alignment object)
	for (int to_idx = 0; to_idx < textLineVEC.size(); to_idx++) {

		if (textLineVEC[to_idx].is_text_alignment_object()) {

			//found one - must now replace any copyThis pointers in their shared tao vector with this
			for (int tao_idx = 0; tao_idx < textLineVEC[to_idx].aligned_objects->size(); tao_idx++) {

				if ((*(textLineVEC[to_idx].aligned_objects))[tao_idx].first == &copyThis) {

					(*(textLineVEC[to_idx].aligned_objects))[tao_idx].first = this;
				}
			}
		}
	}

	return *this;
}

//check for tao_line_in pointers in tao shared vector and replace them with pointer to tao_line_out
void TextLine::FixTaoLine(TextLine* tao_line)
{
	//check all text objects on this line to see if they are tao (text alignment object)
	for (int to_idx = 0; to_idx < textLineVEC.size(); to_idx++) {

		if (textLineVEC[to_idx].is_text_alignment_object()) {

			//found one - must now replace any tao_line pointers in their shared tao vector with this
			for (int tao_idx = 0; tao_idx < textLineVEC[to_idx].aligned_objects->size(); tao_idx++) {

				if ((*(textLineVEC[to_idx].aligned_objects))[tao_idx].first == tao_line) {

					(*(textLineVEC[to_idx].aligned_objects))[tao_idx].first = this;
				}
			}
		}
	}
}

void TextLine::Draw(float dX, float dY, float widthLimit) 
{
	for (int i = 0; i <= LastElem(); i++) {

		textLineVEC[i].Draw(dX, dY, widthLimit);
	}
}

INT2 TextLine::IndextoINT2Index(int idx) 
{
	INT2 i2idx;
	int totalLength1 = 0, totalLength2 = 0;

	//first find TextObject for which idx indexes its text
	for (int i = 0; i <= LastElem(); i++) {

		totalLength2 += textLineVEC[i].displayedtext_length();

		if (idx < totalLength2) break;
		else { i2idx.i++; totalLength1 = totalLength2; }
	}

	//finally find index inside the i2idx.i TextObject
	i2idx.j = idx - totalLength1;

	if (i2idx.i > LastElem()) {

		i2idx.i = LastElem();
		i2idx.j = textLineVEC[i2idx.i].displayedtext_length();
	}

	return i2idx;
}

void TextLine::RecalculateRectangles(int idx) 
{
	if (idx <= 0) {

		//rectangles must start from 0 X coordinate
		textLineVEC[0].set_top_left(0, textLineVEC[0].top());
		idx = 1;
	}

	for (int i = idx; i <= LastElem(); i++) {

		float X = 0, Y = 0;
		textLineVEC[i - 1].top_right(X, Y);
		textLineVEC[i].set_top_left(X, Y);
	}
}

int TextLine::get_char_index(INT2 relMouse) 
{
	int length = 0;

	for (int i = 0; i < textLineVEC.size(); i++) {

		if (relMouse.i <= (int)textLineVEC[i].right()) {

			//found clicked TextObject with index i
			return (length + textLineVEC[i].get_char_index(relMouse));
		}

		length += textLineVEC[i].displayedtext_length();
	}

	return -1;
}

int TextLine::get_textobject_index(INT2 relMouse) 
{
	for (int i = 0; i < textLineVEC.size(); i++) {

		if (relMouse.i <= (int)textLineVEC[i].right()) {

			//found clicked TextObject with index i
			return i;
		}
	}

	return -1;
}

INT2 TextLine::get_char_and_textobject_index(INT2 relMouse) 
{
	int length = 0;

	for (int i = 0; i < textLineVEC.size(); i++) {

		if (relMouse.i <= (int)textLineVEC[i].right()) {

			//found clicked TextObject with index i
			return INT2(i, length + textLineVEC[i].get_char_index(relMouse));
		}

		length += textLineVEC[i].displayedtext_length();
	}

	return INT2(-1, -1);
}

//fix the rectangle of the last TextObject on this line (i.e. make it snap to the text) - used on lines which have been split, since their last rectangle is elongated up to the window width.
TextLine& TextLine::FixSplitLineRectangle(void) 
{
	if (textLineVEC[LastElem()].IsBlank() && !textLineVEC[LastElem()].IsActionSet())
		delobject(LastElem());
	else textLineVEC[LastElem()].CalculateRectangle();

	return *this;
}

//elongate rectangle of last TextObject on this line up to the widthLimit value - used for split lines
TextLine& TextLine::SetSplitLineRectangle(float widthLimit) 
{
	//insert a blank text object with same formatting as the last one, except for no outline, and elongate it to width limit.
	FormatSpecifier fs;
	fs.textOutline = TOO_NONE;
	fs.bgrndColor = D2D1::ColorF(0, 0, 0, 0);

	push(TextObject("", fs));
	textLineVEC[LastElem()].toRect.right = widthLimit;

	return *this;
}

//Get text of TextLine object into a std::string
std::string& operator<<(std::string &lhs, const TextLine &rhs) 
{
	lhs = "";

	for (int i = 0; i < rhs.textLineVEC.size(); i++) {

		lhs += rhs.textLineVEC[i];
	}

	return lhs;
}

//add two text lines into one (e.g. textLine += addThisLine;)
TextLine& operator+=(TextLine &lhs, TextLine &rhs) 
{
	//append elements in textLineVEC from rhs to lhs
	//At the point of append, check if the adjacent TextObjects have identical formatting. If so, join them into a single object.

	if (lhs[lhs.LastElem()].SameType(rhs[0])) {

		lhs[lhs.LastElem()] += rhs[0].GetText();
		lhs.textLineVEC.insert(lhs.textLineVEC.end(), rhs.textLineVEC.begin() + 1, rhs.textLineVEC.end());
	}
	else lhs.textLineVEC.insert(lhs.textLineVEC.end(), rhs.textLineVEC.begin(), rhs.textLineVEC.end());

	//rectangles will not be correct now, so fix this : we want rhs to be a continuation of lhs now
	lhs.RecalculateRectangles();

	return lhs;
}

//Add a TextObject at the end of this line, replacing the last TextObject if blank
void TextLine::push(const TextObject &newTextObject) 
{
	//If last TextObject is blank then replace it with current one (there's no need to keep a blank object inserted)
	if (textLineVEC[LastElem()].IsBlank() && !textLineVEC[LastElem()].IsActionSet())
		textLineVEC.pop_back();

	textLineVEC.push_back(newTextObject);

	RecalculateRectangles(LastElem());
}

template void TextLine::insert<bool>(int charIdx, bool &insertText);
template void TextLine::insert<char>(int charIdx, char &insertText);
template void TextLine::insert<int>(int charIdx, int &insertText);
template void TextLine::insert<size_t>(int charIdx, size_t &insertText);
template void TextLine::insert<float>(int charIdx, float &insertText);
template void TextLine::insert<double>(int charIdx, double &insertText);
template void TextLine::insert<std::string>(int charIdx, std::string &insertText);
template void TextLine::insert<INT2>(int charIdx, INT2 &insertText);
template void TextLine::insert<INT3>(int charIdx, INT3 &insertText);

template <typename VType> void TextLine::insert(int charIdx, VType &insertText)
{
	//allow charIdx = length(), which results in appending text at the end
	if (GoodIdx(length(), charIdx)) {

		INT2 i2idx = IndextoINT2Index(charIdx);
		textLineVEC[i2idx.i].insert(i2idx.j, insertText);
		RecalculateRectangles(i2idx.i + 1);
	}
}

bool TextLine::delchar(int delIdx) 
{
	//delete character at the required index.
	//If deletion results in empty std::string for that TextObject, then erase the TextObject too, unless it is the last one on the line. In the latter case return true (deletion has resulted in empty textLine)
	//Finally recalculate rectangles so they are contiguous.

	if (!GoodIdx(length() - 1, delIdx)) return false;

	INT2 i2idx = IndextoINT2Index(delIdx);

	if (textLineVEC[i2idx.i].delchar(i2idx.j)) {

		//erase TextObject if it's blank now, unless it is the last on this text line
		if (size() > 1) {

			textLineVEC.erase(textLineVEC.begin() + i2idx.i);
			RecalculateRectangles(i2idx.i);
			return false;
		}
		else return true;
	}
	else {

		RecalculateRectangles(i2idx.i + 1);
		return false;
	}
}

bool TextLine::deltext(int startIdx, int numChars) 
{
	bool empty = false;

	if (numChars < 0) numChars = length() - startIdx;

	for (int idx = 0; idx < numChars; idx++)
		empty |= delchar(startIdx);

	return empty;
}

void TextLine::delobject(int delIdx)
{
	if (LastElem() == 0) {

		//do not delete last remaining element, but set its text to blank
		textLineVEC[0].set("");
		return;
	}

	if (!GoodIdx(LastElem(), delIdx)) delIdx = LastElem();

	textLineVEC.erase(textLineVEC.begin() + delIdx);
	RecalculateRectangles(delIdx);
}

int TextLine::CharIndex_of_TextObject(TextObject* pTO)
{
	int length = 0;

	for (int to_idx = 0; to_idx < textLineVEC.size(); to_idx++) {

		if (&(textLineVEC[to_idx]) == pTO) return length;

		length += textLineVEC[to_idx].displayedtext_length();
	}

	//-1 if couldn't find it
	return -1;
}

int TextLine::length(void) 
{
	int textLength = 0;

	for (int i = 0; i <= LastElem(); i++) {

		textLength += textLineVEC[i].displayedtext_length();
	}

	return textLength;
}

D2D1_RECT_F TextLine::get_char_rect(int charIdx) 
{
	INT2 i2idx = IndextoINT2Index(charIdx);

	return textLineVEC[i2idx.i].substr_rect(i2idx.j);
}

void TextLine::set_top_left(float X, float Y) 
{
	textLineVEC[0].set_top_left(X, Y);

	RecalculateRectangles(1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////TEXT LINES

TextLines::TextLines(void) {

	this->pActionHandler = nullptr;

	textLinesVEC.push_back(TextLine());

	this->widthLimit = 0;
}

TextLines::TextLines(SimTOFunct *pActionHandler) {

	this->pActionHandler = pActionHandler;

	textLinesVEC.push_back(TextLine());

	this->widthLimit = 0;
}

void TextLines::RecalculateRectangles(int idx) {

	if (idx <= 0) {

		float X = 0, Y = 0;
		textLinesVEC[0].set_top_left(0, 0);

		idx = 1;
	}

	for (int i = idx; i <= LastLine(); i++) {

		float X = 0, Y = 0;
		textLinesVEC[i - 1].bottom_left(X, Y);
		textLinesVEC[i].set_top_left(X, Y);
	}
}

bool TextLines::CheckLine_InteractiveObjectState(int lineIdx) 
{
	bool stateChanged = false;

	//check TextObjects on this line one by one
	for (int i = 0; i < textLinesVEC[lineIdx].size(); i++) {

		if (i < textLinesVEC[lineIdx].size() && textLinesVEC[lineIdx][i].IsActionSet()) {

			//Action object - check it
			InteractiveObjectStateChange IOS_change = textLinesVEC[lineIdx][i].CheckStateChange();
			if (IOS_change.stateChanged) {

				//This object has changed. Check special state values and make changes as required.
				stateChanged = true;

				switch (textLinesVEC[lineIdx][i].get_iop_state()) {

				case IOS_DELETING:
					//must delete this TextObject
					textLinesVEC[lineIdx].delobject(i);
					break;

				case IOS_DELETINGPARAGRAPH:
					del_paragraph(lineIdx);
					//must return as this line has been deleted : nothing left to do here
					return stateChanged;
					break;

				case IOS_REPLACINGPARAGRAPH:
				{
					replace_paragraph(lineIdx, BuildFormattedTextLine(IOS_change.textMessage, lineIdx));
					//must return as this line has been replaced : nothing left to do here
					return stateChanged;
				}
				break;

				case IOS_WASLASTINLIST:
				{
					//this object was last in list, but is no longer last - it is asking for a new element to be inserted after it (i.e. the next element after it in the list - this could be the last in list now, but the state handler will verify this on next call)
					textLinesVEC[lineIdx][i].set_iop_state(IOS_ON);

					//next line down - note, in a list the entries occupy one paragraph each (but there may be multiple interactive objects per list entry - it could be that insertFormattedText will spawn multiple interactive objects)
					int insertIdx = GetSplitEndIndex(lineIdx) + 1;
					if (insertIdx <= LastLine()) {

						//The message may consist of multiple lines - insert them one by one
						//Check for newline and carriage return separators (\n, \r)
						std::vector<std::string> messageLines = split(IOS_change.textMessage, "\r", "\n");

						for (int idx = 0; idx < messageLines.size(); idx++) {

							insert_paragraph(idx + insertIdx, BuildFormattedTextLine(messageLines[idx], idx + insertIdx));
						}
					}
					else {
						//The message may consist of multiple lines - insert them one by one
						//Check for newline and carriage return separators (\n, \r)
						std::vector<std::string> messageLines = split(IOS_change.textMessage, "\r", "\n");

						for (int idx = 0; idx < messageLines.size(); idx++) {

							push(BuildFormattedTextLine(messageLines[idx], idx + insertIdx));
						}
					}
				}
				break;

				case IOS_MISALIGNEDTAO:
				{
					textLinesVEC[lineIdx][i].set_iop_state(IOS_ON);

					//go through all taos in the shared vector of our tao at [lineIdx][i]

					TextLine* tl = (*(textLinesVEC[lineIdx][i].aligned_objects))[0].first;
					TextObject* to = (*(textLinesVEC[lineIdx][i].aligned_objects))[0].second;

					int max_startchar_idx = tl->CharIndex_of_TextObject(to);

					//pass 1 : find maximum starting index
					for (int to_idx = 1; to_idx < textLinesVEC[lineIdx][i].aligned_objects->size(); to_idx++) {

						tl = (*(textLinesVEC[lineIdx][i].aligned_objects))[to_idx].first;
						to = (*(textLinesVEC[lineIdx][i].aligned_objects))[to_idx].second;
						
						int start_char_idx = tl->CharIndex_of_TextObject(to);
						max_startchar_idx = (max_startchar_idx > start_char_idx ? max_startchar_idx : start_char_idx);
					}
					
					//pass 2 : adjust taos so they they end at max_startchar_idx : this means at least one of them will have length zero.
					for (int to_idx = 0; to_idx < textLinesVEC[lineIdx][i].aligned_objects->size(); to_idx++) {

						tl = (*(textLinesVEC[lineIdx][i].aligned_objects))[to_idx].first;
						to = (*(textLinesVEC[lineIdx][i].aligned_objects))[to_idx].second;

						int start_char_idx = tl->CharIndex_of_TextObject(to);
						int length = max_startchar_idx - start_char_idx;
						to->set(std::string(length, ' '));

						//the line also needs the rectangles recalculated
						tl->RecalculateRectangles();
					}
				}
				break;

				default:
					break;
				}

				//recalculate rectangles on this line if not empty
				if (textLinesVEC[lineIdx].length()) {

					textLinesVEC[lineIdx].RecalculateRectangles();
				}
				else {

					//delete empty line - note, using delline not del_paragraph, as this line might be part of a paragraph (split line sequence)
					delline(lineIdx);
					break;
				}
			}
		}
	}

	//if any TextObject state has changed then must also recalculate rectangles from this line upwards.
	if (stateChanged) RecalculateRectangles(lineIdx);

	return stateChanged;
}

std::vector<TextLine> TextLines::SplitLinetoFitWidth(TextLine& textLine)
{
	std::vector<TextLine> split_textLine;

	if (widthLimit > 0 && !textLine.IsBlank()) {

		TextLine nonOverflowingTextLine;
		TextLine textRemaining = textLine;

		float appendWidth = 0;

		//go through the TextObjects on the line one by one : always work on the 0th element by deleting objects from the start after processing them
		while (!textRemaining.IsBlank()) {

			if (appendWidth + textRemaining[0].width() <= widthLimit) {

				appendWidth += textRemaining[0].width();

				//TextObject does not exceed width limit: add it and remove it from textRemaining so we can look at the next object
				nonOverflowingTextLine.push(textRemaining[0]);
				textRemaining.delobject(0);
			}
			else if (!textRemaining[0].IsActionSet()) {
				
				//TextObject exceeds text limit: split it, unless it is an interactive object (action handler set)
				TextObject nonOverflowingTextObject;
				textRemaining[0].SplitOnOverflow(nonOverflowingTextObject, widthLimit - appendWidth, ' ');

				//append the first part of the split TextObject to finish the non-overflowing line
				if (!nonOverflowingTextObject.IsBlank())
					nonOverflowingTextLine.push(nonOverflowingTextObject);

				//if nonOverflowingTextLine is blank then no object could be pushed on this line, either whole or split - block of text too long.
				//In this case just add the next word as is : line will overflow, but cannot be helped.
				if (nonOverflowingTextLine.IsBlank()) {

					//get next word from textRemaining
					textRemaining[0].ShiftFirstWord(nonOverflowingTextObject);

					//add it to the line being built, but do not delete the first object from textRemaining - there may be further words in it. (i.e. only delete it if blank)
					nonOverflowingTextLine.push(nonOverflowingTextObject);
				}

				//push the finished line in the vector
				split_textLine.push_back(nonOverflowingTextLine);

				//start again for any remaining text
				nonOverflowingTextLine.clear();
				appendWidth = 0;
			}
			else if (nonOverflowingTextLine.IsBlank()) {

				//textRemaining[0] is an interactive type and couldn't be added to the line because it's too long - interactive objects should be fairly short!
				//In this case just add it as is : line will overflow, but cannot be helped.
				nonOverflowingTextLine.push(textRemaining[0]);
				textRemaining.delobject(0);

				//push the finished line in the vector
				split_textLine.push_back(nonOverflowingTextLine);

				//start again for any remaining text
				nonOverflowingTextLine.clear();
				appendWidth = 0;
			}
			else {

				//push the finished line in the vector
				split_textLine.push_back(nonOverflowingTextLine);

				//start again for any remaining text
				nonOverflowingTextLine.clear();
				appendWidth = 0;
			}
		}

		//add last text line obtaining after processing last TextObject
		if (!nonOverflowingTextLine.IsBlank()) {

			split_textLine.push_back(nonOverflowingTextLine);
		}

		//used push_back to build nonOverflowingTextLine which was then used to build split_textLine : the TextLine pointers in the shared tao vector sitll has &textRemaining : replace these as required
		for (int idx = 0; idx < split_textLine.size(); idx++) {
			
			split_textLine[idx].FixTaoLine(&textRemaining);
		}
	}
	else {

		split_textLine.push_back(textLine);
	}

	return split_textLine;
}

void TextLines::FitLinetoWidth(int &lineIdx)
{
	//If this is part of a previously split line, get the line index of the SPLIT_START line.
	lineIdx = GetSplitStartIndex(lineIdx);
	
	//if we have a split line sequence, recombine it into a single line before splitting it again to the new width limit
	if (textLinesVEC[lineIdx].splitLine == SPLIT_START) {

		//append first line with SPLIT_START
		TextLine recombinedLine = textLinesVEC[lineIdx].FixSplitLineRectangle();
		recombinedLine.splitLine = SPLIT_NO;

		int endIdx = GetSplitEndIndex(lineIdx);

		for (int idx = lineIdx + 1; idx < endIdx; idx++) {

			recombinedLine += textLinesVEC[idx].FixSplitLineRectangle();
			recombinedLine.FixTaoLine(&textLinesVEC[idx]);
		}

		recombinedLine += textLinesVEC[endIdx];
		recombinedLine.FixTaoLine(&textLinesVEC[endIdx]);

		textLinesVEC.erase(textLinesVEC.begin() + lineIdx, textLinesVEC.begin() + endIdx);

		//now replace the line at lineIdx with the original unsplit line
		textLinesVEC[lineIdx] = recombinedLine;

		RecalculateRectangles(lineIdx);
	}
	
	//for the line at lineIdx, split it into multiple lines which fit to widthLimit
	std::vector<TextLine> split_textLine = SplitLinetoFitWidth(textLinesVEC[lineIdx]);
	
	if (split_textLine.size() > 1) {

		//1. Mark start, end, and in-between lines so they can be recognized later as having been part of a single line.
		//Extend rectangle of last object on each line up to exactly the current window width - any mismatch will be used to trigger recalculation of previously split lines.

		split_textLine.front().SetSplitLineRectangle(widthLimit).splitLine = SPLIT_START;

		for (int i = 1; i < (int)split_textLine.size() - 1; i++)
			split_textLine[i].SetSplitLineRectangle(widthLimit).splitLine = SPLIT_MIDDLE;

		//do not elongate rectangle as with the other lines
		split_textLine.back().splitLine = SPLIT_END;

		//lineIdx is now set at SPLIT_START. Be careful to avoid infinite loops when a split text line cannot be fitted to widthLimit (e.g. window too small for word length)
		//Also note, these have not been drawn yet. Draw them and continue from line just after SPLIT_END line.
	}

	//replace original line with the split line sequence
	textLinesVEC.erase(textLinesVEC.begin() + lineIdx);
	textLinesVEC.insert(textLinesVEC.begin() + lineIdx, split_textLine.begin(), split_textLine.end());
	
	RecalculateRectangles(lineIdx);
}

bool TextLines::LineNeedsRecalculation(int lineIdx) 
{
	//either line overflows width, or a split line needs to be recalculated (its width doesn't match the current width)
	return ((textLinesVEC[lineIdx].right() > widthLimit && IsNZ(widthLimit)) ||
		((textLinesVEC[lineIdx].splitLine == SPLIT_START || textLinesVEC[lineIdx].splitLine == SPLIT_MIDDLE) && IsNZ(textLinesVEC[lineIdx].right() - widthLimit)));
}

int TextLines::GetSplitStartIndex(int lineIdx) 
{
	while (lineIdx > 0 && (textLinesVEC[lineIdx].splitLine == SPLIT_MIDDLE || textLinesVEC[lineIdx].splitLine == SPLIT_END)) lineIdx--;

	return lineIdx;
}

int TextLines::GetSplitEndIndex(int lineIdx) 
{
	while (lineIdx < LastLine() && (textLinesVEC[lineIdx].splitLine == SPLIT_START || textLinesVEC[lineIdx].splitLine == SPLIT_MIDDLE)) lineIdx++;

	return lineIdx;
}

//from lineIdx and paragraph character index given, find paragraph which contains lineIdx and return index of line and character index in that line containing the indexed character
INT2 TextLines::GetParaLine_CharIdx(INT2 i2Idx) 
{
	int idx = GetSplitStartIndex(i2Idx.i);

	for (; idx <= GetSplitEndIndex(i2Idx.i); idx++) {

		i2Idx.j -= textLinesVEC[idx].length();

		if (i2Idx.j < 0) {

			return INT2(idx, i2Idx.j + textLinesVEC[idx].length());
		}
	}

	//not a character index, this is the index on last line just after the last character - e.g. this could be the prompter position at the end of an entry line (paragraph)
	return INT2(idx - 1, i2Idx.j + textLinesVEC[idx - 1].length());
}

int TextLines::GetParaCharIdx(INT2 i2Idx) 
{
	int length = 0;

	for (int idx = GetSplitStartIndex(i2Idx.i); idx < i2Idx.i; idx++) {

		length += textLinesVEC[idx].length();
	}

	return length + i2Idx.j;
}

bool TextLines::Draw(int &topLine, float height, float dX, float dY, bool doDraw) 
{
	bool linesRecalculated = false;

	if (!GoodIdx(LastLine(), topLine)) return linesRecalculated;

	float relY = textLinesVEC[topLine].top();

	/*
	//OLD METHOD
	//first go through all the lines on screen to check for any changes needed - e.g. interactive object state changes and to fit lines in console width
	for (int i = topLine; ; i++) {

		//LastLine() can change value during the loop
		if (i > LastLine()) break;

		//Check visible interactive object states on this line before drawing - if any of them have changed then recalculate lines before drawing. 
		//Also delete TextObject if IOS_DELETING state is set. Insert objects if IOS_WASLASTINLIST is set. All this is done in the routine below.
		//If a change did occur then start again - a change in an object can affect text above.
		
		if (CheckLine_InteractiveObjectState(i)) { i = topLine; continue; }
		
		//now check this line actually fits within the width limit - if not, split it into multiple lines which do fit (this will mean LastLine() value increases)
		if (LineNeedsRecalculation(i)) {

			//the topLine index could also change (sometimes alot). e.g. window has small size so text is wrapped in short lines, then window is enlarged. When scrolling the sudden change in number of lines can be significant.
			//measure the change in topLine index by look ast LastLine() index before and after.
			int before = LastLine();

			FitLinetoWidth(i);
			linesRecalculated = true;

			int after = LastLine();

			//now adjust topLine index (was passed through reference). This seems to work very well, i.e. text displays where you think it should.
			topLine += (after - before);
			if (topLine < 0) topLine = 0;
			if (topLine > LastLine()) topLine = after;
			
			//advance index i to end of freshly split line and check if still in view
			i = GetSplitEndIndex(i);
			if (textLinesVEC[i].bottom() - relY > height) break;	//only need lines in view
		}
		else if (textLinesVEC[i].bottom() - relY > height) break;	//only need lines in view
	}
	*/
	if (doDraw) {

		//first go through all the lines on screen to check for any changes needed to interactive objects
		for (int i = topLine; ; i++) {

			//LastLine() can change value during the loop now since we could be adding or deleting objects.
			if (i > LastLine()) break;

			//Check visible interactive object states on this line before drawing - if any of them have changed then recalculate lines before drawing. 
			//Also delete TextObject if IOS_DELETING state is set. Insert objects if IOS_WASLASTINLIST is set. All this is done in the routine below.

			CheckLine_InteractiveObjectState(i);

			//only need lines in view
			if (textLinesVEC[i].bottom() - relY > height) break;
		}
	}

	//next fit lines in console width
	for (int i = topLine; ; i++) {

		//LastLine() can change value during the loop now since we could be splitting or recombining lines.
		if (i > LastLine()) break;

		//now check this line actually fits within the width limit - if not, split it into multiple lines which do fit (this will mean LastLine() value increases)
		if (LineNeedsRecalculation(i)) {

			//the topLine index could also change (sometimes alot). e.g. window has small size so text is wrapped in short lines, then window is enlarged. When scrolling the sudden change in number of lines can be significant.
			//measure the change in topLine index by look ast LastLine() index before and after.
			int before = LastLine();

			FitLinetoWidth(i);
			linesRecalculated = true;

			int after = LastLine();

			//now adjust topLine index (was passed through reference). This seems to work very well from a user point of view, i.e. text displays where you think it should (POLA!)
			topLine += (after - before);
			if (topLine < 0) topLine = 0;
			if (topLine > LastLine()) topLine = after;

			//advance index i to end of freshly split line and check if still in view
			i = GetSplitEndIndex(i);
			if (textLinesVEC[i].bottom() - relY > height) break;	//only need lines in view
		}
		else if (textLinesVEC[i].bottom() - relY > height) break;	//only need lines in view
	}

	if (doDraw) {

		//now all the lines have been adjusted we can just draw them
		for (int i = topLine; i <= LastLine(); i++) {

			if (textLinesVEC[i].bottom() - relY <= height) {

				//in some cases we don't want to graphically update the display, but just recalculate the text lines - e.g. if called from somewhere outside of the required { pBG->BeginD3DDraw(); ... pBG->EndD3DDraw(); }
				textLinesVEC[i].Draw(dX, dY - relY, widthLimit);
			}
			else return linesRecalculated;		//only draw the lines in view
		}
	}

	return linesRecalculated;
}

//Faster version of Draw where we don't refresh any interactive objects but draw them in their current state. Also doesn't re-check for re-alignment of text objects, assumes everything is correct.
//Thus this method purely draws the screen from current state, which is assumed to be in the correct.
//use this whenever you know there cannot be any changes to interactive objects or any other settings (e.g. window dimensions)
//when there are alot of interactive objects on screen the refresh rate can drop significantly making the interface sluggish, so use the full Draw method sparingly)
void TextLines::Draw_Quick(int &topLine, float height, float dX, float dY)
{
	if (!GoodIdx(LastLine(), topLine)) return;

	float relY = textLinesVEC[topLine].top();

	//now all the lines have been adjusted we can just draw them
	for (int i = topLine; i <= LastLine(); i++) {

		if (textLinesVEC[i].bottom() - relY <= height) {

			//in some cases we don't want to graphically update the display, but just recalculate the text lines - e.g. if called from somewhere outside of the required { pBG->BeginD3DDraw(); ... pBG->EndD3DDraw(); }
			textLinesVEC[i].Draw(dX, dY - relY, widthLimit);
		}
		else return;		//only draw the lines in view
	}
}

//Set a single text line
void TextLines::set(const TextLine &setThisTextLine) 
{
	textLinesVEC.clear();
	textLinesVEC.push_back(setThisTextLine);
	RecalculateRectangles();

	int lineIdx = LastLine();
	if (LineNeedsRecalculation(lineIdx)) FitLinetoWidth(lineIdx);
}

//append text to last TextObject on last TextLine
void TextLines::operator+=(const std::string &rhs)
{
	textLinesVEC[LastLine()] += rhs;

	int lineIdx = LastLine();
	if (LineNeedsRecalculation(lineIdx)) FitLinetoWidth(lineIdx);
}


//add a text line at the end, replacing the last line if it is blank
void TextLines::push(const TextLine &newTextLine)
{
	float X = 0, Y = 0;

	if (!textLinesVEC[LastLine()].length()) {

		textLinesVEC[LastLine()].top_left(X, Y);
		textLinesVEC.pop_back();
	}
	else textLinesVEC[LastLine()].bottom_left(X, Y);

	textLinesVEC.push_back(newTextLine);

	textLinesVEC[LastLine()].set_top_left(X, Y);

	int lineIdx = LastLine();
	if (LineNeedsRecalculation(lineIdx)) FitLinetoWidth(lineIdx);
}

///insert text line at given index, snapped to the start of a line sequence
void TextLines::insert_paragraph(int lineIdx, const TextLine &insertTextLine)
{
	if (!GoodIdx(LastLine(), lineIdx)) return;

	//line must be inserted before the start of a split sequence
	lineIdx = GetSplitStartIndex(lineIdx);

	textLinesVEC.insert(textLinesVEC.begin() + lineIdx, insertTextLine);
	//do not use this method to insert part of a split line
	textLinesVEC[lineIdx].splitLine = SPLIT_NO;

	RecalculateRectangles(lineIdx);

	if (LineNeedsRecalculation(lineIdx)) FitLinetoWidth(lineIdx);
}

//replace paragraph at given index (i.e. replace entire split line sequence which includes the given lineIdx)
void TextLines::replace_paragraph(int lineIdx, const TextLine &replaceTextLine)
{
	if (!GoodIdx(LastLine(), lineIdx)) return;

	//entire split line sequence must be replaced
	lineIdx = GetSplitStartIndex(lineIdx);

	//erase all but the lineIdx line
	textLinesVEC.erase(textLinesVEC.begin() + lineIdx + 1, textLinesVEC.begin() + GetSplitEndIndex(lineIdx) + 1);

	//now replace it
	textLinesVEC[lineIdx] = replaceTextLine;
	textLinesVEC[lineIdx].splitLine = SPLIT_NO;

	RecalculateRectangles(lineIdx);

	if (LineNeedsRecalculation(lineIdx)) FitLinetoWidth(lineIdx);
}

//delete paragraph at given index (i.e. delete entire split line sequence which includes the given lineIdx)
void TextLines::del_paragraph(int lineIdx) 
{
	if (!GoodIdx(LastLine(), lineIdx)) return;

	//entire split line sequence must be deleted
	lineIdx = GetSplitStartIndex(lineIdx);

	textLinesVEC.erase(textLinesVEC.begin() + lineIdx, textLinesVEC.begin() + GetSplitEndIndex(lineIdx) + 1);

	RecalculateRectangles(lineIdx);
}

void TextLines::delline(int lineIdx) {

	//this allows deletion of a line in a split line sequence

	if (!GoodIdx(LastLine(), lineIdx)) return;

	//need to consider special cases of SPLIT_START and SPLIT_END : move these flags up or down in the split line sequence
	if (textLinesVEC[lineIdx].splitLine == SPLIT_START && lineIdx < LastLine()) {

		//move SPLIT_START indicator up one as this line is about to be deleted (unless the next line ends the split sequence)
		if (textLinesVEC[lineIdx + 1].splitLine != SPLIT_END)
			textLinesVEC[lineIdx + 1].splitLine = SPLIT_START;
		else textLinesVEC[lineIdx + 1].splitLine = SPLIT_NO;
	}

	if (textLinesVEC[lineIdx].splitLine == SPLIT_END && lineIdx > 0) {

		//move SPLIT_END indicator down one as this line is about to be deleted (unless line below starts the split sequence)
		if (textLinesVEC[lineIdx - 1].splitLine != SPLIT_START)
			textLinesVEC[lineIdx - 1].splitLine = SPLIT_END;
		else textLinesVEC[lineIdx - 1].splitLine = SPLIT_NO;
	}

	//now delete it
	textLinesVEC.erase(textLinesVEC.begin() + lineIdx);

	RecalculateRectangles(lineIdx);
}

std::string TextLines::get_paragraph(int lineIdx) {

	std::string text;

	for (int idx = GetSplitStartIndex(lineIdx); idx <= GetSplitEndIndex(lineIdx); idx++)
		text += textLinesVEC[idx];

	return text;
}

int TextLines::GetLineIdx(INT2 relCoord, int topLine) {

	if (!GoodIdx(LastLine(), topLine)) return -1;

	float Y = (float)relCoord.j + textLinesVEC[topLine].top();

	for (int i = topLine; i < textLinesVEC.size(); i++) {

		//found line
		if (Y <= textLinesVEC[i].bottom()) return i;
	}

	return -1;
}
INT2 TextLines::GetTextObjectIdx(INT2 relCoord, int topLine) {

	if (!GoodIdx(LastLine(), topLine)) return INT2(-1, -1);

	float Y = (float)relCoord.j + textLinesVEC[topLine].top();

	for (int i = topLine; i < textLinesVEC.size(); i++) {

		//found line
		if (Y <= textLinesVEC[i].bottom()) {

			//find text object index on this line for these coordinates (if any)
			int j = textLinesVEC[i].get_textobject_index(relCoord);

			return INT2(i, j);
		}
	}

	return INT2(-1, -1);
}

INT2 TextLines::GetCharacterIdx(INT2 relCoord, int topLine) {

	if (!GoodIdx(LastLine(), topLine)) return INT2(-1, -1);

	float Y = (float)relCoord.j + textLinesVEC[topLine].top();

	for (int i = topLine; i < textLinesVEC.size(); i++) {

		//found line
		if (Y <= textLinesVEC[i].bottom()) {

			//find character index on this line for these coordinates (if any)
			int j = textLinesVEC[i].get_char_index(relCoord);

			return INT2(i, j);
		}
	}

	return INT2(-1, -1);
}

INT3 TextLines::GetFullIdx(INT2 relCoord, int topLine) {

	if (!GoodIdx(LastLine(), topLine)) return INT3(-1, -1, -1);

	float Y = (float)relCoord.j + textLinesVEC[topLine].top();

	for (int i = topLine; i < textLinesVEC.size(); i++) {

		//found line
		if (Y <= textLinesVEC[i].bottom()) {

			//find character index on this line for these coordinates (if any)
			INT2 jk = textLinesVEC[i].get_char_and_textobject_index(relCoord);

			//line, textobject, character
			return INT3(i, jk.i, jk.j);
		}
	}

	return INT3(-1, -1, -1);
}

int TextLines::para_length(int lineIdx) {

	int length = 0;

	for (int idx = GetSplitStartIndex(lineIdx); idx <= GetSplitEndIndex(lineIdx); idx++)
		length += textLinesVEC[idx].length();

	return length;
}

//get smallest rectangle encompassing all text lines
D2D1_RECT_F TextLines::get_rect(void)
{
	if (!textLinesVEC.size()) return  D2D1::RectF(0, 0, 0, 0);

	//first line rect
	D2D1_RECT_F lines_rect = textLinesVEC[0].get_rect();

	//check all other lines
	for (int idx = 1; idx < textLinesVEC.size(); idx++) {

		D2D1_RECT_F line_rect = textLinesVEC[idx].get_rect();

		if (line_rect.left < lines_rect.left) lines_rect.left = line_rect.left;
		if (line_rect.top < lines_rect.top) lines_rect.top = line_rect.top;
		if (line_rect.right > lines_rect.right) lines_rect.right = line_rect.right;
		if (line_rect.bottom > lines_rect.bottom) lines_rect.bottom = line_rect.bottom;
	}

	return lines_rect;
}

//Display new line using text formatting: 
//<b> </b> : bold

//<i> </i> : italic

//[on], [os], [or] : text outline (TOO_NONE, TOO_SQUARE, TOO_ROUND)

//[tcr,g,b,a/tc] : text color specified using r,g,b,a

//[bcr,g,b,a/bc] : background color

//[iomajorId,minorId,auxId,textId/io] : interactive object settings: first two numbers identify the major and minor TextObject type (used to identify the type of action generated when interacting with this TextObject). Next is an auxiliary id number. Last field contains a text identifier.
//</io> use this to specify interactive object end

//[anumber/a]	: text alignment object. This is a special type of interactive text object, whose state is checked internally (in CheckLine_InteractiveObjectsState) - not by the state handler method, and has no action assigned to it. 
//It is used to align text across multiple neighboring lines. This is the first text aligner, see next. The number must be unique on this text line. Number is the majorId in the interactive object.
//[sanumber/sa] : text alignment object as above. This object will attempt to synchronise with another text aligner object of same number on the line above, and if that fails on the line below it.
//Text alignment objects hold a vector of smart pointers to all other text alignment objects to which it is synchronised. 
//They consist simply of spaces and the number of spaces is modified in the internal state handler (at display update time) so that they all finish at the same character number on each of their respective lines, whilst keeping the number of spaces to a minimum. 
//This means the text that comes after a text alignment object always starts at the same character index across the text lines holding these synchronised objects (so looks aligned in the console if using a monospaced font). 
//Can use multiple text alignment objects on the same line, just make sure you number them differently.

//</c> clear formatting back to default

//in addition to text also pass the intended lineIndex where this must be placed. This is used by synchronising text alignment objects - if a new stao is being made, it will search for taos to synchronise to at lineIndex - 1.
TextLine TextLines::BuildFormattedTextLine(std::string text, int lineIndex)
{
	TextLine newLine;

	//Defined in TextFormatting.h
	std::vector<std::string>& formattingVector = TF::textformattingVector;

	//Default settings : modified by formatting
	FormatSpecifier fs;

	SimTOFunct *pActionHandler_ = nullptr;

	std::string remainingText;

	while (true) {

		if (!text.length()) break;

		//search for the first occurence of a format specifier
		int formCode = SplitOnFormatSpecifier(text, remainingText, formattingVector);

		if (text.length()) newLine.push(TextObject(text, fs, pActionHandler_));

		//If no remaining text, then there's nothing left to do
		if (!remainingText.length()) break;

		//change formatting depending on the format specifier found. There may be multiple specifiers one after another.
		while (true) {

			//make changes to formatting if any required
			switch (formCode) {

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
				std::string save_text = remainingText;	//save remaining text in case the passed text does not have the correct formatting sequence
				text = remainingText;
				formCode = SplitOnFormatSpecifier(text, remainingText, formattingVector);
				if (formCode == TF_TCEND) {

					std::vector<std::string> rgba = split(text, ",");
					if (rgba.size() == 4) {
						D3DCOLORVALUE newColor = { ToNum(rgba[0]),ToNum(rgba[1]),ToNum(rgba[2]),ToNum(rgba[3]) };
						fs.textColor = newColor;
					}
					else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
				}
				else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
			}
			break;

			case TF_BCSTART:
			{
				std::string save_text = remainingText;	//save remaining text in case the passed text does not have the correct formatting sequence
				text = remainingText;
				formCode = SplitOnFormatSpecifier(text, remainingText, formattingVector);
				if (formCode == TF_BCEND) {

					std::vector<std::string> rgba = split(text, ",");
					if (rgba.size() == 4) {
						D3DCOLORVALUE newColor = { ToNum(rgba[0]),ToNum(rgba[1]),ToNum(rgba[2]),ToNum(rgba[3]) };
						fs.bgrndColor = newColor;
					}
					else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
				}
				else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
			}
			break;

			case TF_SETACTIONSPECSTART:
			{
				std::string save_text = remainingText;	//save remaining text in case the passed text does not have the correct formatting sequence
				text = remainingText;
				formCode = SplitOnFormatSpecifier(text, remainingText, formattingVector);
				if (formCode == TF_SETACTIONSPECEND) {

					std::vector<std::string> ids = split(text, ",");
					if (ids.size() >= 4) {

						pActionHandler_ = pActionHandler;
						if (pActionHandler_) {

							pActionHandler_->set_iop(InteractiveObjectProperties(ToNum(ids[0]), ToNum(ids[1]), ToNum(ids[2]), combine(subvec(ids, 3), ",")));
						}
					}
					else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
				}
				else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
			}
			break;

			case TF_TAOSIMPLESTART:
			{
				std::string save_text = remainingText;	//save remaining text in case the passed text does not have the correct formatting sequence
				text = remainingText;
				formCode = SplitOnFormatSpecifier(text, remainingText, formattingVector);
				if (formCode == TF_TAOSIMPLEEND) {

					std::vector<std::string> ids = split(text, ",");
					if (ids.size() >= 1) {

						pActionHandler_ = pActionHandler;
						if (pActionHandler_) {

							int tao_number = ToNum(ids[0]);
							pActionHandler_->set_iop(InteractiveObjectProperties(tao_number));

							//Now immediately store this text alignment object with blank text and reset pActionHandler_. After this continue - there will probably be another specifier right after.
							newLine.push(TextObject("", fs, pActionHandler_, &newLine));
							pActionHandler_ = nullptr;
						}
					}
					else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
				}
				else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
			}
			break;

			case TF_TAOSYNCHROSTART:
			{
				std::string save_text = remainingText;	//save remaining text in case the passed text does not have the correct formatting sequence
				text = remainingText;
				formCode = SplitOnFormatSpecifier(text, remainingText, formattingVector);
				if (formCode == TF_TAOSYNCHROEND) {

					std::vector<std::string> ids = split(text, ",");
					if (ids.size() >= 1) {

						pActionHandler_ = pActionHandler;
						if (pActionHandler_) {

							int tao_number = ToNum(ids[0]);
							bool synchronised = false;

							pActionHandler_->set_iop(InteractiveObjectProperties(tao_number));

							//search for tao to synchronise to on line above
							if (GoodIdx(LastLine(), lineIndex - 1)) {

								for (int to_idx = 0; to_idx < textLinesVEC[lineIndex - 1].size(); to_idx++) {

									if (textLinesVEC[lineIndex - 1][to_idx].is_text_alignment_object() &&
										textLinesVEC[lineIndex - 1][to_idx].get_text_alignment_object_number() == tao_number) {

										//found a tao to synchronise to - pass it to the TextObject constructor
										newLine.push(TextObject("", fs, pActionHandler_, &newLine, &(textLinesVEC[lineIndex - 1][to_idx])));
										pActionHandler_ = nullptr;
										synchronised = true;
										break;		//stop at first one found
									}
								}
							}

							if (!synchronised) {

								//search for tao to synchronise to on line below if the above failed
								if (GoodIdx(LastLine(), lineIndex + 1)) {

									for (int to_idx = 0; to_idx < textLinesVEC[lineIndex + 1].size(); to_idx++) {

										if (textLinesVEC[lineIndex + 1][to_idx].is_text_alignment_object() &&
											textLinesVEC[lineIndex + 1][to_idx].get_text_alignment_object_number() == tao_number) {

											//found a tao to synchronise to - pass it to the TextObject constructor
											newLine.push(TextObject("", fs, pActionHandler_, &newLine, &(textLinesVEC[lineIndex + 1][to_idx])));
											pActionHandler_ = nullptr;
											synchronised = true;
											break;		//stop at first one found
										}
									}
								}
							}

							//couldn't synchronise to anything so just make a simple tao
							if (!synchronised) {

								//If we are here then couldn't find a tao to synchronise to - start a new tao list
								//Now immediately store this text alignment object with blank text and reset pActionHandler_. After this continue - there will probably be another specifier right after.
								newLine.push(TextObject("", fs, pActionHandler_, &newLine));
								pActionHandler_ = nullptr;
							}
						}
					}
					else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
				}
				else { remainingText = save_text; }	//incorrect formatting sequence - restore and try next format specifier
			}
			break;

			case TF_SETACTIONEND:
				pActionHandler_ = nullptr;
				break;

			case TF_CLEARFORMATTING:
				fs = FormatSpecifier();
				break;

			default:
				break;
			}

			//check to see if other formatting specifiers follow immediately.
			text = remainingText;
			if (!text.length()) break;
			std::string remainingText_;
			formCode = SplitOnFormatSpecifier(text, remainingText_, formattingVector);

			//if text contains something then no immediate specifiers found, otherwise take another iteration to process formCode and check for further specifiers.
			if (text.length() > 0) break;
			else remainingText = remainingText_;
		}

		//more text to process
		text = remainingText;
	}

	return newLine;
}

#endif