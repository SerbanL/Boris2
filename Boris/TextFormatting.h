#pragma once

#include <vector>
#include <string>

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

//Text format specifiers - 1 to 1 correspondence with the format specifiers in formattingVector below
enum TF_ {
	TF_BOLDSTART = 0, TF_BOLDEND, TF_ITALICSTART, TF_ITALICEND,
	TF_ONONE, TF_OSQUARE, TF_OROUND, TF_TCSTART, TF_TCEND, TF_BCSTART, TF_BCEND,
	TF_SETACTIONSPECSTART, TF_SETACTIONSPECEND, TF_SETACTIONEND,
	TF_TAOSIMPLESTART, TF_TAOSIMPLEEND, TF_TAOSYNCHROSTART, TF_TAOSYNCHROEND,
	TF_CLEARFORMATTING
};

namespace TF {

	//all formatting strings in 1-2-1 correspondence with the TF_ enum
	static std::vector<std::string> textformattingVector = { 
		"<b>", "</b>", "<i>", "</i>", "[on]", "[os]", "[or]", 
		"[tc", "/tc]", "[bc", "/bc]", "[io", "/io]", 
		"</io>", 
		"[a", "/a]", 
		"[sa", "/sa]", 
		"</c>" 
	};

	//list all text formatting which have a start, end, and contain further information between them. 
	//e.g. "[tc", "/tc]" -> these contain the text color info between them
	//these block specifiers are useful e.g. if we want to strip all formatting from formatted text
	static std::vector<std::pair<std::string, std::string>> textblockSpecifiers = {

		{ "[tc", "/tc]" }, { "[bc", "/bc]"}, { "[io", "/io]" }, { "[a", "/a]"}, { "[sa", "/sa]"}
	};

	//all other formatting specifiers
	static std::vector<std::string> simpleFormattingSpecifiers = {

		"<b>", "</b>", "<i>", "</i>", "[on]", "[os]", "[or]", "</io>", "</c>"
	};
};
