#include "stdafx.h"
#include "BorisDisplayNonGraphical.h"

#if GRAPHICS == 0

void BorisDisplay::DisplayFormattedConsoleMessage(std::string text)
{
	//first trim all formmatting specifiers in block form (i.e. with extra text info between start and end text formatting specifiers)
	for (int idx = 0; idx < TF::textblockSpecifiers.size(); idx++) {

		text = trimblock(text, TF::textblockSpecifiers[idx].first, TF::textblockSpecifiers[idx].second);
	}

	//now trim all other formatting specifiers
	for (int idx = 0; idx < TF::simpleFormattingSpecifiers.size(); idx++) {

		text = trim(text, TF::simpleFormattingSpecifiers[idx]);
	}

	//finally display the unformatted text
	std::cout << text << std::endl;
}

#endif