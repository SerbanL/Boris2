#pragma once

#include <string>
#include <vector>
#include <cctype>
#include <algorithm>

#include "Funcs_Vectors.h"

///////////////////////////////////////////////////////////////////////////////
//
//	Case change

//from text return a new std::string with lowercase only
inline std::string lowercase(const std::string& text)
{
	std::string text_lowercase = text;

	std::transform(text_lowercase.begin(), text_lowercase.end(), text_lowercase.begin(), [](unsigned char c)
	{
		return std::tolower(c);
	});

	return text_lowercase;
}

//from text return a new std::string with uppercase only
inline std::string uppercase(const std::string& text)
{
	std::string text_uppercase = text;

	std::transform(text_uppercase.begin(), text_uppercase.end(), text_uppercase.begin(), [](unsigned char c)
	{
		return std::toupper(c);
	});

	return text_uppercase;
}

///////////////////////////////////////////////////////////////////////////////
//
//	Trimming

//remove all match substring - variants for a single match std::string (and in particular a dedicated function for spaces)
//return trimmed std::string (but leave original intact)
inline std::string trim(std::string text, const std::string& match)
{
	//trim can be called only as trim(text), where match = "" by default - this is used in the parameter pack version when recursing.
	if (!match.length()) return text;

	size_t pos = 0;
	while ((pos = text.find(match, pos)) != std::string::npos) {

		text.erase(pos, match.length());
	}

	return text;
}

//parameter pack version to allow multiple match strings
template <typename... T> std::string trim(std::string text, const std::string& match, const T&... match_further)
{
	//call the original trim function (compiler prioritises calling overloaded function rather than recursing)
	text = trim(text, match);
	//if match_further... expands to 2 or more arguments then code will recurse on this function, so the next match std::string will be operated on
	//if match_further... expands to fewer than 2 arguments the overloaded version is called
	text = trim(text, match_further...);
	return text;
}

//remove any text from input found between the start and end strings, including these matches (both matches must be found)
inline std::string trimblock(std::string text, const std::string& start, const std::string& end)
{
	if (!start.length() || !end.length()) return text;

	size_t pos_start = 0;
	size_t pos_end = 0;

	while (pos_start < text.length()) {

		pos_start = text.find(start, pos_start);

		if (pos_start == std::string::npos) break;

		//found start match, now get end match

		pos_end = text.find(end, pos_start);

		if (pos_end == std::string::npos) break;

		//found both start and end matches : erase
		text.erase(pos_start, pos_end - pos_start + end.length());
	}

	return text;
}

inline std::string trimspaces(std::string text) { return trim(text, " "); }

inline std::string trimendspaces(std::string text)
{
	if (text.length()) {

		int s = 0, e = text.length();

		if (text[0] == ' ') s++;
		if (text.back() == ' ') e--;

		if (e - s > 0) return text.substr(s, e - s);
		else return "";
	}
	else return text;
}

//remove all leading spaces from text
inline std::string trim_leading_spaces(const std::string& text)
{
	int start_idx = 0;

	for (int idx = 0; idx < text.length(); idx++) {

		if (text[idx] == ' ') start_idx++;
		else break;
	}

	return text.substr(start_idx);
}

//remove all trailing spaces from text
inline std::string trim_trailing_spaces(const std::string& text)
{
	int end_idx = text.length() - 1;

	for (int idx = text.length() - 1; idx >= 0; idx--) {

		if (text[idx] == ' ') end_idx--;
		else break;
	}

	return text.substr(0, end_idx + 1);
}

///////////////////////////////////////////////////////////////////////////////
//
//	Splitting

//split input std::string by using separators stored in the std::vector. Split on first separator found and continue with the rest of the std::string. Do not include separators in split text.
inline std::vector<std::string> split(const std::string& textstring, const std::vector<std::string>& separators)
{
	std::vector<std::string> stringsVector;

	if (!separators.size()) { stringsVector.push_back(textstring); return stringsVector; }

	if (!textstring.length()) return stringsVector;

	size_t pos1 = 0, pos2 = 0;

	while (pos2 < (int)textstring.length()) {

		pos2 = std::string::npos;

		//find next separator from list
		int vecIdx = (int)separators.size();
		for (int i = 0; i < (int)separators.size(); i++) {

			size_t pos = textstring.find(separators[i], pos1);

			if (pos2 == std::string::npos) {

				pos2 = pos;
				if (pos != std::string::npos) vecIdx = i;
			}
			//if a separator was found then pos2 now either has this separator position or a previous one. copy this position to pos2 if it is closer to it.
			if (pos != std::string::npos && pos < pos2) { pos2 = pos; vecIdx = i; }
		}

		//if no separator found then get the remaining part of the std::string and finish
		if (vecIdx == (int)separators.size()) {

			stringsVector.push_back(textstring.substr(pos1, std::string::npos));
			break;
		}
		else {

			//found a separator, store substring so far
			stringsVector.push_back(textstring.substr(pos1, pos2 - pos1));
			//prepare for next substring - skip over separator char
			pos2 += separators[vecIdx].length();
			pos1 = pos2;
		}
	}

	return stringsVector;
}

//split input std::string by using separators in a parameter pack. Split on first separator found and continue with the rest of the std::string. Do not include separators in split text.
template <typename... T>
std::vector<std::string> split(const std::string& textstring, const T& ... separators)
{
	return split(textstring, make_vector((std::string)separators...));
}

//split using space
inline std::vector<std::string> split(const std::string& textstring)
{
	return split(textstring, {" "});
}

//split input std::string into substrings using space as a separator but only when the spaces delimit a maximally numeric substring.
//A numeric substring is a std::string which contains only digits and spaces. A maximally numeric substring is a numeric substring which cannot be extended by adding characters either to the left or to the right.
//e.g. "c:\file path\file name 2.txt 0 1 2 text 3 4 5" will be split into "c:\file path\file name 2.txt", "0 1 2", "text", "3 4 5"
inline std::vector<std::string> split_numeric(const std::string& textstring)
{
	std::vector<std::string> stringsVector;

	if (!textstring.length()) return stringsVector;

	std::string numeric_character = "0123456789eE-+.";

	int start_text_idx = 0;
	bool isnumeric = false;

	for (int idx = 0; idx < textstring.length(); idx++) {

		isnumeric = false;

		//check for start of a potentially numeric substring : must have a space followed by a number, or be the first index
		if (textstring[idx] == ' ' || idx == 0) {

			//if first idx then check for start of a potentially numeric substring, else increment idx (was space at idx) and check if we finished textstring
			if (idx != 0 && ++idx >= textstring.length()) break;

			if (numeric_character.find(textstring[idx]) != std::string::npos) {

				//found start of a potentially numeric substring
				int start_idx = idx;
				int end_idx = start_idx + 1;

				//reached end : the last std::string was numeric
				if (end_idx == textstring.length()) isnumeric = true;

				//is it a numeric substring? There must be a space after a digit to be a real numeric substring
				for (end_idx = start_idx + 1; end_idx < textstring.length(); end_idx++) {

					//is it the last character?
					if (end_idx == textstring.length() - 1) {

						//if a number then we have a numeric substring
						if (numeric_character.find(textstring[end_idx]) != std::string::npos) isnumeric = true;

						//if a space then we also have a numeric substring but end_idx must not include it
						if (textstring[end_idx] == ' ') { isnumeric = true; break; }
					}

					//skip over digits
					if (numeric_character.find(textstring[end_idx]) != std::string::npos) continue;

					//found a space? then it's a numeric substring
					if (textstring[end_idx] == ' ') {

						isnumeric = true;
						continue;
					}

					//if we get here then the character is not part of a numeric substring : finish current search
					break;
				}

				//did the search identify a numeric substring? If so separate it out...
				if (isnumeric) {

					//substring starting at start_text_idx with length start_idx - start_text_idx is a normal text substring
					if (start_idx != start_text_idx)
						stringsVector.push_back(trim_trailing_spaces(textstring.substr(start_text_idx, start_idx - start_text_idx)));

					//substring starting at start_idx with length end_idx - start_idx is a numeric substring
					stringsVector.push_back(trim_trailing_spaces(textstring.substr(start_idx, end_idx - start_idx)));

					//start of new normal text substring
					start_text_idx = end_idx;
				}

				//and continue search for numeric substrings
				idx = end_idx;
			}
		}
	}

	if (!isnumeric) stringsVector.push_back(trim_trailing_spaces(textstring.substr(start_text_idx)));

	return stringsVector;
}

//split input std::string in two if it contains a substring from formatSpecifiers, keeping the part to the left of the match in the input std::string, and the part to the right of the match in the output std::string (the matched std::string is not kept). Return the index to the matched std::string.
inline int SplitOnFormatSpecifier(std::string &stringIN, std::string &stringOUT, const std::vector<std::string> &formatSpecifiers)
{
	size_t pos = std::string::npos, first = std::string::npos;

	int fsIdx = -1;

	for (int i = 0; i < formatSpecifiers.size(); i++) {

		pos = stringIN.find(formatSpecifiers[i]);
		if (first == std::string::npos) { first = pos; fsIdx = i; }
		if (pos != std::string::npos && pos < first) { first = pos; fsIdx = i; }
	}

	if (first != std::string::npos) {

		stringOUT = stringIN.substr(first + formatSpecifiers[fsIdx].length());
		stringIN = stringIN.substr(0, first);
	}
	else stringOUT = "";

	return fsIdx;
}

///////////////////////////////////////////////////////////////////////////////
//
//	Joining

//combine strings in std::vector into a single std::string, separating them using the separator
inline std::string combine(const std::vector<std::string> &textComponents, const std::string& separator)
{
	std::string text;

	for (int idx = 0; idx < (int)textComponents.size(); idx++) {

		text += textComponents[idx];

		if (idx != (int)textComponents.size() - 1) text += separator;
	}

	return text;
}

///////////////////////////////////////////////////////////////////////////////
//
//	Substitutions / Erasures

//replace all match substrings in text by replace std::string - return true if any replacements made
inline bool replaceall(std::string &text, const std::string& match, const std::string& replace)
{
	bool replaced = false;

	size_t pos = 0;

	while (pos < text.length()) {

		pos = text.find(match, pos);

		if (pos == std::string::npos) break;
		else replaced = true;

		text.replace(pos, match.length(), replace);
		pos += replace.length();
	}

	return replaced;
}

//get first substring contained between start and end substrings, also deleting it from the original std::string (including the start and end)
inline std::string remove_first_contained_substring(std::string& text, const std::string& start, const std::string& end)
{
	size_t pos_start = text.find(start) + start.length();
	size_t pos_end = text.find(end);

	if (pos_start != std::string::npos && pos_end != std::string::npos && pos_end >= pos_start) {

		std::string match = text.substr(pos_start, pos_end - pos_start);
		text.erase(pos_start - start.length(), pos_end - pos_start + start.length() + end.length());
		return match;
	}

	return "";
}

//get last substring contained between start and end substrings, also deleting it from the original std::string (including the start and end)
inline std::string remove_last_contained_substring(std::string& text, const std::string& start, const std::string& end)
{
	std::string text_reverse = text;
	std::string start_reverse = start;
	std::string end_reverse = end;

	std::reverse(text_reverse.begin(), text_reverse.end());
	std::reverse(start_reverse.begin(), start_reverse.end());
	std::reverse(end_reverse.begin(), end_reverse.end());

	size_t pos_end = text_reverse.find(end_reverse) + end.length();
	size_t pos_start = text_reverse.find(start_reverse);

	if (pos_start != std::string::npos && pos_end != std::string::npos && pos_start >= pos_end) {

		pos_end = text.length() - pos_end;
		pos_start = text.length() - pos_start;

		std::string match = text.substr(pos_start, pos_end - pos_start);
		text.erase(pos_start - start.length(), pos_end - pos_start + start.length() + end.length());
		return match;
	}

	return "";
}

///////////////////////////////////////////////////////////////////////////////
//
//	Special

//from the right std::string shift the first substring up to and including the separator character (if any) in leftString.
//Return true if a shift occured, false otherwise (no separator found, not even at the end -> no shifts occur in this case)
inline bool ShiftSubstring_R2L(std::string &leftString, std::string &rightString, char separator)
{
	size_t pos = rightString.find_first_of(separator);

	if (pos != std::string::npos) {

		//found separator	
		leftString += rightString.substr(0, pos + 1); //shift substring, including separator
		rightString = rightString.substr(pos + 1);

		return true;
	}

	return false;
}

//As above but left to right : if there is a separator at the end of leftString, shift substring from next-to-last separator (but not including it) to rightString.
//If there is no separator at the end then shift from last separator (but not including it) to rightString. Return false if no separator found at all -> no shifts made in this case
//NOTE: Typically ShiftSubstring_L2R is meant to be used after a sequence of ShiftSubstring_R2L calls, so there should be a separator at the end.
inline bool ShiftSubstring_L2R(std::string &leftString, std::string &rightString, char separator)
{
	size_t pos = leftString.find_last_of(separator);

	if (pos != std::string::npos) {

		//if separator is last, then search for next-to-last separator
		if (pos == leftString.length() - 1) pos = leftString.substr(0, pos).find_last_of(separator);
		//if that was the only separator found will shift everything
		if (pos == std::string::npos) pos = 0;
		else pos++;

		//shift substring, but not including the separator at pos
		rightString = leftString.substr(pos) + rightString;

		//remove substring from left side, leaving the separator at pos at the end
		leftString = leftString.substr(0, pos);

		return true;
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////
//
//	Getters

//return std::string before match std::string, excluding the match (return empty if not match)
inline std::string get_before_match(const std::string& text, const std::string& match)
{
	size_t pos_match = text.find(match);

	if (pos_match != std::string::npos) {

		return text.substr(0, pos_match);
	}
	else return "";
}

//return std::string after match std::string, excluding the match (return empty if not match)
inline std::string get_after_match(const std::string& text, const std::string& match)
{
	size_t pos_match = text.find(match);

	if (pos_match != std::string::npos && pos_match + match.length() <= text.length()) {

		return text.substr(pos_match + match.length());
	}
	else return "";
}

//get first substring contained between start and end substrings
inline std::string get_first_contained_substring(const std::string& text, const std::string& start, const std::string& end)
{
	size_t pos_start = text.find(start) + start.length();
	size_t pos_end = text.find(end);

	if (pos_start != std::string::npos && pos_end != std::string::npos && pos_end >= pos_start)
		return text.substr(pos_start, pos_end - pos_start);

	return "";
}

//get last substring contained between start and end substrings
inline std::string get_last_contained_substring(const std::string& text, std::string& start, const std::string& end)
{
	std::string text_reverse = text;
	std::string start_reverse = start;
	std::string end_reverse = end;

	std::reverse(text_reverse.begin(), text_reverse.end());
	std::reverse(start_reverse.begin(), start_reverse.end());
	std::reverse(end_reverse.begin(), end_reverse.end());

	size_t pos_end = text_reverse.find(end_reverse) + end.length();
	size_t pos_start = text_reverse.find(start_reverse);

	if (pos_start != std::string::npos && pos_end != std::string::npos && pos_start >= pos_end) {

		pos_end = text.length() - pos_end;
		pos_start = text.length() - pos_start;

		return text.substr(pos_start, pos_end - pos_start);
	}

	return "";
}

//get the wordIdx-th word from given std::string (if any)
inline std::string get_word(const std::string& text, int wordIdx = 0)
{
	size_t pos1 = 0, pos2 = 0;

	//skip any leading spaces
	while (pos2 < text.length() && text[pos2] == ' ')
		pos2++;

	for (int idx = 0; idx <= wordIdx; idx++) {

		if (pos2 >= text.length()) return "";

		//pos1 and pos2 are at the start of idx-th word
		pos1 = pos2;
		pos2 = text.find(" ", pos1);

		//if this is the last idx value then get this word (it could also be the last word in text, in which case pos2 is std::string::npos)
		if (idx == wordIdx) {

			return text.substr(pos1, pos2 - pos1);
		}

		//if space not found then index-th word could not be found - return empty
		if (pos2 == std::string::npos) return "";

		//skip over any spaces to look for next word
		while (pos2 < text.length() && text[pos2] == ' ')
			pos2++;
	}

	return "";
}

//get start and end indexes for word containing the given character index
inline std::pair<int, int> get_word_indexes(const std::string& text, int charIdx)
{
	std::pair<int, int> wordIndexes;

	if (charIdx < 0 || charIdx >(int)text.length()) return wordIndexes;

	wordIndexes = std::pair<int, int>(charIdx, 0);

	//first find right index
	size_t pos = text.find(" ", (size_t)charIdx);

	if (pos == std::string::npos) wordIndexes.second = (int)text.length();
	else wordIndexes.second = (int)pos - 1;

	//now find left index
	for (wordIndexes.first = charIdx; wordIndexes.first >= 0; wordIndexes.first--) {

		if (text[wordIndexes.first] == ' ') { wordIndexes.first++; break; }
	}

	if (wordIndexes.first < 0) wordIndexes.first = 0;

	return wordIndexes;
}

///////////////////////////////////////////////////////////////////////////////
//
//	Properties

//check the input std::string contains only numbers : non-empty and, apart from separators list of characters, must contain only alphanumeric characters, e, -, +, .
inline bool has_numbers_only(const std::string& text, std::string separators)
{
	if (!text.length()) return false;

	std::string allowed_characters = separators + "0123456789eE-+.";

	return (std::find_if(text.begin(), text.end(), [&](const char& c) { return (allowed_characters.find(c) == std::string::npos); }) == text.end());
}

//check the input std::string contains only digits
inline bool has_digits_only(const std::string& text)
{
	if (!text.length()) return false;

	return (std::find_if(text.begin(), text.end(), [&](const char& c) { return !isdigit(c); }) == text.end());
}