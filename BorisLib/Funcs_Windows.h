#pragma once

// Windows Header Files:
#include <Windows.h>

#include <string>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline std::string GetClipboardText(void) 
{	
	if(!OpenClipboard(nullptr)) return "";

	HANDLE hData = GetClipboardData(CF_TEXT);
	if(!hData) {
		
		CloseClipboard();
		return "";
	}

	LPVOID lpvptr = GlobalLock(hData);
	if (!lpvptr) {

		GlobalUnlock(hData);
		CloseClipboard();
		return "";
	}

	std::string text(static_cast<char*>(lpvptr));

	GlobalUnlock(hData);		
	CloseClipboard();

	return text;
}

inline void SetClipboardText(const std::string& text) {

	if(!OpenClipboard(nullptr) || !text.length()) return;

	EmptyClipboard();

	HGLOBAL hg = GlobalAlloc(GMEM_MOVEABLE, text.size() + 1);
	
	if (!hg) {

		CloseClipboard();
		return;
	}

	memcpy(GlobalLock(hg), text.c_str(), text.size() + 1);
	
	GlobalUnlock(hg);
	
	SetClipboardData(CF_TEXT, hg);
	
	CloseClipboard();
	GlobalFree(hg);
}