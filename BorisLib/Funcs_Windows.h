#pragma once

// Windows Header Files:
#include <Windows.h>

#include <shlobj.h>
#include <atlstr.h>
#include <strsafe.h>
#include <shellapi.h>

#include <string>

#include "Funcs_Files.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline std::string GetUserDocumentsPath(void)
{
	wchar_t my_documents[MAX_PATH];

	HRESULT result = SHGetFolderPath(NULL, CSIDL_PERSONAL, NULL, SHGFP_TYPE_CURRENT, my_documents);

	if (result != S_OK) {

		return GetDirectory();
	}
	else {
		
		return WideStringtoString(std::wstring(my_documents)) + "\\";
	}
}

inline bool MakeDirectory(std::string directory)
{
	if (SHCreateDirectoryEx(NULL, StringtoWCHARPointer(directory), NULL) == ERROR_SUCCESS || ERROR_ALREADY_EXISTS == GetLastError()) return true;
	else return false;
}

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