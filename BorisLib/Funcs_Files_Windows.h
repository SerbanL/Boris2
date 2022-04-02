#pragma once

#include "BorisLib_Config.h"
#if OPERATING_SYSTEM == OS_WIN

//Windows Header Files:
#include <atlstr.h>
#include <strsafe.h>
#include <shellapi.h>
#include <Windows.h>
#include <shlobj.h>

//Others
#include <string>
#include <algorithm>

#include "Funcs_Files.h"
#include "Funcs_Conv_Windows.h"

#ifndef FILEROWCHARS
#define FILEROWCHARS	50000	//maximum number of characters per input file row
#endif

//Get current directory (initially this will be the executable directory but can change if you load/save files through a windows dialog)
inline std::string GetDirectory(void)
{
	TCHAR path[FILEROWCHARS];
	GetCurrentDirectory(FILEROWCHARS, path);
	std::string directory = std::string(CW2A(path)) + "/";

	return directory;
}

//Get executable file directory
inline std::string GetExeDirectory(void)
{
	TCHAR path[FILEROWCHARS];
	GetModuleFileName(NULL, path, FILEROWCHARS);
	std::string fullPath = std::string(CW2A(path));
	std::string directory = fullPath.substr(0, fullPath.find_last_of("\\/")) + "/";

	return directory;
}

//Get executable file name
inline std::string GetExeFilename(void)
{
	TCHAR path[FILEROWCHARS];
	GetModuleFileName(NULL, path, FILEROWCHARS);
	std::string fullPath = std::string(CW2A(path));
	
	std::string fileName;
	
	size_t pos = fullPath.find_last_of("\\/");
	if (pos != std::string::npos) fileName = fullPath.substr(pos + 1);

	return fileName;
}

//return all files sharing the given base (in specified directory) and termination. Return them ordered (including full path) by creation time
inline std::vector<std::string> GetFilesInDirectory(std::string directory, std::string baseFileName, std::string termination)
{
	std::vector<std::pair<std::string, double>> fileNames_creationTimes;

	WIN32_FIND_DATA ffd;
	TCHAR szDir[MAX_PATH];
	HANDLE hFind = INVALID_HANDLE_VALUE;

	std::wstring wstr_directory = StringtoWideString(directory);

	StringCchCopy(szDir, MAX_PATH, wstr_directory.c_str());
	StringCchCat(szDir, MAX_PATH, TEXT("\\*"));

	// Find the first file in the directory.
	hFind = FindFirstFile(szDir, &ffd);

	if (INVALID_HANDLE_VALUE == hFind) return {};

	// Get all the files in the directory.
	do {

		//only get file names, not directories
		if (!(ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {

			//creation time : form a double as uli.QuadPart
			ULARGE_INTEGER uli;
			uli.LowPart = ffd.ftCreationTime.dwLowDateTime;
			uli.HighPart = ffd.ftCreationTime.dwHighDateTime;

			//get the file name
			std::string fileName = WideStringtoString(std::wstring(ffd.cFileName));

			//if filename including termination is too short then skip it
			if (fileName.length() < termination.length()) continue;

			if (!baseFileName.length()) {

				if (!termination.length()) {

					//if no basefile name or termination is specified then just get the file
					fileNames_creationTimes.push_back({ directory + fileName, (double)uli.QuadPart });
				}
				else if (fileName.substr(fileName.length() - termination.length()).compare(termination) == 0) {

					//if no basefile name is specified but the termination is, then it must match that of the file (the ending of the file name)
					fileNames_creationTimes.push_back({ directory + fileName, (double)uli.QuadPart });
				}
			}
			else if (fileName.substr(0, baseFileName.length()).compare(baseFileName) == 0) {

				//if basefile name is specified then it must match that of the file (the beggining of the file name)
				if (!termination.length()) {

					fileNames_creationTimes.push_back({ directory + fileName, (double)uli.QuadPart });
				}
				else if (fileName.substr(fileName.length() - termination.length()).compare(termination) == 0) {

					fileNames_creationTimes.push_back({ directory + fileName, (double)uli.QuadPart });
				}
			}
		}
	} while (FindNextFile(hFind, &ffd) != 0);			//get next file in directory

	//finish
	FindClose(hFind);

	auto compare = [&](const std::pair<std::string, double>& first, const std::pair<std::string, double>& second) -> bool { return first.second < second.second; };

	//sort by creation time order
	std::sort(fileNames_creationTimes.begin(), fileNames_creationTimes.end(), compare);

	//extract and return fileNames only
	std::vector<std::string> fileNames(fileNames_creationTimes.size());
	std::transform(fileNames_creationTimes.begin(), fileNames_creationTimes.end(), fileNames.begin(), [](auto const& pair) { return pair.first; });
	return fileNames;
}

//open a file programatically : the file name includes a path
inline void open_file(std::string fileName, std::string parameters = "")
{
	std::string command = "open";

	if (!parameters.length())
		ShellExecute(GetDesktopWindow(), StringtoWCHARPointer(command), StringtoWCHARPointer(fileName), NULL, NULL, SW_SHOWNORMAL);
	else
		ShellExecute(GetDesktopWindow(), StringtoWCHARPointer(command), StringtoWCHARPointer(fileName), StringtoWCHARPointer(parameters), NULL, SW_SHOWNORMAL);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline std::string GetUserDocumentsPath(void)
{
	wchar_t my_documents[MAX_PATH];

	HRESULT result = SHGetFolderPath(NULL, CSIDL_PERSONAL, NULL, SHGFP_TYPE_CURRENT, my_documents);

	if (result != S_OK) {

		return GetDirectory();
	}
	else {

		std::string path = WideStringtoString(std::wstring(my_documents)) + "/";
		std::for_each(path.begin(), path.end(), [](char& c){ if (c == '\\') c = '/'; });

		return path;
	}
}

inline bool MakeDirectory(std::string directory)
{
	if (SHCreateDirectoryEx(NULL, StringtoWCHARPointer(directory), NULL) == ERROR_SUCCESS || ERROR_ALREADY_EXISTS == GetLastError()) return true;
	else return false;
}

inline std::string GetClipboardText(void)
{
	if (!OpenClipboard(nullptr)) return "";

	HANDLE hData = GetClipboardData(CF_TEXT);
	if (!hData) {

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

	if (!OpenClipboard(nullptr) || !text.length()) return;

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

#endif