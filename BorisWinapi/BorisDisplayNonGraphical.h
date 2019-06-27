#pragma once

//Contains BorisDisplay when running in non-graphics mode (basic text console only)

#include "CompileFlags.h"
#if GRAPHICS == 0

#include <string>
#include <vector>
#include <iostream>

#include "TextFormatting.h"
#include "BorisLib.h"

#include <SFML/System.hpp>
#include <SFML/Graphics/Texture.hpp>

#pragma comment(lib, "sfml-system.lib")
#pragma comment(lib, "sfml-graphics.lib")

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Basic functionality in console-only mode. Portable. Still need to shape meshes using mask files in console-only mode.

class BorisGraphics
{

public:

	BorisGraphics(void) {}

	//from image file extract a raw bitmap (BYTE array with 4 BYTE-sized entries as B-G-R-A for each pixel) to specified pixel size (so image file is rescaled to specified size)
	void GetBitmapFromImage(string fileName, vector<unsigned char>& bitmap, INT2 n_plane);

	bool MakeVideoFromFileSequence(string directory, vector<string>& fileNames, unsigned int fps, double scaling, int quality)
	{
		cout << "Not available in console-only mode" << endl;

		return false;
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Boris Display container and director object. All access to display elements and drawing done through here.

class BorisDisplay
{
private: //Private data

	BorisGraphics *pBG = nullptr;

public:  //Public data

private: //Private methods	

public:

	BorisDisplay() 
	{
		pBG = new BorisGraphics();
	}

	~BorisDisplay() 
	{
		delete pBG;
	}

	//-----------------------------------------------Thread-safe access points to BorisDisplay

	void ClearDataBox(void) {}

	//display various message types in the console
	void DisplayConsoleMessage(string text) { cout << text << endl; }
	void DisplayConsoleError(string text) { cout << text << endl; }
	void DisplayConsoleListing(string text) { cout << text << endl; }

	//display message after stripping out the formatting
	void DisplayFormattedConsoleMessage(string text);

	//add a new entry in the data box as an interative object (use formatted text)
	void NewDataBoxField(string formattedText) {}
	//update the data box at given lineIdx to the given value, passed as a string
	void UpdateDataBoxField(int lineIdx, string value_string) {}

	//set the console line entry text
	void SetConsoleEntryLineText(string text) { cout << text; }

	bool SaveMeshImage(string fileName) 
	{ 
		cout << "Not available in console-only mode" << endl; 
		return false;
	}

	//allows access to Boris Graphics public methods.
	BorisGraphics* BGMethods(void) { return pBG; }
};

#endif