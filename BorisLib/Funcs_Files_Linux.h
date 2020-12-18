#pragma once

#include "BorisLib_Config.h"
#if OPERATING_SYSTEM == OS_LIN

#include <libgen.h>
#include <unistd.h>
#include <pwd.h>
#include <filesystem>
#include <limits.h>
#include <X11/Xlib.h>
#include "Funcs_Files.h"

#ifndef FILEROWCHARS
#define FILEROWCHARS	50000	//maximum number of characters per input file row
#endif

//Get current directory (initially this will be the executable directory but can change if you load/save files through a windows dialog : on Linux this is the same as GetExeDirectory)
inline std::string GetDirectory(void)
{
	char result[ FILEROWCHARS ];
  	ssize_t count = readlink( "/proc/self/exe", result, FILEROWCHARS );

	if (count != -1) return std::string(dirname(result)) + "/";
	else return "";
}

//Get executable file directory
inline std::string GetExeDirectory(void)
{
	char result[ FILEROWCHARS ];
  	ssize_t count = readlink( "/proc/self/exe", result, FILEROWCHARS );

	if (count != -1) return std::string(dirname(result)) + "/";
	else return "";
}

//return all files sharing the given base (in specified directory) and termination. Return them ordered (including full path) by creation time
inline std::vector<std::string> GetFilesInDirectory(std::string directory, std::string baseFileName, std::string termination)
{
	std::vector<std::string> fileNames;
	std::vector<double> creationTimes;

	//TO DO

	return fileNames;
}

//open a file programatically : the file name includes a path
inline void open_file(std::string fileName, std::string parameters = "")
{
	pid_t pid = fork();
    if (pid == 0) {

		if (!parameters.length())
			execl("/usr/bin/xdg-open", "xdg-open", fileName.c_str(), (char *)0);
		else
			execl("/usr/bin/xdg-open", "xdg-open", fileName.c_str(), parameters.c_str());
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline std::string GetUserDocumentsPath(void)
{
	struct passwd *pw;
    uid_t uid;
    int c;

    uid = geteuid();
    pw = getpwuid(uid);
	if (pw) {
		
		return std::string("/home/") + std::string(pw->pw_name) + std::string("/Documents/");
	}
	else return "";
}

inline bool MakeDirectory(std::string directory)
{
	//available in C++17
	std::filesystem::create_directories(directory);

	//TO DO: check if directory already exists and return true; for all other errors return false
	return true;
}

//This function has been adapted from : https://stackoverflow.com/questions/27378318/c-get-std::string-from-clipboard-on-linux
//Needs "-lX11" for linking
inline std::string GetClipboardText(void)
{
	auto PrintSelection = [](Display *display, Window window, const char *bufname, const char *fmtname, std::string& text) -> bool
	{
		char *result;
    	unsigned long ressize, restail;
    	int resbits;
    	Atom bufid = XInternAtom(display, bufname, false);
    	Atom fmtid = XInternAtom(display, fmtname, false);
    	Atom propid = XInternAtom(display, "XSEL_DATA", false);
    	Atom incrid = XInternAtom(display, "INCR", false);
       
    	XEvent event;
    
    	XSelectInput (display, window, PropertyChangeMask);
    	XConvertSelection(display, bufid, fmtid, propid, window, CurrentTime);
    
    	do {

        	XNextEvent(display, &event);
        
    	} while (event.type != SelectionNotify || event.xselection.selection != bufid);

  		if (event.xselection.property) {
      
      		XGetWindowProperty(display, window, propid, 0, LONG_MAX/4, True, AnyPropertyType, &fmtid, &resbits, &ressize, &restail, (unsigned char**)&result);
      
      		if (fmtid != incrid) {
          
          		text = std::string(result);
        		XFree(result);
        		return true;
      		}
      
	  		//can't get text in one go : need to do it incrementally
      		if (fmtid == incrid) {	
			  
        		do {
              
            		do {
            
                		XNextEvent(display, &event);
        
            		} while (event.type != PropertyNotify || event.xproperty.atom != propid || event.xproperty.state != PropertyNewValue);

            		XGetWindowProperty(display, window, propid, 0, LONG_MAX/4, True, AnyPropertyType, &fmtid, &resbits, &ressize, &restail, (unsigned char**)&result);
        
            		text += std::string(result);
            		XFree(result);
        
        		} while (ressize > 0);

        		return true;
        	}
    
    	}
    	else return false;
	};

	std::string text;
    Display *display = XOpenDisplay(NULL);
    unsigned long color = BlackPixel(display, DefaultScreen(display));
    Window window = XCreateSimpleWindow(display, DefaultRootWindow(display), 0,0, 1,1, 0, color, color);
    
    bool success = PrintSelection(display, window, "CLIPBOARD", "UTF8_STRING", text);
    if (!success) success = PrintSelection(display, window, "CLIPBOARD", "STRING", text);

    XDestroyWindow(display, window);
    XCloseDisplay(display);

    if (success) return text;
    else return "";
}

inline void SetClipboardText(const std::string& text) 
{
	//TO DO
}

#endif