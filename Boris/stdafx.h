// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "CompileFlags.h"
#if OPERATING_SYSTEM == OS_WIN

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
#define _CRT_SECURE_NO_WARNINGS			// Don't show deprecation warning

// Windows Header Files:
#include <Windows.h>
#include <Windowsx.h>
#include <process.h>
#include <Shellapi.h>

// C RunTime Header Files
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <tchar.h>

// C++ Header Files
#include <string>

// Additional header files

#endif
