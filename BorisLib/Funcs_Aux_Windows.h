#pragma once

#include "BorisLib_Config.h"
#if OPERATING_SYSTEM == OS_WIN

#define _CRT_SECURE_NO_WARNINGS

#include <d2d1_1.h>
#include <windows.h>

#include "Funcs_Aux_base.h"
#include "Types_VAL.h"

///////////////////////////////////////////////////////////////////////////////
//MEMORY

//Return total free memory in MB
inline size_t MemGetFreeMB(void)
{
	MEMORYSTATUSEX statex;
	statex.dwLength = sizeof(statex);
	GlobalMemoryStatusEx(&statex);

	return statex.ullAvailPhys;
}

//Return total memory in MB
inline size_t MemGetTotalMB(void)
{
	MEMORYSTATUSEX statex;
	statex.dwLength = sizeof(statex);
	GlobalMemoryStatusEx(&statex);

	return statex.ullTotalPhys;
}

///////////////////////////////////////////////////////////////////////////////
//RECTS

inline void ShiftRect(D2D1_RECT_F &rect, float dX, float dY) { rect.left += dX; rect.right += dX; rect.top += dY; rect.bottom += dY; }
inline void ShiftRect(D2D1_RECT_F &rect, FLT2 delta) { rect.left += delta.x; rect.right += delta.x; rect.top += delta.y; rect.bottom += delta.y; }

inline D2D1_RECT_F GetShiftedRect(D2D1_RECT_F rect, float dX, float dY) { rect.left += dX; rect.right += dX; rect.top += dY; rect.bottom += dY; return rect; }
inline D2D1_RECT_F GetShiftedRect(D2D1_RECT_F rect, FLT2 delta) { rect.left += delta.x; rect.right += delta.x; rect.top += delta.y; rect.bottom += delta.y; return rect; }

inline D2D1_RECT_F GetNormalisedRect(D2D1_RECT_F rect, float width, float height) { return D2D1::RectF(rect.left / width, rect.top / height, rect.right / width, rect.bottom / height); }

//Shrink or enlarge rectangle uniformly : change > 0 means shrinking
inline void UniformResizeRect(D2D1_RECT_F &rect, float change) { rect.left += change; rect.right -= change; rect.top += change; rect.bottom -= change; }

inline bool IsInside(D2D1_RECT_F rect, INT2 mouse) { return (mouse.i >= rect.left && mouse.i <= rect.right && mouse.j >= rect.top && mouse.j <= rect.bottom); }

///////////////////////////////////////////////////////////////////////////////
//Others

inline bool SameColor(D2D1_COLOR_F color1, D2D1_COLOR_F color2) { return (IsE(color1.r, color2.r) && IsE(color1.g, color2.g) && IsE(color1.b, color2.b) && IsE(color1.a, color2.a)); }

inline unsigned GetSystemTickCount() { return GetTickCount(); }

#endif