#pragma once

#include "Funcs_Aux_base.h"
#include "Types.h"

///////////////////////////////////////////////////////////////////////////////
//RECTS

inline void ShiftRect(D2D1_RECT_F &rect, float dX, float dY) { rect.left += dX; rect.right += dX; rect.top += dY; rect.bottom += dY; }
inline void ShiftRect(D2D1_RECT_F &rect, FLT2 delta) { rect.left += delta.x; rect.right += delta.x; rect.top += delta.y; rect.bottom += delta.y; }

inline D2D1_RECT_F GetShiftedRect(D2D1_RECT_F rect, float dX, float dY) { rect.left += dX; rect.right += dX; rect.top += dY; rect.bottom += dY; return rect;}
inline D2D1_RECT_F GetShiftedRect(D2D1_RECT_F rect, FLT2 delta) { rect.left += delta.x; rect.right += delta.x; rect.top += delta.y; rect.bottom += delta.y; return rect; }

inline D2D1_RECT_F GetNormalisedRect(D2D1_RECT_F rect, float width, float height) { return D2D1::RectF(rect.left/width, rect.top/height, rect.right/width, rect.bottom/height); }

//Shrink or enlarge rectangle uniformly : change > 0 means shrinking
inline void UniformResizeRect(D2D1_RECT_F &rect, float change) { rect.left += change; rect.right -= change; rect.top += change; rect.bottom -= change; }

inline bool IsInside(D2D1_RECT_F rect, INT2 mouse) { return (mouse.i >= rect.left && mouse.i <= rect.right && mouse.j >= rect.top && mouse.j <= rect.bottom); }