#pragma once

#include "cuTypes.h"
#include "cuFuncs_Aux.h"

//--------------------------------------------MULTIPLE ENTRIES SETTERS

//set value in rectangle (i.e. in cells intersecting the rectangle), where the rectangle is relative to this cuVEC's rectangle.
template <typename VType> 
__host__ void cuVEC<VType>::setrect(const cuRect& rectangle, VType value)
{
	if (!get_gpu_value(rect).intersects(rectangle + get_gpu_value(rect).s)) return;

	cuBox cells_box = box_from_rect_max_cpu(rectangle + get_gpu_value(rect).s);

	setbox(cells_box, value);
}
