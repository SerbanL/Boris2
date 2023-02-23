#pragma once

#include "VEC.h"

//--------------------------------------------GETTERS

//get index of cell which contains position (absolute value, not relative to start of rectangle), capped to mesh size
template <typename VType>
INT3 VEC<VType>::cellidx_from_position(const DBL3& absolute_position) const
{
	INT3 ijk = INT3((int)floor_epsilon((absolute_position.x - rect.s.x) / h.x),
		(int)floor_epsilon((absolute_position.y - rect.s.y) / h.y),
		(int)floor_epsilon((absolute_position.z - rect.s.z) / h.z));

	if (ijk.i < 0) ijk.i = 0;
	if (ijk.j < 0) ijk.j = 0;
	if (ijk.k < 0) ijk.k = 0;

	if (ijk.i > n.x) ijk.i = n.x;
	if (ijk.j > n.y) ijk.j = n.y;
	if (ijk.k > n.z) ijk.k = n.z;

	return ijk;
}

//extract box of cells intersecting with the given rectangle (rectangle is in absolute coordinates). Cells in box : from and including start, up to but not including end; Limited to VEC sizes.
template <typename VType>
Box VEC<VType>::box_from_rect_max(const Rect& rectangle) const
{
	if (!rectangle.intersects(rect)) return Box();

	//get start point. this will be limited to 0 to n (inclusive)
	INT3 start = cellidx_from_position(rectangle.s);

	//the Rect could be a plane rectangle on one of the surfaces of this mesh, so adjust start point for this
	if (start.x >= n.x) start.x = n.x - 1;
	if (start.y >= n.y) start.y = n.y - 1;
	if (start.z >= n.z) start.z = n.z - 1;

	//get end point. this will be limited to 0 to n (inclusive)
	INT3 end = cellidx_from_position(rectangle.e);

	DBL3 snap = (h & end) + rect.s;

	//add 1 since end point must be included in the box, unless the rectangle end point is already at the end of a cell
	if ((IsNZ(snap.x - rectangle.e.x) || start.x == end.x) && end.x < n.x) end.x++;
	if ((IsNZ(snap.y - rectangle.e.y) || start.y == end.y) && end.y < n.y) end.y++;
	if ((IsNZ(snap.z - rectangle.e.z) || start.z == end.z) && end.z < n.z) end.z++;

	return Box(start, end);
}

//extract box of cells completely included in the given rectangle
template <typename VType>
Box VEC<VType>::box_from_rect_min(const Rect& rectangle) const
{
	if (!rectangle.intersects(rect)) return Box();

	//get i,j,k indexes of cells containing the start and end points of mesh_intersection
	INT3 start = cellidx_from_position(rectangle.s);
	INT3 end = cellidx_from_position(rectangle.e);

	//adjust start so that Box(start, end) represents the set of cells completely included in mesh_intersection
	Rect cell_rect = get_cellrect(start);
	if (!rectangle.contains(cell_rect)) start += INT3(1, 0, 0);

	cell_rect = get_cellrect(start);
	if (!rectangle.contains(cell_rect)) start += INT3(-1, 1, 0);

	cell_rect = get_cellrect(start);
	if (!rectangle.contains(cell_rect)) start += INT3(0, -1, 1);

	cell_rect = get_cellrect(start);
	if (!rectangle.contains(cell_rect)) start += INT3(1, 1, 0);

	if (end.x < start.x || end.y < start.y || end.z < start.z) return Box();

	return Box(start, end);
}

//count cells which don't have a null value set : i.e. non-empty.
template <typename VType>
int VEC<VType>::get_nonempty_cells(void) const
{
	int non_empty = 0;

#pragma omp parallel for reduction(+:non_empty)
	for (int idx = 0; idx < n.dim(); idx++) {

		if (quantity[idx] != VType()) non_empty = non_empty + 1;
	}

	return non_empty;
}

//check if all cells intersecting the rectangle (absolute coordinates) are empty
template <typename VType>
bool VEC<VType>::is_empty(const Rect& rectangle) const
{
	Box cells = box_from_rect_max(rectangle);

	for (int i = cells.s.x; i < cells.e.x; i++) {
		for (int j = cells.s.y; j < cells.e.y; j++) {
			for (int k = cells.s.z; k < cells.e.z; k++) {

				if (is_not_empty(INT3(i, j, k))) return false;
			}
		}
	}

	return true;
}

//check if all cells intersecting the rectangle (absolute coordinates) are not empty
template <typename VType>
bool VEC<VType>::is_not_empty(const Rect& rectangle) const
{
	Box cells = box_from_rect_max(rectangle);

	for (int i = cells.s.x; i < cells.e.x; i++) {
		for (int j = cells.s.y; j < cells.e.y; j++) {
			for (int k = cells.s.z; k < cells.e.z; k++) {

				if (is_empty(INT3(i, j, k))) return false;
			}
		}
	}

	return true;
}