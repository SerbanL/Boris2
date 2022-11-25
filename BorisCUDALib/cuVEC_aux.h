#pragma once

#include "cuVEC.h"

//--------------------------------------------GETTERS

//from cell index return cell center coordinates (relative to start of rectangle)
template <typename VType>
__device__ cuReal3 cuVEC<VType>::cellidx_to_position(int idx) const
{
	cuReal3 ijk_pos = cuReal3((idx % n.x) + 0.5, ((idx / n.x) % n.y) + 0.5, (idx / (n.x*n.y)) + 0.5);
	return (h & ijk_pos);
}

template <typename VType>
__host__ cuReal3 cuVEC<VType>::cellidx_to_position_cpu(int idx)
{
	cuSZ3 n_ = get_gpu_value(n);
	cuReal3 ijk_pos = cuReal3((idx % n_.x) + 0.5, ((idx / n_.x) % n_.y) + 0.5, (idx / (n_.x*n_.y)) + 0.5);
	return (get_gpu_value(h) & ijk_pos);
}

//from cell index return cell center coordinates (relative to start of rectangle)
template <typename VType>
__device__ cuReal3 cuVEC<VType>::cellidx_to_position(const cuINT3& ijk) const
{
	cuReal3 ijk_pos = cuReal3(ijk.i + 0.5, ijk.j + 0.5, ijk.k + 0.5);
	return (h & ijk_pos);
}

template <typename VType>
__host__ cuReal3 cuVEC<VType>::cellidx_to_position_cpu(cuINT3 ijk)
{
	cuReal3 ijk_pos = cuReal3(ijk.i + 0.5, ijk.j + 0.5, ijk.k + 0.5);
	return (get_gpu_value(h) & ijk_pos);
}

template <typename VType>
__host__ int cuVEC<VType>::position_to_cellidx_cpu(const cuReal3& position)
{
	cuReal3 h_cpu = cellsize_cpu();
	cuSZ3 n_cpu = size_cpu();
	return (int)cu_floor_epsilon(position.x / h_cpu.x) + (int)cu_floor_epsilon(position.y / h_cpu.y) * n_cpu.x + (int)cu_floor_epsilon(position.z / h_cpu.z) * n_cpu.x * n_cpu.y;
}

//get index of cell which contains position (absolute value, not relative to start of rectangle), capped to mesh size
template <typename VType>
__device__ cuINT3 cuVEC<VType>::cellidx_from_position(const cuReal3& absolute_position) const
{
	cuINT3 ijk = cuINT3(
		(int)cu_floor_epsilon((absolute_position.x - rect.s.x) / h.x),
		(int)cu_floor_epsilon((absolute_position.y - rect.s.y) / h.y),
		(int)cu_floor_epsilon((absolute_position.z - rect.s.z) / h.z));

	if (ijk.i < 0) ijk.i = 0;
	if (ijk.j < 0) ijk.j = 0;
	if (ijk.k < 0) ijk.k = 0;

	if (ijk.i > n.x) ijk.i = n.x;
	if (ijk.j > n.y) ijk.j = n.y;
	if (ijk.k > n.z) ijk.k = n.z;

	return ijk;
}

template <typename VType>
__host__ cuINT3 cuVEC<VType>::cellidx_from_position_cpu(cuReal3 absolute_position)
{
	cuSZ3 n_ = get_gpu_value(n);
	cuRect rect_ = get_gpu_value(rect);
	cuReal3 h_ = get_gpu_value(h);

	cuINT3 ijk = cuINT3(
		(int)cu_floor_epsilon((absolute_position.x - rect_.s.x) / h_.x),
		(int)cu_floor_epsilon((absolute_position.y - rect_.s.y) / h_.y),
		(int)cu_floor_epsilon((absolute_position.z - rect_.s.z) / h_.z));

	if (ijk.i < 0) ijk.i = 0;
	if (ijk.j < 0) ijk.j = 0;
	if (ijk.k < 0) ijk.k = 0;

	if (ijk.i > n_.x) ijk.i = n_.x;
	if (ijk.j > n_.y) ijk.j = n_.y;
	if (ijk.k > n_.z) ijk.k = n_.z;

	return ijk;
}

//get cell rectangle (absolute values, not relative to start of mesh rectangle) for cell with index ijk
template <typename VType>
__device__ cuRect cuVEC<VType>::get_cellrect(const cuINT3& ijk) const
{
	return cuRect(rect.s + (h & ijk), rect.s + (h & ijk) + h);
}

template <typename VType>
__host__ cuRect cuVEC<VType>::get_cellrect_cpu(cuINT3 ijk)
{
	cuRect rect_ = get_gpu_value(rect);
	cuReal3 h_ = get_gpu_value(h);

	return cuRect(rect_.s + (h_ & ijk), rect_.s + (h_ & ijk) + h_);
}

//get_cellrect using single index.
template <typename VType>
__device__ cuRect cuVEC<VType>::get_cellrect(int idx) const
{
	cuINT3 ijk = cuINT3((idx % n.x), (idx / n.x) % n.y, idx / (n.x*n.y));
	return cuRect(rect.s + (h & ijk), rect.s + (h & ijk) + h);
}

template <typename VType>
__host__ cuRect cuVEC<VType>::get_cellrect_cpu(int idx)
{
	cuSZ3 n_ = get_gpu_value(n);
	cuReal3 h_ = get_gpu_value(h);
	cuRect rect_ = get_gpu_value(rect);

	cuINT3 ijk = cuINT3((idx % n_.x), (idx / n_.x) % n_.y, idx / (n_.x*n_.y));
	return cuRect(rect_.s + (h_ & ijk), rect_.s + (h_ & ijk) + h_);
}

//extract box of cells intersecting with the given rectangle (rectangle is in absolute coordinates). Cells in box : from and including start, up to but not including end; Limited to cuVEC sizes.
template <typename VType>
__device__ cuBox cuVEC<VType>::box_from_rect_max(const cuRect& rectangle) const
{
	if (!rectangle.intersects(rect)) return cuBox();

	//get start point. this will be limited to 0 to n (inclusive)
	cuINT3 start = cellidx_from_position(rectangle.s);

	//the Rect could be a plane rectangle on one of the surfaces of this mesh, so adjust start point for this
	if (start.x >= n.x) start.x = n.x - 1;
	if (start.y >= n.y) start.y = n.y - 1;
	if (start.z >= n.z) start.z = n.z - 1;

	//get end point. this will be limited to 0 to n (inclusive)
	cuINT3 end = cellidx_from_position(rectangle.e);

	cuReal3 snap = (h & end) + rect.s;

	//add 1 since end point must be included in the box, unless the rectangle end point is already at the end of a cell
	if ((cuIsNZ(snap.x - rectangle.e.x) || start.x == end.x) && end.x < n.x) end.x++;
	if ((cuIsNZ(snap.y - rectangle.e.y) || start.y == end.y) && end.y < n.y) end.y++;
	if ((cuIsNZ(snap.z - rectangle.e.z) || start.z == end.z) && end.z < n.z) end.z++;

	return cuBox(start, end);
}

template <typename VType>
__host__ cuBox cuVEC<VType>::box_from_rect_max_cpu(cuRect rectangle)
{
	cuSZ3 n_ = get_gpu_value(n);
	cuReal3 h_ = get_gpu_value(h);
	cuRect rect_ = get_gpu_value(rect);

	if (!rectangle.intersects(rect_)) return cuBox();

	//get start point. this will be limited to 0 to n (inclusive)
	cuINT3 start = cellidx_from_position_cpu(rectangle.s);

	//the Rect could be a plane rectangle on one of the surfaces of this mesh, so adjust start point for this
	if (start.x >= int(n_.x)) start.x = int(n_.x) - 1;
	if (start.y >= int(n_.y)) start.y = int(n_.y) - 1;
	if (start.z >= int(n_.z)) start.z = int(n_.z) - 1;

	//get end point. this will be limited to 0 to n (inclusive)
	cuINT3 end = cellidx_from_position_cpu(rectangle.e);

	cuReal3 snap = (h_ & end) + rect_.s;

	//add 1 since end point must be included in the box, unless the rectangle end point is already at the end of a cell
	if ((cuIsNZ(snap.x - rectangle.e.x) || start.x == end.x) && end.x < int(n_.x)) end.x++;
	if ((cuIsNZ(snap.y - rectangle.e.y) || start.y == end.y) && end.y < int(n_.y)) end.y++;
	if ((cuIsNZ(snap.z - rectangle.e.z) || start.z == end.z) && end.z < int(n_.z)) end.z++;

	return cuBox(start, end);
}

//extract box of cells completely included in the given rectangle (rectangle is in absolute coordinates).
template <typename VType>
__device__ cuBox cuVEC<VType>::box_from_rect_min(const cuRect& rectangle) const
{
	if (!rectangle.intersects(rect)) return cuBox();

	//get i,j,k indexes of cells containing the start and end points of mesh_intersection
	cuINT3 start = cellidx_from_position(rectangle.s);
	cuINT3 end = cellidx_from_position(rectangle.e);

	//adjust start so that Box(start, end) represents the set of cells completely included in mesh_intersection
	cuRect cell_rect = get_cellrect(start);
	if (!rectangle.contains(cell_rect)) start += cuINT3(1, 0, 0);

	cell_rect = get_cellrect(start);
	if (!rectangle.contains(cell_rect)) start += cuINT3(-1, 1, 0);

	cell_rect = get_cellrect(start);
	if (!rectangle.contains(cell_rect)) start += cuINT3(0, -1, 1);

	cell_rect = get_cellrect(start);
	if (!rectangle.contains(cell_rect)) start += cuINT3(1, 1, 0);

	return cuBox(start, end);
}

template <typename VType>
__host__ cuBox cuVEC<VType>::box_from_rect_min_cpu(cuRect rectangle)
{
	if (!rectangle.intersects(get_gpu_value(rect))) return cuBox();

	//get i,j,k indexes of cells containing the start and end points of mesh_intersection
	cuINT3 start = cellidx_from_position_cpu(rectangle.s);
	cuINT3 end = cellidx_from_position_cpu(rectangle.e);

	//adjust start so that Box(start, end) represents the set of cells completely included in mesh_intersection
	cuRect cell_rect = get_cellrect_cpu(start);
	if (!rectangle.contains(cell_rect)) start += cuINT3(1, 0, 0);

	cell_rect = get_cellrect_cpu(start);
	if (!rectangle.contains(cell_rect)) start += cuINT3(-1, 1, 0);

	cell_rect = get_cellrect_cpu(start);
	if (!rectangle.contains(cell_rect)) start += cuINT3(0, -1, 1);

	cell_rect = get_cellrect_cpu(start);
	if (!rectangle.contains(cell_rect)) start += cuINT3(1, 1, 0);

	return cuBox(start, end);
}

template <typename VType>
__host__ int cuVEC<VType>::get_nonempty_cells_cpu(size_t arr_size)
{
	if (!arr_size) count_nonempty_cells(linear_size_cpu());
	else count_nonempty_cells(arr_size);
	
	return get_gpu_value(aux_integer);
}