#pragma once

#include "mcuVEC.h"

//clear all allocated memory
template <typename VType, typename MType>
void mcuVEC<VType, MType>::clear_memory_aux(void)
{
	if (pn_d) {

		delete[] pn_d;
		pn_d = nullptr;
	}

	if (prect_d) {

		delete[] prect_d;
		prect_d = nullptr;
	}

	if (pbox_d) {

		delete[] pbox_d;
		pbox_d = nullptr;
	}

	//delete any memory allocated for communication between devices
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		if (phalo_ngbr_n[mGPU] != nullptr) {

			delete phalo_ngbr_n[mGPU];
			phalo_ngbr_n[mGPU] = nullptr;
		}

		if (phalo_ngbr_p[mGPU] != nullptr) {

			delete phalo_ngbr_p[mGPU];
			phalo_ngbr_p[mGPU] = nullptr;
		}
	}
}

//for given h and rect values, find n
template <typename VType, typename MType>
cuSZ3 mcuVEC<VType, MType>::get_n_from_h_and_rect(const cuReal3& h_, const cuRect& rect_)
{
	//calculate new n from rect and current h
	cuSZ3 new_n = cu_round(rect_ / h_);
	if (new_n.x < 1) new_n.x = 1;
	if (new_n.y < 1) new_n.y = 1;
	if (new_n.z < 1) new_n.z = 1;

	return new_n;
}

//for configured number of devices and given n value, calculate the number of cells (as cuSZ3) each device should be assigned
//this is returned as first: n value for all devices but last. second: n value for last device (this can differ if devices cannot be assigned all same n value)
template <typename VType, typename MType>
std::pair<cuSZ3, cuSZ3> mcuVEC<VType, MType>::get_devices_n_values(cuSZ3 new_n)
{
	cuSZ3 nd, ndl;

	//break up space along largest dimension, with z dimension as default
	if (new_n.x > new_n.y && new_n.x > new_n.z) {

		nd = cuSZ3((new_n.x - (new_n.x % mGPU.get_num_devices())) / mGPU.get_num_devices(), new_n.y, new_n.z);
		ndl = cuSZ3(new_n.x - nd.x * (mGPU.get_num_devices() - 1), new_n.y, new_n.z);
	}
	else if (new_n.y > new_n.z) {

		nd = cuSZ3(new_n.x, (new_n.y - (new_n.y % mGPU.get_num_devices())) / mGPU.get_num_devices(), new_n.z);
		ndl = cuSZ3(new_n.x, new_n.y - nd.y * (mGPU.get_num_devices() - 1), new_n.z);
	}
	else {

		nd = cuSZ3(new_n.x, new_n.y, (new_n.z - (new_n.z % mGPU.get_num_devices())) / mGPU.get_num_devices());
		ndl = cuSZ3(new_n.x, new_n.y, new_n.z - nd.z * (mGPU.get_num_devices() - 1));
	}

	return std::pair<cuSZ3, cuSZ3>(nd, ndl);
}

//--------------------------------------------GETTERS

//from cell index return cell center coordinates (relative to start of rectangle)
template <typename VType, typename MType>
cuReal3 mcuVEC<VType, MType>::cellidx_to_position_cpu(int idx)
{
	cuReal3 ijk_pos = cuReal3((idx % n.x) + 0.5, ((idx / n.x) % n.y) + 0.5, (idx / (n.x*n.y)) + 0.5);
	return (h & ijk_pos);
}

//from cell index return cell center coordinates (relative to start of rectangle)
template <typename VType, typename MType>
cuReal3 mcuVEC<VType, MType>::cellidx_to_position_cpu(cuINT3 ijk)
{
	cuReal3 ijk_pos = cuReal3(ijk.i + 0.5, ijk.j + 0.5, ijk.k + 0.5);
	return (h & ijk_pos);
}

//return cell index from relative position : the inverse of cellidx_to_position
template <typename VType, typename MType>
int mcuVEC<VType, MType>::position_to_cellidx_cpu(const cuReal3& position)
{
	return (int)cu_floor_epsilon(position.x / h.x) + (int)cu_floor_epsilon(position.y / h.y) * n.x + (int)cu_floor_epsilon(position.z / h.z) * n.x * n.y;
}

//get index of cell which contains position (absolute value, not relative to start of rectangle), capped to mesh size
template <typename VType, typename MType>
cuINT3 mcuVEC<VType, MType>::cellidx_from_position_cpu(cuReal3 absolute_position)
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

//get cell rectangle (absolute values, not relative to start of mesh rectangle) for cell with index ijk
template <typename VType, typename MType>
cuRect mcuVEC<VType, MType>::get_cellrect_cpu(cuINT3 ijk)
{
	return cuRect(rect.s + (h & ijk), rect.s + (h & ijk) + h);
}

//get_cellrect using single index.
template <typename VType, typename MType>
cuRect mcuVEC<VType, MType>::get_cellrect_cpu(int idx)
{
	cuINT3 ijk = cuINT3((idx % n.x), (idx / n.x) % n.y, idx / (n.x*n.y));
	return cuRect(rect.s + (h & ijk), rect.s + (h & ijk) + h);
}

//extract box of cells intersecting with the given rectangle (rectangle is in absolute coordinates). Cells in box : from and including start, up to but not including end; Limited to cuVEC sizes.
template <typename VType, typename MType>
cuBox mcuVEC<VType, MType>::box_from_rect_max_cpu(cuRect rectangle)
{
	if (!rectangle.intersects(rect)) return cuBox();

	//get start point. this will be limited to 0 to n (inclusive)
	cuINT3 start = cellidx_from_position_cpu(rectangle.s);

	//the Rect could be a plane rectangle on one of the surfaces of this mesh, so adjust start point for this
	if (start.x >= int(n.x)) start.x = int(n.x) - 1;
	if (start.y >= int(n.y)) start.y = int(n.y) - 1;
	if (start.z >= int(n.z)) start.z = int(n.z) - 1;

	//get end point. this will be limited to 0 to n (inclusive)
	cuINT3 end = cellidx_from_position_cpu(rectangle.e);

	cuReal3 snap = (h & end) + rect.s;

	//add 1 since end point must be included in the box, unless the rectangle end point is already at the end of a cell
	if ((cuIsNZ(snap.x - rectangle.e.x) || start.x == end.x) && end.x < int(n.x)) end.x++;
	if ((cuIsNZ(snap.y - rectangle.e.y) || start.y == end.y) && end.y < int(n.y)) end.y++;
	if ((cuIsNZ(snap.z - rectangle.e.z) || start.z == end.z) && end.z < int(n.z)) end.z++;

	return cuBox(start, end);
}

//extract box of cells completely included in the given rectangle (rectangle is in absolute coordinates).
template <typename VType, typename MType>
cuBox mcuVEC<VType, MType>::box_from_rect_min_cpu(cuRect rectangle)
{
	if (!rectangle.intersects(rect)) return cuBox();

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

//count cells which don't have a null value set : i.e. non-empty; set result in aux_integer
template <typename VType, typename MType>
void mcuVEC<VType, MType>::count_nonempty_cells(void)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->count_nonempty_cells(pn_d[mGPU].dim());
	}
}
