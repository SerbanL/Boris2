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

//--------------------------------------------OPERATIONS

template <typename VType> 
__device__ VType cuVEC<VType>::weighted_average(const cuReal3& coord, const cuReal3& stencil)
{
	//smoothedValue must be initialized to zero
	VType smoothedValue = VType();

	//1. find indexes of lower-left and uper-right cells intersecting with the stencil (i.e. find cells which contain its lower-left corner and upper-right corner)
	//Note, it could be a single cell containing both corners : either the stencil is strictly included in it or coincides with it

	//positions of lower-left and uper-right corners of cell relative to containing mesh
	cuReal3 pos_ll = coord - stencil / 2;
	cuReal3 pos_ur = pos_ll + stencil;

	cuINT3 idx_ll = cu_floor(pos_ll / h);
	cuINT3 idx_ur = cu_ceil(pos_ur / h);

	//2. obtain weighted average value from cells intersecting with the stencil (the larger the distance of the cell centre from coord, the lower the weight)
	cuReal d_max = cu_GetMagnitude(stencil / 2 + h / 2);						//this is the maximum possible distance (a cell at this distance will get a weight of zero)	
	cuReal d_total = 0;														//total reciprocal distance
	
	for (int ii = (idx_ll.i >= 0 ? idx_ll.i : 0); ii < (idx_ur.i < n.x ? idx_ur.i : n.x); ii++) {
		for (int jj = (idx_ll.j >= 0 ? idx_ll.j : 0); jj < (idx_ur.j < n.y ? idx_ur.j : n.y); jj++) {
			for (int kk = (idx_ll.k >= 0 ? idx_ll.k : 0); kk < (idx_ur.k < n.z ? idx_ur.k : n.z); kk++) {

				//find reciprocal distance for each cell : this is its weight * total reciprocal distance (to divide by d_total at the end)
				cuReal d_recip = d_max - cu_get_distance(coord, cuReal3(((cuReal)ii + 0.5)*h.x, ((cuReal)jj + 0.5)*h.y, ((cuReal)kk + 0.5)*h.z));
				d_total += d_recip;

				smoothedValue += d_recip * quantity[ii + jj * n.x + kk * n.x*n.y];
			}
		}
	}

	//finally divide by total reciprocal distance to finish off the averaging
	if(d_total) smoothedValue /= d_total;

	return smoothedValue;
}

template <typename VType>
__device__ VType cuVEC<VType>::weighted_average(const cuRect& rectangle)
{
	//smoothedValue must be initialized to zero
	VType smoothedValue = VType();

	//1. find indexes of lower-left and uper-right cells intersecting with the stencil (i.e. find cells which contain its lower-left corner and upper-right corner)
	//Note, it could be a single cell containing both corners : either the stencil is strictly included in it or coincides with it

	//indexes of lower-left and uper-right cells
	cuINT3 idx_ll = cu_floor((rectangle.s - rect.s) / h);
	cuINT3 idx_ur = cu_ceil((rectangle.e - rect.s) / h);

	cuReal3 stencil = rectangle.e - rectangle.s;
	cuReal3 center = ((rectangle.e + rectangle.s) / 2) - rect.s;

	//2. obtain weighted average value from cells intersecting with the stencil (the larger the distance of the cell centre from coord, the lower the weight)
	cuReal d_max = cu_GetMagnitude(stencil / 2 + h / 2);						//this is the maximum possible distance (a cell at this distance will get a weight of zero)	
	cuReal d_total = 0;														//total reciprocal distance

	for (int ii = (idx_ll.i >= 0 ? idx_ll.i : 0); ii < (idx_ur.i < n.x ? idx_ur.i : n.x); ii++) {
		for (int jj = (idx_ll.j >= 0 ? idx_ll.j : 0); jj < (idx_ur.j < n.y ? idx_ur.j : n.y); jj++) {
			for (int kk = (idx_ll.k >= 0 ? idx_ll.k : 0); kk < (idx_ur.k < n.z ? idx_ur.k : n.z); kk++) {

				//find reciprocal distance for each cell : this is its weight * total reciprocal distance (to divide by d_total at the end)
				cuReal d_recip = d_max - cu_get_distance(center, cuReal3(((cuReal)ii + 0.5)*h.x, ((cuReal)jj + 0.5)*h.y, ((cuReal)kk + 0.5)*h.z));
				d_total += d_recip;

				smoothedValue += d_recip * quantity[ii + jj * n.x + kk * n.x*n.y];
			}
		}
	}

	//finally divide by total reciprocal distance to finish off the averaging
	if (d_total) smoothedValue /= d_total;

	return smoothedValue;
}

template <typename VType>
__device__ VType cuVEC<VType>::weighted_average(const cuINT3& ijk, const cuReal3& cs)
{
	if (cs == h) return quantity[ijk.i + ijk.j * n.x + ijk.k * n.x * n.y];
	
	//smoothedValue must be initialized to zero
	VType smoothedValue = VType();

	//1. find indexes of lower-left and uper-right cells intersecting with the stencil (i.e. find cells which contain its lower-left corner and upper-right corner)
	//Note, it could be a single cell containing both corners : either the stencil is strictly included in it or coincides with it
	
	//positions of lower-left and uper-right corners of cell relative to containing mesh
	cuReal3 pos_ll = cs & ijk;
	cuReal3 pos_ur = pos_ll + cs;
	
	cuINT3 idx_ll = cu_floor(pos_ll / h);
	cuINT3 idx_ur = cu_ceil(pos_ur / h);

	cuReal3 coord = pos_ll + (cs / 2);
	
	//2. obtain weighted average value from cells intersecting with the stencil (the larger the distance of the cell centre from coord, the lower the weight)
	cuReal d_max = cu_GetMagnitude(cs / 2 + h / 2);						//this is the maximum possible distance (a cell at this distance will get a weight of zero)	
	cuReal d_total = 0;														//total reciprocal distance
	
	for (int ii = (idx_ll.i >= 0 ? idx_ll.i : 0); ii < (idx_ur.i < n.x ? idx_ur.i : n.x); ii++) {
		for (int jj = (idx_ll.j >= 0 ? idx_ll.j : 0); jj < (idx_ur.j < n.y ? idx_ur.j : n.y); jj++) {
			for (int kk = (idx_ll.k >= 0 ? idx_ll.k : 0); kk < (idx_ur.k < n.z ? idx_ur.k : n.z); kk++) {

				//find reciprocal distance for each cell : this is its weight * total reciprocal distance (to divide by d_total at the end)
				cuReal d_recip = d_max - cu_get_distance(coord, cuReal3(((cuReal)ii + 0.5)*h.x, ((cuReal)jj + 0.5)*h.y, ((cuReal)kk + 0.5)*h.z));
				d_total += d_recip;

				smoothedValue += d_recip * quantity[ii + jj * n.x + kk * n.x*n.y];
			}
		}
	}
	
	//finally divide by total reciprocal distance to finish off the averaging
	if (d_total) smoothedValue /= d_total;
	
	return smoothedValue;
}