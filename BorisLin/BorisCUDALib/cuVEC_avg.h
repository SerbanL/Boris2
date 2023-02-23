#pragma once

#include "cuTypes.h"
#include "cuFuncs_Aux.h"

//--------------------------------------------AVERAGING OPERATIONS

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
	cuBReal d_max = cu_GetMagnitude(stencil / 2 + h / 2);						//this is the maximum possible distance (a cell at this distance will get a weight of zero)	
	cuBReal d_total = 0;														//total reciprocal distance
	
	for (int ii = (idx_ll.i >= 0 ? idx_ll.i : 0); ii < (idx_ur.i < n.x ? idx_ur.i : n.x); ii++) {
		for (int jj = (idx_ll.j >= 0 ? idx_ll.j : 0); jj < (idx_ur.j < n.y ? idx_ur.j : n.y); jj++) {
			for (int kk = (idx_ll.k >= 0 ? idx_ll.k : 0); kk < (idx_ur.k < n.z ? idx_ur.k : n.z); kk++) {

				//find reciprocal distance for each cell : this is its weight * total reciprocal distance (to divide by d_total at the end)
				cuBReal d_recip = d_max - cu_get_distance(coord, cuReal3(((cuBReal)ii + 0.5)*h.x, ((cuBReal)jj + 0.5)*h.y, ((cuBReal)kk + 0.5)*h.z));
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
	cuBReal d_max = cu_GetMagnitude(cs / 2 + h / 2);						//this is the maximum possible distance (a cell at this distance will get a weight of zero)	
	cuBReal d_total = 0;														//total reciprocal distance
	
	for (int ii = (idx_ll.i >= 0 ? idx_ll.i : 0); ii < (idx_ur.i < n.x ? idx_ur.i : n.x); ii++) {
		for (int jj = (idx_ll.j >= 0 ? idx_ll.j : 0); jj < (idx_ur.j < n.y ? idx_ur.j : n.y); jj++) {
			for (int kk = (idx_ll.k >= 0 ? idx_ll.k : 0); kk < (idx_ur.k < n.z ? idx_ur.k : n.z); kk++) {

				//find reciprocal distance for each cell : this is its weight * total reciprocal distance (to divide by d_total at the end)
				cuBReal d_recip = d_max - cu_get_distance(coord, cuReal3(((cuBReal)ii + 0.5)*h.x, ((cuBReal)jj + 0.5)*h.y, ((cuBReal)kk + 0.5)*h.z));
				d_total += d_recip;

				smoothedValue += d_recip * quantity[ii + jj * n.x + kk * n.x*n.y];
			}
		}
	}
	
	//finally divide by total reciprocal distance to finish off the averaging
	if (d_total) smoothedValue /= d_total;
	
	return smoothedValue;
}

//full average in given rectangle (relative coordinates).
template <typename VType>
__device__ VType cuVEC<VType>::average(const cuRect& rectangle)
{
	cuBox box = box_from_rect_max(rectangle + rect.s);

	VType av = VType();
	int count = 0;

	for (int i = (box.s.x >= 0 ? box.s.x : 0); i < (box.e.x <= n.x ? box.e.x : n.x); i++) {
		for (int j = (box.s.y >= 0 ? box.s.y : 0); j < (box.e.y <= n.y ? box.e.y : n.y); j++) {
			for (int k = (box.s.z >= 0 ? box.s.z : 0); k < (box.e.z <= n.z ? box.e.z : n.z); k++) {

				//get running average
				VType value = quantity[i + j * n.x + k * n.x*n.y];

				count++;
				av = (av * (count - 1) + value) / count;
			}
		}
	}

	return av;
}

//average in given rectangle (relative coordinates), excluding zero points (assumed empty).
template <typename VType>
__device__ VType cuVEC<VType>::average_nonempty(const cuRect& rectangle)
{
	cuBox box = box_from_rect_max(rectangle + rect.s);

	VType av = VType();
	int count = 0;

	for (int i = (box.s.x >= 0 ? box.s.x : 0); i < (box.e.x <= n.x ? box.e.x : n.x); i++) {
		for (int j = (box.s.y >= 0 ? box.s.y : 0); j < (box.e.y <= n.y ? box.e.y : n.y); j++) {
			for (int k = (box.s.z >= 0 ? box.s.z : 0); k < (box.e.z <= n.z ? box.e.z : n.z); k++) {

				//get running average
				VType value = quantity[i + j * n.x + k * n.x*n.y];

				if (cuIsNZ(cu_GetMagnitude(value))) {

					count++;
					av = (av * (count - 1) + value) / count;
				}
			}
		}
	}

	return av;
}