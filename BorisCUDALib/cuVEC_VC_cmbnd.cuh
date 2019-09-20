#pragma once

#include "cuVEC_VC.h"
#include "cuVEC_VC_cmbnd.h"
#include "cuFuncs_Math.h"
#include "launchers.h"
#include "cuVEC_aux.cuh"

#include "CUDAError.h"

//--------------------------------------------CALCULATE COMPOSITE MEDIA BOUNDARY VALUES : CONTINUOUS

//-------------------- GLOBAL

template <typename VType, typename Class_CMBND>
__global__ void set_cmbnd_continuous_kernel(
	cuVEC_VC<VType>& V_sec, cuVEC_VC<VType>& V_pri, 
	Class_CMBND& cmbndFuncs_sec, Class_CMBND& cmbndFuncs_pri, 
	CMBNDInfoCUDA& contact)
{
	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 box_sizes = contact.cells_box.size();

	if (box_idx < box_sizes.dim()) {

		int i = (box_idx % box_sizes.x) + contact.cells_box.s.i;
		int j = ((box_idx / box_sizes.x) % box_sizes.y) + contact.cells_box.s.j;
		int k = (box_idx / (box_sizes.x * box_sizes.y)) + contact.cells_box.s.k;

		//cellsizes perpendicular to interface
		cuBReal hL = contact.hshift_secondary.norm();
		cuBReal hR = contact.hshift_primary.norm();
		cuBReal hmax = (hL > hR ? hL : hR);

		int cell1_idx = i + j * V_pri.n.x + k * V_pri.n.x*V_pri.n.y;

		if (V_pri.is_empty(cell1_idx) || V_pri.is_not_cmbnd(cell1_idx)) return;

		//calculate second primary cell index
		int cell2_idx = (i + contact.cell_shift.i) + (j + contact.cell_shift.j) * V_pri.n.x + (k + contact.cell_shift.k) * V_pri.n.x*V_pri.n.y;

		//cell values either side of the boundary: V_m2 V_m1 | V_1 V_2; positions as : -2 -1 | 1 2
		//V_m2 and V_2 are known. We need to set values V_m1 and V_1. Here we only set V_1. For this primary, secondary mesh pair there will be another one with order reversed, and there our V_m1 value will be set - so don't worry about it here!
		//NOTE : meshes at composite media boundaries must always be at least 2 non-empty cells thick in directions perpendicular to interface !!! -> this is actually checked when the list of contacts is made

		//relative position of cell -1 in secondary mesh
		cuReal3 relpos_m1 = V_pri.rect.s - V_sec.rect.s + ((cuReal3(i, j, k) + cuReal3(0.5)) & V_pri.h) + (contact.hshift_primary + contact.hshift_secondary) / 2;

		//stencil is used for weighted_average to obtain values in the secondary mesh : has size equal to primary cellsize area on interface with thickness set by secondary cellsize thickness
		cuReal3 stencil = V_pri.h - cu_mod(contact.hshift_primary) + cu_mod(contact.hshift_secondary);

		//potential values at cells -2 and 2
		VType V_2 = V_pri[cell2_idx];
		VType V_m2 = V_sec.weighted_average(relpos_m1 + contact.hshift_secondary, stencil);

		//obtain a and b values used to define the flux as f(V) = a + b V', both on primary and secondary

		//a values
		VType a_val_sec = cmbndFuncs_sec.a_func_sec(relpos_m1, contact.hshift_secondary, stencil);
		VType a_val_pri = cmbndFuncs_pri.a_func_pri(cell1_idx, cell2_idx, contact.hshift_secondary);

		//b values adjusted with weights
		cuBReal b_val_sec = cmbndFuncs_sec.b_func_sec(relpos_m1, contact.hshift_secondary, stencil) * contact.weights.i;
		cuBReal b_val_pri = cmbndFuncs_pri.b_func_pri(cell1_idx, cell2_idx) * contact.weights.j;

		//V'' values at cell positions -1 and 1
		VType Vdiff2_sec = cmbndFuncs_sec.diff2_sec(relpos_m1, stencil);
		VType Vdiff2_pri = cmbndFuncs_pri.diff2_pri(cell1_idx);

		//Formula for V1
		
		V_pri[cell1_idx] = (V_m2 * 2 * b_val_sec / 3 + V_2 * (b_val_pri + b_val_sec / 3)
			- Vdiff2_sec * b_val_sec * hL * hL - Vdiff2_pri * b_val_pri * hR * hR
			+ (a_val_pri - a_val_sec) * hmax) / (b_val_sec + b_val_pri);
	}
}

//-------------------- LAUNCHER

template <typename VType>
template <typename Class_CMBND>
__host__ void cuVEC_VC<VType>::set_cmbnd_continuous(size_t size, cuVEC_VC<VType>& V_sec, Class_CMBND& cmbndFuncs_sec, Class_CMBND& cmbndFuncs_pri, CMBNDInfoCUDA& contact)
{	
	set_cmbnd_continuous_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (V_sec, *this, cmbndFuncs_sec, cmbndFuncs_pri, contact);
}

//--------------------------------------------CALCULATE COMPOSITE MEDIA BOUNDARY VALUES : CONTINUOUS FLUX ONLY

//-------------------- GLOBAL

template <typename VType, typename Class_CMBND, typename Class_CMBND_S>
__global__ void set_cmbnd_continuousflux_kernel(
	cuVEC_VC<VType>& V_sec, cuVEC_VC<VType>& V_pri,
	Class_CMBND& cmbndFuncs_sec, Class_CMBND& cmbndFuncs_pri,
	Class_CMBND_S& cmbndFuncs_s,
	CMBNDInfoCUDA& contact)
{
	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 box_sizes = contact.cells_box.size();

	if (box_idx < box_sizes.dim()) {

		int i = (box_idx % box_sizes.x) + contact.cells_box.s.i;
		int j = ((box_idx / box_sizes.x) % box_sizes.y) + contact.cells_box.s.j;
		int k = (box_idx / (box_sizes.x * box_sizes.y)) + contact.cells_box.s.k;

		//cellsizes perpendicular to interface
		cuBReal hL = contact.hshift_secondary.norm();
		cuBReal hR = contact.hshift_primary.norm();
		cuBReal hmax = (hL > hR ? hL : hR);

		int cell1_idx = i + j * V_pri.n.x + k * V_pri.n.x*V_pri.n.y;

		if (V_pri.is_empty(cell1_idx) || V_pri.is_not_cmbnd(cell1_idx)) return;

		//calculate second primary cell index
		int cell2_idx = (i + contact.cell_shift.i) + (j + contact.cell_shift.j) * V_pri.n.x + (k + contact.cell_shift.k) * V_pri.n.x*V_pri.n.y;

		//cell values either side of the boundary: V_m2 V_m1 | V_1 V_2; positions as : -2 -1 | 1 2
		//V_m2 and V_2 are known. We need to set values V_m1 and V_1. Here we only set V_1. For this primary, secondary mesh pair there will be another one with order reversed, and there our V_m1 value will be set - so don't worry about it here!
		//NOTE : meshes at composite media boundaries must always be at least 2 non-empty cells thick in directions perpendicular to interface !!! -> this is actually checked when the list of contacts is made

		//relative position of cell -1 in secondary mesh
		cuReal3 relpos_m1 = V_pri.rect.s - V_sec.rect.s + ((cuReal3(i, j, k) + cuReal3(0.5)) & V_pri.h) + (contact.hshift_primary + contact.hshift_secondary) / 2;

		//stencil is used for weighted_average to obtain values in the secondary mesh : has size equal to primary cellsize area on interface with thickness set by secondary cellsize thickness
		cuReal3 stencil = V_pri.h - cu_mod(contact.hshift_primary) + cu_mod(contact.hshift_secondary);

		//potential values at cells -2 and 2
		VType V_2 = V_pri[cell2_idx];
		VType V_m2 = V_sec.weighted_average(relpos_m1 + contact.hshift_secondary, stencil);

		//obtain a and b values used to define the flux as f(V) = a + b V', both on primary and secondary

		//a values
		VType a_val_sec = cmbndFuncs_sec.a_func_sec(relpos_m1, contact.hshift_secondary, stencil);
		VType a_val_pri = cmbndFuncs_pri.a_func_pri(cell1_idx, cell2_idx, contact.hshift_secondary);

		//b values adjusted with weights
		cuBReal b_val_sec = cmbndFuncs_sec.b_func_sec(relpos_m1, contact.hshift_secondary, stencil) * contact.weights.i;
		cuBReal b_val_pri = cmbndFuncs_pri.b_func_pri(cell1_idx, cell2_idx) * contact.weights.j;

		//V'' values at cell positions -1 and 1
		VType Vdiff2_sec = cmbndFuncs_sec.diff2_sec(relpos_m1, stencil);
		VType Vdiff2_pri = cmbndFuncs_pri.diff2_pri(cell1_idx);

		//A value
		VType A_val = cmbndFuncs_s.A_func(cell1_idx, cell2_idx, relpos_m1, contact.hshift_secondary, stencil, cmbndFuncs_sec, cmbndFuncs_pri);

		//B value, then use it to find G value
		cuBReal B_val = cmbndFuncs_s.B_func(cell1_idx, cell2_idx, relpos_m1, contact.hshift_secondary, stencil, cmbndFuncs_sec, cmbndFuncs_pri);

		cuBReal G_val = 0.0;
		if (B_val) G_val = 2 * b_val_pri * b_val_sec / (3 * B_val * hmax);

		//Formula for V1 - note, this reduces to the continuous case for G_val = 0 (or B_val tends to infinity) and A_val = VType(0)
		V_pri[cell1_idx] = (V_m2 * 2 * b_val_sec / 3 + V_2 * (b_val_pri + b_val_sec / 3 + G_val)
			- Vdiff2_sec * b_val_sec * hL * hL - Vdiff2_pri * (b_val_pri + G_val) * hR * hR
			+ (a_val_pri * (1 + G_val / b_val_pri) - a_val_sec) * hmax
			- (G_val / b_val_pri) * A_val * hmax) / (b_val_sec + b_val_pri + G_val);
	}
}

//-------------------- LAUNCHER

template <typename VType>
template <typename Class_CMBND, typename Class_CMBND_S>
__host__ void cuVEC_VC<VType>::set_cmbnd_continuousflux(size_t size, cuVEC_VC<VType>& V_sec, Class_CMBND& cmbndFuncs_sec, Class_CMBND& cmbndFuncs_pri, Class_CMBND_S& cmbndFuncs_s, CMBNDInfoCUDA& contact)
{
	set_cmbnd_continuousflux_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (V_sec, *this, cmbndFuncs_sec, cmbndFuncs_pri, cmbndFuncs_s, contact);
}

//--------------------------------------------CALCULATE COMPOSITE MEDIA BOUNDARY VALUES : DISCONTINUOUS

//-------------------- GLOBAL

template <typename VType, typename Class_CMBND, typename Class_CMBND_S>
__global__ void set_cmbnd_discontinuous_kernel(
	cuVEC_VC<VType>& V_sec, cuVEC_VC<VType>& V_pri,
	Class_CMBND& cmbndFuncs_sec, Class_CMBND& cmbndFuncs_pri,
	Class_CMBND_S& cmbndFuncs_s,
	CMBNDInfoCUDA& contact)
{
	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 box_sizes = contact.cells_box.size();

	if (box_idx < box_sizes.dim()) {

		int i = (box_idx % box_sizes.x) + contact.cells_box.s.i;
		int j = ((box_idx / box_sizes.x) % box_sizes.y) + contact.cells_box.s.j;
		int k = (box_idx / (box_sizes.x * box_sizes.y)) + contact.cells_box.s.k;

		//cellsizes perpendicular to interface
		cuBReal hL = contact.hshift_secondary.norm();
		cuBReal hR = contact.hshift_primary.norm();
		cuBReal hmax = (hL > hR ? hL : hR);

		int cell1_idx = i + j * V_pri.n.x + k * V_pri.n.x*V_pri.n.y;

		if (V_pri.is_empty(cell1_idx) || V_pri.is_not_cmbnd(cell1_idx)) return;

		//calculate second primary cell index
		int cell2_idx = (i + contact.cell_shift.i) + (j + contact.cell_shift.j) * V_pri.n.x + (k + contact.cell_shift.k) * V_pri.n.x*V_pri.n.y;

		//cell values either side of the boundary: V_m2 V_m1 | V_1 V_2; positions as : -2 -1 | 1 2
		//V_m2 and V_2 are known. We need to set values V_m1 and V_1. Here we only set V_1. For this primary, secondary mesh pair there will be another one with order reversed, and there our V_m1 value will be set - so don't worry about it here!
		//NOTE : meshes at composite media boundaries must always be at least 2 non-empty cells thick in directions perpendicular to interface !!! -> this is actually checked when the list of contacts is made

		//relative position of cell -1 in secondary mesh
		cuReal3 relpos_m1 = V_pri.rect.s - V_sec.rect.s + ((cuReal3(i, j, k) + cuReal3(0.5)) & V_pri.h) + (contact.hshift_primary + contact.hshift_secondary) / 2;

		//stencil is used for weighted_average to obtain values in the secondary mesh : has size equal to primary cellsize area on interface with thickness set by secondary cellsize thickness
		cuReal3 stencil = V_pri.h - cu_mod(contact.hshift_primary) + cu_mod(contact.hshift_secondary);

		//potential values at cells -2 and 2
		VType V_2 = V_pri[cell2_idx];
		VType V_m2 = V_sec.weighted_average(relpos_m1 + contact.hshift_secondary, stencil);

		//obtain a and b values used to define the flux as f(V) = a + b V', both on primary and secondary

		//a values
		VType a_val_sec = cmbndFuncs_sec.a_func_sec(relpos_m1, contact.hshift_secondary, stencil);
		VType a_val_pri = cmbndFuncs_pri.a_func_pri(cell1_idx, cell2_idx, contact.hshift_secondary);

		//b values adjusted with weights
		cuBReal b_val_sec = cmbndFuncs_sec.b_func_sec(relpos_m1, contact.hshift_secondary, stencil) * contact.weights.i;
		cuBReal b_val_pri = cmbndFuncs_pri.b_func_pri(cell1_idx, cell2_idx) * contact.weights.j;

		//V'' values at cell positions -1 and 1
		VType Vdiff2_sec = cmbndFuncs_sec.diff2_sec(relpos_m1, stencil);
		VType Vdiff2_pri = cmbndFuncs_pri.diff2_pri(cell1_idx);

		//A values
		VType A_val_sec = cmbndFuncs_s.A_func_sec(cell1_idx, cell2_idx, relpos_m1, contact.hshift_secondary, stencil, cmbndFuncs_sec, cmbndFuncs_pri);
		VType A_val_pri = cmbndFuncs_s.A_func_pri(cell1_idx, cell2_idx, relpos_m1, contact.hshift_secondary, stencil, cmbndFuncs_sec, cmbndFuncs_pri);

		//B values
		auto B_val_sec = cmbndFuncs_s.B_func_sec(cell1_idx, cell2_idx, relpos_m1, contact.hshift_secondary, stencil, cmbndFuncs_sec, cmbndFuncs_pri);
		auto B_val_pri = cmbndFuncs_s.B_func_pri(cell1_idx, cell2_idx, relpos_m1, contact.hshift_secondary, stencil, cmbndFuncs_sec, cmbndFuncs_pri);

		//c values
		cuBReal c_m1 = cmbndFuncs_sec.c_func_sec(relpos_m1, stencil);
		cuBReal c_m2 = cmbndFuncs_sec.c_func_sec(relpos_m1 + contact.hshift_secondary, stencil);
		cuBReal c_1 = cmbndFuncs_pri.c_func_pri(cell1_idx);
		cuBReal c_2 = cmbndFuncs_pri.c_func_pri(cell2_idx);

		//Form N and N3 values:
		decltype(B_val_sec) N = hmax * c_m2 * B_val_sec + 2.0 * b_val_sec * cu_ident<decltype(B_val_sec)>();
		decltype(B_val_sec) N3 = 3.0 * hmax * c_m1 * B_val_sec + 2.0 * b_val_sec * cu_ident<decltype(B_val_sec)>();

		//Form G value:
		decltype(B_val_pri) G = (B_val_pri * cu_inverse<decltype(N3)>(N3)) * hmax;

		//Form inverse D value:
		decltype(B_val_pri) D_inv = cu_inverse<decltype(B_val_pri)>(9.0 * hmax * c_m1 * c_1 * G * B_val_sec - 3.0 * hmax * c_1 * B_val_pri - 2.0 * b_val_pri * cu_ident<decltype(B_val_pri)>());

		//Formula for V1
		V_pri[cell1_idx] = D_inv *
			((hmax * c_m2 * B_val_pri - 3.0 * c_m1 * G * N) * V_m2
			+ (3.0 * hmax * c_m1 * c_2 * G * B_val_sec - hmax * c_2 * B_val_pri - 2.0 * b_val_pri * cu_ident<decltype(B_val_pri)>()) * V_2
			+ 6.0 * c_m1 * b_val_sec * hL * hL * G * Vdiff2_sec + 2.0 * b_val_pri * hR * hR * Vdiff2_pri
			+ (2.0 * (A_val_pri - a_val_pri) - 6.0 * c_m1 * G * (A_val_sec - a_val_sec)) * hmax);
	}
}

//-------------------- LAUNCHER

template <typename VType>
template <typename Class_CMBND, typename Class_CMBND_S>
__host__ void cuVEC_VC<VType>::set_cmbnd_discontinuous(size_t size, cuVEC_VC<VType>& V_sec, Class_CMBND& cmbndFuncs_sec, Class_CMBND& cmbndFuncs_pri, Class_CMBND_S& cmbndFuncs_s, CMBNDInfoCUDA& contact)
{
	set_cmbnd_discontinuous_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (V_sec, *this, cmbndFuncs_sec, cmbndFuncs_pri, cmbndFuncs_s, contact);
}