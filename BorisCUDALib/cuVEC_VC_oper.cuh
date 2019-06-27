#pragma once

#include "cuVEC_VC.h"
#include "cuFuncs_Math.h"
#include "launchers.h"
#include "Reduction.cuh"

//--------------------------------------------MULTIPLE ENTRIES SETTERS - OTHERS

//------------------------------------------------------------------- SETNONEMPTY

template <typename VType>
__global__ void setnonempty_kernel(cuSZ3& n, int*& ngbrFlags, VType*& quantity, VType value)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < n.dim()) {

		if (ngbrFlags[idx] & NF_NOTEMPTY) {

			quantity[idx] = value;
		}
	}
}

template void cuVEC_VC<float>::setnonempty(float value);
template void cuVEC_VC<double>::setnonempty(double value);

template void cuVEC_VC<cuFLT3>::setnonempty(cuFLT3 value);
template void cuVEC_VC<cuDBL3>::setnonempty(cuDBL3 value);

template void cuVEC_VC<cuFLT33>::setnonempty(cuFLT33 value);
template void cuVEC_VC<cuDBL33>::setnonempty(cuDBL33 value);

template <typename VType>
__host__ void cuVEC_VC<VType>::setnonempty(VType value)
{
	setnonempty_kernel <<< (get_gpu_value(n).dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (n, ngbrFlags, quantity, value);
}

//------------------------------------------------------------------- SETRECTNONEMPTY

template <typename VType>
__global__ void setrectnonempty_kernel(cuSZ3& n, int*& ngbrFlags, VType*& quantity, cuBox box, VType value)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	if (idx < n.dim() && box.Contains(ijk)) {

		if (ngbrFlags[idx] & NF_NOTEMPTY) {

			quantity[idx] = value;
		}
	}
}

template void cuVEC_VC<float>::setrectnonempty(const cuRect& rectangle, float value);
template void cuVEC_VC<double>::setrectnonempty(const cuRect& rectangle, double value);

template void cuVEC_VC<cuFLT3>::setrectnonempty(const cuRect& rectangle, cuFLT3 value);
template void cuVEC_VC<cuDBL3>::setrectnonempty(const cuRect& rectangle, cuDBL3 value);

template void cuVEC_VC<cuFLT33>::setrectnonempty(const cuRect& rectangle, cuFLT33 value);
template void cuVEC_VC<cuDBL33>::setrectnonempty(const cuRect& rectangle, cuDBL33 value);

//set value in non-empty cells only in given rectangle (relative coordinates)
template <typename VType>
__host__ void cuVEC_VC<VType>::setrectnonempty(const cuRect& rectangle, VType value)
{
	cuBox box = box_from_rect_max_cpu(rectangle + get_gpu_value(rect).s);

	setrectnonempty_kernel <<< (get_gpu_value(n).dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, ngbrFlags, quantity, box, value);
}

//------------------------------------------------------------------- SCALE VALUES

template <typename VType>
__global__ void scale_values_kernel(cuSZ3& n, int*& ngbrFlags, VType*& quantity, cuReal constant)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < n.dim()) {

		if (ngbrFlags[idx] & NF_NOTEMPTY) {

			quantity[idx] *= constant;
		}
	}
}

template void cuVEC_VC<float>::scale_values(size_t size, cuReal constant);
template void cuVEC_VC<double>::scale_values(size_t size, cuReal constant);

template void cuVEC_VC<cuFLT3>::scale_values(size_t size, cuReal constant);
template void cuVEC_VC<cuDBL3>::scale_values(size_t size, cuReal constant);

template void cuVEC_VC<cuFLT33>::scale_values(size_t size, cuReal constant);
template void cuVEC_VC<cuDBL33>::scale_values(size_t size, cuReal constant);

//scale all stored values by the given constant
template <typename VType>
void cuVEC_VC<VType>::scale_values(size_t size, cuReal constant)
{
	scale_values_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, ngbrFlags, quantity, constant);
}

//------------------------------------------------------------------- RENORMALIZE (cuVEC_VC)

template <typename VType, typename PType>
__global__ void cuvec_vc_renormalize_kernel(cuSZ3& n, int*& ngbrFlags, VType*& quantity, PType new_norm)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < n.dim()) {

		PType curr_norm = cu_GetMagnitude(quantity[idx]);

		if ((ngbrFlags[idx] & NF_NOTEMPTY) && cuIsNZ(curr_norm)) {

			quantity[idx] *= new_norm / curr_norm;
		}
	}
}

template void cuVEC_VC<float>::renormalize(size_t arr_size, float new_norm);
template void cuVEC_VC<double>::renormalize(size_t arr_size, double new_norm);

template void cuVEC_VC<cuFLT3>::renormalize(size_t arr_size, float new_norm);
template void cuVEC_VC<cuDBL3>::renormalize(size_t arr_size, double new_norm);

//re-normalize all non-zero values to have the new magnitude (multiply by new_norm and divide by current magnitude)
template <typename VType>
template <typename PType>
__host__ void cuVEC_VC<VType>::renormalize(size_t arr_size, PType new_norm)
{
	cuvec_vc_renormalize_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, ngbrFlags, quantity, new_norm);
}

//------------------------------------------------------------------- AVERAGE NONEMPTY (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_average_nonempty_kernel(cuSZ3& n, int*& ngbrFlags, VType*& quantity, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	reduction_avg(idx, n.dim(), quantity, average_value, points_count, ngbrFlags[idx] & NF_NOTEMPTY);
}

template <typename VType>
__global__ void cuvec_vc_average_nonempty_kernel(cuSZ3& n, cuBox box, int*& ngbrFlags, VType*& quantity, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	reduction_avg(idx, n.dim(), quantity, average_value, points_count, box.Contains(ijk) && (ngbrFlags[idx] & NF_NOTEMPTY));
}

template float cuVEC_VC<float>::average_nonempty(size_t arr_size, cuBox box);
template double cuVEC_VC<double>::average_nonempty(size_t arr_size, cuBox box);

template cuFLT3 cuVEC_VC<cuFLT3>::average_nonempty(size_t arr_size, cuBox box);
template cuDBL3 cuVEC_VC<cuDBL3>::average_nonempty(size_t arr_size, cuBox box);

template <typename VType>
__host__ VType cuVEC_VC<VType>::average_nonempty(size_t arr_size, cuBox box)
{
	zero_aux_values <<<1, CUDATHREADS >>> (aux_value, aux_value2, aux_value3, aux_real, aux_integer);

	if (box.IsNull()) {

		cuvec_vc_average_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, ngbrFlags, quantity, aux_value, aux_integer);
	}
	else cuvec_vc_average_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, box, ngbrFlags, quantity, aux_value, aux_integer);

	VType av = get_gpu_value(aux_value);
	size_t points = get_gpu_value(aux_integer);

	if (points) return av / points;
	else return VType();
}

//------------------------------------------------------------------- AVERAGE NONEMPTY using a passed Rect

template <typename VType>
__global__ void average_nonempty_kernel(cuSZ3& n, cuRect rectangle, cuVEC_VC<VType>& cuvec, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	reduction_avg(idx, n.dim(), cuvec.data(), average_value, points_count, cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(ijk));
}

template float cuVEC_VC<float>::average_nonempty(size_t arr_size, cuRect rectangle);
template double cuVEC_VC<double>::average_nonempty(size_t arr_size, cuRect rectangle);

template cuFLT3 cuVEC_VC<cuFLT3>::average_nonempty(size_t arr_size, cuRect rectangle);
template cuDBL3 cuVEC_VC<cuDBL3>::average_nonempty(size_t arr_size, cuRect rectangle);

template <typename VType>
__host__ VType cuVEC_VC<VType>::average_nonempty(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average_nonempty(arr_size, cuBox());

	zero_aux_values <<<1, CUDATHREADS >>> (aux_value, aux_value2, aux_value3, aux_real, aux_integer);

	average_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, rectangle, *this, aux_value, aux_integer);

	VType av = get_gpu_value(aux_value);
	size_t points = get_gpu_value(aux_integer);

	if (points) return av / points;
	else return VType();
}

//------------------------------------------------------------------- SHIFT : x

template <typename VType>
__global__ void shift_x_left1_kernel(cuRect shift_rect, cuVEC_VC<VType>& cuVEC, VType*& aux_block_values)
{
	__shared__ VType shared_memory[CUDATHREADS];

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = cuVEC.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBox shift_box = cuVEC.box_from_rect_min(shift_rect);
	//go from shift_box.s.x to shift_box.e.x - 2 inclusive : shift_box.e.x - 2 takes value from shift_box.e.x - 1
	shift_box.e.x -= 2;

	if (idx < n.dim()) {

		shared_memory[threadIdx.x] = cuVEC[idx];
	}
	else shared_memory[threadIdx.x] = VType();

	//all values in this block must be transferred to shared memory before proceeding
	__syncthreads();

	//shift values within the box using shared_memory (both destination and source must not be empty).
	//We cannot shift the first element in this block since it must go to the block before this - store it in aux_block_values for later.
	//Similarly we cannot write the last element in this block since it needs a value from the next block.
	if (threadIdx.x == 0) aux_block_values[blockIdx.x] = shared_memory[0];

	if (shift_box.Contains(ijk) && cuVEC.is_not_empty(idx) && cuVEC.is_not_empty(idx + 1)) {

		if (threadIdx.x < CUDATHREADS - 1) {

			cuVEC[idx] = shared_memory[threadIdx.x + 1];
		}
	}
}

template <typename VType>
__global__ void shift_x_left1_stitch_kernel(cuRect shift_rect, cuVEC_VC<VType>& cuVEC, VType*& aux_block_values)
{
	//index in aux_block_values
	int aux_blocks_idx = blockIdx.x * blockDim.x + threadIdx.x;

	//index in cuVEC : aux_block_values stored block beginning values which must be shifted to the cell to the left
	int cell_idx = aux_blocks_idx * CUDATHREADS - 1;

	cuSZ3 n = cuVEC.n;
	cuINT3 ijk = cuINT3(cell_idx % n.x, (cell_idx / n.x) % n.y, cell_idx / (n.x*n.y));

	cuBox shift_box = cuVEC.box_from_rect_min(shift_rect);
	shift_box.e.x -= 2;
	
	if (shift_box.Contains(ijk) && cuVEC.is_not_empty(cell_idx) && cuVEC.is_not_empty(cell_idx + 1)) {

		cuVEC[cell_idx] = aux_block_values[aux_blocks_idx];
	}
}

template <typename VType>
__global__ void shift_x_right1_kernel(cuRect shift_rect, cuVEC_VC<VType>& cuVEC, VType*& aux_block_values)
{
	__shared__ VType shared_memory[CUDATHREADS];

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = cuVEC.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBox shift_box = cuVEC.box_from_rect_min(shift_rect);
	//go from shift_box.s.x + 1 to shift_box.e.x - 1 inclusive : shift_box.s.x + 1 takes value from shift_box.s.x
	shift_box.s.x++;
	shift_box.e.x--;

	if (idx < n.dim()) {

		shared_memory[threadIdx.x] = cuVEC[idx];
	}
	else shared_memory[threadIdx.x] = VType();

	//all values in this block must be transferred to shared memory before proceeding
	__syncthreads();

	//shift values within the box using shared_memory (both destination and source must not be empty).
	//We cannot shift the last element in this block since it must go to the block after this - store it in aux_block_values for later.
	//Similarly we cannot write the first element in this block since it needs a value from the previous block.
	if (threadIdx.x == CUDATHREADS - 1) aux_block_values[blockIdx.x] = shared_memory[CUDATHREADS - 1];

	if (shift_box.Contains(ijk) && cuVEC.is_not_empty(idx) && cuVEC.is_not_empty(idx - 1)) {

		if (threadIdx.x > 0) {

			cuVEC[idx] = shared_memory[threadIdx.x - 1];
		}
	}
}

template <typename VType>
__global__ void shift_x_right1_stitch_kernel(cuRect shift_rect, cuVEC_VC<VType>& cuVEC, VType*& aux_block_values)
{
	//index in aux_block_values
	int aux_blocks_idx = blockIdx.x * blockDim.x + threadIdx.x;

	//index in cuVEC : aux_block_values stored block ending values which must be shifted to the cell to the right
	int cell_idx = aux_blocks_idx * CUDATHREADS + CUDATHREADS;

	cuSZ3 n = cuVEC.n;
	cuINT3 ijk = cuINT3(cell_idx % n.x, (cell_idx / n.x) % n.y, cell_idx / (n.x*n.y));

	cuBox shift_box = cuVEC.box_from_rect_min(shift_rect);
	shift_box.s.x++;
	shift_box.e.x--;

	if (shift_box.Contains(ijk) && cuVEC.is_not_empty(cell_idx) && cuVEC.is_not_empty(cell_idx - 1)) {

		cuVEC[cell_idx] = aux_block_values[aux_blocks_idx];
	}
}

template void cuVEC_VC<float>::shift_x(size_t size, cuReal delta, cuRect shift_rect);
template void cuVEC_VC<double>::shift_x(size_t size, cuReal delta, cuRect shift_rect);

template void cuVEC_VC<cuFLT3>::shift_x(size_t size, cuReal delta, cuRect shift_rect);
template void cuVEC_VC<cuDBL3>::shift_x(size_t size, cuReal delta, cuRect shift_rect);

//shift all the values in this VEC by the given delta (units same as h). Shift values in given shift_rect (absolute coordinates).
//Also keep magnitude in each cell (e.g. use for vectorial quantities, such as magnetization, to shift only the direction).
template <typename VType>
__host__ void cuVEC_VC<VType>::shift_x(size_t size, cuReal delta, cuRect shift_rect)
{
	cuReal3 shift_debt_cpu = get_gpu_value(shift_debt);
	cuReal3 h_cpu = get_gpu_value(h);

	if (fabs(shift_debt_cpu.x + delta) < h_cpu.x) {

		//total shift not enough : bank it and return
		shift_debt_cpu.x += delta;
		set_gpu_value(shift_debt, shift_debt_cpu);
		return;
	}

	//only shift an integer number of cells : there might be a sub-cellsize remainder so just bank it to be used next time
	int cells_shift = (int)((shift_debt_cpu.x + delta) / h_cpu.x);
	shift_debt_cpu.x -= h_cpu.x * cells_shift - delta;
	set_gpu_value(shift_debt, shift_debt_cpu);
	
	if (cells_shift < 0) {

		//only shift one cell at a time - for a moving mesh algorithm it would be very unusual to have to shift by more than one cell at a time if configured properly (mesh trigger from finest mesh)
		//one-call shift routines for cells_shift > 1 are not straight-forward so not worth implementing for now
		while (cells_shift < 0) {

			shift_x_left1_kernel << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (shift_rect, *this, aux_block_values);

			size_t stitch_size = (size + CUDATHREADS) / CUDATHREADS;
			shift_x_left1_stitch_kernel << < (stitch_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (shift_rect, *this, aux_block_values);
			
			cells_shift++;
		}
	}
	else {

		while (cells_shift > 0) {

			shift_x_right1_kernel << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (shift_rect, *this, aux_block_values);

			size_t stitch_size = (size + CUDATHREADS) / CUDATHREADS;
			shift_x_right1_stitch_kernel << < (stitch_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (shift_rect, *this, aux_block_values);

			cells_shift--;
		}
	}
}

//------------------------------------------------------------------- SHIFT : y