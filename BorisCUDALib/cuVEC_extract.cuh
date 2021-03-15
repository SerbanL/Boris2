#pragma once

#include "cuVEC.h"
#include "launchers.h"

//1. Profile values only, without stencil operation, and with mesh wrap-around

//--------------------------------------------EXTRACT A LINE PROFILE - values only without stencil

template <typename VType>
__global__ void extract_profilevalues_kernel(size_t size, cuVEC<VType>& cuvec, VType* profile_gpu, cuReal3 start, cuReal3 step)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cuReal3 position = start + (cuBReal)idx * step;

		//if position is outside mesh then wrap around
		cuReal3 meshDim = cuReal3(cuvec.rect.e.x - cuvec.rect.s.x, cuvec.rect.e.y - cuvec.rect.s.y, cuvec.rect.e.z - cuvec.rect.s.z);
		position.x -= cu_floor_epsilon(position.x / meshDim.x) * meshDim.x;
		position.y -= cu_floor_epsilon(position.y / meshDim.y) * meshDim.y;
		position.z -= cu_floor_epsilon(position.z / meshDim.z) * meshDim.z;

		VType value = cuvec.weighted_average(position, cuvec.h);
		profile_gpu[idx] = value;
	}
}

template void cuVEC<char>::extract_profilevalues(size_t size, cu_arr<char>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<int>::extract_profilevalues(size_t size, cu_arr<int>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<unsigned>::extract_profilevalues(size_t size, cu_arr<unsigned>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<long>::extract_profilevalues(size_t size, cu_arr<long>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<size_t>::extract_profilevalues(size_t size, cu_arr<size_t>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<float>::extract_profilevalues(size_t size, cu_arr<float>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<double>::extract_profilevalues(size_t size, cu_arr<double>& profile_gpu, cuReal3 start, cuReal3 step);

template void cuVEC<cuINT3>::extract_profilevalues(size_t size, cu_arr<cuINT3>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ3>::extract_profilevalues(size_t size, cu_arr<cuSZ3>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT3>::extract_profilevalues(size_t size, cu_arr<cuFLT3>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL3>::extract_profilevalues(size_t size, cu_arr<cuDBL3>& profile_gpu, cuReal3 start, cuReal3 step);

template void cuVEC<cuINT4>::extract_profilevalues(size_t size, cu_arr<cuINT4>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ4>::extract_profilevalues(size_t size, cu_arr<cuSZ4>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT4>::extract_profilevalues(size_t size, cu_arr<cuFLT4>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL4>::extract_profilevalues(size_t size, cu_arr<cuDBL4>& profile_gpu, cuReal3 start, cuReal3 step);

//extract profile to a cu_arr : extract size points starting at start in the direction step; use weighted average to extract profile
template <typename VType>
__host__ void cuVEC<VType>::extract_profilevalues(size_t size, cu_arr<VType>& profile_gpu, cuReal3 start, cuReal3 step)
{
	extract_profilevalues_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(size, *this, (VType*)profile_gpu, start, step);
}

//-----------------------

template <typename VType>
__global__ void extract_profilevalues_component_x_kernel(size_t size, cuVEC<VType>& cuvec, cuBReal* profile_gpu, cuReal3 start, cuReal3 step)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cuReal3 position = start + (cuBReal)idx * step;

		//if position is outside mesh then wrap around
		cuReal3 meshDim = cuReal3(cuvec.rect.e.x - cuvec.rect.s.x, cuvec.rect.e.y - cuvec.rect.s.y, cuvec.rect.e.z - cuvec.rect.s.z);
		position.x -= cu_floor_epsilon(position.x / meshDim.x) * meshDim.x;
		position.y -= cu_floor_epsilon(position.y / meshDim.y) * meshDim.y;
		position.z -= cu_floor_epsilon(position.z / meshDim.z) * meshDim.z;

		VType value = cuvec.weighted_average(position, cuvec.h);
		profile_gpu[idx] = value.x;
	}
}

template void cuVEC<cuINT3>::extract_profilevalues_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ3>::extract_profilevalues_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT3>::extract_profilevalues_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL3>::extract_profilevalues_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);

template void cuVEC<cuINT4>::extract_profilevalues_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ4>::extract_profilevalues_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT4>::extract_profilevalues_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL4>::extract_profilevalues_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);

//these specifically apply for VType == cuReal3, allowing extraction of the x, y, z components separately
template <typename VType>
__host__ void cuVEC<VType>::extract_profilevalues_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step)
{
	extract_profilevalues_component_x_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(size, *this, profile_gpu, start, step);
}

//-----------------------

template <typename VType>
__global__ void extract_profilevalues_component_y_kernel(size_t size, cuVEC<VType>& cuvec, cuBReal* profile_gpu, cuReal3 start, cuReal3 step)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cuReal3 position = start + (cuBReal)idx * step;

		//if position is outside mesh then wrap around
		cuReal3 meshDim = cuReal3(cuvec.rect.e.x - cuvec.rect.s.x, cuvec.rect.e.y - cuvec.rect.s.y, cuvec.rect.e.z - cuvec.rect.s.z);
		position.x -= cu_floor_epsilon(position.x / meshDim.x) * meshDim.x;
		position.y -= cu_floor_epsilon(position.y / meshDim.y) * meshDim.y;
		position.z -= cu_floor_epsilon(position.z / meshDim.z) * meshDim.z;

		VType value = cuvec.weighted_average(position, cuvec.h);
		profile_gpu[idx] = value.y;
	}
}

template void cuVEC<cuINT3>::extract_profilevalues_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ3>::extract_profilevalues_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT3>::extract_profilevalues_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL3>::extract_profilevalues_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);

template void cuVEC<cuINT4>::extract_profilevalues_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ4>::extract_profilevalues_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT4>::extract_profilevalues_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL4>::extract_profilevalues_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);

template <typename VType>
__host__ void cuVEC<VType>::extract_profilevalues_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step)
{
	extract_profilevalues_component_y_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(size, *this, profile_gpu, start, step);
}

//-----------------------

template <typename VType>
__global__ void extract_profilevalues_component_z_kernel(size_t size, cuVEC<VType>& cuvec, cuBReal* profile_gpu, cuReal3 start, cuReal3 step)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cuReal3 position = start + (cuBReal)idx * step;

		//if position is outside mesh then wrap around
		cuReal3 meshDim = cuReal3(cuvec.rect.e.x - cuvec.rect.s.x, cuvec.rect.e.y - cuvec.rect.s.y, cuvec.rect.e.z - cuvec.rect.s.z);
		position.x -= cu_floor_epsilon(position.x / meshDim.x) * meshDim.x;
		position.y -= cu_floor_epsilon(position.y / meshDim.y) * meshDim.y;
		position.z -= cu_floor_epsilon(position.z / meshDim.z) * meshDim.z;

		VType value = cuvec.weighted_average(position, cuvec.h);
		profile_gpu[idx] = value.z;
	}
}

template void cuVEC<cuINT3>::extract_profilevalues_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ3>::extract_profilevalues_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT3>::extract_profilevalues_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL3>::extract_profilevalues_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);

template void cuVEC<cuINT4>::extract_profilevalues_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ4>::extract_profilevalues_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT4>::extract_profilevalues_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL4>::extract_profilevalues_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);

template <typename VType>
__host__ void cuVEC<VType>::extract_profilevalues_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step)
{
	extract_profilevalues_component_z_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(size, *this, profile_gpu, start, step);
}

//--------------------------------------------EXTRACT A LINE PROFILE - position and values with stencil operation

////////////////////////////////// SETUP KERNELS

//setup internal storage profiles and zero y values for reduction
template <typename VType>
__global__ void setup_profile_components_kernel(cuVEC<VType>& cuvec, cuBReal& step,
	cuReal2*& line_profile_component_x, cuReal2*& line_profile_component_y, cuReal2*& line_profile_component_z, size_t& line_profile_component_size, size_t*& line_profile_avpoints)
{
	//launch with 4 * line_profile_component_size size kernel
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < line_profile_component_size) {

		line_profile_component_x[idx].i = (cuBReal)idx * step;
		line_profile_component_x[idx].j = 0.0;
	}
	else if (idx < 2 * line_profile_component_size) {

		line_profile_component_y[idx - line_profile_component_size].i = (cuBReal)(idx - line_profile_component_size) * step;
		line_profile_component_y[idx - line_profile_component_size].j = 0.0;
	}
	else if (idx < 3 * line_profile_component_size) {

		line_profile_component_z[idx - 2 * line_profile_component_size].i = (cuBReal)(idx - 2 * line_profile_component_size) * step;
		line_profile_component_z[idx - 2 * line_profile_component_size].j = 0.0;
	}
	else if (idx < 4 * line_profile_component_size) {

		line_profile_avpoints[idx - 3 * line_profile_component_size] = 0;
	}
}

//setup profile and zero y values for a cu_arr (so we have cuReal2* here, not cuReal2*&)
template <typename VType>
__global__ void setup_profile_component_cuarr_kernel(cuVEC<VType>& cuvec, cuBReal& step, cuReal2* line_profile, size_t line_profile_size, size_t*& line_profile_avpoints)
{
	//launch with line_profile_size size kernel
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < line_profile_size) {

		line_profile[idx].i = (cuBReal)idx * step;
		line_profile[idx].j = 0.0;

		line_profile_avpoints[idx] = 0;
	}
}

//setup internal storage profiles and external cu_arr storage and zero y values for reduction
template <typename VType>
__global__ void setup_profile_components_combined_kernel(cuVEC<VType>& cuvec, cuBReal& step,
	cuReal2*& line_profile_component_x, cuReal2*& line_profile_component_y, cuReal2*& line_profile_component_z, cuReal2* line_profile, size_t& line_profile_component_size, size_t*& line_profile_avpoints)
{
	//launch with 5 * line_profile_component_size size kernel
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < line_profile_component_size) {

		line_profile_component_x[idx].i = (cuBReal)idx * step;
		line_profile_component_x[idx].j = 0.0;
	}
	else if (idx < 2 * line_profile_component_size) {

		line_profile_component_y[idx - line_profile_component_size].i = (cuBReal)(idx - line_profile_component_size) * step;
		line_profile_component_y[idx - line_profile_component_size].j = 0.0;
	}
	else if (idx < 3 * line_profile_component_size) {

		line_profile_component_z[idx - 2 * line_profile_component_size].i = (cuBReal)(idx - 2 * line_profile_component_size) * step;
		line_profile_component_z[idx - 2 * line_profile_component_size].j = 0.0;
	}
	else if (idx < 4 * line_profile_component_size) {

		line_profile[idx - 3 * line_profile_component_size].i = (cuBReal)(idx - 3 * line_profile_component_size) * step;
		line_profile[idx - 3 * line_profile_component_size].j = 0.0;
	}
	else if (idx < 5 * line_profile_component_size) {

		line_profile_avpoints[idx - 4 * line_profile_component_size] = 0;
	}
}

//zero cu_arr external storage, and line_profile_avpoints internal storage before reduction
template <typename VType>
__global__ void zero_profilevalues_cuarr_kernel(VType* profile_gpu, size_t*& line_profile_avpoints, size_t& line_profile_component_size)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < line_profile_component_size) {

		line_profile_avpoints[idx] = 0;
		profile_gpu[idx] = VType();
	}
}

//zero internal storage, and line_profile_avpoints internal storage before reduction
template <typename VType>
__global__ void zero_profilevalues_kernel(VType*& line_profile, size_t*& line_profile_avpoints, size_t& line_profile_component_size)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < line_profile_component_size) {

		line_profile_avpoints[idx] = 0;
		line_profile[idx] = VType();
	}
}

////////////////////////////////// COPY KERNELS

//copy line profile component from cuvec in line_profile depending on aux_integer value (0: x, 1: y, 2: z)
template <typename VType>
__global__ void copy_component_to_cuarr_kernel(cuVEC<VType>& cuvec, cuReal2* line_profile, size_t line_profile_size, size_t& aux_integer)
{
	//launch with line_profile_size size kernel
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < line_profile_size) {

		if (aux_integer == 0) {

			line_profile[idx] = cuvec.get_line_profile_component_x()[idx];
		}
		else if (aux_integer == 1) {

			line_profile[idx] = cuvec.get_line_profile_component_y()[idx];
		}
		else {

			line_profile[idx] = cuvec.get_line_profile_component_z()[idx];
		}
	}
}

////////////////////////////////// AVERAGE FINISHING KERNELS

//divide by number of averaging points to finish off : internal storage version
template <typename VType>
__global__ void finish_profileaverage_kernel(cuVEC<VType>& cuvec, size_t*& line_profile_avpoints, size_t& line_profile_avpoints_size)
{
	//launch with 3 * line_profile_avpoints_size size kernel
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < line_profile_avpoints_size) {

		size_t points = line_profile_avpoints[idx];
		if (points) cuvec.get_line_profile_component_x()[idx].j /= points;
	}
	else if (idx < 2 * line_profile_avpoints_size) {

		size_t points = line_profile_avpoints[idx - line_profile_avpoints_size];
		if (points) cuvec.get_line_profile_component_y()[idx - line_profile_avpoints_size].j /= points;
	}
	else if (idx < 3 * line_profile_avpoints_size) {

		size_t points = line_profile_avpoints[idx - 2 * line_profile_avpoints_size];
		if (points) cuvec.get_line_profile_component_z()[idx - 2 * line_profile_avpoints_size].j /= points;
	}
}

//divide by number of averaging points to finish off : cu_arr version
template <typename VType>
__global__ void finish_profileaverage_cuarr_kernel(cuVEC<VType>& cuvec, cuReal2* line_profile, size_t*& line_profile_avpoints, size_t& line_profile_avpoints_size)
{
	//launch with line_profile_avpoints_size size kernel
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < line_profile_avpoints_size) {

		size_t points = line_profile_avpoints[idx];
		if (points) line_profile[idx].j /= points;
	}
}

template <typename VType>
__global__ void finish_profileaveragevalues_cuarr_kernel(VType* profile_gpu, size_t*& line_profile_avpoints, size_t& line_profile_component_size)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < line_profile_component_size) {

		size_t points = line_profile_avpoints[idx];
		if (points) profile_gpu[idx] /= line_profile_avpoints[idx];
	}
}

template <typename VType>
__global__ void finish_profileaveragevalues_kernel(VType*& line_profile, size_t*& line_profile_avpoints, size_t& line_profile_component_size)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < line_profile_component_size) {

		size_t points = line_profile_avpoints[idx];
		if (points) line_profile[idx] /= line_profile_avpoints[idx];
	}
}

////////////////////////////////// X REDUCTION KERNEL : 1) internal storage, 2) external storage

//reduction extraction for internal storage array : x component
template <typename VType>
__global__ void extract_profile_component_x_reduction_kernel(
	cuVEC<VType>& cuvec,
	cuReal3& start, cuReal3& end, cuBReal& step, cuReal3& stencil, int& num_points, size_t*& line_profile_avpoints, int profile_idx)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value = 0.0;
	bool include_in_reduction = false;

	if (idx < num_points) {

		//profile point position
		cuReal3 position = start + (cuBReal)profile_idx * step * (end - start).normalized();

		//dimensions in number of cells, of averaging box
		cuINT3 n = stencil / cuvec.h;
		//i, j, k indexes in averaging box
		int i = idx % n.x;
		int j = (idx / n.x) % n.y;
		int k = idx / (n.x * n.y);

		//position in mesh for this kernel thread - this is where we have to get value from for reduction
		cuReal3 ker_position = position - stencil / 2 + (cuINT3(i, j, k) & cuvec.h);

		//ker_position must be inside mesh
		cuReal3 meshDim = cuvec.rect.size();
		if (ker_position >= cuReal3() && ker_position <= meshDim) {

			//form cell index in cuvec
			int cell_idx = cuvec.position_to_cellidx(ker_position);

			if (cuIsNZ(cu_GetMagnitude(cuvec[cell_idx]))) {

				value = cuvec[cell_idx].x;
				include_in_reduction = true;
			}
		}
	}

	reduction_avg(0, 1, &value, cuvec.get_line_profile_component_x()[profile_idx].j, line_profile_avpoints[profile_idx], include_in_reduction);
}

//reduction extraction for cuarr (so passed in as cuReal2*) : x component
template <typename VType>
__global__ void extract_profile_component_x_reduction_kernel(
	cuVEC<VType>& cuvec,
	cuReal3& start, cuReal3& end, cuBReal& step, cuReal3& stencil, int& num_points, size_t*& line_profile_avpoints,
	cuReal2* line_profile_component_x, int profile_idx)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value = 0.0;
	bool include_in_reduction = false;

	if (idx < num_points) {

		//profile point position
		cuReal3 position = start + (cuBReal)profile_idx * step * (end - start).normalized();

		//dimensions in number of cells, of averaging box
		cuINT3 n = stencil / cuvec.h;
		//i, j, k indexes in averaging box
		int i = idx % n.x;
		int j = (idx / n.x) % n.y;
		int k = idx / (n.x * n.y);

		//position in mesh for this kernel thread - this is where we have to get value from for reduction
		cuReal3 ker_position = position - stencil / 2 + (cuINT3(i, j, k) & cuvec.h);

		//ker_position must be inside mesh
		cuReal3 meshDim = cuvec.rect.size();
		if (ker_position >= cuReal3() && ker_position <= meshDim) {

			//form cell index in cuvec
			int cell_idx = cuvec.position_to_cellidx(ker_position);

			if (cuIsNZ(cu_GetMagnitude(cuvec[cell_idx]))) {

				value = cuvec[cell_idx].x;
				include_in_reduction = true;
			}
		}
	}

	reduction_avg(0, 1, &value, line_profile_component_x[profile_idx].j, line_profile_avpoints[profile_idx], include_in_reduction);
}

////////////////////////////////// Y REDUCTION KERNEL : 1) internal storage, 2) external storage

//reduction extraction for internal storage array : y component
template <typename VType>
__global__ void extract_profile_component_y_reduction_kernel(
	cuVEC<VType>& cuvec,
	cuReal3& start, cuReal3& end, cuBReal& step, cuReal3& stencil, int& num_points, size_t*& line_profile_avpoints, int profile_idx)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value = 0.0;
	bool include_in_reduction = false;

	if (idx < num_points) {

		//profile point position
		cuReal3 position = start + (cuBReal)profile_idx * step * (end - start).normalized();

		//dimensions in number of cells, of averaging box
		cuINT3 n = stencil / cuvec.h;
		//i, j, k indexes in averaging box
		int i = idx % n.x;
		int j = (idx / n.x) % n.y;
		int k = idx / (n.x * n.y);

		//position in mesh for this kernel thread - this is where we have to get value from for reduction
		cuReal3 ker_position = position - stencil / 2 + (cuINT3(i, j, k) & cuvec.h);

		//ker_position must be inside mesh
		cuReal3 meshDim = cuvec.rect.size();
		if (ker_position >= cuReal3() && ker_position <= meshDim) {

			//form cell index in cuvec
			int cell_idx = cuvec.position_to_cellidx(ker_position);

			if (cuIsNZ(cu_GetMagnitude(cuvec[cell_idx]))) {

				value = cuvec[cell_idx].y;
				include_in_reduction = true;
			}
		}
	}

	reduction_avg(0, 1, &value, cuvec.get_line_profile_component_y()[profile_idx].j, line_profile_avpoints[profile_idx], include_in_reduction);
}

//reduction extraction for cuarr (so passed in as cuReal2*) : y component
template <typename VType>
__global__ void extract_profile_component_y_reduction_kernel(
	cuVEC<VType>& cuvec,
	cuReal3& start, cuReal3& end, cuBReal& step, cuReal3& stencil, int& num_points, size_t*& line_profile_avpoints,
	cuReal2* line_profile_component_y, int profile_idx)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value = 0.0;
	bool include_in_reduction = false;

	if (idx < num_points) {

		//profile point position
		cuReal3 position = start + (cuBReal)profile_idx * step * (end - start).normalized();

		//dimensions in number of cells, of averaging box
		cuINT3 n = stencil / cuvec.h;
		//i, j, k indexes in averaging box
		int i = idx % n.x;
		int j = (idx / n.x) % n.y;
		int k = idx / (n.x * n.y);

		//position in mesh for this kernel thread - this is where we have to get value from for reduction
		cuReal3 ker_position = position - stencil / 2 + (cuINT3(i, j, k) & cuvec.h);

		//ker_position must be inside mesh
		cuReal3 meshDim = cuvec.rect.size();
		if (ker_position >= cuReal3() && ker_position <= meshDim) {

			//form cell index in cuvec
			int cell_idx = cuvec.position_to_cellidx(ker_position);

			if (cuIsNZ(cu_GetMagnitude(cuvec[cell_idx]))) {

				value = cuvec[cell_idx].y;
				include_in_reduction = true;
			}
		}
	}

	reduction_avg(0, 1, &value, line_profile_component_y[profile_idx].j, line_profile_avpoints[profile_idx], include_in_reduction);
}

////////////////////////////////// Z REDUCTION KERNEL : 1) internal storage, 2) external storage

//reduction extraction for internal storage array : z component
template <typename VType>
__global__ void extract_profile_component_z_reduction_kernel(
	cuVEC<VType>& cuvec,
	cuReal3& start, cuReal3& end, cuBReal& step, cuReal3& stencil, int& num_points, size_t*& line_profile_avpoints, int profile_idx)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value = 0.0;
	bool include_in_reduction = false;

	if (idx < num_points) {

		//profile point position
		cuReal3 position = start + (cuBReal)profile_idx * step * (end - start).normalized();

		//dimensions in number of cells, of averaging box
		cuINT3 n = stencil / cuvec.h;
		//i, j, k indexes in averaging box
		int i = idx % n.x;
		int j = (idx / n.x) % n.y;
		int k = idx / (n.x * n.y);

		//position in mesh for this kernel thread - this is where we have to get value from for reduction
		cuReal3 ker_position = position - stencil / 2 + (cuINT3(i, j, k) & cuvec.h);

		//ker_position must be inside mesh
		cuReal3 meshDim = cuvec.rect.size();
		if (ker_position >= cuReal3() && ker_position <= meshDim) {

			//form cell index in cuvec
			int cell_idx = cuvec.position_to_cellidx(ker_position);

			if (cuIsNZ(cu_GetMagnitude(cuvec[cell_idx]))) {

				value = cuvec[cell_idx].z;
				include_in_reduction = true;
			}
		}
	}

	reduction_avg(0, 1, &value, cuvec.get_line_profile_component_z()[profile_idx].j, line_profile_avpoints[profile_idx], include_in_reduction);
}

//reduction extraction for cuarr (so passed in as cuReal2*) : z component
template <typename VType>
__global__ void extract_profile_component_z_reduction_kernel(
	cuVEC<VType>& cuvec,
	cuReal3& start, cuReal3& end, cuBReal& step, cuReal3& stencil, int& num_points, size_t*& line_profile_avpoints,
	cuReal2* line_profile_component_z, int profile_idx)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value = 0.0;
	bool include_in_reduction = false;

	if (idx < num_points) {

		//profile point position
		cuReal3 position = start + (cuBReal)profile_idx * step * (end - start).normalized();

		//dimensions in number of cells, of averaging box
		cuINT3 n = stencil / cuvec.h;
		//i, j, k indexes in averaging box
		int i = idx % n.x;
		int j = (idx / n.x) % n.y;
		int k = idx / (n.x * n.y);

		//position in mesh for this kernel thread - this is where we have to get value from for reduction
		cuReal3 ker_position = position - stencil / 2 + (cuINT3(i, j, k) & cuvec.h);

		//ker_position must be inside mesh
		cuReal3 meshDim = cuvec.rect.size();
		if (ker_position >= cuReal3() && ker_position <= meshDim) {

			//form cell index in cuvec
			//int cell_idx = cuvec.position_to_cellidx(position - stencil / 2) + i + j * cuvec.n.x + k * cuvec.n.x*cuvec.n.y;
			int cell_idx = cuvec.position_to_cellidx(ker_position);

			if (cuIsNZ(cu_GetMagnitude(cuvec[cell_idx]))) {

				value = cuvec[cell_idx].z;
				include_in_reduction = true;
			}
		}
	}
	
	reduction_avg(0, 1, &value, line_profile_component_z[profile_idx].j, line_profile_avpoints[profile_idx], include_in_reduction);
}

////////////////////////////////// MAX REDUCTION KERNEL : x, y, and z in internal storage + detect absolute maximum at start

template <typename VType>
__global__ void extract_profile_component_max_kernel(
	cuVEC<VType>& cuvec, 
	cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, 
	size_t& line_profile_component_size,
	size_t& aux_integer)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal2*& line_profile_component_x = cuvec.get_line_profile_component_x();
	cuReal2*& line_profile_component_y = cuvec.get_line_profile_component_y();
	cuReal2*& line_profile_component_z = cuvec.get_line_profile_component_z();

	if (idx < line_profile_component_size) {

		cuReal3 position = start + (cuBReal)idx * step * (end - start).normalized();

		line_profile_component_x[idx].i = (cuBReal)idx * step;
		line_profile_component_y[idx].i = (cuBReal)idx * step;
		line_profile_component_z[idx].i = (cuBReal)idx * step;

		VType value = VType();

		//position must be inside mesh
		cuReal3 meshDim = cuvec.rect.size();
		if (position >= cuReal3() && position <= meshDim) {

			value = cuvec.average_nonempty(cuRect(position - stencil / 2, position + stencil / 2));
		}
		
		line_profile_component_x[idx].j = value.x;
		line_profile_component_y[idx].j = value.y;
		line_profile_component_z[idx].j = value.z;

		if (idx == 0) {

			if (fabs(value.x) > fabs(value.y) && fabs(value.x) > fabs(value.z)) aux_integer = 0;
			else if (fabs(value.y) > fabs(value.x) && fabs(value.y) > fabs(value.z)) aux_integer = 1;
			else aux_integer = 2;
		}
	}
}

////////////////////////////////// FULL REDUCTION KERNEL : 1) internal storage, 2) external storage

template <typename VType>
__global__ void extract_profilevalues_reduction_kernel(
	cuVEC<VType>& cuvec,
	cuReal3& start, cuReal3& end, cuBReal& step, cuReal3& stencil, int& num_points, size_t*& line_profile_avpoints,
	VType* line_profile, int profile_idx)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	VType value = 0.0;
	bool include_in_reduction = false;

	if (idx < num_points) {

		//profile point position
		cuReal3 position = start + (cuBReal)profile_idx * step * (end - start).normalized();

		//dimensions in number of cells, of averaging box
		cuINT3 n = stencil / cuvec.h;
		//i, j, k indexes in averaging box
		int i = idx % n.x;
		int j = (idx / n.x) % n.y;
		int k = idx / (n.x * n.y);

		//position in mesh for this kernel thread - this is where we have to get value from for reduction
		cuReal3 ker_position = position - stencil / 2 + (cuINT3(i, j, k) & cuvec.h);

		//ker_position must be inside mesh
		cuReal3 meshDim = cuvec.rect.size();
		if (ker_position >= cuReal3() && ker_position <= meshDim) {

			//form cell index in cuvec
			//int cell_idx = cuvec.position_to_cellidx(position - stencil / 2) + i + j * cuvec.n.x + k * cuvec.n.x*cuvec.n.y;
			int cell_idx = cuvec.position_to_cellidx(ker_position);

			if (cuIsNZ(cu_GetMagnitude(cuvec[cell_idx]))) {

				value = cuvec[cell_idx];
				include_in_reduction = true;
			}
		}
	}

	reduction_avg(0, 1, &value, line_profile[profile_idx], line_profile_avpoints[profile_idx], include_in_reduction);
}

template <typename VType>
__global__ void extract_profilevalues_reduction_kernel(
	cuVEC<VType>& cuvec, cuReal3& start, cuReal3& end, cuBReal& step, cuReal3& stencil, int& num_points, size_t*& line_profile_avpoints, int profile_idx)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	VType value = 0.0;
	bool include_in_reduction = false;

	if (idx < num_points) {

		//profile point position
		cuReal3 position = start + (cuBReal)profile_idx * step * (end - start).normalized();

		//dimensions in number of cells, of averaging box
		cuINT3 n = stencil / cuvec.h;
		//i, j, k indexes in averaging box
		int i = idx % n.x;
		int j = (idx / n.x) % n.y;
		int k = idx / (n.x * n.y);

		//position in mesh for this kernel thread - this is where we have to get value from for reduction
		cuReal3 ker_position = position - stencil / 2 + (cuINT3(i, j, k) & cuvec.h);

		//ker_position must be inside mesh
		cuReal3 meshDim = cuvec.rect.size();
		if (ker_position >= cuReal3() && ker_position <= meshDim) {

			//form cell index in cuvec
			//int cell_idx = cuvec.position_to_cellidx(position - stencil / 2) + i + j * cuvec.n.x + k * cuvec.n.x*cuvec.n.y;
			int cell_idx = cuvec.position_to_cellidx(ker_position);

			if (cuIsNZ(cu_GetMagnitude(cuvec[cell_idx]))) {

				value = cuvec[cell_idx];
				include_in_reduction = true;
			}
		}
	}

	reduction_avg(0, 1, &value, cuvec.get_line_profile()[profile_idx], line_profile_avpoints[profile_idx], include_in_reduction);
}

////////////////////////////////// FULL PROFILE KERNEL : 1) internal storage, 2) external storage

template <typename VType>
__global__ void extract_profilevalues_cuarr_kernel(
	cuVEC<VType>& cuvec,
	cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, 
	VType* profile, size_t profile_size)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < profile_size) {

		cuReal3 position = start + (cuBReal)idx * step * (end - start).normalized();

		//position must be inside mesh
		cuReal3 meshDim = cuvec.rect.size();
		if (position >= cuReal3() && position <= meshDim) {

			profile[idx] = cuvec.average_nonempty(cuRect(position - stencil / 2, position + stencil / 2));
		}
		else profile[idx] = VType();
	}
}

template <typename VType>
__global__ void extract_profilevalues_kernel(
	cuVEC<VType>& cuvec, cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, size_t& line_profile_component_size)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < line_profile_component_size) {

		cuReal3 position = start + (cuBReal)idx * step * (end - start).normalized();

		//position must be inside mesh
		cuReal3 meshDim = cuvec.rect.size();
		if (position >= cuReal3() && position <= meshDim) {

			cuvec.get_line_profile()[idx] = cuvec.average_nonempty(cuRect(position - stencil / 2, position + stencil / 2));
		}
		else {

			cuvec.get_line_profile()[idx] = VType();
		}
	}
}

////////////////////////////////// FUNCTION : X COMPONENT TO EXTERNAL GPU STORAGE

template bool cuVEC<cuFLT3>::extract_profile_component_x(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu);
template bool cuVEC<cuDBL3>::extract_profile_component_x(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu);

template <typename VType>
__host__ bool cuVEC<VType>::extract_profile_component_x(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu)
{
	if (step > 0) {

		size_t size = round((end - start).norm() / step) + 1;

		//make sure memory is allocated correctly for auxiliary arrays
		if (!allocate_profile_component_memory(size)) return false;

		cuReal3 h_cpu = get_gpu_value(h);
		cuINT3 nstencil = stencil / h_cpu;
		int num_stencil_points = nstencil.dim();

		//if stencil has more points than the profile, then better to launch multiple reduction kernels : one per each profile point
		if (num_stencil_points > size) {

			//set data in gpu memory once so we don't have to do it every iteration
			cu_obj<cuReal3> start_gpu;
			start_gpu.from_cpu(start);
			cu_obj<cuReal3> end_gpu;
			end_gpu.from_cpu(end);
			cu_obj<cuBReal> step_gpu;
			step_gpu.from_cpu(step);
			cu_obj<cuReal3> stencil_gpu;
			stencil_gpu.from_cpu(stencil);
			cu_obj<int> num_stencil_points_gpu;
			num_stencil_points_gpu.from_cpu(num_stencil_points);

			//setup the profile and zero y values for reduction
			setup_profile_component_cuarr_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, step_gpu, profile_gpu, profile_gpu.size(), line_profile_avpoints);

			//extract
			for (int idx = 0; idx < size; idx++) {

				extract_profile_component_x_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, profile_gpu, idx);
			}

			//divide by number of averaging points to finish off
			finish_profileaverage_cuarr_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, profile_gpu, line_profile_avpoints, line_profile_component_size);
		}
		else {

			//zero the profile storage
			cu_obj<cuBReal> step_gpu;
			step_gpu.from_cpu(step);
			setup_profile_components_kernel <<< (4 * size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, step_gpu, line_profile_component_x, line_profile_component_y, line_profile_component_z, line_profile_component_size, line_profile_avpoints);

			//extract all 3 components
			extract_profile_component_max_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start, end, step, stencil, line_profile_component_size, aux_integer);

			//select x component
			set_gpu_value(aux_integer, (size_t)0);

			//copy over the correct component to profile_gpu
			copy_component_to_cuarr_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, profile_gpu, profile_gpu.size(), aux_integer);
		}

		//all done
		return true;
	}
	else return false;
}

////////////////////////////////// FUNCTION : X COMPONENT TO EXTERNAL CPU STORAGE

template bool cuVEC<cuFLT3>::extract_profile_component_x(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);
template bool cuVEC<cuDBL3>::extract_profile_component_x(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);

template <typename VType>
__host__ bool cuVEC<VType>::extract_profile_component_x(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu)
{
	if (step > 0) {

		size_t size = round((end - start).norm() / step) + 1;

		//make sure memory is allocated correctly
		if (!allocate_profile_component_memory(size)) return false;

		cuReal3 h_cpu = get_gpu_value(h);
		cuINT3 nstencil = stencil / h_cpu;
		int num_stencil_points = nstencil.dim();

		//if stencil has more points than the profile, then better to launch multiple reduction kernels : one per each profile point

		if (num_stencil_points > size) {

			//set data in gpu memory once so we don't have to do it every iteration
			cu_obj<cuReal3> start_gpu;
			start_gpu.from_cpu(start);
			cu_obj<cuReal3> end_gpu;
			end_gpu.from_cpu(end);
			cu_obj<cuBReal> step_gpu;
			step_gpu.from_cpu(step);
			cu_obj<cuReal3> stencil_gpu;
			stencil_gpu.from_cpu(stencil);
			cu_obj<int> num_stencil_points_gpu;
			num_stencil_points_gpu.from_cpu(num_stencil_points);

			//zero the profile storage
			setup_profile_components_kernel <<< (4 * size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, step_gpu, line_profile_component_x, line_profile_component_y, line_profile_component_z, line_profile_component_size, line_profile_avpoints);

			//extract
			for (int idx = 0; idx < size; idx++) {

				extract_profile_component_x_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, idx);
			}

			//divide by number of averaging points to finish off
			finish_profileaverage_kernel <<< (3 * size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, line_profile_avpoints, line_profile_component_size);
		}
		else {

			extract_profile_component_max_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start, end, step, stencil, line_profile_component_size, aux_integer);
		}

		//copy over to profile_cpu
		if (profile_cpu.size() != size) if (!malloc_vector(profile_cpu, size)) return false;

		cudaError_t error = gpu_to_cpu_managed(profile_cpu.data(), line_profile_component_x, size);
		if (error != cudaSuccess) return false;

		return true;
	}
	else return false;
}

////////////////////////////////// FUNCTION : Y COMPONENT TO EXTERNAL GPU STORAGE

template bool cuVEC<cuFLT3>::extract_profile_component_y(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu);
template bool cuVEC<cuDBL3>::extract_profile_component_y(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu);

template <typename VType>
__host__ bool cuVEC<VType>::extract_profile_component_y(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu)
{
	if (step > 0) {

		size_t size = round((end - start).norm() / step) + 1;

		//make sure memory is allocated correctly for auxiliary arrays
		if (!allocate_profile_component_memory(size)) return false;

		cuReal3 h_cpu = get_gpu_value(h);
		cuINT3 nstencil = stencil / h_cpu;
		int num_stencil_points = nstencil.dim();

		//if stencil has more points than the profile, then better to launch multiple reduction kernels : one per each profile point
		if (num_stencil_points > size) {

			//set data in gpu memory once so we don't have to do it every iteration
			cu_obj<cuReal3> start_gpu;
			start_gpu.from_cpu(start);
			cu_obj<cuReal3> end_gpu;
			end_gpu.from_cpu(end);
			cu_obj<cuBReal> step_gpu;
			step_gpu.from_cpu(step);
			cu_obj<cuReal3> stencil_gpu;
			stencil_gpu.from_cpu(stencil);
			cu_obj<int> num_stencil_points_gpu;
			num_stencil_points_gpu.from_cpu(num_stencil_points);

			//setup the profile and zero y values for reduction
			setup_profile_component_cuarr_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, step_gpu, profile_gpu, profile_gpu.size(), line_profile_avpoints);

			//extract
			for (int idx = 0; idx < size; idx++) {

				extract_profile_component_y_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, profile_gpu, idx);
			}

			//divide by number of averaging points to finish off
			finish_profileaverage_cuarr_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, profile_gpu, line_profile_avpoints, line_profile_component_size);
		}
		else {

			//zero the profile storage
			cu_obj<cuBReal> step_gpu;
			step_gpu.from_cpu(step);
			setup_profile_components_kernel <<< (4 * size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, step_gpu, line_profile_component_x, line_profile_component_y, line_profile_component_z, line_profile_component_size, line_profile_avpoints);

			//extract all 3 components
			extract_profile_component_max_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start, end, step, stencil, line_profile_component_size, aux_integer);

			//select y component
			set_gpu_value(aux_integer, (size_t)1);

			//copy over the correct component to profile_gpu
			copy_component_to_cuarr_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, profile_gpu, profile_gpu.size(), aux_integer);
		}

		//all done
		return true;
	}
	else return false;
}

////////////////////////////////// FUNCTION : Y COMPONENT TO EXTERNAL CPU STORAGE

template bool cuVEC<cuFLT3>::extract_profile_component_y(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);
template bool cuVEC<cuDBL3>::extract_profile_component_y(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);

template <typename VType>
__host__ bool cuVEC<VType>::extract_profile_component_y(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu)
{
	if (step > 0) {

		size_t size = round((end - start).norm() / step) + 1;

		//make sure memory is allocated correctly
		if (!allocate_profile_component_memory(size)) return false;

		cuReal3 h_cpu = get_gpu_value(h);
		cuINT3 nstencil = stencil / h_cpu;
		int num_stencil_points = nstencil.dim();

		//if stencil has more points than the profile, then better to launch multiple reduction kernels : one per each profile point

		if (num_stencil_points > size) {

			//set data in gpu memory once so we don't have to do it every iteration
			cu_obj<cuReal3> start_gpu;
			start_gpu.from_cpu(start);
			cu_obj<cuReal3> end_gpu;
			end_gpu.from_cpu(end);
			cu_obj<cuBReal> step_gpu;
			step_gpu.from_cpu(step);
			cu_obj<cuReal3> stencil_gpu;
			stencil_gpu.from_cpu(stencil);
			cu_obj<int> num_stencil_points_gpu;
			num_stencil_points_gpu.from_cpu(num_stencil_points);

			//zero the profile storage
			setup_profile_components_kernel <<< (4 * size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, step_gpu, line_profile_component_x, line_profile_component_y, line_profile_component_z, line_profile_component_size, line_profile_avpoints);

			//extract
			for (int idx = 0; idx < size; idx++) {

				extract_profile_component_y_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, idx);
			}

			//divide by number of averaging points to finish off
			finish_profileaverage_kernel <<< (3 * size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, line_profile_avpoints, line_profile_component_size);
		}
		else {

			extract_profile_component_max_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start, end, step, stencil, line_profile_component_size, aux_integer);
		}

		//copy over to profile_cpu
		if (profile_cpu.size() != size) if (!malloc_vector(profile_cpu, size)) return false;

		cudaError_t error = gpu_to_cpu_managed(profile_cpu.data(), line_profile_component_y, size);
		if (error != cudaSuccess) return false;

		return true;
	}
	else return false;
}

////////////////////////////////// FUNCTION : Z COMPONENT TO EXTERNAL GPU STORAGE

template bool cuVEC<cuFLT3>::extract_profile_component_z(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu);
template bool cuVEC<cuDBL3>::extract_profile_component_z(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu);

template <typename VType>
__host__ bool cuVEC<VType>::extract_profile_component_z(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu)
{
	if (step > 0) {

		size_t size = round((end - start).norm() / step) + 1;

		//make sure memory is allocated correctly for auxiliary arrays
		if (!allocate_profile_component_memory(size)) return false;

		cuReal3 h_cpu = get_gpu_value(h);
		cuINT3 nstencil = stencil / h_cpu;
		int num_stencil_points = nstencil.dim();

		//if stencil has more points than the profile, then better to launch multiple reduction kernels : one per each profile point
		if (num_stencil_points > size) {

			//set data in gpu memory once so we don't have to do it every iteration
			cu_obj<cuReal3> start_gpu;
			start_gpu.from_cpu(start);
			cu_obj<cuReal3> end_gpu;
			end_gpu.from_cpu(end);
			cu_obj<cuBReal> step_gpu;
			step_gpu.from_cpu(step);
			cu_obj<cuReal3> stencil_gpu;
			stencil_gpu.from_cpu(stencil);
			cu_obj<int> num_stencil_points_gpu;
			num_stencil_points_gpu.from_cpu(num_stencil_points);

			//setup the profile and zero y values for reduction
			setup_profile_component_cuarr_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, step_gpu, profile_gpu, profile_gpu.size(), line_profile_avpoints);

			//extract
			for (int idx = 0; idx < size; idx++) {

				extract_profile_component_z_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, profile_gpu, idx);
			}

			//divide by number of averaging points to finish off
			finish_profileaverage_cuarr_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, profile_gpu, line_profile_avpoints, line_profile_component_size);
		}
		else {

			//zero the profile storage
			cu_obj<cuBReal> step_gpu;
			step_gpu.from_cpu(step);
			setup_profile_components_kernel <<< (4 * size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, step_gpu, line_profile_component_x, line_profile_component_y, line_profile_component_z, line_profile_component_size, line_profile_avpoints);

			//extract all 3 components
			extract_profile_component_max_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start, end, step, stencil, line_profile_component_size, aux_integer);

			//select z component
			set_gpu_value(aux_integer, (size_t)2);

			//copy over the correct component to profile_gpu
			copy_component_to_cuarr_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, profile_gpu, profile_gpu.size(), aux_integer);
		}

		//all done
		return true;
	}
	else return false;
}

////////////////////////////////// FUNCTION : Z COMPONENT TO EXTERNAL CPU STORAGE

template bool cuVEC<cuFLT3>::extract_profile_component_z(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);
template bool cuVEC<cuDBL3>::extract_profile_component_z(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);

template <typename VType>
__host__ bool cuVEC<VType>::extract_profile_component_z(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu)
{
	if (step > 0) {

		size_t size = round((end - start).norm() / step) + 1;

		//make sure memory is allocated correctly
		if (!allocate_profile_component_memory(size)) return false;

		cuReal3 h_cpu = get_gpu_value(h);
		cuINT3 nstencil = stencil / h_cpu;
		int num_stencil_points = nstencil.dim();

		//if stencil has more points than the profile, then better to launch multiple reduction kernels : one per each profile point

		if (num_stencil_points > size) {

			//set data in gpu memory once so we don't have to do it every iteration
			cu_obj<cuReal3> start_gpu;
			start_gpu.from_cpu(start);
			cu_obj<cuReal3> end_gpu;
			end_gpu.from_cpu(end);
			cu_obj<cuBReal> step_gpu;
			step_gpu.from_cpu(step);
			cu_obj<cuReal3> stencil_gpu;
			stencil_gpu.from_cpu(stencil);
			cu_obj<int> num_stencil_points_gpu;
			num_stencil_points_gpu.from_cpu(num_stencil_points);

			//zero the profile storage
			setup_profile_components_kernel <<< (4 * size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, step_gpu, line_profile_component_x, line_profile_component_y, line_profile_component_z, line_profile_component_size, line_profile_avpoints);

			//extract
			for (int idx = 0; idx < size; idx++) {

				extract_profile_component_z_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, idx);
			}

			//divide by number of averaging points to finish off
			finish_profileaverage_kernel <<< (3 * size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, line_profile_avpoints, line_profile_component_size);
		}
		else {

			extract_profile_component_max_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start, end, step, stencil, line_profile_component_size, aux_integer);
		}

		//copy over to profile_cpu
		if (profile_cpu.size() != size) if (!malloc_vector(profile_cpu, size)) return false;

		cudaError_t error = gpu_to_cpu_managed(profile_cpu.data(), line_profile_component_z, size);
		if (error != cudaSuccess) return false;

		return true;
	}
	else return false;
}

////////////////////////////////// FUNCTION : MAX COMPONENT TO EXTERNAL GPU STORAGE

template bool cuVEC<cuFLT3>::extract_profile_component_max(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu);
template bool cuVEC<cuDBL3>::extract_profile_component_max(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu);

template <typename VType>
__host__ bool cuVEC<VType>::extract_profile_component_max(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu)
{
	if (step > 0) {

		size_t size = round((end - start).norm() / step) + 1;

		//make sure memory is allocated correctly
		if (!allocate_profile_component_memory(size)) return false;

		cuReal3 h_cpu = get_gpu_value(h);
		cuINT3 nstencil = stencil / h_cpu;
		int num_stencil_points = nstencil.dim();

		//if stencil has more points than the profile, then better to launch multiple reduction kernels : one per each profile point

		if (num_stencil_points > size) {

			//set data in gpu memory once so we don't have to do it every iteration
			cu_obj<cuReal3> start_gpu;
			start_gpu.from_cpu(start);
			cu_obj<cuReal3> end_gpu;
			end_gpu.from_cpu(end);
			cu_obj<cuBReal> step_gpu;
			step_gpu.from_cpu(step);
			cu_obj<cuReal3> stencil_gpu;
			stencil_gpu.from_cpu(stencil);
			cu_obj<int> num_stencil_points_gpu;
			num_stencil_points_gpu.from_cpu(num_stencil_points);

			//zero the profile storage
			setup_profile_components_combined_kernel <<< (5 * size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, step_gpu, line_profile_component_x, line_profile_component_y, line_profile_component_z, profile_gpu, line_profile_component_size, line_profile_avpoints);

			//reduce first point for all 3 components so we can find maximum
			extract_profile_component_x_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, 0);

			//get first point number of averaging points count - this will be added to below, so will need to reset after
			size_t avpoints = 0;
			gpu_to_cpu_managed(&avpoints, line_profile_avpoints, 1);

			//reduce first point for all 3 components so we can find maximum
			extract_profile_component_y_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, 0);

			//reduce first point for all 3 components so we can find maximum
			extract_profile_component_z_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, 0);

			cuReal2 value_x, value_y, value_z;
			gpu_to_cpu_managed(&value_x, line_profile_component_x, 1);
			gpu_to_cpu_managed(&value_y, line_profile_component_y, 1);
			gpu_to_cpu_managed(&value_z, line_profile_component_z, 1);

			//restore line_profile_avpoints first value to correct value : will need it when finishing off profile_gpu data
			cpu_to_gpu_managed(line_profile_avpoints, &avpoints, 1);

			int component;
			if (value_x.j > value_y.j && value_x.j > value_z.j) component = 0;
			else if (value_y.j > value_x.j && value_y.j > value_z.j) component = 1;
			else component = 2;

			if (component == 0) {

				//copy over 1st point we already have, then continue
				gpu_to_gpu_managed2nd(profile_gpu.get_array(), line_profile_component_x, 1);

				for (int idx = 1; idx < size; idx++) {

					extract_profile_component_x_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
						(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, profile_gpu, idx);
				}
			}
			else if (component == 1) {

				//copy over 1st point we already have, then continue
				gpu_to_gpu_managed2nd(profile_gpu.get_array(), line_profile_component_y, 1);

				for (int idx = 1; idx < size; idx++) {

					extract_profile_component_y_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
						(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, profile_gpu, idx);
				}
			}
			else {

				//copy over 1st point we already have, then continue
				gpu_to_gpu_managed2nd(profile_gpu.get_array(), line_profile_component_z, 1);

				for (int idx = 1; idx < size; idx++) {

					extract_profile_component_z_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
						(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, profile_gpu, idx);
				}
			}

			//divide by number of averaging points to finish off
			finish_profileaverage_cuarr_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, profile_gpu, line_profile_avpoints, line_profile_component_size);
		}
		else {

			//zero the profile storage
			cu_obj<cuBReal> step_gpu;
			step_gpu.from_cpu(step);
			setup_profile_components_combined_kernel <<< (5 * size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, step_gpu, line_profile_component_x, line_profile_component_y, line_profile_component_z, profile_gpu, line_profile_component_size, line_profile_avpoints);

			//extract all 3 components
			extract_profile_component_max_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start, end, step, stencil, line_profile_component_size, aux_integer);

			//copy over the correct component to profile_gpu
			copy_component_to_cuarr_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, profile_gpu, profile_gpu.size(), aux_integer);
		}

		//all done
		return true;
	}
	else return false;
}

////////////////////////////////// FUNCTION : MAX COMPONENT TO EXTERNAL CPU STORAGE

template bool cuVEC<cuFLT3>::extract_profile_component_max(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);
template bool cuVEC<cuDBL3>::extract_profile_component_max(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);

template <typename VType>
__host__ bool cuVEC<VType>::extract_profile_component_max(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu)
{
	if (step > 0) {

		size_t size = round((end - start).norm() / step) + 1;

		//make sure profile_cpu has correct size
		if (profile_cpu.size() != size) if (!malloc_vector(profile_cpu, size)) return false;

		//make sure memory is allocated correctly
		if (!allocate_profile_component_memory(size)) return false;

		cuReal3 h_cpu = get_gpu_value(h);
		cuINT3 nstencil = stencil / h_cpu;
		int num_stencil_points = nstencil.dim();

		//if stencil has more points than the profile, then better to launch multiple reduction kernels : one per each profile point

		if (num_stencil_points > size) {

			//set data in gpu memory once so we don't have to do it every iteration
			cu_obj<cuReal3> start_gpu;
			start_gpu.from_cpu(start);
			cu_obj<cuReal3> end_gpu;
			end_gpu.from_cpu(end);
			cu_obj<cuBReal> step_gpu;
			step_gpu.from_cpu(step);
			cu_obj<cuReal3> stencil_gpu;
			stencil_gpu.from_cpu(stencil);
			cu_obj<int> num_stencil_points_gpu;
			num_stencil_points_gpu.from_cpu(num_stencil_points);

			//zero the profile storage
			setup_profile_components_kernel <<< (4 * size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, step_gpu, line_profile_component_x, line_profile_component_y, line_profile_component_z, line_profile_component_size, line_profile_avpoints);

			//reduce first point for all 3 components so we can find maximum
			extract_profile_component_x_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, 0);
			
			//get first point number of averaging points count - this will be added to below, so will need to reset after
			size_t avpoints = 0;
			gpu_to_cpu_managed(&avpoints, line_profile_avpoints, 1);

			//reduce first point for all 3 components so we can find maximum
			extract_profile_component_y_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, 0);

			//reduce first point for all 3 components so we can find maximum
			extract_profile_component_z_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, 0);

			cuReal2 value_x, value_y, value_z;
			gpu_to_cpu_managed(&value_x, line_profile_component_x, 1);
			gpu_to_cpu_managed(&value_y, line_profile_component_y, 1);
			gpu_to_cpu_managed(&value_z, line_profile_component_z, 1);
			
			//restore line_profile_avpoints first value to correct value : will need it when finishing off profile data
			cpu_to_gpu_managed(line_profile_avpoints, &avpoints, 1);

			int component;
			if (value_x.j > value_y.j && value_x.j > value_z.j) component = 0;
			else if (value_y.j > value_x.j && value_y.j > value_z.j) component = 1;
			else component = 2;

			if (component == 0) {

				//extract
				for (int idx = 1; idx < size; idx++) {

					extract_profile_component_x_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
						(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, idx);
				}

				//divide by number of averaging points to finish off
				finish_profileaverage_kernel <<< (3 * size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(*this, line_profile_avpoints, line_profile_component_size);

				//copy over to profile_cpu
				cudaError_t error = gpu_to_cpu_managed(profile_cpu.data(), line_profile_component_x, size);
				if (error != cudaSuccess) return false;
			}
			else if (component == 1) {

				//extract
				for (int idx = 1; idx < size; idx++) {

					extract_profile_component_y_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
						(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, idx);
				}

				//divide by number of averaging points to finish off
				finish_profileaverage_kernel <<< (3 * size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(*this, line_profile_avpoints, line_profile_component_size);

				//copy over to profile_cpu
				cudaError_t error = gpu_to_cpu_managed(profile_cpu.data(), line_profile_component_y, size);
				if (error != cudaSuccess) return false;
			}
			else {

				//extract
				for (int idx = 1; idx < size; idx++) {

					extract_profile_component_z_reduction_kernel << < (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >
						(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, idx);
				}

				//divide by number of averaging points to finish off
				finish_profileaverage_kernel <<< (3 * size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(*this, line_profile_avpoints, line_profile_component_size);

				//copy over to profile_cpu
				cudaError_t error = gpu_to_cpu_managed(profile_cpu.data(), line_profile_component_z, size);
				if (error != cudaSuccess) return false;
			}

			return true;
		}
		else {

			extract_profile_component_max_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start, end, step, stencil, line_profile_component_size, aux_integer);

			size_t component_max = get_gpu_value(aux_integer);

			if (component_max == 0) {

				//copy over to profile_cpu
				cudaError_t error = gpu_to_cpu_managed(profile_cpu.data(), line_profile_component_x, size);
				if (error != cudaSuccess) return false;
			}
			else if (component_max == 1) {

				//copy over to profile_cpu
				cudaError_t error = gpu_to_cpu_managed(profile_cpu.data(), line_profile_component_y, size);
				if (error != cudaSuccess) return false;
			}
			else {

				//copy over to profile_cpu
				cudaError_t error = gpu_to_cpu_managed(profile_cpu.data(), line_profile_component_z, size);
				if (error != cudaSuccess) return false;
			}

			return true;
		}
	}
	else return false;
}

template <typename VType>
__host__ bool cuVEC<VType>::extract_profile(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<VType>& profile_gpu)
{
	if (step > 0) {

		size_t size = round((end - start).norm() / step) + 1;

		//make sure memory is allocated correctly for auxiliary arrays
		if (!allocate_profile_component_memory(size)) return false;

		cuReal3 h_cpu = get_gpu_value(h);
		cuINT3 nstencil = stencil / h_cpu;
		int num_stencil_points = nstencil.dim();

		//if stencil has more points than the profile, then better to launch multiple reduction kernels : one per each profile point
		if (num_stencil_points > size) {

			//set data in gpu memory once so we don't have to do it every iteration
			cu_obj<cuReal3> start_gpu;
			start_gpu.from_cpu(start);
			cu_obj<cuReal3> end_gpu;
			end_gpu.from_cpu(end);
			cu_obj<cuBReal> step_gpu;
			step_gpu.from_cpu(step);
			cu_obj<cuReal3> stencil_gpu;
			stencil_gpu.from_cpu(stencil);
			cu_obj<int> num_stencil_points_gpu;
			num_stencil_points_gpu.from_cpu(num_stencil_points);

			//zero values for reduction
			zero_profilevalues_cuarr_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(profile_gpu, line_profile_avpoints, line_profile_component_size);

			//extract
			for (int idx = 0; idx < size; idx++) {

				extract_profilevalues_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, profile_gpu, idx);
			}

			//divide by number of averaging points to finish off
			finish_profileaveragevalues_cuarr_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(profile_gpu, line_profile_avpoints, line_profile_component_size);
		}
		else {

			//extract profile
			extract_profilevalues_cuarr_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start, end, step, stencil, profile_gpu, profile_gpu.size());
		}

		//all done
		return true;
	}
	else return false;
}

template <typename VType>
template <typename SType>
__host__ bool cuVEC<VType>::extract_profile(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<SType>& profile_cpu)
{
	if (step > 0) {

		size_t size = round((end - start).norm() / step) + 1;

		//make sure memory is allocated correctly for auxiliary arrays
		if (!allocate_profile_component_memory(size)) return false;

		cuReal3 h_cpu = get_gpu_value(h);
		cuINT3 nstencil = stencil / h_cpu;
		int num_stencil_points = nstencil.dim();

		//if stencil has more points than the profile, then better to launch multiple reduction kernels : one per each profile point
		if (num_stencil_points > size) {

			//set data in gpu memory once so we don't have to do it every iteration
			cu_obj<cuReal3> start_gpu;
			start_gpu.from_cpu(start);
			cu_obj<cuReal3> end_gpu;
			end_gpu.from_cpu(end);
			cu_obj<cuBReal> step_gpu;
			step_gpu.from_cpu(step);
			cu_obj<cuReal3> stencil_gpu;
			stencil_gpu.from_cpu(stencil);
			cu_obj<int> num_stencil_points_gpu;
			num_stencil_points_gpu.from_cpu(num_stencil_points);

			//zero values for reduction
			zero_profilevalues_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(line_profile, line_profile_avpoints, line_profile_component_size);

			//extract
			for (int idx = 0; idx < size; idx++) {

				extract_profilevalues_reduction_kernel <<< (num_stencil_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(*this, start_gpu, end_gpu, step_gpu, stencil_gpu, num_stencil_points_gpu, line_profile_avpoints, idx);
			}

			//divide by number of averaging points to finish off
			finish_profileaveragevalues_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(line_profile, line_profile_avpoints, line_profile_component_size);
		}
		else {

			//extract profile
			extract_profilevalues_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(*this, start, end, step, stencil, line_profile_component_size);
		}

		//copy over to profile_cpu
		if (profile_cpu.size() != size) if (!malloc_vector(profile_cpu, size)) return false;

		cudaError_t error = gpu_to_cpu_managed(profile_cpu.data(), line_profile, size);
		if (error != cudaSuccess) return false;

		//all done
		return true;
	}
	else return false;
}