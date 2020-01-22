#pragma once

#include "cuVEC.h"
#include "launchers.h"

//--------------------------------------------EXTRACT A LINE PROFILE

template <typename VType>
__global__ void extract_profile_kernel(size_t size, cuVEC<VType>& cuvec, VType* profile_gpu, cuReal3 start, cuReal3 step)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cuReal3 position = start + ((cuBReal)idx + 0.5) * step;

		VType value = cuvec.weighted_average(position, cuvec.h);
		profile_gpu[idx] = value;
	}
}

template void cuVEC<char>::extract_profile(size_t size, cu_arr<char>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<int>::extract_profile(size_t size, cu_arr<int>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<unsigned>::extract_profile(size_t size, cu_arr<unsigned>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<long>::extract_profile(size_t size, cu_arr<long>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<size_t>::extract_profile(size_t size, cu_arr<size_t>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<float>::extract_profile(size_t size, cu_arr<float>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<double>::extract_profile(size_t size, cu_arr<double>& profile_gpu, cuReal3 start, cuReal3 step);

template void cuVEC<cuINT3>::extract_profile(size_t size, cu_arr<cuINT3>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ3>::extract_profile(size_t size, cu_arr<cuSZ3>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT3>::extract_profile(size_t size, cu_arr<cuFLT3>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL3>::extract_profile(size_t size, cu_arr<cuDBL3>& profile_gpu, cuReal3 start, cuReal3 step);

//extract profile to a cu_arr : extract size points starting at start in the direction step; use weighted average to extract profile
template <typename VType>
__host__ void cuVEC<VType>::extract_profile(size_t size, cu_arr<VType>& profile_gpu, cuReal3 start, cuReal3 step)
{
	extract_profile_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (size, *this, (VType*)profile_gpu, start, step);
}

//-----------------------

template <typename VType>
__global__ void extract_profile_component_x_kernel(size_t size, cuVEC<VType>& cuvec, cuBReal* profile_gpu, cuReal3 start, cuReal3 step)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cuReal3 position = start + ((cuBReal)idx + 0.5) * step;

		VType value = cuvec.weighted_average(position, cuvec.h);
		profile_gpu[idx] = value.x;
	}
}

template void cuVEC<cuINT3>::extract_profile_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ3>::extract_profile_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT3>::extract_profile_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL3>::extract_profile_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);

//these specifically apply for VType == cuReal3, allowing extraction of the x, y, z components separately
template <typename VType>
__host__ void cuVEC<VType>::extract_profile_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step)
{
	extract_profile_component_x_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (size, *this, profile_gpu, start, step);
}

//-----------------------

template <typename VType>
__global__ void extract_profile_component_y_kernel(size_t size, cuVEC<VType>& cuvec, cuBReal* profile_gpu, cuReal3 start, cuReal3 step)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cuReal3 position = start + ((cuBReal)idx + 0.5) * step;

		VType value = cuvec.weighted_average(position, cuvec.h);
		profile_gpu[idx] = value.y;
	}
}

template void cuVEC<cuINT3>::extract_profile_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ3>::extract_profile_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT3>::extract_profile_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL3>::extract_profile_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);

template <typename VType>
__host__ void cuVEC<VType>::extract_profile_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step)
{
	extract_profile_component_y_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (size, *this, profile_gpu, start, step);
}

//-----------------------

template <typename VType>
__global__ void extract_profile_component_z_kernel(size_t size, cuVEC<VType>& cuvec, cuBReal* profile_gpu, cuReal3 start, cuReal3 step)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cuReal3 position = start + ((cuBReal)idx + 0.5) * step;

		VType value = cuvec.weighted_average(position, cuvec.h);
		profile_gpu[idx] = value.z;
	}
}

template void cuVEC<cuINT3>::extract_profile_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ3>::extract_profile_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT3>::extract_profile_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL3>::extract_profile_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);

template <typename VType>
__host__ void cuVEC<VType>::extract_profile_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step)
{
	extract_profile_component_z_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (size, *this, profile_gpu, start, step);
}

//-----------------------

//--------------------------------------------EXTRACT A LINE PROFILE USING POINTS (no weighted average)

template <typename VType>
__global__ void extract_profilepoints_kernel(size_t size, cuVEC<VType>& cuvec, VType* profile_gpu, cuReal3 start, cuReal3 step)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cuReal3 position = start + ((cuBReal)idx + 0.5) * step;

		VType value = cuvec[position];
		profile_gpu[idx] = value;
	}
}

template void cuVEC<char>::extract_profilepoints(size_t size, cu_arr<char>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<int>::extract_profilepoints(size_t size, cu_arr<int>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<unsigned>::extract_profilepoints(size_t size, cu_arr<unsigned>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<long>::extract_profilepoints(size_t size, cu_arr<long>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<size_t>::extract_profilepoints(size_t size, cu_arr<size_t>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<float>::extract_profilepoints(size_t size, cu_arr<float>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<double>::extract_profilepoints(size_t size, cu_arr<double>& profile_gpu, cuReal3 start, cuReal3 step);

template void cuVEC<cuINT3>::extract_profilepoints(size_t size, cu_arr<cuINT3>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ3>::extract_profilepoints(size_t size, cu_arr<cuSZ3>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT3>::extract_profilepoints(size_t size, cu_arr<cuFLT3>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL3>::extract_profilepoints(size_t size, cu_arr<cuDBL3>& profile_gpu, cuReal3 start, cuReal3 step);

//extract profile to a cu_arr : extract size points starting at start in the direction step; use weighted average to extract profile
template <typename VType>
__host__ void cuVEC<VType>::extract_profilepoints(size_t size, cu_arr<VType>& profile_gpu, cuReal3 start, cuReal3 step)
{
	extract_profilepoints_kernel << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (size, *this, (VType*)profile_gpu, start, step);
}

//-----------------------

template <typename VType>
__global__ void extract_profilepoints_component_x_kernel(size_t size, cuVEC<VType>& cuvec, cuBReal* profile_gpu, cuReal3 start, cuReal3 step)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cuReal3 position = start + ((cuBReal)idx + 0.5) * step;

		VType value = cuvec[position];
		profile_gpu[idx] = value.x;
	}
}

template void cuVEC<cuINT3>::extract_profilepoints_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ3>::extract_profilepoints_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT3>::extract_profilepoints_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL3>::extract_profilepoints_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);

//these specifically apply for VType == cuReal3, allowing extraction of the x, y, z components separately
template <typename VType>
__host__ void cuVEC<VType>::extract_profilepoints_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step)
{
	extract_profilepoints_component_x_kernel << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (size, *this, profile_gpu, start, step);
}

//-----------------------

template <typename VType>
__global__ void extract_profilepoints_component_y_kernel(size_t size, cuVEC<VType>& cuvec, cuBReal* profile_gpu, cuReal3 start, cuReal3 step)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cuReal3 position = start + ((cuBReal)idx + 0.5) * step;

		VType value = cuvec[position];
		profile_gpu[idx] = value.y;
	}
}

template void cuVEC<cuINT3>::extract_profilepoints_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ3>::extract_profilepoints_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT3>::extract_profilepoints_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL3>::extract_profilepoints_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);

template <typename VType>
__host__ void cuVEC<VType>::extract_profilepoints_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step)
{
	extract_profilepoints_component_y_kernel << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (size, *this, profile_gpu, start, step);
}

//-----------------------

template <typename VType>
__global__ void extract_profilepoints_component_z_kernel(size_t size, cuVEC<VType>& cuvec, cuBReal* profile_gpu, cuReal3 start, cuReal3 step)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cuReal3 position = start + ((cuBReal)idx + 0.5) * step;

		VType value = cuvec[position];
		profile_gpu[idx] = value.z;
	}
}

template void cuVEC<cuINT3>::extract_profilepoints_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuSZ3>::extract_profilepoints_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuFLT3>::extract_profilepoints_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
template void cuVEC<cuDBL3>::extract_profilepoints_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);

template <typename VType>
__host__ void cuVEC<VType>::extract_profilepoints_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step)
{
	extract_profilepoints_component_z_kernel << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (size, *this, profile_gpu, start, step);
}

//-----------------------