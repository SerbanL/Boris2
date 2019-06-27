#pragma once

#include "launchers.h"

#include "cuTypes.h"
#include "cuArray.h"

///////////////////////////////////////////////////////////////////////////////

template <typename Type>
__global__ void cu_transpose_xy_kernel(Type* input, Type* output, cuSZ3& inembedN, cuSZ3& inN, cuSZ3& inNtransp)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < inembedN.dim()) {

		//ijk index in the input array recovered from the linear index idx which spans only a size of inembedN.dim(), and is meant to span the entire inembedN region (which start from 0,0,0)
		cuINT3 ijk = cuINT3(idx % inembedN.x, (idx / inembedN.x) % inembedN.y, idx / (inembedN.x*inembedN.y));

		//get linear index in the input array
		int idx_in = ijk.i + ijk.j * inN.x + ijk.k * inN.x * inN.y;

		//output index in the transposed array
		int idx_out = ijk.j + ijk.i * inNtransp.y + ijk.k * inNtransp.x * inNtransp.y;

		output[idx_out] = input[idx_in];
	}
}

template <typename Type>
__global__ void cu_transpose3_xy_kernel(
	Type* input_x, Type* input_y, Type* input_z, 
	Type* output_x, Type* output_y, Type* output_z,
	cuSZ3& inembedN, cuSZ3& inN, cuSZ3& inNtransp)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < inembedN.dim()) {

		//ijk index in the input array recovered from the linear index idx which spans only a size of inembedN.dim(), and is meant to span the entire inembedN region (which start from 0,0,0)
		cuINT3 ijk = cuINT3(idx % inembedN.x, (idx / inembedN.x) % inembedN.y, idx / (inembedN.x*inembedN.y));

		//get linear index in the input array
		int idx_in = ijk.i + ijk.j * inN.x + ijk.k * inN.x * inN.y;

		//output index in the transposed array
		int idx_out = ijk.j + ijk.i * inNtransp.y + ijk.k * inNtransp.x * inNtransp.y;

		output_x[idx_out] = input_x[idx_in];
		output_y[idx_out] = input_y[idx_in];
		output_z[idx_out] = input_z[idx_in];
	}
}

//out-of-place transpose in the xy planes.
//We have:
//1. Input logical array size cuinN
//2. Input embedded array size cuinembedN
//3. The output transposed array size is obtained from cuinNtransp - this is untransposed, thus the transposed array size is (cuinNtransp.y, cuinNtransp.x, cuinNtransp.z)
//
//We are only transposing the embedded array, thus size = cuinembedN.dim(), into the output array. The output array must be large enough to accept the transposed points (and certainly of linear size <= cuinNtransp.dim()).
//Normally cuinN = cuinNtransp, but there might be special cases where this is not true, e.g. cuinembedN = cuinN transposed into a larger array cuinNtransp.
template <typename Type>
__host__ void cu_transpose_xy(size_t size, cu_arr<Type>& input, cu_arr<Type>& output, cu_obj<cuSZ3>& cuinembedN, cu_obj<cuSZ3>& cuinN, cu_obj<cuSZ3>& cuinNtransp)
{
	cu_transpose_xy_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> ((Type*)input, (Type*)output, cuinembedN, cuinN, cuinNtransp);
}

//as above but 3-component version
template <typename Type>
__host__ void cu_transpose_xy(
	size_t size, 
	cu_arr<Type>& input_x, cu_arr<Type>& input_y, cu_arr<Type>& input_z, 
	cu_arr<Type>& output_x, cu_arr<Type>& output_y, cu_arr<Type>& output_z,
	cu_obj<cuSZ3>& cuinembedN, cu_obj<cuSZ3>& cuinN, cu_obj<cuSZ3>& cuinNtransp)
{
	cu_transpose3_xy_kernel << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (
		(Type*)input_x, (Type*)input_y, (Type*)input_z,
		(Type*)output_x, (Type*)output_y, (Type*)output_z,
		cuinembedN, cuinN, cuinNtransp);
}

///////////////////////////////////////////////////////////////////////////////

template <typename Type>
__global__ void cu_zeropad_kernel(Type* input, cuSZ3& N, cuRect& region)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//size of the region to zero
	cuSZ3 regionN = region.size();

	if (idx < regionN.dim()) {

		//ijk index in the input array
		cuINT3 ijk = cuINT3(idx % regionN.x, (idx / regionN.x) % regionN.y, idx / (regionN.x*regionN.y)) + region.s;

		input[ijk.i + ijk.j * N.x + ijk.k * N.x*N.y] = Type();
	}
}

template <typename Type>
__global__ void cu_zeropad3_kernel(Type* input_x, Type* input_y, Type* input_z, cuSZ3& N, cuRect& region)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//size of the region to zero
	cuSZ3 regionN = region.size();

	if (idx < regionN.dim()) {

		//ijk index in the input array
		cuINT3 ijk = cuINT3(idx % regionN.x, (idx / regionN.x) % regionN.y, idx / (regionN.x*regionN.y)) + region.s;

		input_x[ijk.i + ijk.j * N.x + ijk.k * N.x*N.y] = Type();
		input_y[ijk.i + ijk.j * N.x + ijk.k * N.x*N.y] = Type();
		input_z[ijk.i + ijk.j * N.x + ijk.k * N.x*N.y] = Type();
	}
}

//In the input array zero the specified region
//the input array has xyz dimensions specified by cuN
//size must be the linear size of the region i.e. size = region.size().dim()
template <typename Type>
__host__ void cu_zeropad(size_t size, cu_arr<Type>& input, cu_obj<cuSZ3>& cuN, cu_obj<cuRect>& region)
{
	cu_zeropad_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> ((Type*)input, cuN, region);
}

//as above but 3-component version
template <typename Type>
__host__ void cu_zeropad(size_t size, cu_arr<Type>& input_x, cu_arr<Type>& input_y, cu_arr<Type>& input_z, cu_obj<cuSZ3>& cuN, cu_obj<cuRect>& region)
{
	cu_zeropad3_kernel << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > ((Type*)input_x, (Type*)input_y, (Type*)input_z, cuN, region);
}

///////////////////////////////////////////////////////////////////////////////

template <typename Type>
__global__ void cu_copy_kernel(Type* input, Type* output, cuSZ3& Nin, cuSZ3& Nout)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Nin.dim()) {

		//ijk index in the input array
		cuINT3 ijk_in = cuINT3(idx % Nin.x, (idx / Nin.x) % Nin.y, idx / (Nin.x*Nin.y));
		
		//calculate output index to embed value in output array
		int idx_out = ijk_in.i + ijk_in.j * Nout.x + ijk_in.k * Nout.x * Nout.y;

		output[idx_out] = input[idx];
	}
}

template <typename Type>
__global__ void cu_copy3_kernel(
	Type* input_x, Type* input_y, Type* input_z, 
	Type* output_x, Type* output_y, Type* output_z,
	cuSZ3& Nin, cuSZ3& Nout)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Nin.dim()) {

		//ijk index in the input array
		cuINT3 ijk_in = cuINT3(idx % Nin.x, (idx / Nin.x) % Nin.y, idx / (Nin.x*Nin.y));

		//calculate output index to embed value in output array
		int idx_out = ijk_in.i + ijk_in.j * Nout.x + ijk_in.k * Nout.x * Nout.y;

		output_x[idx_out] = input_x[idx];
		output_y[idx_out] = input_y[idx];
		output_z[idx_out] = input_z[idx];
	}
}


//copy input array by embedding in output array, where input array has logical dimensions cuNin and output array has logical dimensions cuNout, such that cuNin <= cuNout
//size = cuNin.dim()
template <typename Type>
__host__ void cu_copy(size_t size, cu_arr<Type>& input, cu_arr<Type>& output, cu_obj<cuSZ3>& cuNin, cu_obj<cuSZ3>& cuNout)
{
	cu_copy_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> ((Type*)input, (Type*)output, cuNin, cuNout);
}

//as above but 3-component version
template <typename Type>
__host__ void cu_copy3(
	size_t size, 
	cu_arr<Type>& input_x, cu_arr<Type>& input_y, cu_arr<Type>& input_z,
	cu_arr<Type>& output_x, cu_arr<Type>& output_y, cu_arr<Type>& output_z, 
	cu_obj<cuSZ3>& cuNin, cu_obj<cuSZ3>& cuNout)
{
	cu_copy_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (
		(Type*)input_x, (Type*)input_y, (Type*)input_z, 
		(Type*)output_x, (Type*)output_y, (Type*)output_z,
		cuNin, cuNout);
}

///////////////////////////////////////////////////////////////////////////////
