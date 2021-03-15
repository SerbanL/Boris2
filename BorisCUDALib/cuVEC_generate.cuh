#pragma once

#include "cuVEC.h"

//--------------------------------------------VEC GENERATORS : cuVEC_generate.cuh

// GENERATE LINEAR

template bool cuVEC<float>::generate_linear(cuReal3 new_h, cuRect new_rect, cuReal3 position1, float value1, cuReal3 position2, float value2);
template bool cuVEC<double>::generate_linear(cuReal3 new_h, cuRect new_rect, cuReal3 position1, double value1, cuReal3 position2, double value2);

template bool cuVEC<cuFLT3>::generate_linear(cuReal3 new_h, cuRect new_rect, cuReal3 position1, cuFLT3 value1, cuReal3 position2, cuFLT3 value2);
template bool cuVEC<cuDBL3>::generate_linear(cuReal3 new_h, cuRect new_rect, cuReal3 position1, cuDBL3 value1, cuReal3 position2, cuDBL3 value2);

//linear : use interpolation to set values in this VEC based on projected distance between position1 and position2 and given fixed end values.
template <typename VType>
__host__ bool cuVEC<VType>::generate_linear(cuReal3 new_h, cuRect new_rect, cuReal3 position1, VType value1, cuReal3 position2, VType value2)
{
	if (!resize(new_h, new_rect)) return false;

	set_linear(position1, value1, position2, value2);

	return true;
}

template <typename VType>
__global__ void set_linear_kernel(cuVEC<VType>& vec, cuReal3 position1, VType value1, cuReal3 position2, VType value2)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < vec.n.dim()) {

		cuReal3 del_p = position2 - position1;
		cuBReal del_p_sq = del_p * del_p;

		//the absolute position of this cell center for which we are setting the value
		cuReal3 P = vec.cellidx_to_position(idx) + vec.rect.s;

		//p is the point on the line containing position1 and position2, such that the segment (P - p) is orthogonal to it
		cuReal3 p = (del_p * (P - position1)) * del_p / del_p_sq + position1;

		//use linear interpolation along the line containing position1 and position2 to set value for this cell
		vec[idx] = (value1 * (position2 - p).norm() + value2 * (p - position1).norm()) / sqrt(del_p_sq);
	}
}

// SET LINEAR

template void cuVEC<float>::set_linear(cuReal3 position1, float value1, cuReal3 position2, float value2);
template void cuVEC<double>::set_linear(cuReal3 position1, double value1, cuReal3 position2, double value2);

template void cuVEC<cuFLT3>::set_linear(cuReal3 position1, cuFLT3 value1, cuReal3 position2, cuFLT3 value2);
template void cuVEC<cuDBL3>::set_linear(cuReal3 position1, cuDBL3 value1, cuReal3 position2, cuDBL3 value2);

template void cuVEC<cuFLT4>::set_linear(cuReal3 position1, cuFLT4 value1, cuReal3 position2, cuFLT4 value2);
template void cuVEC<cuDBL4>::set_linear(cuReal3 position1, cuDBL4 value1, cuReal3 position2, cuDBL4 value2);

//similar to generate_linear except new dimensions not set
template <typename VType>
__host__ void cuVEC<VType>::set_linear(cuReal3 position1, VType value1, cuReal3 position2, VType value2)
{
	if (position1 == position2) return;

	set_linear_kernel <<< (get_gpu_value(n).dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, position1, value1, position2, value2);
}