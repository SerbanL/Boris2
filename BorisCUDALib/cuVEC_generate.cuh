#pragma once

#include "cuVEC.h"

//--------------------------------------------VEC GENERATORS : cuVEC_generate.cuh

// GENERATE LINEAR

template bool cuVEC<float>::generate_linear(cuReal3 new_h, cuRect new_rect, cuRect contact1, float value1, cuRect contact2, float value2);
template bool cuVEC<double>::generate_linear(cuReal3 new_h, cuRect new_rect, cuRect contact1, double value1, cuRect contact2, double value2);

template bool cuVEC<cuFLT3>::generate_linear(cuReal3 new_h, cuRect new_rect, cuRect contact1, cuFLT3 value1, cuRect contact2, cuFLT3 value2);
template bool cuVEC<cuDBL3>::generate_linear(cuReal3 new_h, cuRect new_rect, cuRect contact1, cuDBL3 value1, cuRect contact2, cuDBL3 value2);

//linear : use interpolation to set values in this VEC based on projected distance between position1 and position2 and given fixed end values.
template <typename VType>
__host__ bool cuVEC<VType>::generate_linear(cuReal3 new_h, cuRect new_rect, cuRect contact1, VType value1, cuRect contact2, VType value2)
{
	if (!resize(new_h, new_rect)) return false;

	set_linear(contact1, value1, contact2, value2);

	return true;
}

template <typename VType>
__global__ void set_linear_kernel(cuVEC<VType>& vec, cuRect contact1, VType value1, cuRect contact2, VType value2, cuReal2 degeneracy)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < vec.n.dim()) {

		//the absolute position of this cell center for which we are setting the value
		cuReal3 P = vec.cellidx_to_position(idx) + vec.rect.s;

		//distance to position 1
		cuBReal d1 = contact1.get_closest_distance(P);
		//distance to position 2
		cuBReal d2 = contact2.get_closest_distance(P);
		//total distance
		cuBReal d = d1 + d2;

		//value to set
		VType value = (value1 + value2) / 2;
		if (d) value = (value1 * (d - d1) + value2 * (d - d2)) / d;

		if (degeneracy.second > 1) {
			//degeneracy set : unweighted average, and first index (0) sets value, the rest add
			if (degeneracy.first) vec[idx] += value / degeneracy.second;
			//Caller must make sure index 0 goes first
			else vec[idx] = value / degeneracy.second;
		}
		else {
			//no degeneracy : just set value
			vec[idx] = value;
		}
	}
}

// SET LINEAR

template void cuVEC<float>::set_linear(cuRect contact1, float value1, cuRect contact2, float value2, cuReal2 degeneracy);
template void cuVEC<double>::set_linear(cuRect contact1, double value1, cuRect contact2, double value2, cuReal2 degeneracy);

template void cuVEC<cuFLT3>::set_linear(cuRect contact1, cuFLT3 value1, cuRect contact2, cuFLT3 value2, cuReal2 degeneracy);
template void cuVEC<cuDBL3>::set_linear(cuRect contact1, cuDBL3 value1, cuRect contact2, cuDBL3 value2, cuReal2 degeneracy);

//similar to generate_linear except new dimensions not set
template <typename VType>
__host__ void cuVEC<VType>::set_linear(cuRect contact1, VType value1, cuRect contact2, VType value2, cuReal2 degeneracy)
{
	set_linear_kernel <<< (get_gpu_value(n).dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, contact1, value1, contact2, value2, degeneracy);
}