#pragma once

#include "mcuVEC.h"
#include "launchers.h"

//--------------------------------------------AUXILIARY : mcuVEC_histo.cuh

template <typename VType>
__global__ void add_histograms_kernel(VType** arr_col, size_t num_arrays, size_t array_size)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < array_size) {

		for (int arr_idx = 1; arr_idx < num_arrays; arr_idx++) {

			arr_col[0][idx] += arr_col[arr_idx][idx];
		}
	}
}

template void mcuVEC<int, cuVEC<int>>::add_histograms(void);
template void mcuVEC<float, cuVEC<float>>::add_histograms(void);
template void mcuVEC<double, cuVEC<double>>::add_histograms(void);

template void mcuVEC<cuINT3, cuVEC<cuINT3>>::add_histograms(void);
template void mcuVEC<cuFLT3, cuVEC<cuFLT3>>::add_histograms(void);
template void mcuVEC<cuDBL3, cuVEC<cuDBL3>>::add_histograms(void);

template void mcuVEC<cuINT4, cuVEC<cuINT4>>::add_histograms(void);
template void mcuVEC<cuFLT4, cuVEC<cuFLT4>>::add_histograms(void);
template void mcuVEC<cuDBL4, cuVEC<cuDBL4>>::add_histograms(void);

template void mcuVEC<int, cuVEC_VC<int>>::add_histograms(void);
template void mcuVEC<float, cuVEC_VC<float>>::add_histograms(void);
template void mcuVEC<double, cuVEC_VC<double>>::add_histograms(void);

template void mcuVEC<cuINT3, cuVEC_VC<cuINT3>>::add_histograms(void);
template void mcuVEC<cuFLT3, cuVEC_VC<cuFLT3>>::add_histograms(void);
template void mcuVEC<cuDBL3, cuVEC_VC<cuDBL3>>::add_histograms(void);

template void mcuVEC<cuINT4, cuVEC_VC<cuINT4>>::add_histograms(void);
template void mcuVEC<cuFLT4, cuVEC_VC<cuFLT4>>::add_histograms(void);
template void mcuVEC<cuDBL4, cuVEC_VC<cuDBL4>>::add_histograms(void);

//add arrays in histogram_base_aux_col together, setting result in first one
template <typename VType, typename MType>
void mcuVEC<VType, MType>::add_histograms(void)
{
	add_histograms_kernel<cuBReal> <<< (histogram_base_aux[0].size() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (histogram_base_aux_col, histogram_base_aux_col.size(), histogram_base_aux[0].size());
}