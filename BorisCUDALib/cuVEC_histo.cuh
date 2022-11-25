#pragma once

#include "cuVEC.h"
#include "launchers.h"

//-------------------------------------------- Auxiliary

//zero cu_arr external storage, and line_profile_avpoints internal storage before reduction
inline __global__ void zero_histogramvalues_kernel(cuBReal*& histogram, size_t& histogram_size, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < histogram_size) {

		histogram[idx] = 0.0;

		if (idx == 0) min = 1e38;
		if (idx == 1) max = -1e38;
	}
}

template <typename VType>
__global__ void zero_histogram_preaverage_values_kernel(VType*& histogram_preaverage, size_t& histogram_preaverage_size, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < histogram_preaverage_size) {

		histogram_preaverage[idx] = VType();

		if (idx == 0) min = 1e38;
		if (idx == 1) max = -1e38;
	}
}

template <typename VType>
inline __global__ void get_minmax_angles_kernel(cuVEC<VType>& vec, VType ndir, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value = 0.0;
	bool include_in_reduction = false;

	if (idx < vec.linear_size()) {

		if (vec.is_not_empty(idx)) {

			value = acos(vec[idx].normalized() * ndir);
			include_in_reduction = true;
		}
	}

	reduction_minmax(0, 1, &value, min, max, include_in_reduction);
}

//-------------------------------------------- Pre-averaging using segment reduction

template <typename VType>
__global__ void preaverage_kernel(
	cuVEC<VType>& cuvec, VType*& histogram_preaverage, cuINT3 num_av_cells, cuINT3 av_cell_dims)
{
	//launched with num_av_cells.dim(), (1024 or CUDATHREADS) kernel dimensions
	//i.e there are num_av_cells.dim() segments, each of size of av_cell_dims.dim()

	//segment size
	size_t K = av_cell_dims.dim();

	//linear index in this segment, starting at threadIdx.x value
	int linear_idx = threadIdx.x;

	//partial segment sum in this thread
	VType sum = 0;

	//each segment receives up to 1024 worker threads. first use them to load all input data in current segment.
	while (linear_idx < K) {

		//segment ijk values
		int i_seg = blockIdx.x % num_av_cells.x;
		int j_seg = (blockIdx.x / num_av_cells.x) % num_av_cells.y;
		int k_seg = blockIdx.x / (num_av_cells.x * num_av_cells.y);

		//convert linear segment index to cuvec ijk index for this segment
		int i = linear_idx % av_cell_dims.x;
		int j = (linear_idx / av_cell_dims.x) % av_cell_dims.y;
		int k = linear_idx / (av_cell_dims.x * av_cell_dims.y);

		//finally required ijk index in cuvec
		cuINT3 ijk = cuINT3(i_seg * av_cell_dims.i + i, j_seg * av_cell_dims.j + j, k_seg * av_cell_dims.k + k);

		if (ijk < cuvec.n && cuvec.is_not_empty(ijk)) sum += cuvec[ijk] / K;

		linear_idx += blockDim.x;
	}

	//now reduced all partial segment sums in this block
	if (blockDim.x == 1024) reduction_sum_blocksize1024(0, 1, &sum, histogram_preaverage[blockIdx.x]);
	else reduction_sum(0, 1, &sum, histogram_preaverage[blockIdx.x]);
}

template <typename VType>
__global__ void average_preaverage_kernel(VType*& histogram_preaverage, size_t& histogram_preaverage_size, cuBReal& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value = 0.0;
	bool include_in_reduction = false;

	if (idx < histogram_preaverage_size) {

		value = cu_GetMagnitude(histogram_preaverage[idx]);
		include_in_reduction = true;
	}

	reduction_avg(0, 1, &value, average_value, points_count, include_in_reduction);
}

//-------------------------------------------- Min-max reduction on pre-averaged data

template <typename VType>
__global__ void preaveraged_minmax(
	VType*& histogram_preaverage, size_t& histogram_preaverage_size, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value = 0.0;
	bool include_in_reduction = false;

	if (idx < histogram_preaverage_size) {

		value = cu_GetMagnitude(histogram_preaverage[idx]);
		include_in_reduction = true;
	}

	reduction_minmax(0, 1, &value, min, max, include_in_reduction);
}

template <typename VType>
__global__ void preaveraged_minmax_angles(
	VType*& histogram_preaverage, size_t& histogram_preaverage_size, VType ndir, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value = 0.0;
	bool include_in_reduction = false;

	if (idx < histogram_preaverage_size) {

		if (!histogram_preaverage[idx].IsNull()) value = acos(histogram_preaverage[idx].normalized() * ndir);
		include_in_reduction = true;
	}

	reduction_minmax(0, 1, &value, min, max, include_in_reduction);
}

//-------------------------------------------- Histogram kernels (magnitude)

//extract histogram directly from mesh without pre-averaging
template <typename VType>
__global__ void get_mag_histogram_kernel(
	cuVEC<VType>& vec, int num_bins, cuBReal bin, cuBReal min, cuBReal max, cuBReal*& histogram, size_t num_nonempty_cells)
{
	//launched with one block per bin, and each block has blockDim.x threads (set it to 1024 if needed)

	int idx = threadIdx.x;

	cuBReal histogram_value = 0.0;
	bool include_in_reduction = false;
	int bin_idx = 0;

	while (idx < vec.linear_size()) {

		if (vec.is_not_empty(idx)) {

			cuBReal mag_value = cu_GetMagnitude(vec[idx]);

			bin_idx = floor((mag_value - (min - bin / 2)) / bin);

			//here we only want to collect values for the bin corresponding to block number
			if (bin_idx == blockIdx.x) {

				histogram_value += 1.0 / num_nonempty_cells;
				include_in_reduction = true;
			}
		}

		idx += blockDim.x;
	}

	//do a reduction in this block to get total probability value in this bin
	if (blockDim.x == 1024) reduction_sum_blocksize1024(0, 1, &histogram_value, histogram[blockIdx.x], include_in_reduction);
	else reduction_sum(0, 1, &histogram_value, histogram[blockIdx.x], include_in_reduction);
}

//extract histogram after pre-averaging
template <typename VType>
__global__ void get_mag_histogram_kernel(
	VType*& histogram_preaverage, size_t& histogram_preaverage_size, int num_bins, cuBReal bin, cuBReal min, cuBReal max, cuBReal*& histogram, size_t num_nonempty_cells)
{
	//launched with one block per bin, and each block has blockDim.x threads (set it to 1024 if needed)

	int idx = threadIdx.x;

	cuBReal histogram_value = 0.0;
	bool include_in_reduction = false;
	int bin_idx = 0;

	while (idx < histogram_preaverage_size) {

		cuBReal mag_value = cu_GetMagnitude(histogram_preaverage[idx]);

		bin_idx = floor((mag_value - (min - bin / 2)) / bin);

		//here we only want to collect values for the bin corresponding to block number
		if (bin_idx == blockIdx.x) {

			histogram_value += 1.0 / num_nonempty_cells;
			include_in_reduction = true;
		}

		idx += blockDim.x;
	}

	//do a reduction in this block to get total probability value in this bin
	if (blockDim.x == 1024) reduction_sum_blocksize1024(0, 1, &histogram_value, histogram[blockIdx.x], include_in_reduction);
	else reduction_sum(0, 1, &histogram_value, histogram[blockIdx.x], include_in_reduction);
}

//-------------------------------------------- Histogram kernels (angle)

//extract histogram directly from mesh without pre-averaging

template <typename VType>
__global__ void get_ang_histogram_kernel(
	cuVEC<VType>& vec, VType ndir, int num_bins, cuBReal bin, cuBReal min, cuBReal max, cuBReal*& histogram, size_t num_nonempty_cells)
{
	//launched with one block per bin, and each block has blockDim.x threads (set it to 1024 if needed)

	int idx = threadIdx.x;

	cuBReal histogram_value = 0.0;
	bool include_in_reduction = false;
	int bin_idx = 0;

	while (idx < vec.linear_size()) {

		if (vec.is_not_empty(idx)) {

			cuBReal ang_value = acos(vec[idx].normalized() * ndir);

			bin_idx = floor((ang_value - (min - bin / 2)) / bin);

			//here we only want to collect values for the bin corresponding to block number
			if (bin_idx == blockIdx.x) {

				histogram_value += 1.0 / num_nonempty_cells;
				include_in_reduction = true;
			}
		}

		idx += blockDim.x;
	}

	//do a reduction in this block to get total probability value in this bin
	if (blockDim.x == 1024) reduction_sum_blocksize1024(0, 1, &histogram_value, histogram[blockIdx.x], include_in_reduction);
	else reduction_sum(0, 1, &histogram_value, histogram[blockIdx.x], include_in_reduction);
}

//extract histogram after pre-averaging
template <typename VType>
__global__ void get_ang_histogram_kernel(
	VType*& histogram_preaverage, size_t& histogram_preaverage_size, VType ndir, int num_bins, cuBReal bin, cuBReal min, cuBReal max, cuBReal*& histogram, size_t num_nonempty_cells)
{
	//launched with one block per bin, and each block has blockDim.x threads (set it to 1024 if needed)

	int idx = threadIdx.x;

	cuBReal histogram_value = 0.0;
	bool include_in_reduction = false;
	int bin_idx = 0;

	while (idx < histogram_preaverage_size) {

		if (!histogram_preaverage[idx].IsNull()) {

			cuBReal ang_value = acos(histogram_preaverage[idx].normalized() * ndir);

			bin_idx = floor((ang_value - (min - bin / 2)) / bin);

			//here we only want to collect values for the bin corresponding to block number
			if (bin_idx == blockIdx.x) {

				histogram_value += 1.0 / num_nonempty_cells;
				include_in_reduction = true;
			}
		}

		idx += blockDim.x;
	}

	//do a reduction in this block to get total probability value in this bin
	if (blockDim.x == 1024) reduction_sum_blocksize1024(0, 1, &histogram_value, histogram[blockIdx.x], include_in_reduction);
	else reduction_sum(0, 1, &histogram_value, histogram[blockIdx.x], include_in_reduction);
}

//-------------------------------------------- HISTOGRAMS : cuVEC_histo.cuh

// Magnitude

template bool cuVEC<float>::get_mag_histogram(std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, int num_bins, double& min, double& max, size_t num_nonempty_cells, cuINT3 macrocell_dims, bool set_output, int stage_control);
template bool cuVEC<double>::get_mag_histogram(std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, int num_bins, double& min, double& max, size_t num_nonempty_cells, cuINT3 macrocell_dims, bool set_output, int stage_control);

template bool cuVEC<cuFLT3>::get_mag_histogram(std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, int num_bins, double& min, double& max, size_t num_nonempty_cells, cuINT3 macrocell_dims, bool set_output, int stage_control);
template bool cuVEC<cuDBL3>::get_mag_histogram(std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, int num_bins, double& min, double& max, size_t num_nonempty_cells, cuINT3 macrocell_dims, bool set_output, int stage_control);

//compute magnitude histogram data
//extract histogram between magnitudes min and max with given number of bins. if min max not given (set them to zero) then determine them first.
//if num_bins not given then use default value of 100
//if macrocell_dims greater than 1 in any dimension then first average mesh data in macrocells of given size
//without macrocell then pass in num_nonempty_cells (number of nonempty cells); if this is not passed it it's counted first (costs another kernel launch)
//output transferred to cpu in histogram_cpu
template <typename VType>
__host__ bool cuVEC<VType>::get_mag_histogram(
	std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu,
	int num_bins, double& min, double& max,
	size_t num_nonempty_cells, cuINT3 macrocell_dims, 
	bool set_output, int stage_control)
{
	cuSZ3 n_cpu = size_cpu();

	if (num_bins < 2) num_bins = 100;

	//number of cells used to perform pre-averaging
	cuSZ3 num_av_cells = cuSZ3();
	if (macrocell_dims != cuINT3(1)) num_av_cells = cu_round((cuReal3)n_cpu / macrocell_dims);

	//make sure there's enough memory allocated
	if (!allocate_histogram_memory(num_bins, num_av_cells.dim())) return false;
	if (histogram_x_cpu.size() != num_bins) if (!malloc_vector(histogram_x_cpu, num_bins)) return false;
	if (histogram_p_cpu.size() != num_bins) if (!malloc_vector(histogram_p_cpu, num_bins)) return false;

	//bin size : will be adjusted below
	double bin = 0.0;

	if (num_av_cells.IsNull()) {

		//build histogram from all cells separately

		//first determine minimum and maximum values if we need to
		if ((min <= 0.0 && max <= 0.0) || max <= 0.0) {

			cuReal2 minmax = get_minmax(n_cpu.dim(), cuBox());
			min = minmax.i;
			max = minmax.j;
		}

		//size of bin
		bin = (max - min) / (num_bins - 1);

		//zero reduction data
		zero_histogramvalues_kernel <<< (num_bins + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (histogram, histogram_size, aux_real, aux_real2);

		//count number of non-empty cells if we have to
		if (!num_nonempty_cells) {

			count_nonempty_cells(n_cpu.dim());
			num_nonempty_cells = get_gpu_value(aux_integer);
		}

		//now get histogram
		get_mag_histogram_kernel <<< num_bins, 1024 >>>
			(*this, num_bins, (cuBReal)bin, (cuBReal)min, (cuBReal)max, histogram, num_nonempty_cells);
	}
	else {

		//build histogram by first averaging in macrocells

		//if stage_control is 2, then skip initialization (method called with stage_control 1 previously)
		if (stage_control != 2) {

			//zero pre-average data
			zero_histogram_preaverage_values_kernel <<< (num_av_cells.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (histogram_preaverage, histogram_preaverage_size, aux_real, aux_real2);

			//zero reduction data
			zero_histogramvalues_kernel <<< (num_bins + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (histogram, histogram_size, aux_real, aux_real2);

			//do pre-averaging
			if (macrocell_dims.dim() < 1024)
				preaverage_kernel <<< num_av_cells.dim(), CUDATHREADS >>> (*this, histogram_preaverage, num_av_cells, macrocell_dims);
			else
				preaverage_kernel <<< num_av_cells.dim(), 1024 >>> (*this, histogram_preaverage, num_av_cells, macrocell_dims);

			//first determine minimum and maximum values if we need to
			if ((min <= 0.0 && max <= 0.0) || max <= 0.0) {

				preaveraged_minmax <<< (num_av_cells.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(histogram_preaverage, histogram_preaverage_size, aux_real, aux_real2);

				min = get_gpu_value(aux_real);
				max = get_gpu_value(aux_real2);
			}
		}

		//for stage_control 1 skip final calculation
		if (stage_control == 1) return true;

		//size of bin
		bin = (max - min) / (num_bins - 1);

		//in pre-averaging mode total number of cells equals total number of macrocells
		size_t num_cells = num_av_cells.dim();

		//build histogram from averaged data
		get_mag_histogram_kernel <<< num_bins, 1024 >>>
			(histogram_preaverage, histogram_preaverage_size, num_bins, (cuBReal)bin, (cuBReal)min, (cuBReal)max, histogram, num_cells);
	}

	if (set_output) {

		//copy over to histogram_cpu
#pragma omp parallel for
		for (int idx = 0; idx < num_bins; idx++) {

			histogram_x_cpu[idx] = min + idx * bin;
		}

		cudaError_t error = gpu_to_cpu_managed(histogram_p_cpu.data(), histogram, num_bins);
		if (error != cudaSuccess) return false;
	}

	return true;
}

//-------------------------------------------- HISTOGRAMS : cuVEC_histo.cuh

//Angular

template bool cuVEC<cuFLT3>::get_ang_histogram(std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, int num_bins, double& min, double& max, size_t num_nonempty_cells, cuINT3 macrocell_dims, cuFLT3 ndir, bool set_output, int stage_control);
template bool cuVEC<cuDBL3>::get_ang_histogram(std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, int num_bins, double& min, double& max, size_t num_nonempty_cells, cuINT3 macrocell_dims, cuDBL3 ndir, bool set_output, int stage_control);

//compute magnitude histogram data
//extract histogram between magnitudes min and max with given number of bins. if min max not given (set them to zero) then determine them first.
//if num_bins not given then use default value of 100
//if macrocell_dims greater than 1 in any dimension then first average mesh data in macrocells of given size
//without macrocell then pass in num_nonempty_cells (number of nonempty cells); if this is not passed it it's counted first (costs another kernel launch)
//output transferred to cpu in histogram_cpu
template <typename VType>
__host__ bool cuVEC<VType>::get_ang_histogram(
	std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu,
	int num_bins, double& min, double& max,
	size_t num_nonempty_cells, cuINT3 macrocell_dims, VType ndir,
	bool set_output, int stage_control)
{
	cuSZ3 n_cpu = size_cpu();

	if (num_bins < 2) num_bins = 100;
	//ndir should be normalized already, but just in case
	if (!ndir.IsNull()) ndir = ndir.normalized();

	//number of cells used to perform pre-averaging
	cuSZ3 num_av_cells = cuSZ3();
	if (macrocell_dims != cuINT3(1)) num_av_cells = cu_round((cuReal3)n_cpu / macrocell_dims);

	//make sure there's enough memory allocated
	if (!allocate_histogram_memory(num_bins, num_av_cells.dim())) return false;
	if (histogram_x_cpu.size() != num_bins) if (!malloc_vector(histogram_x_cpu, num_bins)) return false;
	if (histogram_p_cpu.size() != num_bins) if (!malloc_vector(histogram_p_cpu, num_bins)) return false;

	//bin size : will be adjusted below
	double bin = 0.0;

	//determine average direction if we have to
	if (ndir.IsNull()) {

		VType average = average_nonempty(n_cpu.dim());

		if (!average.IsNull()) ndir = average.normalized();
		else ndir = VType(1, 0, 0);
	}

	if (num_av_cells.IsNull()) {

		//build histogram from all cells separately

		//if stage_control is 2, then skip initialization (method called with stage_control 1 previously)
		if (stage_control != 2) {

			//zero reduction data
			zero_histogramvalues_kernel <<< (num_bins + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (histogram, histogram_size, aux_real, aux_real2);

			//determine minimum and maximum values if we need to
			if ((min <= 0.0 && max <= 0.0) || max <= 0.0) {

				get_minmax_angles_kernel << < (n_cpu.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (*this, ndir, aux_real, aux_real2);

				min = get_gpu_value(aux_real);
				max = get_gpu_value(aux_real2);
			}
		}

		//for stage_control 1 skip final calculation
		if (stage_control == 1) return true;

		//size of bin
		bin = (max - min) / (num_bins - 1);

		//count number of non-empty cells if we have to
		if (!num_nonempty_cells) {

			count_nonempty_cells(n_cpu.dim());
			num_nonempty_cells = get_gpu_value(aux_integer);
		}

		//now get histogram
		get_ang_histogram_kernel <<< num_bins, 1024 >>>
			(*this, ndir, num_bins, (cuBReal)bin, (cuBReal)min, (cuBReal)max, histogram, num_nonempty_cells);
	}
	else {

		//if stage_control is 2, then skip initialization (method called with stage_control 1 previously)
		if (stage_control != 2) {

			//build histogram by first averaging in macrocells

			//zero pre-average data
			zero_histogram_preaverage_values_kernel <<< (num_av_cells.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (histogram_preaverage, histogram_preaverage_size, aux_real, aux_real2);

			//zero reduction data
			zero_histogramvalues_kernel << < (num_bins + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (histogram, histogram_size, aux_real, aux_real2);

			//do pre-averaging
			if (macrocell_dims.dim() < 1024)
				preaverage_kernel <<< num_av_cells.dim(), CUDATHREADS >>> (*this, histogram_preaverage, num_av_cells, macrocell_dims);
			else
				preaverage_kernel <<< num_av_cells.dim(), 1024 >>> (*this, histogram_preaverage, num_av_cells, macrocell_dims);

			//determine minimum and maximum values if we need to
			if ((min <= 0.0 && max <= 0.0) || max <= 0.0) {

				preaveraged_minmax_angles <<< (num_av_cells.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
					(histogram_preaverage, histogram_preaverage_size, ndir, aux_real, aux_real2);

				min = get_gpu_value(aux_real);
				max = get_gpu_value(aux_real2);
			}
		}

		//for stage_control 1 skip final calculation
		if (stage_control == 1) return true;

		//size of bin
		bin = (max - min) / (num_bins - 1);

		//in pre-averaging mode total number of cells equals total number of macrocells
		size_t num_cells = num_av_cells.dim();

		//build histogram from averaged data
		get_ang_histogram_kernel <<< num_bins, 1024 >>>
			(histogram_preaverage, histogram_preaverage_size, ndir, num_bins, (cuBReal)bin, (cuBReal)min, (cuBReal)max, histogram, num_cells);
	}

	if (set_output) {

		//copy over to histogram_cpu
#pragma omp parallel for
		for (int idx = 0; idx < num_bins; idx++) {

			histogram_x_cpu[idx] = min + idx * bin;
		}

		cudaError_t error = gpu_to_cpu_managed(histogram_p_cpu.data(), histogram, num_bins);
		if (error != cudaSuccess) return false;
	}

	return true;
}