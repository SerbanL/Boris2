#pragma once

#include "VEC.h"

template bool VEC<float>::get_mag_histogram(std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, int num_bins, double& min, double& max, size_t num_nonempty_cells, INT3 macrocell_dims);
template bool VEC<double>::get_mag_histogram(std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, int num_bins, double& min, double& max, size_t num_nonempty_cells, INT3 macrocell_dims);

template bool VEC<FLT3>::get_mag_histogram(std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, int num_bins, double& min, double& max, size_t num_nonempty_cells, INT3 macrocell_dims);
template bool VEC<DBL3>::get_mag_histogram(std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, int num_bins, double& min, double& max, size_t num_nonempty_cells, INT3 macrocell_dims);

//compute magnitude histogram data
//extract histogram between magnitudes min and max with given number of bins. if min max not given (set them to zero) then determine them first.
//if num_bins not given then use default value of 100
//if macrocell_dims greater than 1 in any dimension then first average mesh data in macrocells of given size
//without macrocell then pass in num_nonempty_cells (number of nonempty cells); if this is not passed it it's counted first (costs another kernel launch)
//output transferred to cpu in histogram_cpu
template <typename VType>
bool VEC<VType>::get_mag_histogram(std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, int num_bins, double& min, double& max, size_t num_nonempty_cells, INT3 macrocell_dims)
{
	int OmpThreads = omp_get_num_procs();

	if (num_bins < 2) num_bins = 100;

	//number of cells used to perform pre-averaging
	SZ3 num_av_cells = SZ3();
	if (macrocell_dims != INT3(1)) num_av_cells = round((DBL3)n / macrocell_dims);

	//make sure there's enough memory allocated
	if (histogram.size() != num_bins) {

		if (!malloc_vector(histogram, num_bins)) return false;

		for (int bin_idx = 0; bin_idx < num_bins; bin_idx++) {

			if (!malloc_vector(histogram[bin_idx], OmpThreads)) return false;
		}
	}

	//bin size : will be adjusted below
	double bin = 0.0;

	if (num_av_cells.IsNull()) {

		//build histogram from all cells separately

		//first determine minimum and maximum values if we need to
		if ((min <= 0.0 && max <= 0.0) || max <= 0.0) {

			DBL2 minmax = get_minmax();
			min = minmax.i;
			max = minmax.j;
		}

		//size of bin
		bin = (max - min) / (num_bins - 1);

		//zero reduction data
#pragma omp parallel for
		for (int bin_idx = 0; bin_idx < num_bins; bin_idx++) {

			for (int tn = 0; tn < OmpThreads; tn++) histogram[bin_idx][tn] = 0.0;
		}

		//count number of non-empty cells if we have to
		if (!num_nonempty_cells) {

			num_nonempty_cells = get_nonempty_cells();
		}

		//now get histogram
#pragma omp parallel for
		for (int idx = 0; idx < n.dim(); idx++) {

			int tn = omp_get_thread_num();

			if (is_not_empty(idx)) {

				double value = GetMagnitude(quantity[idx]);

				int bin_idx = floor((value - (min - bin / 2)) / bin);

				if (bin_idx >= 0 && bin_idx < num_bins) {

					histogram[bin_idx][tn] += 1.0 / num_nonempty_cells;
				}
			}
		}
	}
	else {

		//build histogram by first averaging in macrocells

		//first determine minimum and maximum values if we need to
		if ((min <= 0.0 && max <= 0.0) || max <= 0.0) {

			magnitude_reduction.new_minmax_reduction();

			for (int k_seg = 0; k_seg < num_av_cells.k; k_seg++) {
#pragma omp parallel for
				for (int j_seg = 0; j_seg < num_av_cells.j; j_seg++) {
					for (int i_seg = 0; i_seg < num_av_cells.i; i_seg++) {

						VType average = VType();

						for (int k = 0; k < macrocell_dims.k; k++) {
							for (int j = 0; j < macrocell_dims.j; j++) {
								for (int i = 0; i < macrocell_dims.i; i++) {

									INT3 ijk = INT3(i_seg * macrocell_dims.i + i, j_seg * macrocell_dims.j + j, k_seg * macrocell_dims.k + k);

									if (ijk < n && is_not_empty(ijk)) average += (*this)[ijk] / macrocell_dims.dim();
								}
							}
						}

						magnitude_reduction.reduce_minmax(GetMagnitude(average));
					}
				}
			}

			DBL2 minmax = magnitude_reduction.minmax();
			min = minmax.i;
			max = minmax.j;
		}

		//size of bin
		bin = (max - min) / (num_bins - 1);

		//zero reduction data
#pragma omp parallel for
		for (int bin_idx = 0; bin_idx < num_bins; bin_idx++) {

			for (int tn = 0; tn < OmpThreads; tn++) histogram[bin_idx][tn] = 0.0;
		}

		//in pre-averaging mode total number of cells equals total number of macrocells
		size_t num_cells = num_av_cells.dim();

		//now get histogram
		for (int k_seg = 0; k_seg < num_av_cells.k; k_seg++) {
#pragma omp parallel for
			for (int j_seg = 0; j_seg < num_av_cells.j; j_seg++) {
				for (int i_seg = 0; i_seg < num_av_cells.i; i_seg++) {

					VType average = VType();

					for (int k = 0; k < macrocell_dims.k; k++) {
						for (int j = 0; j < macrocell_dims.j; j++) {
							for (int i = 0; i < macrocell_dims.i; i++) {

								INT3 ijk = INT3(i_seg * macrocell_dims.i + i, j_seg * macrocell_dims.j + j, k_seg * macrocell_dims.k + k);

								if (ijk < n && is_not_empty(ijk)) average += (*this)[ijk] / macrocell_dims.dim();
							}
						}
					}

					int tn = omp_get_thread_num();

					double value = GetMagnitude(average);

					int bin_idx = floor((value - (min - bin / 2)) / bin);

					if (bin_idx >= 0 && bin_idx < num_bins) {

						histogram[bin_idx][tn] += 1.0 / num_cells;
					}
				}
			}
		}
	}

	//copy over to histogram_cpu
	if (histogram_x_cpu.size() != num_bins) if (!malloc_vector(histogram_x_cpu, num_bins)) return false;
	if (histogram_p_cpu.size() != num_bins) if (!malloc_vector(histogram_p_cpu, num_bins)) return false;

#pragma omp parallel for
	for (int bin_idx = 0; bin_idx < num_bins; bin_idx++) {

		histogram_x_cpu[bin_idx] = min + bin_idx * bin;

		histogram_p_cpu[bin_idx] = 0.0;
		for (int tn = 0; tn < OmpThreads; tn++) {

			histogram_p_cpu[bin_idx] += histogram[bin_idx][tn];
		}
	}

	return true;
}

//-------------------------------------------

template bool VEC<FLT3>::get_ang_histogram(std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, int num_bins, double& min, double& max, size_t num_nonempty_cells, INT3 macrocell_dims, FLT3 ndir);
template bool VEC<DBL3>::get_ang_histogram(std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, int num_bins, double& min, double& max, size_t num_nonempty_cells, INT3 macrocell_dims, DBL3 ndir);

//compute magnitude histogram data
//extract histogram between magnitudes min and max with given number of bins. if min max not given (set them to zero) then determine them first.
//if num_bins not given then use default value of 100
//if macrocell_dims greater than 1 in any dimension then first average mesh data in macrocells of given size
//without macrocell then pass in num_nonempty_cells (number of nonempty cells); if this is not passed it it's counted first (costs another kernel launch)
//output transferred to cpu in histogram_cpu
template <typename VType>
bool VEC<VType>::get_ang_histogram(std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, int num_bins, double& min, double& max, size_t num_nonempty_cells, INT3 macrocell_dims, VType ndir)
{
	int OmpThreads = omp_get_num_procs();

	if (num_bins < 2) num_bins = 100;

	//number of cells used to perform pre-averaging
	SZ3 num_av_cells = SZ3();
	if (macrocell_dims != INT3(1)) num_av_cells = round((DBL3)n / macrocell_dims);

	//make sure there's enough memory allocated
	if (histogram.size() != num_bins) {

		if (!malloc_vector(histogram, num_bins)) return false;

		for (int bin_idx = 0; bin_idx < num_bins; bin_idx++) {

			if (!malloc_vector(histogram[bin_idx], OmpThreads)) return false;
		}
	}

	//bin size : will be adjusted below
	double bin = 0.0;

	//first determine average direction if we have to
	if (ndir.IsNull()) {

		VType average = average_nonempty_omp();
		if (!average.IsNull()) ndir = average.normalized();
		else ndir = DBL3(1, 0, 0);
	}

	if (num_av_cells.IsNull()) {

		//build histogram from all cells separately

		//determine minimum and maximum values if we need to
		if ((min <= 0.0 && max <= 0.0) || max <= 0.0) {

			magnitude_reduction.new_minmax_reduction();

#pragma omp parallel for
			for (int idx = 0; idx < n.dim(); idx++) {

				if (is_not_empty(idx)) {

					magnitude_reduction.reduce_minmax(acos(quantity[idx].normalized() * ndir));
				}
			}

			DBL2 minmax = magnitude_reduction.minmax();
			min = minmax.i;
			max = minmax.j;
		}

		//size of bin
		bin = (max - min) / (num_bins - 1);

		//zero reduction data
#pragma omp parallel for
		for (int bin_idx = 0; bin_idx < num_bins; bin_idx++) {

			for (int tn = 0; tn < OmpThreads; tn++) histogram[bin_idx][tn] = 0.0;
		}

		//count number of non-empty cells if we have to
		if (!num_nonempty_cells) {

			num_nonempty_cells = get_nonempty_cells();
		}

		//now get histogram
#pragma omp parallel for
		for (int idx = 0; idx < n.dim(); idx++) {

			int tn = omp_get_thread_num();

			if (is_not_empty(idx)) {

				double value = acos(quantity[idx].normalized() * ndir);

				int bin_idx = floor((value - (min - bin / 2)) / bin);

				if (bin_idx >= 0 && bin_idx < num_bins) {

					histogram[bin_idx][tn] += 1.0 / num_nonempty_cells;
				}
			}
		}
	}
	else {

		//build histogram by first averaging in macrocells

		//determine minimum and maximum values if we need to
		if ((min <= 0.0 && max <= 0.0) || max <= 0.0) {

			magnitude_reduction.new_minmax_reduction();

			for (int k_seg = 0; k_seg < num_av_cells.k; k_seg++) {
#pragma omp parallel for
				for (int j_seg = 0; j_seg < num_av_cells.j; j_seg++) {
					for (int i_seg = 0; i_seg < num_av_cells.i; i_seg++) {

						VType average = VType();

						for (int k = 0; k < macrocell_dims.k; k++) {
							for (int j = 0; j < macrocell_dims.j; j++) {
								for (int i = 0; i < macrocell_dims.i; i++) {

									INT3 ijk = INT3(i_seg * macrocell_dims.i + i, j_seg * macrocell_dims.j + j, k_seg * macrocell_dims.k + k);

									if (ijk < n && is_not_empty(ijk)) average += (*this)[ijk] / macrocell_dims.dim();
								}
							}
						}
						
						if (!average.IsNull()) magnitude_reduction.reduce_minmax(acos(average.normalized() * ndir));
					}
				}
			}

			DBL2 minmax = magnitude_reduction.minmax();
			min = minmax.i;
			max = minmax.j;
		}

		//size of bin
		bin = (max - min) / (num_bins - 1);

		//zero reduction data
#pragma omp parallel for
		for (int bin_idx = 0; bin_idx < num_bins; bin_idx++) {

			for (int tn = 0; tn < OmpThreads; tn++) histogram[bin_idx][tn] = 0.0;
		}

		//in pre-averaging mode total number of cells equals total number of macrocells
		size_t num_cells = num_av_cells.dim();

		//now get histogram
		for (int k_seg = 0; k_seg < num_av_cells.k; k_seg++) {
#pragma omp parallel for
			for (int j_seg = 0; j_seg < num_av_cells.j; j_seg++) {
				for (int i_seg = 0; i_seg < num_av_cells.i; i_seg++) {

					VType average = VType();

					for (int k = 0; k < macrocell_dims.k; k++) {
						for (int j = 0; j < macrocell_dims.j; j++) {
							for (int i = 0; i < macrocell_dims.i; i++) {

								INT3 ijk = INT3(i_seg * macrocell_dims.i + i, j_seg * macrocell_dims.j + j, k_seg * macrocell_dims.k + k);

								if (ijk < n && is_not_empty(ijk)) average += (*this)[ijk] / macrocell_dims.dim();
							}
						}
					}

					if (!average.IsNull()) {

						int tn = omp_get_thread_num();

						double value = acos(average.normalized() * ndir);

						int bin_idx = floor((value - (min - bin / 2)) / bin);

						if (bin_idx >= 0 && bin_idx < num_bins) {

							histogram[bin_idx][tn] += 1.0 / num_cells;
						}
					}
				}
			}
		}
	}

	//copy over to histogram_cpu
	if (histogram_x_cpu.size() != num_bins) if (!malloc_vector(histogram_x_cpu, num_bins)) return false;
	if (histogram_p_cpu.size() != num_bins) if (!malloc_vector(histogram_p_cpu, num_bins)) return false;

#pragma omp parallel for
	for (int bin_idx = 0; bin_idx < num_bins; bin_idx++) {

		histogram_x_cpu[bin_idx] = min + bin_idx * bin;

		histogram_p_cpu[bin_idx] = 0.0;
		for (int tn = 0; tn < OmpThreads; tn++) {

			histogram_p_cpu[bin_idx] += histogram[bin_idx][tn];
		}
	}

	return true;
}