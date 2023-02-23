#pragma once

#include "mcuVEC.h"

//--------------------------------------------HISTOGRAMS : mcuVEC_histo.h

template <typename VType, typename MType>
void mcuVEC<VType, MType>::setup_histogram_transfers(void)
{
	//setup histogram transfers for each device : main handles to histogram on each device
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		histogram_transf.setup_device_memory_handle_managed(mGPU, mng(mGPU)->histogram_ref());
	}
}

template <typename VType, typename MType>
bool mcuVEC<VType, MType>::get_mag_histogram(
	std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu,
	int num_bins, double& min, double& max, size_t num_nonempty_cells, cuINT3 macrocell_dims)
{
	int stage_control = 0;

	//1a. if min max not given then determine them from the entire VEC
	if ((min <= 0.0 && max <= 0.0) || max <= 0.0) {

		if (macrocell_dims == cuINT3(1)) {

			//no macrocell : just find min max from entire cuVEC
			cuReal2 minmax = get_minmax(cuBox());
			min = minmax.i;
			max = minmax.j;
		}
		else {

			//in macrocell mode, need to do pre-averages in each subVEC and get min max values for each
			stage_control = 1;
			
			min = 1e38;
			max = -1e38;

			for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

				double lmin = 0, lmax = 0;
				if (!mng(mGPU)->get_mag_histogram(histogram_x_cpu, histogram_p_cpu, num_bins, lmin, lmax, num_nonempty_cells, macrocell_dims, false, stage_control)) return false;
				if (lmin < min) min = lmin;
				if (lmax > max) max = lmax;
			}

			//next time we call histogram methods, pre-averages will be skipped as already available
			stage_control = 2;
		}
	}

	//1b. if num_nonempty_cells not given find it
	if (!num_nonempty_cells) num_nonempty_cells = get_nonempty_cells_cpu();

	//2. make sure enough memory is allocated for histogram_x_cpu, histogram_p_cpu
	if (num_bins < 2) num_bins = 100;
	if (histogram_x_cpu.size() != num_bins) if (!malloc_vector(histogram_x_cpu, num_bins)) return false;
	if (histogram_p_cpu.size() != num_bins) if (!malloc_vector(histogram_p_cpu, num_bins)) return false;

	//3. Generate data in histogram_x_cpu
	//size of bin
	double bin = (max - min) / (num_bins - 1);

#pragma omp parallel for
	for (int idx = 0; idx < num_bins; idx++) {

		histogram_x_cpu[idx] = min + idx * bin;
	}

	//4. Compute sub-histograms on each device, but do not set output yet
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		if (!mng(mGPU)->get_mag_histogram(histogram_x_cpu, histogram_p_cpu, num_bins, min, max, num_nonempty_cells, macrocell_dims, false, stage_control)) return false;
	}

	//5. Add together all sub-histograms to obtain an overall histogram, stil in gpu memory

	//SETUP TRANSFERS

	//allocate histogram arrays on base gpu if needed
	mGPU.select_device(0);
	size_t histogram_size = mng(0)->get_histogram_size_cpu();
	bool histogram_size_changed = false;
	for (int idx = 0; idx < histogram_base_aux.size(); idx++) {

		if (histogram_size != histogram_base_aux[idx].size()) {

			//allocate...
			histogram_base_aux[idx].resize(histogram_size);

			//...and refresh extra histogram handles on base device
			histogram_transf.setup_extra_device_memory_handle(0, idx + 1, histogram_base_aux[idx].get_array());

			//set flag true : we'll need to refresh histogram_base_aux_col handles too
			histogram_size_changed = true;
		}
	}
	
	//refresh histogram_base_aux_col if needed
	if (histogram_size_changed) {

		histogram_base_aux_col.clear();

		for (int idx = 0; idx < histogram_base_aux.size(); idx++) {

			histogram_base_aux_col.push_back(histogram_base_aux[idx].get_managed_array());
		}

		//if histogram size has changed, then also need to refresh base handles as these will also have been re-allocated
		setup_histogram_transfers();
	}

	//setup transfer size
	histogram_transf.set_transfer_size(num_bins);

	//DO TRANSFERS TO BASE DEVICE (from histogram on each device to auxiliary spaces on base device)
	for (int idx = 0; idx < histogram_base_aux.size(); idx++) {

		histogram_transf.transfer(0, idx + 1, idx, 0);
	}
	
	//ADD HISTOGRAMS ON BASE DEVICE
	add_histograms();

	//6. Finally set output in histogram_p_cpu
	histogram_base_aux[0].copy_to_vector(histogram_p_cpu);

	return true;
}

template <typename VType, typename MType>
bool mcuVEC<VType, MType>::get_ang_histogram(
	std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu,
	int num_bins, double& min, double& max, size_t num_nonempty_cells, cuINT3 macrocell_dims, VType ndir)
{
	int stage_control = 0;

	//1a. If ndir not given, find it
	if (ndir.IsNull()) {

		VType average = average_nonempty(cuBox());

		if (!average.IsNull()) ndir = average.normalized();
		else ndir = VType(1, 0, 0);
	}
	//must make sure ndir is normalized
	else ndir = ndir.normalized();

	//1b. if min max (angle deviation from ndir) not given then determine them from the entire VEC
	if ((min <= 0.0 && max <= 0.0) || max <= 0.0) {

		//in macrocell mode, need to do pre-averages in each subVEC and get min max values for each
		stage_control = 1;

		min = 1e38;
		max = -1e38;

		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			double lmin = 0, lmax = 0;
			if (!mng(mGPU)->get_ang_histogram(histogram_x_cpu, histogram_p_cpu, num_bins, lmin, lmax, num_nonempty_cells, macrocell_dims, ndir, false, stage_control)) return false;
			if (lmin < min) min = lmin;
			if (lmax > max) max = lmax;
		}

		//next time we call histogram methods, pre-averages will be skipped as already available
		stage_control = 2;
	}

	//1c. if num_nonempty_cells not given find it
	if (!num_nonempty_cells) num_nonempty_cells = get_nonempty_cells_cpu();

	//2. make sure enough memory is allocated for histogram_x_cpu, histogram_p_cpu
	if (num_bins < 2) num_bins = 100;
	if (histogram_x_cpu.size() != num_bins) if (!malloc_vector(histogram_x_cpu, num_bins)) return false;
	if (histogram_p_cpu.size() != num_bins) if (!malloc_vector(histogram_p_cpu, num_bins)) return false;

	//3. Generate data in histogram_x_cpu
	//size of bin
	double bin = (max - min) / (num_bins - 1);

#pragma omp parallel for
	for (int idx = 0; idx < num_bins; idx++) {

		histogram_x_cpu[idx] = min + idx * bin;
	}

	//4. Compute sub-histograms on each device, but do not set output yet
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		if (!mng(mGPU)->get_ang_histogram(histogram_x_cpu, histogram_p_cpu, num_bins, min, max, num_nonempty_cells, macrocell_dims, ndir, false, stage_control)) return false;
	}

	//5. Add together all sub-histograms to obtain an overall histogram, stil in gpu memory

	//SETUP TRANSFERS

	//allocate histogram arrays on base gpu if needed
	mGPU.select_device(0);
	size_t histogram_size = mng(0)->get_histogram_size_cpu();
	bool histogram_size_changed = false;
	for (int idx = 0; idx < histogram_base_aux.size(); idx++) {

		if (histogram_size != histogram_base_aux[idx].size()) {

			//allocate...
			histogram_base_aux[idx].resize(histogram_size);

			//...and refresh extra histogram handles on base device
			histogram_transf.setup_extra_device_memory_handle(0, idx, histogram_base_aux[idx].get_array());

			//set flag true : we'll need to refresh histogram_base_aux_col handles too
			histogram_size_changed = true;
		}
	}

	//refresh histogram_base_aux_col if needed
	if (histogram_size_changed) {

		histogram_base_aux_col.clear();

		for (int idx = 0; idx < histogram_base_aux.size(); idx++) {

			histogram_base_aux_col.push_back(histogram_base_aux[idx].get_managed_array());
		}

		//if histogram size has changed, then also need to refresh base handles as these will also have been re-allocated
		setup_histogram_transfers();
	}

	//setup transfer size
	histogram_transf.set_transfer_size(num_bins);

	//DO TRANSFERS TO BASE DEVICE (from histogram on each device to auxiliary spaces on base device)
	for (int idx = 0; idx < histogram_base_aux.size(); idx++) {

		histogram_transf.transfer(0, idx + 1, idx, 0);
	}

	//ADD HISTOGRAMS ON BASE DEVICE
	add_histograms();

	//6. Finally set output in histogram_p_cpu
	histogram_base_aux[0].copy_to_vector(histogram_p_cpu);

	return true;
}
