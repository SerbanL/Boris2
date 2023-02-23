#pragma once

#include "mcuVEC.h"

//------------------------------------------------------------------- ALLOCATION

//allocate memory for halos
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::allocate_halos(void)
{
	//only use if halos applicable
	if (halo_depth == 0 || mGPU.get_num_devices() == 1) return;

	int num_halo_cells = 0;

	//now determine halo direction and number of cells we need to allocate for halos
	//x direction
	if (n.x > n.y && n.x > n.z) {

		halo_flag_n = NF2_HALONX; halo_flag_p = NF2_HALOPX; halo_flag = NF2_HALOX;
		num_halo_cells = halo_depth * n.y * n.z;
	}
	//y direction
	else if (n.y > n.z) {

		halo_flag_n = NF2_HALONY; halo_flag_p = NF2_HALOPY; halo_flag = NF2_HALOY;
		num_halo_cells = halo_depth * n.x * n.z;
	}
	//z direction
	else {

		halo_flag_n = NF2_HALONZ; halo_flag_p = NF2_HALOPZ; halo_flag = NF2_HALOZ;
		num_halo_cells = halo_depth * n.x * n.y;
	}

	bool reallocated = false;

	//allocate memory (if needed)
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		if (phalo_ngbr_p[mGPU] == nullptr || phalo_ngbr_p[mGPU]->size() != num_halo_cells) {

			if (phalo_ngbr_p[mGPU]) {

				delete phalo_ngbr_p[mGPU];
				phalo_ngbr_p[mGPU] = nullptr;
			}

			phalo_ngbr_p[mGPU] = new cu_arr<int>(num_halo_cells);
			reallocated = true;
		}

		if (phalo_ngbr_n[mGPU] == nullptr || phalo_ngbr_n[mGPU]->size() != num_halo_cells) {

			if (phalo_ngbr_n[mGPU]) {

				delete phalo_ngbr_n[mGPU];
				phalo_ngbr_n[mGPU] = nullptr;
			}

			phalo_ngbr_n[mGPU] = new cu_arr<int>(num_halo_cells);
			reallocated = true;
		}
	}

	if (reallocated) {

		//setup memory transfer objects
		halo_ngbr_n_transf.set_transfer_size(num_halo_cells);
		halo_ngbr_p_transf.set_transfer_size(num_halo_cells);
		halo_quant_n_transf.set_transfer_size(num_halo_cells);
		halo_quant_p_transf.set_transfer_size(num_halo_cells);

		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			halo_ngbr_n_transf.setup_device_memory_handle(mGPU, phalo_ngbr_n[mGPU]->get_array());
			halo_ngbr_p_transf.setup_device_memory_handle(mGPU, phalo_ngbr_p[mGPU]->get_array());
		}
	}
}

//------------------------------------------------------------------- SETUP

//coordinate ngbr flag exchanges to set halo conditions in managed devices
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::set_halo_conditions(void)
{
	//only use if halos applicable
	if (halo_depth == 0 || mGPU.get_num_devices() == 1) return;

	//for each device copy ngbr flags region from adjacent device, then set halo flags
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//cells box on n and p sides; halo flags on n and p sides.
		cuBox box_n, box_p;

		switch (halo_flag) {

		case NF2_HALOX:
			if (mGPU + 1 < mGPU.get_num_devices()) box_p = cuBox(1, pn_d[mGPU + 1].y, pn_d[mGPU + 1].z);
			if (mGPU > 0) box_n = cuBox(pn_d[mGPU - 1].x - 1, 0, 0, pn_d[mGPU - 1].x, pn_d[mGPU - 1].y, pn_d[mGPU - 1].z);
			break;
		case NF2_HALOY:
			if (mGPU + 1 < mGPU.get_num_devices()) box_p = cuBox(pn_d[mGPU + 1].x, 1, pn_d[mGPU + 1].z);
			if (mGPU > 0) box_n = cuBox(0, pn_d[mGPU - 1].y - 1, 0, pn_d[mGPU - 1].x, pn_d[mGPU - 1].y, pn_d[mGPU - 1].z);
			break;
		case NF2_HALOZ:
			if (mGPU + 1 < mGPU.get_num_devices()) box_p = cuBox(pn_d[mGPU + 1].x, pn_d[mGPU + 1].y, 1);
			if (mGPU > 0) box_n = cuBox(0, 0, pn_d[mGPU - 1].z - 1, pn_d[mGPU - 1].x, pn_d[mGPU - 1].y, pn_d[mGPU - 1].z);
			break;
		}

		//get ngbr flags from device on p side (if possible)
		//to do this: 1) select next device, 2) copy flags region from that device in linear storage space here, 3) transfer to current device in one operation
		if (mGPU.select_next_device() && mGPU + 1 < mGPU.get_num_devices()) {

			//2)
			mng(mGPU + 1)->extract_ngbrFlags(box_p, *phalo_ngbr_p[mGPU + 1]);

			//3) transfer phalo_ngbr_p[mGPU + 1] to phalo_ngbr_p[mGPU]
			halo_ngbr_p_transf.transfer(mGPU, mGPU + 1);
		}

		//get ngbr flags from device on n side (if possible)
		if (mGPU.select_previous_device() && mGPU > 0) {

			//2)
			mng(mGPU - 1)->extract_ngbrFlags(box_n, *phalo_ngbr_n[mGPU - 1]);

			//3) transfer phalo_ngbr_n[mGPU - 1] to phalo_ngbr_n[mGPU]
			halo_ngbr_n_transf.transfer(mGPU, mGPU - 1);
		}

		//set flags in current device
		mGPU.select_current_device();
		if (mGPU < mGPU.get_num_devices() - 1) mng(mGPU)->set_halo_conditions(*phalo_ngbr_p[mGPU], halo_flag_p, halo_depth);
		if (mGPU > 0) mng(mGPU)->set_halo_conditions(*phalo_ngbr_n[mGPU], halo_flag_n, halo_depth);

		//now that required halos have been allocated, configure transfer object
		switch (halo_flag) {

		case NF2_HALOX:
			halo_quant_n_transf.setup_device_memory_handle_managed(mGPU, mng(mGPU)->halo_nx_ref());
			halo_quant_p_transf.setup_device_memory_handle_managed(mGPU, mng(mGPU)->halo_px_ref());
			break;

		case NF2_HALOY:
			halo_quant_n_transf.setup_device_memory_handle_managed(mGPU, mng(mGPU)->halo_ny_ref());
			halo_quant_p_transf.setup_device_memory_handle_managed(mGPU, mng(mGPU)->halo_py_ref());
			break;

		case NF2_HALOZ:
			halo_quant_n_transf.setup_device_memory_handle_managed(mGPU, mng(mGPU)->halo_nz_ref());
			halo_quant_p_transf.setup_device_memory_handle_managed(mGPU, mng(mGPU)->halo_pz_ref());
			break;
		}
	}
}

//------------------------------------------------------------------- RUNTIME : EXCHANGE

//exchange values in all halos
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::exchange_halos(void)
{
	//only use if halos applicable
	if (halo_depth == 0 || mGPU.get_num_devices() == 1) return;

	//1. in each device extract values to halos

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		if (mGPU > 0 && mGPU < mGPU.get_num_devices() - 1) {

			mng(mGPU)->extract_halo(halo_flag_n, halo_quant_n_transf.get_transfer_size());
			mng(mGPU)->extract_halo(halo_flag_p, halo_quant_p_transf.get_transfer_size());
		}
		else if (mGPU > 0) mng(mGPU)->extract_halo(halo_flag_p, halo_quant_p_transf.get_transfer_size());
		else mng(mGPU)->extract_halo(halo_flag_n, halo_quant_n_transf.get_transfer_size());
	}

	//2. transfer n and p halos between devices:
	//for p halos transfer in increasing index order from idx to idx - 1
	//for n halos transfer in decreasing index order from idx - 1 to idx

	for (int idx = 1; idx < mGPU.get_num_devices(); idx++) {

		mGPU.select_device(idx - 1);
		halo_quant_p_transf.transfer(idx - 1, idx);
	}

	for (int idx = mGPU.get_num_devices() - 1; idx > 0; idx--) {

		mGPU.select_device(idx);
		halo_quant_n_transf.transfer(idx, idx - 1);
	}
}