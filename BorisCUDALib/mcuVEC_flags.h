#pragma once

#include "mcuVEC.h"

//--------------------------------------------FLAG CHECKING : mcuVEC_flags.h

template <typename VType, typename MType>
int mcuVEC<VType, MType>::get_nonempty_cells_cpu(void)
{
	int nonempty_cells = 0;

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		nonempty_cells += mng(mGPU)->get_nonempty_cells_cpu();
	}

	return nonempty_cells;
}

//--------------------------------------------SET CELL FLAGS - EXTERNAL USE : mcuVEC_flags.h

//set dirichlet boundary conditions from surface_rect (must be a rectangle intersecting with one of the surfaces of this mesh) and value
//return false on memory allocation failure only, otherwise return true even if surface_rect was not valid
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
bool mcuVEC<VType, MType>::set_dirichlet_conditions(cuRect surface_rect, VType value)
{
	bool success = true;

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		success &= mng(mGPU)->set_dirichlet_conditions(surface_rect, value);
	}

	return success;
}

//clear all dirichlet flags and vectors
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::clear_dirichlet_flags(void)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->clear_dirichlet_flags();
	}
}

//set pbc conditions : setting any to false clears flags. Applicable to cuVEC_VC only
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::set_pbc(int pbc_x_, int pbc_y_, int pbc_z_)
{
	pbc_x = pbc_x_;
	pbc_y = pbc_y_;
	pbc_z = pbc_z_;

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->set_pbc(pbc_x, pbc_y, pbc_z);
	}

	//if pbc is along halo direction, need to adjust halo conditions
	if ((pbc_x && halo_flag == NF2_HALOX) || (pbc_y && halo_flag == NF2_HALOY) || (pbc_z && halo_flag == NF2_HALOZ)) set_halo_conditions();
}

template <typename VType, typename MType>
template <typename cpuVEC_VC, typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::set_pbc_from(cpuVEC_VC& vec_vc)
{
	set_pbc(vec_vc.is_pbc_x(), vec_vc.is_pbc_y(), vec_vc.is_pbc_z());
}

template <typename VType, typename MType>
template <typename cpuVEC_VC, typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::set_pbc_to(cpuVEC_VC& vec_vc)
{
	vec_vc.set_pbc(pbc_x, pbc_y, pbc_z);
}

//clear all pbc flags : can also be achieved setting all flags to false in set_pbc but this one is more readable. Applicable to cuVEC_VC only
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::clear_pbc(void)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->clear_pbc();
	}
}

//clear all composite media boundary flags. Applicable to cuVEC_VC only
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::clear_cmbnd_flags(void)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->clear_cmbnd_flags();
	}
}

//mark cells included in this rectangle (absolute coordinates) to be skipped during some computations (if status true, else clear the skip cells flags in this rectangle). Applicable to cuVEC_VC only
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::set_skipcells(cuRect rectangle, bool status)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->set_skipcells(rectangle, status);
	}
}

//clear all skip cell flags. Applicable to cuVEC_VC only
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::clear_skipcells(void)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->clear_skipcells();
	}
}

//Robin conditions. Applicable to cuVEC_VC only
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::set_robin_conditions(cuReal2 robin_v_, cuReal2 robin_px_, cuReal2 robin_nx_, cuReal2 robin_py_, cuReal2 robin_ny_, cuReal2 robin_pz_, cuReal2 robin_nz_)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->set_robin_conditions(robin_v_, robin_px_, robin_nx_, robin_py_, robin_ny_, robin_pz_, robin_nz_);
	}
}

//clear all Robin boundary conditions and values. Applicable to cuVEC_VC only
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::clear_robin_conditions(void)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->clear_robin_conditions();
	}
}