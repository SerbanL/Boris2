#pragma once

#include "cuVEC_VC.h"

//-------------------------------- CURL OPERATOR

//curl operator. Use Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
//can only be applied if VType is a VAL3
template <typename VType>
__device__ VType cuVEC_VC<VType>::curl_neu(int idx) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
		curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOX)) {

		//halo index
		int hidx = ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y;

		if (ngbrFlags2[idx] & NF2_HALOPX) {

			if (ngbrFlags[idx] & NF_NNX) {

				curl.y -= (halo_px[hidx].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (halo_px[hidx].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (halo_px[hidx].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (halo_px[hidx].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.x);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPX) {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - halo_nx[hidx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - halo_nx[hidx].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (cuVEC<VType>::quantity[idx].z - halo_nx[hidx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx].y - halo_nx[hidx].y) / (2 * cuVEC<VType>::h.x);
			}
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.x;
			curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / cuVEC<VType>::h.x;
			curl.z += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / cuVEC<VType>::h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.x);
			}
		}

		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
		}
	}

	//y direction differentials
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
		curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOY)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPY) {

			if (ngbrFlags[idx] & NF_NNY) {

				curl.x += (halo_py[hidx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (halo_py[hidx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (halo_py[hidx].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (halo_py[hidx].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.y);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPY) {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - halo_ny[hidx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - halo_ny[hidx].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (cuVEC<VType>::quantity[idx].z - halo_ny[hidx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx].x - halo_ny[hidx].x) / (2 * cuVEC<VType>::h.y);
			}
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.y;
			curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / cuVEC<VType>::h.y;
			curl.z -= (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cuVEC<VType>::h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.y);
			}
		}

		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
		}
	}

	//z direction differentials
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
		curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOZ)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPZ) {

			if (ngbrFlags[idx] & NF_NNZ) {

				curl.x -= (halo_pz[hidx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (halo_pz[hidx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (halo_pz[hidx].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (halo_pz[hidx].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.z);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPZ) {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - halo_nz[hidx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - halo_nz[hidx].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (cuVEC<VType>::quantity[idx].y - halo_nz[hidx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx].x - halo_nz[hidx].x) / (2 * cuVEC<VType>::h.z);
			}
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.z;
			curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cuVEC<VType>::h.z;
			curl.y += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cuVEC<VType>::h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.z);
			}
		}

		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
		}
	}

	return curl;
}


//curl operator. Use non-homogeneous Neumann boundary conditions.
//Can be used at composite media boundaries where sided differentials will be used instead.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
//can only be applied if VType is a VAL3
template <typename VType>
template <typename Class_BDiff>
__device__ VType cuVEC_VC<VType>::curl_nneu(int idx, const Class_BDiff& bdiff_class) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
		curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOX)) {

		//halo index
		int hidx = ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y;

		if (ngbrFlags2[idx] & NF2_HALOPX) {

			if (ngbrFlags[idx] & NF_NNX) {

				curl.y -= (halo_px[hidx].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (halo_px[hidx].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (halo_px[hidx].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (halo_px[hidx].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.x);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPX) {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - halo_nx[hidx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - halo_nx[hidx].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (cuVEC<VType>::quantity[idx].z - halo_nx[hidx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx].y - halo_nx[hidx].y) / (2 * cuVEC<VType>::h.x);
			}
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.x;
			curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / cuVEC<VType>::h.x;
			curl.z += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / cuVEC<VType>::h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= ((cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.x + bdiff_val.x.z) / 2;
				curl.z += ((cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.x + bdiff_val.x.y) / 2;
			}
		}

		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= ((cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / cuVEC<VType>::h.x + bdiff_val.x.z) / 2;
				curl.z += ((cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / cuVEC<VType>::h.x + bdiff_val.x.y) / 2;
			}
		}
	}

	//y direction differentials
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
		curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOY)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPY) {

			if (ngbrFlags[idx] & NF_NNY) {

				curl.x += (halo_py[hidx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (halo_py[hidx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (halo_py[hidx].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (halo_py[hidx].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.y);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPY) {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - halo_ny[hidx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - halo_ny[hidx].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (cuVEC<VType>::quantity[idx].z - halo_ny[hidx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx].x - halo_ny[hidx].x) / (2 * cuVEC<VType>::h.y);
			}
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.y;
			curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / cuVEC<VType>::h.y;
			curl.z -= (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cuVEC<VType>::h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.y + bdiff_val.y.z) / 2;
				curl.z -= ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.y + bdiff_val.y.x) / 2;
			}
		}

		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += ((cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / cuVEC<VType>::h.y + bdiff_val.y.z) / 2;
				curl.z -= ((cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cuVEC<VType>::h.y + bdiff_val.y.x) / 2;
			}
		}
	}

	//z direction differentials
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
		curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOZ)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPZ) {

			if (ngbrFlags[idx] & NF_NNZ) {

				curl.x -= (halo_pz[hidx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (halo_pz[hidx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (halo_pz[hidx].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (halo_pz[hidx].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.z);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPZ) {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - halo_nz[hidx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - halo_nz[hidx].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (cuVEC<VType>::quantity[idx].y - halo_nz[hidx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx].x - halo_nz[hidx].x) / (2 * cuVEC<VType>::h.z);
			}
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.z;
			curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cuVEC<VType>::h.z;
			curl.y += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cuVEC<VType>::h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.z + bdiff_val.z.y) / 2;
				curl.y += ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.z + bdiff_val.z.x) / 2;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= ((cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cuVEC<VType>::h.z + bdiff_val.z.y) / 2;
				curl.y += ((cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cuVEC<VType>::h.z + bdiff_val.z.x) / 2;
			}
		}
	}

	return curl;
}

//Same as above but boundary conditions specified using a constant
template <typename VType>
__device__ VType cuVEC_VC<VType>::curl_nneu(int idx, const cuVAL3<VType>& bdiff) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
		curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOX)) {

		//halo index
		int hidx = ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y;

		if (ngbrFlags2[idx] & NF2_HALOPX) {

			if (ngbrFlags[idx] & NF_NNX) {

				curl.y -= (halo_px[hidx].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (halo_px[hidx].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (halo_px[hidx].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (halo_px[hidx].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.x);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPX) {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - halo_nx[hidx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - halo_nx[hidx].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (cuVEC<VType>::quantity[idx].z - halo_nx[hidx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx].y - halo_nx[hidx].y) / (2 * cuVEC<VType>::h.x);
			}
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.x;
			curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / cuVEC<VType>::h.x;
			curl.z += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / cuVEC<VType>::h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= ((cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.x + bdiff.x.z) / 2;
				curl.z += ((cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.x + bdiff.x.y) / 2;
			}
		}

		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= ((cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / cuVEC<VType>::h.x + bdiff.x.z) / 2;
				curl.z += ((cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / cuVEC<VType>::h.x + bdiff.x.y) / 2;
			}
		}
	}

	//y direction differentials
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
		curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOY)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPY) {

			if (ngbrFlags[idx] & NF_NNY) {

				curl.x += (halo_py[hidx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (halo_py[hidx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (halo_py[hidx].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (halo_py[hidx].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.y);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPY) {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - halo_ny[hidx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - halo_ny[hidx].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (cuVEC<VType>::quantity[idx].z - halo_ny[hidx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx].x - halo_ny[hidx].x) / (2 * cuVEC<VType>::h.y);
			}
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.y;
			curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / cuVEC<VType>::h.y;
			curl.z -= (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cuVEC<VType>::h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.y + bdiff.y.z) / 2;
				curl.z -= ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.y + bdiff.y.x) / 2;
			}
		}

		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += ((cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / cuVEC<VType>::h.y + bdiff.y.z) / 2;
				curl.z -= ((cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cuVEC<VType>::h.y + bdiff.y.x) / 2;
			}
		}
	}

	//z direction differentials
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
		curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOZ)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPZ) {

			if (ngbrFlags[idx] & NF_NNZ) {

				curl.x -= (halo_pz[hidx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (halo_pz[hidx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (halo_pz[hidx].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (halo_pz[hidx].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.z);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPZ) {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - halo_nz[hidx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - halo_nz[hidx].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (cuVEC<VType>::quantity[idx].y - halo_nz[hidx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx].x - halo_nz[hidx].x) / (2 * cuVEC<VType>::h.z);
			}
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.z;
			curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cuVEC<VType>::h.z;
			curl.y += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cuVEC<VType>::h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.z + bdiff.z.y) / 2;
				curl.y += ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.z + bdiff.z.x) / 2;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= ((cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cuVEC<VType>::h.z + bdiff.z.y) / 2;
				curl.y += ((cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cuVEC<VType>::h.z + bdiff.z.x) / 2;
			}
		}
	}

	return curl;
}

//curl operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
//can only be applied if VType is a VAL3
template <typename VType>
__device__ VType cuVEC_VC<VType>::curl_diri(int idx) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
		curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOX)) {

		//halo index
		int hidx = ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y;

		if (ngbrFlags2[idx] & NF2_HALOPX) {

			if (ngbrFlags[idx] & NF_NNX) {

				curl.y -= (halo_px[hidx].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (halo_px[hidx].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (halo_px[hidx].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (halo_px[hidx].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.x);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPX) {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - halo_nx[hidx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - halo_nx[hidx].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (cuVEC<VType>::quantity[idx].z - halo_nx[hidx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx].y - halo_nx[hidx].y) / (2 * cuVEC<VType>::h.x);
			}
		}
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

			curl.y -= (cuVEC<VType>::quantity[idx + 1].z + cuVEC<VType>::quantity[idx].z - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx).z) / (2 * cuVEC<VType>::h.x);
			curl.z += (cuVEC<VType>::quantity[idx + 1].y + cuVEC<VType>::quantity[idx].y - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx).y) / (2 * cuVEC<VType>::h.x);
		}
		else {

			curl.y -= (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx).z - cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
			curl.z += (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx).y - cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.x;
			curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / cuVEC<VType>::h.x;
			curl.z += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / cuVEC<VType>::h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.x);
			}
		}

		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
		}
	}

	//y direction differentials
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
		curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOY)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPY) {

			if (ngbrFlags[idx] & NF_NNY) {

				curl.x += (halo_py[hidx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (halo_py[hidx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (halo_py[hidx].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (halo_py[hidx].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.y);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPY) {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - halo_ny[hidx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - halo_ny[hidx].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (cuVEC<VType>::quantity[idx].z - halo_ny[hidx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx].x - halo_ny[hidx].x) / (2 * cuVEC<VType>::h.y);
			}
		}
	}
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

			curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z + cuVEC<VType>::quantity[idx].z - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx).z) / (2 * cuVEC<VType>::h.y);
			curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x + cuVEC<VType>::quantity[idx].x - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx).x) / (2 * cuVEC<VType>::h.y);
		}
		else {

			curl.x += (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx).z - cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
			curl.z -= (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx).x - cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.y;
			curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / cuVEC<VType>::h.y;
			curl.z -= (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cuVEC<VType>::h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.y);
			}
		}

		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
		}
	}

	//z direction differentials
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
		curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOZ)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPZ) {

			if (ngbrFlags[idx] & NF_NNZ) {

				curl.x -= (halo_pz[hidx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (halo_pz[hidx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (halo_pz[hidx].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (halo_pz[hidx].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.z);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPZ) {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - halo_nz[hidx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - halo_nz[hidx].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (cuVEC<VType>::quantity[idx].y - halo_nz[hidx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx].x - halo_nz[hidx].x) / (2 * cuVEC<VType>::h.z);
			}
		}
	}
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

			curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y + cuVEC<VType>::quantity[idx].y - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx).y) / (2 * cuVEC<VType>::h.z);
			curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x + cuVEC<VType>::quantity[idx].x - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx).x) / (2 * cuVEC<VType>::h.z);
		}
		else {

			curl.x -= (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx).y - cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
			curl.y += (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx).x - cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.z;
			curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cuVEC<VType>::h.z;
			curl.y += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cuVEC<VType>::h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.z);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
		}
	}

	return curl;
}

//curl operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
//Can be used at composite media boundaries where sided differentials will be used instead.
//can only be applied if VType is a VAL3
template <typename VType>
template <typename Class_BDiff>
__device__ VType cuVEC_VC<VType>::curl_diri_nneu(int idx, const Class_BDiff& bdiff_class) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
		curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOX)) {

		//halo index
		int hidx = ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y;

		if (ngbrFlags2[idx] & NF2_HALOPX) {

			if (ngbrFlags[idx] & NF_NNX) {

				curl.y -= (halo_px[hidx].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (halo_px[hidx].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (halo_px[hidx].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (halo_px[hidx].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.x);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPX) {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - halo_nx[hidx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - halo_nx[hidx].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (cuVEC<VType>::quantity[idx].z - halo_nx[hidx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx].y - halo_nx[hidx].y) / (2 * cuVEC<VType>::h.x);
			}
		}
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

			curl.y -= (cuVEC<VType>::quantity[idx + 1].z + cuVEC<VType>::quantity[idx].z - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx).z) / (2 * cuVEC<VType>::h.x);
			curl.z += (cuVEC<VType>::quantity[idx + 1].y + cuVEC<VType>::quantity[idx].y - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx).y) / (2 * cuVEC<VType>::h.x);
		}
		else {

			curl.y -= (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx).z - cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
			curl.z += (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx).y - cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.x;
			curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / cuVEC<VType>::h.x;
			curl.z += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / cuVEC<VType>::h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= ((cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.x + bdiff_val.x.z) / 2;
				curl.z += ((cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.x + bdiff_val.x.y) / 2;
			}
		}

		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= ((cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / cuVEC<VType>::h.x + bdiff_val.x.z) / 2;
				curl.z += ((cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / cuVEC<VType>::h.x + bdiff_val.x.y) / 2;
			}
		}
	}

	//y direction differentials
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
		curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOY)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPY) {

			if (ngbrFlags[idx] & NF_NNY) {

				curl.x += (halo_py[hidx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (halo_py[hidx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (halo_py[hidx].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (halo_py[hidx].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.y);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPY) {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - halo_ny[hidx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - halo_ny[hidx].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (cuVEC<VType>::quantity[idx].z - halo_ny[hidx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx].x - halo_ny[hidx].x) / (2 * cuVEC<VType>::h.y);
			}
		}
	}
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

			curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z + cuVEC<VType>::quantity[idx].z - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx).z) / (2 * cuVEC<VType>::h.y);
			curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x + cuVEC<VType>::quantity[idx].x - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx).x) / (2 * cuVEC<VType>::h.y);
		}
		else {

			curl.x += (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx).z - cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
			curl.z -= (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx).x - cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.y;
			curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / cuVEC<VType>::h.y;
			curl.z -= (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cuVEC<VType>::h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.y + bdiff_val.y.z) / 2;
				curl.z -= ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.y + bdiff_val.y.x) / 2;
			}
		}

		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += ((cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / cuVEC<VType>::h.y + bdiff_val.y.z) / 2;
				curl.z -= ((cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cuVEC<VType>::h.y + bdiff_val.y.x) / 2;
			}
		}
	}

	//z direction differentials
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
		curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOZ)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPZ) {

			if (ngbrFlags[idx] & NF_NNZ) {

				curl.x -= (halo_pz[hidx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (halo_pz[hidx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (halo_pz[hidx].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (halo_pz[hidx].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.z);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPZ) {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - halo_nz[hidx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - halo_nz[hidx].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (cuVEC<VType>::quantity[idx].y - halo_nz[hidx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx].x - halo_nz[hidx].x) / (2 * cuVEC<VType>::h.z);
			}
		}
	}
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

			curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y + cuVEC<VType>::quantity[idx].y - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx).y) / (2 * cuVEC<VType>::h.z);
			curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x + cuVEC<VType>::quantity[idx].x - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx).x) / (2 * cuVEC<VType>::h.z);
		}
		else {

			curl.x -= (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx).y - cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
			curl.y += (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx).x - cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.z;
			curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cuVEC<VType>::h.z;
			curl.y += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cuVEC<VType>::h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.z + bdiff_val.z.y) / 2;
				curl.y += ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.z + bdiff_val.z.x) / 2;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= ((cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cuVEC<VType>::h.z + bdiff_val.z.y) / 2;
				curl.y += ((cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cuVEC<VType>::h.z + bdiff_val.z.x) / 2;
			}
		}
	}

	return curl;
}

//Same as above but boundary conditions specified using a constant
template <typename VType>
__device__ VType cuVEC_VC<VType>::curl_diri_nneu(int idx, const cuVAL3<VType>& bdiff) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
		curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOX)) {

		//halo index
		int hidx = ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y;

		if (ngbrFlags2[idx] & NF2_HALOPX) {

			if (ngbrFlags[idx] & NF_NNX) {

				curl.y -= (halo_px[hidx].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (halo_px[hidx].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (halo_px[hidx].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (halo_px[hidx].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.x);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPX) {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - halo_nx[hidx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - halo_nx[hidx].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (cuVEC<VType>::quantity[idx].z - halo_nx[hidx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx].y - halo_nx[hidx].y) / (2 * cuVEC<VType>::h.x);
			}
		}
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

			curl.y -= (cuVEC<VType>::quantity[idx + 1].z + cuVEC<VType>::quantity[idx].z - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx).z) / (2 * cuVEC<VType>::h.x);
			curl.z += (cuVEC<VType>::quantity[idx + 1].y + cuVEC<VType>::quantity[idx].y - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx).y) / (2 * cuVEC<VType>::h.x);
		}
		else {

			curl.y -= (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx).z - cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
			curl.z += (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx).y - cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.x;
			curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / cuVEC<VType>::h.x;
			curl.z += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / cuVEC<VType>::h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= ((cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.x + bdiff.x.z) / 2;
				curl.z += ((cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.x + bdiff.x.y) / 2;
			}
		}

		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= ((cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / cuVEC<VType>::h.x + bdiff.x.z) / 2;
				curl.z += ((cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / cuVEC<VType>::h.x + bdiff.x.y) / 2;
			}
		}
	}

	//y direction differentials
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
		curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOY)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPY) {

			if (ngbrFlags[idx] & NF_NNY) {

				curl.x += (halo_py[hidx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (halo_py[hidx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (halo_py[hidx].z - cuVEC<VType>::quantity[idx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (halo_py[hidx].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.y);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPY) {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - halo_ny[hidx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - halo_ny[hidx].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (cuVEC<VType>::quantity[idx].z - halo_ny[hidx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx].x - halo_ny[hidx].x) / (2 * cuVEC<VType>::h.y);
			}
		}
	}
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

			curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z + cuVEC<VType>::quantity[idx].z - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx).z) / (2 * cuVEC<VType>::h.y);
			curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x + cuVEC<VType>::quantity[idx].x - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx).x) / (2 * cuVEC<VType>::h.y);
		}
		else {

			curl.x += (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx).z - cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
			curl.z -= (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx).x - cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.y;
			curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / cuVEC<VType>::h.y;
			curl.z -= (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cuVEC<VType>::h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.y + bdiff.y.z) / 2;
				curl.z -= ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.y + bdiff.y.x) / 2;
			}
		}

		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += ((cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / cuVEC<VType>::h.y + bdiff.y.z) / 2;
				curl.z -= ((cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cuVEC<VType>::h.y + bdiff.y.x) / 2;
			}
		}
	}

	//z direction differentials
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
		curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOZ)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPZ) {

			if (ngbrFlags[idx] & NF_NNZ) {

				curl.x -= (halo_pz[hidx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (halo_pz[hidx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (halo_pz[hidx].y - cuVEC<VType>::quantity[idx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (halo_pz[hidx].x - cuVEC<VType>::quantity[idx].x) / (2 * cuVEC<VType>::h.z);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPZ) {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - halo_nz[hidx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - halo_nz[hidx].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (cuVEC<VType>::quantity[idx].y - halo_nz[hidx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx].x - halo_nz[hidx].x) / (2 * cuVEC<VType>::h.z);
			}
		}
	}
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

			curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y + cuVEC<VType>::quantity[idx].y - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx).y) / (2 * cuVEC<VType>::h.z);
			curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x + cuVEC<VType>::quantity[idx].x - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx).x) / (2 * cuVEC<VType>::h.z);
		}
		else {

			curl.x -= (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx).y - cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
			curl.y += (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx).x - cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.z;
			curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cuVEC<VType>::h.z;
			curl.y += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cuVEC<VType>::h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.z + bdiff.z.y) / 2;
				curl.y += ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.z + bdiff.z.x) / 2;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= ((cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cuVEC<VType>::h.z + bdiff.z.y) / 2;
				curl.y += ((cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cuVEC<VType>::h.z + bdiff.z.x) / 2;
			}
		}
	}

	return curl;
}

//curl operator. Use sided differentials at boundaries (including at composite media boundaries)
//can only be applied if VType is a VAL3
template <typename VType>
__device__ VType cuVEC_VC<VType>::curl_sided(int idx) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
		curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOX)) {

		//halo index
		int hidx = ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y;

		if (ngbrFlags2[idx] & NF2_HALOPX) {

			if (ngbrFlags[idx] & NF_NNX) {

				curl.y -= (halo_px[hidx].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (halo_px[hidx].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (halo_px[hidx].z - cuVEC<VType>::quantity[idx].z) / (cuVEC<VType>::h.x);
				curl.z += (halo_px[hidx].y - cuVEC<VType>::quantity[idx].y) / (cuVEC<VType>::h.x);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPX) {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - halo_nx[hidx].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - halo_nx[hidx].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (cuVEC<VType>::quantity[idx].z - halo_nx[hidx].z) / (cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx].y - halo_nx[hidx].y) / (cuVEC<VType>::h.x);
			}
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (cuVEC<VType>::quantity[idx + 1].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.x;
				curl.z += (cuVEC<VType>::quantity[idx + 1].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.x;
			}
		}

		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].z - cuVEC<VType>::quantity[idx - 1].z) / (2 * cuVEC<VType>::h.x);
				curl.z += (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].y - cuVEC<VType>::quantity[idx - 1].y) / (2 * cuVEC<VType>::h.x);
			}
			else {

				curl.y -= (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - 1].z) / cuVEC<VType>::h.x;
				curl.z += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - 1].y) / cuVEC<VType>::h.x;
			}
		}
	}

	//y direction differentials
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
		curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOY)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPY) {

			if (ngbrFlags[idx] & NF_NNY) {

				curl.x += (halo_py[hidx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (halo_py[hidx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (halo_py[hidx].z - cuVEC<VType>::quantity[idx].z) / (cuVEC<VType>::h.y);
				curl.z -= (halo_py[hidx].x - cuVEC<VType>::quantity[idx].x) / (cuVEC<VType>::h.y);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPY) {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - halo_ny[hidx].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - halo_ny[hidx].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (cuVEC<VType>::quantity[idx].z - halo_ny[hidx].z) / (cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx].x - halo_ny[hidx].x) / (cuVEC<VType>::h.y);
			}
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.y;
				curl.z -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.y;
			}
		}

		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / (2 * cuVEC<VType>::h.y);
				curl.z -= (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1) * cuVEC<VType>::n.x].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / (2 * cuVEC<VType>::h.y);
			}
			else {

				curl.x += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z) / cuVEC<VType>::h.y;
				curl.z -= (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cuVEC<VType>::h.y;
			}
		}
	}

	//z direction differentials
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
		curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOZ)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPZ) {

			if (ngbrFlags[idx] & NF_NNZ) {

				curl.x -= (halo_pz[hidx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (halo_pz[hidx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (halo_pz[hidx].y - cuVEC<VType>::quantity[idx].y) / (cuVEC<VType>::h.z);
				curl.y += (halo_pz[hidx].x - cuVEC<VType>::quantity[idx].x) / (cuVEC<VType>::h.z);
			}
		}
		else {

			if (ngbrFlags[idx] & NF_NPZ) {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - halo_nz[hidx].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - halo_nz[hidx].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (cuVEC<VType>::quantity[idx].y - halo_nz[hidx].y) / (cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx].x - halo_nz[hidx].x) / (cuVEC<VType>::h.z);
			}
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.z;
				curl.y += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / (2 * cuVEC<VType>::h.z);
				curl.y += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / (2 * cuVEC<VType>::h.z);
			}
			else {

				curl.x -= (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cuVEC<VType>::h.z;
				curl.y += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cuVEC<VType>::h.z;
			}
		}
	}

	return curl;
}