#pragma once

#include "cuVEC_VC.h"

//-------------------------------- NEIGHBOR SUM

//calculate 6-point neighbor sum at given index
	//missing neighbors not added, including at boundaries, but taking into account pbc
template <typename VType>
__device__ VType cuVEC_VC<VType>::ngbr_sum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		sum = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx - 1]);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOX)) {

		//halo index
		int hidx = ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y;

		if (ngbrFlags2[idx] & NF2_HALOPX) {

			if (ngbrFlags[idx] & NF_NNX) sum = (halo_px[hidx] + cuVEC<VType>::quantity[idx - 1]);
			else sum = halo_px[hidx];
		}
		else {

			if (ngbrFlags[idx] & NF_NPX) sum = (cuVEC<VType>::quantity[idx + 1] + halo_nx[hidx]);
			else sum = halo_nx[hidx];
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = (cuVEC<VType>::quantity[idx + 1] + cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1]);
			}
			else {

				sum = cuVEC<VType>::quantity[idx + 1];
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)] + cuVEC<VType>::quantity[idx - 1]);
			}
			else {

				sum = cuVEC<VType>::quantity[idx - 1];
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOY)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPY) {

			if (ngbrFlags[idx] & NF_NNY) sum += (halo_py[hidx] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			else sum += (halo_py[hidx]);
		}
		else {

			if (ngbrFlags[idx] & NF_NPY) sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + halo_ny[hidx]);
			else sum += (halo_ny[hidx]);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]);
			}
			else {

				sum += cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x];
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			}
			else {

				sum += cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x];
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOZ)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPZ) {

			if (ngbrFlags[idx] & NF_NNZ) sum += (halo_pz[hidx] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			else sum += (halo_pz[hidx]);
		}
		else {

			if (ngbrFlags[idx] & NF_NPZ) sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + halo_nz[hidx]);
			else sum += (halo_nz[hidx]);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
			else {

				sum += cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y];
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
			else {

				sum += cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y];
			}
		}
	}

	return sum;
}

//same as ngbr_sum but sum normalised values only; for scalar values this is a trivial operation, but for vectors it's not.
template <typename VType>
__device__ VType cuVEC_VC<VType>::ngbr_dirsum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		sum = (cuVEC<VType>::quantity[idx + 1] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]) + cuVEC<VType>::quantity[idx - 1] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]));
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOX)) {

		//halo index
		int hidx = ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y;

		if (ngbrFlags2[idx] & NF2_HALOPX) {

			if (ngbrFlags[idx] & NF_NNX) sum = (halo_px[hidx] / cu_GetMagnitude(halo_px[hidx]) + cuVEC<VType>::quantity[idx - 1] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]));
			else sum = (halo_px[hidx] / cu_GetMagnitude(halo_px[hidx]));
		}
		else {

			if (ngbrFlags[idx] & NF_NPX) sum = (cuVEC<VType>::quantity[idx + 1] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]) + halo_nx[hidx] / cu_GetMagnitude(halo_nx[hidx]));
			else sum = (halo_nx[hidx] / cu_GetMagnitude(halo_nx[hidx]));
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = (cuVEC<VType>::quantity[idx + 1] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]) + cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1]));
			}
			else {

				sum = cuVEC<VType>::quantity[idx + 1] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)]) + cuVEC<VType>::quantity[idx - 1] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]));
			}
			else {

				sum = cuVEC<VType>::quantity[idx - 1] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]);
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]) + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]));
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOY)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPY) {

			if (ngbrFlags[idx] & NF_NNY) sum += (halo_py[hidx] / cu_GetMagnitude(halo_py[hidx]) + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]));
			else sum += (halo_py[hidx] / cu_GetMagnitude(halo_py[hidx]));
		}
		else {

			if (ngbrFlags[idx] & NF_NPY) sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]) + halo_ny[hidx] / cu_GetMagnitude(halo_ny[hidx]));
			else sum += (halo_ny[hidx] / cu_GetMagnitude(halo_ny[hidx]));
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]) + cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]));
			}
			else {

				sum += cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]) + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]));
			}
			else {

				sum += cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]) + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]));
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOZ)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPZ) {

			if (ngbrFlags[idx] & NF_NNZ) sum += (halo_pz[hidx] / cu_GetMagnitude(halo_pz[hidx]) + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]));
			else sum += (halo_pz[hidx] / cu_GetMagnitude(halo_pz[hidx]));
		}
		else {

			if (ngbrFlags[idx] & NF_NPZ) sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]) + halo_nz[hidx] / cu_GetMagnitude(halo_nz[hidx]));
			else sum += (halo_nz[hidx] / cu_GetMagnitude(halo_nz[hidx]));
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]) + cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y]));
			}
			else {

				sum += cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y]) + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]));
			}
			else {

				sum += cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
		}
	}

	return sum;
}

//calculate 6-point anisotropic neighbor sum at given index as rij x Vj over j points neighboring the point i at this index.
//missing neighbors not added, including at boundaries, but taking into account pbc
//only used if VType is a cuVAL3
template <typename VType>
__device__ VType cuVEC_VC<VType>::anisotropic_ngbr_sum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		sum = VType(0.0, -cuVEC<VType>::quantity[idx + 1].z, cuVEC<VType>::quantity[idx + 1].y) - VType(0.0, -cuVEC<VType>::quantity[idx - 1].z, cuVEC<VType>::quantity[idx - 1].y);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOX)) {

		//halo index
		int hidx = ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y;

		if (ngbrFlags2[idx] & NF2_HALOPX) {

			if (ngbrFlags[idx] & NF_NNX) sum = VType(0.0, -halo_px[hidx].z, halo_px[hidx].y) - VType(0.0, -cuVEC<VType>::quantity[idx - 1].z, cuVEC<VType>::quantity[idx - 1].y);
			else sum = VType(0.0, -halo_px[hidx].z, halo_px[hidx].y);
		}
		else {

			if (ngbrFlags[idx] & NF_NPX) sum = VType(0.0, -cuVEC<VType>::quantity[idx + 1].z, cuVEC<VType>::quantity[idx + 1].y) - VType(0.0, -halo_nx[hidx].z, halo_nx[hidx].y);
			else sum = VType(0.0, halo_nx[hidx].z, -halo_nx[hidx].y);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(0.0, -cuVEC<VType>::quantity[idx + 1].z, cuVEC<VType>::quantity[idx + 1].y) - cu_get_sign(pbc_x) * VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].y);
			}
			else {

				sum = VType(0.0, -cuVEC<VType>::quantity[idx + 1].z, cuVEC<VType>::quantity[idx + 1].y);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = cu_get_sign(pbc_x) * VType(0.0, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].z, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].y) - VType(0.0, -cuVEC<VType>::quantity[idx - 1].z, cuVEC<VType>::quantity[idx - 1].y);
			}
			else {

				sum = VType(0.0, cuVEC<VType>::quantity[idx - 1].z, -cuVEC<VType>::quantity[idx - 1].y);
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x) - VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOY)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPY) {

			if (ngbrFlags[idx] & NF_NNY) sum += VType(halo_py[hidx].z, 0.0, -halo_py[hidx].x) - VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x);
			else sum += VType(halo_py[hidx].z, 0.0, -halo_py[hidx].x);
		}
		else {

			if (ngbrFlags[idx] & NF_NPY) sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x) - VType(halo_ny[hidx].z, 0.0, -halo_ny[hidx].x);
			else sum += VType(-halo_ny[hidx].z, 0.0, halo_ny[hidx].x);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x) - cu_get_sign(pbc_y) * VType(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].x);
			}
			else {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += cu_get_sign(pbc_y) * VType(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].x) - VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x);
			}
			else {

				sum += VType(-cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, 0.0, +cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) - VType(-cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOZ)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPZ) {

			if (ngbrFlags[idx] & NF_NNZ) sum += VType(-halo_pz[hidx].y, halo_pz[hidx].x, 0.0) - VType(-cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0);
			else sum += VType(-halo_pz[hidx].y, halo_pz[hidx].x, 0.0);
		}
		else {

			if (ngbrFlags[idx] & NF_NPZ) sum += VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) - VType(-halo_nz[hidx].y, halo_nz[hidx].x, 0.0);
			else sum += VType(halo_nz[hidx].y, -halo_nz[hidx].x, 0.0);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) - cu_get_sign(pbc_z) * VType(-cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0);
			}
			else {

				sum += VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += cu_get_sign(pbc_z) * VType(-cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) - VType(-cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0);
			}
			else {

				sum += VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0);
			}
		}
	}

	return sum;
}

//same as anisotropic_ngbr_sum but sum normalised values only; for scalar values this is a trivial operation, but for vectors it's not.
template <typename VType>
__device__ VType cuVEC_VC<VType>::anisotropic_ngbr_dirsum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		sum = VType(0.0, -cuVEC<VType>::quantity[idx + 1].z, cuVEC<VType>::quantity[idx + 1].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]) - VType(0.0, -cuVEC<VType>::quantity[idx - 1].z, cuVEC<VType>::quantity[idx - 1].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOX)) {

		//halo index
		int hidx = ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y;

		if (ngbrFlags2[idx] & NF2_HALOPX) {

			if (ngbrFlags[idx] & NF_NNX) sum = VType(0.0, -halo_px[hidx].z, halo_px[hidx].y) / cu_GetMagnitude(halo_px[hidx]) - VType(0.0, -cuVEC<VType>::quantity[idx - 1].z, cuVEC<VType>::quantity[idx - 1].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]);
			else sum = VType(0.0, -halo_px[hidx].z, halo_px[hidx].y) / cu_GetMagnitude(halo_px[hidx]);
		}
		else {

			if (ngbrFlags[idx] & NF_NPX) sum = VType(0.0, -cuVEC<VType>::quantity[idx + 1].z, cuVEC<VType>::quantity[idx + 1].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]) - VType(0.0, -halo_nx[hidx].z, halo_nx[hidx].y) / cu_GetMagnitude(halo_nx[hidx]);
			else sum = VType(0.0, halo_nx[hidx].z, -halo_nx[hidx].y) / cu_GetMagnitude(halo_nx[hidx]);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(0.0, -cuVEC<VType>::quantity[idx + 1].z, cuVEC<VType>::quantity[idx + 1].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]) - cu_get_sign(pbc_x) * VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1]);
			}
			else {

				sum = VType(0.0, -cuVEC<VType>::quantity[idx + 1].z, cuVEC<VType>::quantity[idx + 1].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = cu_get_sign(pbc_x) * VType(0.0, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].z, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)]) - VType(0.0, -cuVEC<VType>::quantity[idx - 1].z, cuVEC<VType>::quantity[idx - 1].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]);
			}
			else {

				sum = VType(0.0, cuVEC<VType>::quantity[idx - 1].z, -cuVEC<VType>::quantity[idx - 1].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]);
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]) - VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOY)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPY) {

			if (ngbrFlags[idx] & NF_NNY) sum += VType(halo_py[hidx].z, 0.0, -halo_py[hidx].x) / cu_GetMagnitude(halo_py[hidx]) - VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			else sum += VType(halo_py[hidx].z, 0.0, -halo_py[hidx].x) / cu_GetMagnitude(halo_py[hidx]);
		}
		else {

			if (ngbrFlags[idx] & NF_NPY) sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]) - VType(halo_ny[hidx].z, 0.0, -halo_ny[hidx].x) / cu_GetMagnitude(halo_ny[hidx]);
			else sum += VType(-halo_ny[hidx].z, 0.0, halo_ny[hidx].x) / cu_GetMagnitude(halo_ny[hidx]);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]) - cu_get_sign(pbc_y) * VType(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]);
			}
			else {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += cu_get_sign(pbc_y) * VType(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]) - VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			}
			else {

				sum += VType(-cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, 0.0, +cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]) - VType(-cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOZ)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPZ) {

			if (ngbrFlags[idx] & NF_NNZ) sum += VType(-halo_pz[hidx].y, halo_pz[hidx].x, 0.0) / cu_GetMagnitude(halo_pz[hidx]) - VType(-cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			else sum += VType(-halo_pz[hidx].y, halo_pz[hidx].x, 0.0) / cu_GetMagnitude(halo_pz[hidx]);
		}
		else {

			if (ngbrFlags[idx] & NF_NPZ) sum += VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]) - VType(-halo_nz[hidx].y, halo_nz[hidx].x, 0.0) / cu_GetMagnitude(halo_nz[hidx]);
			else sum += VType(halo_nz[hidx].y, -halo_nz[hidx].x, 0.0) / cu_GetMagnitude(halo_nz[hidx]);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]) - cu_get_sign(pbc_z) * VType(-cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
			else {

				sum += VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += cu_get_sign(pbc_z) * VType(-cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y]) - VType(-cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
			else {

				sum += VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
		}
	}

	return sum;
}

//-------------------------------- Z SYMMETRY ANISOTROPIC SIMPLE NEIGHBOR SUM

//calculate 6-point anisotropic neighbor sum at given index as (rij x z) x Vj over j points neighboring the point i at this index.
//missing neighbors not added, including at boundaries, but taking into account pbc
//only used if VType is a VAL3
template <typename VType>
__device__ VType cuVEC_VC<VType>::zanisotropic_ngbr_sum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		sum = VType(-cuVEC<VType>::quantity[idx + 1].z, 0.0, cuVEC<VType>::quantity[idx + 1].x) - VType(-cuVEC<VType>::quantity[idx - 1].z, 0.0, cuVEC<VType>::quantity[idx - 1].x);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOX)) {

		//halo index
		int hidx = ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y;

		if (ngbrFlags2[idx] & NF2_HALOPX) {

			if (ngbrFlags[idx] & NF_NNX) sum = VType(-halo_px[hidx].z, 0.0, halo_px[hidx].x) - VType(-cuVEC<VType>::quantity[idx - 1].z, 0.0, cuVEC<VType>::quantity[idx - 1].x);
			else sum = VType(-halo_px[hidx].z, 0.0, halo_px[hidx].x);
		}
		else {

			if (ngbrFlags[idx] & NF_NPX) sum = VType(-cuVEC<VType>::quantity[idx + 1].z, 0.0, cuVEC<VType>::quantity[idx + 1].x) - VType(-halo_nx[hidx].z, 0.0, halo_nx[hidx].x);
			else sum = VType(halo_nx[hidx].z, 0.0, -halo_nx[hidx].x);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-cuVEC<VType>::quantity[idx + 1].z, 0.0, cuVEC<VType>::quantity[idx + 1].x) - cu_get_sign(pbc_x) * VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].z, 0.0, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].x);
			}
			else {

				sum = VType(-cuVEC<VType>::quantity[idx + 1].z, 0.0, cuVEC<VType>::quantity[idx + 1].x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = cu_get_sign(pbc_x) * VType(-cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].z, 0.0, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].x) - VType(-cuVEC<VType>::quantity[idx - 1].z, 0.0, cuVEC<VType>::quantity[idx - 1].x);
			}
			else {

				sum = VType(cuVEC<VType>::quantity[idx - 1].z, 0.0, -cuVEC<VType>::quantity[idx - 1].x);
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y) - VType(0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOY)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPY) {

			if (ngbrFlags[idx] & NF_NNY) sum += VType(0.0, -halo_py[hidx].z, halo_py[hidx].y) - VType(0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y);
			else sum += VType(0.0, -halo_py[hidx].z, halo_py[hidx].y);
		}
		else {

			if (ngbrFlags[idx] & NF_NPY) sum += VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y) - VType(0.0, -halo_ny[hidx].z, halo_ny[hidx].y);
			else sum += VType(0.0, halo_ny[hidx].z, -halo_ny[hidx].y);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y) - cu_get_sign(pbc_y) * VType(0.0, -cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y);
			}
			else {

				sum += VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += cu_get_sign(pbc_y) * VType(0.0, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y) - VType(0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y);
			}
			else {

				sum += VType(0.0, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y);
			}
		}
	}

	return sum;
}

//same as zanisotropic_ngbr_sum but sum normalised values only; for scalar values this is a trivial operation, but for vectors it's not.
template <typename VType>
__device__ VType cuVEC_VC<VType>::zanisotropic_ngbr_dirsum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		sum = VType(-cuVEC<VType>::quantity[idx + 1].z, 0.0, cuVEC<VType>::quantity[idx + 1].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]) - VType(-cuVEC<VType>::quantity[idx - 1].z, 0.0, cuVEC<VType>::quantity[idx - 1].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOX)) {

		//halo index
		int hidx = ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y;

		if (ngbrFlags2[idx] & NF2_HALOPX) {

			if (ngbrFlags[idx] & NF_NNX) sum = VType(-halo_px[hidx].z, 0.0, halo_px[hidx].x) / cu_GetMagnitude(halo_px[hidx]) - VType(-cuVEC<VType>::quantity[idx - 1].z, 0.0, cuVEC<VType>::quantity[idx - 1].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]);
			else sum = VType(-halo_px[hidx].z, 0.0, halo_px[hidx].x) / cu_GetMagnitude(halo_px[hidx]);
		}
		else {

			if (ngbrFlags[idx] & NF_NPX) sum = VType(-cuVEC<VType>::quantity[idx + 1].z, 0.0, cuVEC<VType>::quantity[idx + 1].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]) - VType(-halo_nx[hidx].z, 0.0, halo_nx[hidx].x) / cu_GetMagnitude(halo_nx[hidx]);
			else sum = VType(halo_nx[hidx].z, 0.0, -halo_nx[hidx].x) / cu_GetMagnitude(halo_nx[hidx]);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-cuVEC<VType>::quantity[idx + 1].z, 0.0, cuVEC<VType>::quantity[idx + 1].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]) - cu_get_sign(pbc_x) * VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].z, 0.0, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1]);
			}
			else {

				sum = VType(-cuVEC<VType>::quantity[idx + 1].z, 0.0, cuVEC<VType>::quantity[idx + 1].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = cu_get_sign(pbc_x) * VType(-cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].z, 0.0, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)]) - VType(-cuVEC<VType>::quantity[idx - 1].z, 0.0, cuVEC<VType>::quantity[idx - 1].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]);
			}
			else {

				sum = VType(cuVEC<VType>::quantity[idx - 1].z, 0.0, -cuVEC<VType>::quantity[idx - 1].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]);
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]) - VType(0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOY)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPY) {

			if (ngbrFlags[idx] & NF_NNY) sum += VType(0.0, -halo_py[hidx].z, halo_py[hidx].y) / cu_GetMagnitude(halo_py[hidx]) - VType(0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			else sum += VType(0.0, -halo_py[hidx].z, halo_py[hidx].y) / cu_GetMagnitude(halo_py[hidx]);
		}
		else {

			if (ngbrFlags[idx] & NF_NPY) sum += VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]) - VType(0.0, -halo_ny[hidx].z, halo_ny[hidx].y) / cu_GetMagnitude(halo_ny[hidx]);
			else sum += VType(0.0, halo_ny[hidx].z, -halo_ny[hidx].y) / cu_GetMagnitude(halo_ny[hidx]);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]) - cu_get_sign(pbc_y) * VType(0.0, -cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]);
			}
			else {

				sum += VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += cu_get_sign(pbc_y) * VType(0.0, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]) - VType(0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			}
			else {

				sum += VType(0.0, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			}
		}
	}

	return sum;
}

//-------------------------------- Y SYMMETRY ANISOTROPIC SIMPLE NEIGHBOR SUM

//calculate 6-point anisotropic neighbor sum at given index as (rij x y) x Vj over j points neighboring the point i at this index.
//missing neighbors not added, including at boundaries, but taking into account pbc
//only used if VType is a VAL3
template <typename VType>
__device__ VType cuVEC_VC<VType>::yanisotropic_ngbr_sum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		sum = VType(-cuVEC<VType>::quantity[idx + 1].y, cuVEC<VType>::quantity[idx + 1].x, 0.0) - VType(-cuVEC<VType>::quantity[idx - 1].y, cuVEC<VType>::quantity[idx - 1].x, 0.0);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOX)) {

		//halo index
		int hidx = ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y;

		if (ngbrFlags2[idx] & NF2_HALOPX) {

			if (ngbrFlags[idx] & NF_NNX) sum = VType(-halo_px[hidx].y, halo_px[hidx].x, 0.0) - VType(-cuVEC<VType>::quantity[idx - 1].y, cuVEC<VType>::quantity[idx - 1].x, 0.0);
			else sum = VType(-halo_px[hidx].y, halo_px[hidx].x, 0.0);
		}
		else {

			if (ngbrFlags[idx] & NF_NPX) sum = VType(-cuVEC<VType>::quantity[idx + 1].y, cuVEC<VType>::quantity[idx + 1].x, 0.0) - VType(-halo_nx[hidx].y, halo_nx[hidx].x, 0.0);
			else sum = VType(halo_nx[hidx].y, -halo_nx[hidx].x, 0.0);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-cuVEC<VType>::quantity[idx + 1].y, cuVEC<VType>::quantity[idx + 1].x, 0.0) - cu_get_sign(pbc_x) * VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].y, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].x, 0.0);
			}
			else {

				sum = VType(-cuVEC<VType>::quantity[idx + 1].y, cuVEC<VType>::quantity[idx + 1].x, 0.0);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = cu_get_sign(pbc_x) * VType(-cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].y, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].x, 0.0) - VType(-cuVEC<VType>::quantity[idx - 1].y, cuVEC<VType>::quantity[idx - 1].x, 0.0);
			}
			else {

				sum = VType(cuVEC<VType>::quantity[idx - 1].y, -cuVEC<VType>::quantity[idx - 1].x, 0.0);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += VType(0.0, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y)
			- VType(0.0, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOZ)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPZ) {

			if (ngbrFlags[idx] & NF_NNZ)
				sum += VType(0.0, halo_pz[hidx].z, -halo_pz[hidx].y)
				- VType(0.0, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y);
			else sum += VType(0.0, halo_pz[hidx].z, -halo_pz[hidx].y);
		}
		else {

			if (ngbrFlags[idx] & NF_NPZ) 
				sum += VType(0.0, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y)
				- VType(0.0, halo_nz[hidx].z, -halo_nz[hidx].y);
			else sum += VType(0.0, -halo_nz[hidx].z, halo_nz[hidx].y);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(0.0, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y)
					- cu_get_sign(pbc_z) * VType(0.0, cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].y);
			}
			else {

				sum += VType(0.0, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += cu_get_sign(pbc_z) * VType(0.0, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].y)
					- VType(0.0, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y);
			}
			else {

				sum += VType(0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y);
			}
		}
	}

	return sum;
}

//same as zanisotropic_ngbr_sum but sum normalised values only; for scalar values this is a trivial operation, but for vectors it's not.
template <typename VType>
__device__ VType cuVEC_VC<VType>::yanisotropic_ngbr_dirsum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		sum = VType(-cuVEC<VType>::quantity[idx + 1].y, cuVEC<VType>::quantity[idx + 1].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]) - VType(-cuVEC<VType>::quantity[idx - 1].y, cuVEC<VType>::quantity[idx - 1].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOX)) {

		//halo index
		int hidx = ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y;

		if (ngbrFlags2[idx] & NF2_HALOPX) {

			if (ngbrFlags[idx] & NF_NNX) sum = VType(-halo_px[hidx].y, halo_px[hidx].x, 0.0) / cu_GetMagnitude(halo_px[hidx]) - VType(-cuVEC<VType>::quantity[idx - 1].y, cuVEC<VType>::quantity[idx - 1].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]);
			else sum = VType(-halo_px[hidx].y, halo_px[hidx].x, 0.0) / cu_GetMagnitude(halo_px[hidx]);
		}
		else {

			if (ngbrFlags[idx] & NF_NPX) sum = VType(-cuVEC<VType>::quantity[idx + 1].y, cuVEC<VType>::quantity[idx + 1].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]) - VType(-halo_nx[hidx].y, halo_nx[hidx].x, 0.0) / cu_GetMagnitude(halo_nx[hidx]);
			else sum = VType(halo_nx[hidx].y, -halo_nx[hidx].x, 0.0) / cu_GetMagnitude(halo_nx[hidx]);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-cuVEC<VType>::quantity[idx + 1].y, cuVEC<VType>::quantity[idx + 1].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1])
					- cu_get_sign(pbc_x) * VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].y, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1]);
			}
			else {

				sum = VType(-cuVEC<VType>::quantity[idx + 1].y, cuVEC<VType>::quantity[idx + 1].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = cu_get_sign(pbc_x) * VType(-cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].y, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)])
					- VType(-cuVEC<VType>::quantity[idx - 1].y, cuVEC<VType>::quantity[idx - 1].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]);
			}
			else {

				sum = VType(cuVEC<VType>::quantity[idx - 1].y, -cuVEC<VType>::quantity[idx - 1].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += VType(0.0, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y])
			- VType(0.0, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOZ)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPZ) {

			if (ngbrFlags[idx] & NF_NNZ) 
				sum += VType(0.0, halo_pz[hidx].z, -halo_pz[hidx].y) / cu_GetMagnitude(halo_pz[hidx])
				- VType(0.0, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			else sum += VType(0.0, halo_pz[hidx].z, -halo_pz[hidx].y) / cu_GetMagnitude(halo_pz[hidx]);
		}
		else {

			if (ngbrFlags[idx] & NF_NPZ) 
				sum += VType(0.0, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y])
				- VType(0.0, halo_nz[hidx].z, -halo_nz[hidx].y) / cu_GetMagnitude(halo_nz[hidx]);
			else sum += VType(0.0, -halo_nz[hidx].z, halo_nz[hidx].y) / cu_GetMagnitude(halo_nz[hidx]);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(0.0, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y])
					- cu_get_sign(pbc_z) * VType(0.0, cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
			else {

				sum += VType(0.0, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += cu_get_sign(pbc_z) * VType(0.0, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y])
					- VType(0.0, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
			else {

				sum += VType(0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
		}
	}

	return sum;
}

//-------------------------------- X SYMMETRY ANISOTROPIC SIMPLE NEIGHBOR SUM

//calculate 6-point anisotropic neighbor sum at given index as (rij x x) x Vj over j points neighboring the point i at this index.
//missing neighbors not added, including at boundaries, but taking into account pbc
//only used if VType is a VAL3
template <typename VType>
__device__ VType cuVEC_VC<VType>::xanisotropic_ngbr_sum(int idx) const
{
	VType sum = VType();

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x, 0.0) - VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x, 0.0);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOY)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPY) {

			if (ngbrFlags[idx] & NF_NNY) sum += VType(halo_py[hidx].y, -halo_py[hidx].x, 0.0) - VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x, 0.0);
			else sum += VType(halo_py[hidx].y, -halo_py[hidx].x, 0.0);
		}
		else {

			if (ngbrFlags[idx] & NF_NPY) sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x, 0.0) - VType(halo_ny[hidx].y, -halo_ny[hidx].x, 0.0);
			else sum += VType(-halo_ny[hidx].y, halo_ny[hidx].x, 0.0);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x, 0.0) - cu_get_sign(pbc_y) * VType(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].x, 0.0);
			}
			else {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x, 0.0);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += cu_get_sign(pbc_y) * VType(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].x, 0.0) - VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x, 0.0);
			}
			else {

				sum += VType(-cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x, 0.0);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x)
			- VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOZ)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPZ) {

			if (ngbrFlags[idx] & NF_NNZ) 
				sum += VType(halo_pz[hidx].z, 0.0, -halo_pz[hidx].x)
				- VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x);
			else sum += VType(halo_pz[hidx].z, 0.0, -halo_pz[hidx].x);
		}
		else {

			if (ngbrFlags[idx] & NF_NPZ) 
				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x)
				- VType(halo_nz[hidx].z, 0.0, -halo_nz[hidx].x);
			else sum += VType(-halo_nz[hidx].z, 0.0, halo_nz[hidx].x);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x)
					- cu_get_sign(pbc_z) * VType(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].x);
			}
			else {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += cu_get_sign(pbc_z) * VType(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].x)
					- VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x);
			}
			else {

				sum += VType(-cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x);
			}
		}
	}

	return sum;
}

//same as zanisotropic_ngbr_sum but sum normalised values only; for scalar values this is a trivial operation, but for vectors it's not.
template <typename VType>
__device__ VType cuVEC_VC<VType>::xanisotropic_ngbr_dirsum(int idx) const
{
	VType sum = VType();

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x])
			- VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOY)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + (idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPY) {

			if (ngbrFlags[idx] & NF_NNY) 
				sum += VType(halo_py[hidx].y, -halo_py[hidx].x, 0.0) / cu_GetMagnitude(halo_py[hidx])
				- VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			else sum += VType(halo_py[hidx].y, -halo_py[hidx].x, 0.0) / cu_GetMagnitude(halo_py[hidx]);
		}
		else {

			if (ngbrFlags[idx] & NF_NPY) 
				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x])
				- VType(halo_ny[hidx].y, -halo_ny[hidx].x, 0.0) / cu_GetMagnitude(halo_ny[hidx]);
			else sum += VType(-halo_ny[hidx].y, halo_ny[hidx].x, 0.0) / cu_GetMagnitude(halo_ny[hidx]);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x])
					- cu_get_sign(pbc_y) * VType(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]);
			}
			else {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += cu_get_sign(pbc_y) * VType(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x])
					- VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			}
			else {

				sum += VType(-cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y])
			- VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
	}
	//use halo region?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_HALOZ)) {

		//halo index
		int hidx = (idx % cuVEC<VType>::n.x) + ((idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x;

		if (ngbrFlags2[idx] & NF2_HALOPZ) {

			if (ngbrFlags[idx] & NF_NNZ) 
				sum += VType(halo_pz[hidx].z, 0.0, -halo_pz[hidx].x) / cu_GetMagnitude(halo_pz[hidx])
				- VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			else sum += VType(halo_pz[hidx].z, 0.0, -halo_pz[hidx].x) / cu_GetMagnitude(halo_pz[hidx]);
		}
		else {

			if (ngbrFlags[idx] & NF_NPZ) 
				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y])
				- VType(halo_nz[hidx].z, 0.0, -halo_nz[hidx].x) / cu_GetMagnitude(halo_nz[hidx]);
			else sum += VType(-halo_nz[hidx].z, 0.0, halo_nz[hidx].x) / cu_GetMagnitude(halo_nz[hidx]);
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y])
					- cu_get_sign(pbc_z) * VType(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
			else {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += cu_get_sign(pbc_z) * VType(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1)*cuVEC<VType>::n.x*cuVEC<VType>::n.y])
					- VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
			else {

				sum += VType(-cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z, 0.0, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
		}
	}

	return sum;
}
