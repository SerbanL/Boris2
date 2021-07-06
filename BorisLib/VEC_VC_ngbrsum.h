#pragma once

#include "VEC_VC.h"

//-------------------------------- SIMPLE NEIGHBOR SUM

//calculate 6-point neighbor sum at given index
//missing neighbors not added, including at boundaries, but taking into account pbc
template <typename VType>
VType VEC_VC<VType>::ngbr_sum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		sum = (VEC<VType>::quantity[idx + 1] + VEC<VType>::quantity[idx - 1]);
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = (VEC<VType>::quantity[idx + 1] + get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1]);
			}
			else {

				sum = VEC<VType>::quantity[idx + 1];
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] + VEC<VType>::quantity[idx - 1]);
			}
			else {

				sum = VEC<VType>::quantity[idx - 1];
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += (VEC<VType>::quantity[idx + VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x]);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (VEC<VType>::quantity[idx + VEC<VType>::n.x] + get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x]);
			}
			else {

				sum += VEC<VType>::quantity[idx + VEC<VType>::n.x];
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x]);
			}
			else {

				sum += VEC<VType>::quantity[idx - VEC<VType>::n.x];
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y]);
			}
			else {

				sum += VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y];
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
			}
			else {

				sum += VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y];
			}
		}
	}

	return sum;
}

//same as ngbr_sum but sum normalised values only; for scalar values this is a trivial operation, but for vectors it's not.
template <typename VType>
VType VEC_VC<VType>::ngbr_dirsum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {
		
		//both neighbors available so must be an inner point along this direction
		sum = (VEC<VType>::quantity[idx + 1] / GetMagnitude(VEC<VType>::quantity[idx + 1]) + VEC<VType>::quantity[idx - 1] / GetMagnitude(VEC<VType>::quantity[idx - 1]));
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = (VEC<VType>::quantity[idx + 1] / GetMagnitude(VEC<VType>::quantity[idx + 1]) + get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1] / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x - 1]));
			}
			else {

				sum = VEC<VType>::quantity[idx + 1] / GetMagnitude(VEC<VType>::quantity[idx + 1]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				//sum = (VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] / GetMagnitude(VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)]) + VEC<VType>::quantity[idx - 1] / GetMagnitude(VEC<VType>::quantity[idx - 1]));
				sum = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] / GetMagnitude(VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)]) + VEC<VType>::quantity[idx - 1] / GetMagnitude(VEC<VType>::quantity[idx - 1]));
			}
			else {

				sum = VEC<VType>::quantity[idx - 1] / GetMagnitude(VEC<VType>::quantity[idx - 1]);
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += (VEC<VType>::quantity[idx + VEC<VType>::n.x] / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x]) + VEC<VType>::quantity[idx - VEC<VType>::n.x] / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x]));
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (VEC<VType>::quantity[idx + VEC<VType>::n.x] / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x]) + get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x] / GetMagnitude(VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x]));
			}
			else {

				sum += VEC<VType>::quantity[idx + VEC<VType>::n.x] / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x] / GetMagnitude(VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x]) + VEC<VType>::quantity[idx - VEC<VType>::n.x] / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x]));
			}
			else {

				sum += VEC<VType>::quantity[idx - VEC<VType>::n.x] / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x]);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y]) + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]));
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y]) + get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] / GetMagnitude(VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y]));
			}
			else {

				sum += VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y]);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] / GetMagnitude(VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y]) + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]));
			}
			else {

				sum += VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
			}
		}
	}

	return sum;
}

//-------------------------------- ANISOTROPIC SIMPLE NEIGHBOR SUM

//calculate 6-point anisotropic neighbor sum at given index as rij x Vj over j points neighboring the point i at this index.
//missing neighbors not added, including at boundaries, but taking into account pbc
//only used if VType is a VAL3
template <typename VType>
VType VEC_VC<VType>::anisotropic_ngbr_sum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		sum = VType(0.0, -VEC<VType>::quantity[idx + 1].z, VEC<VType>::quantity[idx + 1].y) - VType(0.0, -VEC<VType>::quantity[idx - 1].z, VEC<VType>::quantity[idx - 1].y);
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(0.0, -VEC<VType>::quantity[idx + 1].z, VEC<VType>::quantity[idx + 1].y) - get_sign(pbc_x) * VType(0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].z, VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].y);
			}
			else {

				sum = VType(0.0, -VEC<VType>::quantity[idx + 1].z, VEC<VType>::quantity[idx + 1].y);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = get_sign(pbc_x) * VType(0.0, -VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].z, VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].y) - VType(0.0, -VEC<VType>::quantity[idx - 1].z, VEC<VType>::quantity[idx - 1].y);
			}
			else {

				sum = VType(0.0, VEC<VType>::quantity[idx - 1].z, -VEC<VType>::quantity[idx - 1].y);
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x].z, 0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x].x) - VType(VEC<VType>::quantity[idx - VEC<VType>::n.x].z, 0.0, -VEC<VType>::quantity[idx - VEC<VType>::n.x].x);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x].z, 0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x].x) - get_sign(pbc_y) * VType(VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].z, 0.0, -VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].x);
			}
			else {

				sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x].z, 0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x].x);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += get_sign(pbc_y) * VType(VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].z, 0.0, -VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].x) - VType(VEC<VType>::quantity[idx - VEC<VType>::n.x].z, 0.0, -VEC<VType>::quantity[idx - VEC<VType>::n.x].x);
			}
			else {

				sum += VType(-VEC<VType>::quantity[idx - VEC<VType>::n.x].z, 0.0, +VEC<VType>::quantity[idx - VEC<VType>::n.x].x);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += VType(-VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y, VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0) - VType(-VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y, VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(-VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y, VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0) - get_sign(pbc_z) * VType(-VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].y, VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0);
			}
			else {

				sum += VType(-VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y, VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += get_sign(pbc_z) * VType(-VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].y, VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0) - VType(-VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y, VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0);
			}
			else {

				sum += VType(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y, -VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0);
			}
		}
	}

	return sum;
}

//same as anisotropic_ngbr_sum but sum normalised values only; for scalar values this is a trivial operation, but for vectors it's not.
template <typename VType>
VType VEC_VC<VType>::anisotropic_ngbr_dirsum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		sum = VType(0.0, -VEC<VType>::quantity[idx + 1].z, VEC<VType>::quantity[idx + 1].y) / GetMagnitude(VEC<VType>::quantity[idx + 1]) - VType(0.0, -VEC<VType>::quantity[idx - 1].z, VEC<VType>::quantity[idx - 1].y) / GetMagnitude(VEC<VType>::quantity[idx - 1]);
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(0.0, -VEC<VType>::quantity[idx + 1].z, VEC<VType>::quantity[idx + 1].y) / GetMagnitude(VEC<VType>::quantity[idx + 1]) - get_sign(pbc_x) * VType(0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].z, VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].y) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x - 1]);
			}
			else {

				sum = VType(0.0, -VEC<VType>::quantity[idx + 1].z, VEC<VType>::quantity[idx + 1].y) / GetMagnitude(VEC<VType>::quantity[idx + 1]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = get_sign(pbc_x) * VType(0.0, -VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].z, VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].y) / GetMagnitude(VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)]) - VType(0.0, -VEC<VType>::quantity[idx - 1].z, VEC<VType>::quantity[idx - 1].y) / GetMagnitude(VEC<VType>::quantity[idx - 1]);
			}
			else {

				sum = VType(0.0, VEC<VType>::quantity[idx - 1].z, -VEC<VType>::quantity[idx - 1].y) / GetMagnitude(VEC<VType>::quantity[idx - 1]);
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x].z, 0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x].x) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x]) - VType(VEC<VType>::quantity[idx - VEC<VType>::n.x].z, 0.0, -VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x]);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x].z, 0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x].x) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x]) - get_sign(pbc_y) * VType(VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].z, 0.0, -VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].x) / GetMagnitude(VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x]);
			}
			else {

				sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x].z, 0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x].x) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += get_sign(pbc_y) * VType(VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].z, 0.0, -VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].x) / GetMagnitude(VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x]) - VType(VEC<VType>::quantity[idx - VEC<VType>::n.x].z, 0.0, -VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x]);
			}
			else {

				sum += VType(-VEC<VType>::quantity[idx - VEC<VType>::n.x].z, 0.0, +VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x]);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += VType(-VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y, VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y]) - VType(-VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y, VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(-VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y, VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y]) - get_sign(pbc_z) * VType(-VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].y, VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y]);
			}
			else {

				sum += VType(-VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y, VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y]);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += get_sign(pbc_z) * VType(-VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].y, VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y]) - VType(-VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y, VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
			}
			else {

				sum += VType(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y, -VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
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
VType VEC_VC<VType>::zanisotropic_ngbr_sum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		sum = VType(-VEC<VType>::quantity[idx + 1].z, 0.0, VEC<VType>::quantity[idx + 1].x) - VType(-VEC<VType>::quantity[idx - 1].z, 0.0, VEC<VType>::quantity[idx - 1].x);
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-VEC<VType>::quantity[idx + 1].z, 0.0, VEC<VType>::quantity[idx + 1].x) - get_sign(pbc_x) * VType(-VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].z, 0.0, VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].x);
			}
			else {

				sum = VType(-VEC<VType>::quantity[idx + 1].z, 0.0, VEC<VType>::quantity[idx + 1].x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = get_sign(pbc_x) * VType(-VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].z, 0.0, VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].x) - VType(-VEC<VType>::quantity[idx - 1].z, 0.0, VEC<VType>::quantity[idx - 1].x);
			}
			else {

				sum = VType(VEC<VType>::quantity[idx - 1].z, 0.0, -VEC<VType>::quantity[idx - 1].x);
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x].z, VEC<VType>::quantity[idx + VEC<VType>::n.x].y) - VType(0.0, -VEC<VType>::quantity[idx - VEC<VType>::n.x].z, VEC<VType>::quantity[idx - VEC<VType>::n.x].y);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x].z, VEC<VType>::quantity[idx + VEC<VType>::n.x].y) - get_sign(pbc_y) * VType(0.0, -VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].z, VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y);
			}
			else {

				sum += VType(0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x].z, VEC<VType>::quantity[idx + VEC<VType>::n.x].y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += get_sign(pbc_y) * VType(0.0, -VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].z, VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y) - VType(0.0, -VEC<VType>::quantity[idx - VEC<VType>::n.x].z, VEC<VType>::quantity[idx - VEC<VType>::n.x].y);
			}
			else {

				sum += VType(0.0, VEC<VType>::quantity[idx - VEC<VType>::n.x].z, -VEC<VType>::quantity[idx - VEC<VType>::n.x].y);
			}
		}
	}

	return sum;
}

//same as zanisotropic_ngbr_sum but sum normalised values only; for scalar values this is a trivial operation, but for vectors it's not.
template <typename VType>
VType VEC_VC<VType>::zanisotropic_ngbr_dirsum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		sum = VType(-VEC<VType>::quantity[idx + 1].z, 0.0, VEC<VType>::quantity[idx + 1].x) / GetMagnitude(VEC<VType>::quantity[idx + 1]) - VType(-VEC<VType>::quantity[idx - 1].z, 0.0, VEC<VType>::quantity[idx - 1].x) / GetMagnitude(VEC<VType>::quantity[idx - 1]);
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-VEC<VType>::quantity[idx + 1].z, 0.0, VEC<VType>::quantity[idx + 1].x) / GetMagnitude(VEC<VType>::quantity[idx + 1]) - get_sign(pbc_x) * VType(-VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].z, 0.0, VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].x) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x - 1]);
			}
			else {

				sum = VType(-VEC<VType>::quantity[idx + 1].z, 0.0, VEC<VType>::quantity[idx + 1].x) / GetMagnitude(VEC<VType>::quantity[idx + 1]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = get_sign(pbc_x) * VType(-VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].z, 0.0, VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].x) / GetMagnitude(VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)]) - VType(-VEC<VType>::quantity[idx - 1].z, 0.0, VEC<VType>::quantity[idx - 1].x) / GetMagnitude(VEC<VType>::quantity[idx - 1]);
			}
			else {

				sum = VType(VEC<VType>::quantity[idx - 1].z, 0.0, -VEC<VType>::quantity[idx - 1].x) / GetMagnitude(VEC<VType>::quantity[idx - 1]);
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x].z, VEC<VType>::quantity[idx + VEC<VType>::n.x].y) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x]) - VType(0.0, -VEC<VType>::quantity[idx - VEC<VType>::n.x].z, VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x]);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x].z, VEC<VType>::quantity[idx + VEC<VType>::n.x].y) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x]) - get_sign(pbc_y) * VType(0.0, -VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].z, VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y) / GetMagnitude(VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x]);
			}
			else {

				sum += VType(0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x].z, VEC<VType>::quantity[idx + VEC<VType>::n.x].y) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += get_sign(pbc_y) * VType(0.0, -VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].z, VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y) / GetMagnitude(VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x]) - VType(0.0, -VEC<VType>::quantity[idx - VEC<VType>::n.x].z, VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x]);
			}
			else {

				sum += VType(0.0, VEC<VType>::quantity[idx - VEC<VType>::n.x].z, -VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x]);
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
VType VEC_VC<VType>::yanisotropic_ngbr_sum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		sum = VType(-VEC<VType>::quantity[idx + 1].y, VEC<VType>::quantity[idx + 1].x, 0.0) - VType(-VEC<VType>::quantity[idx - 1].y, VEC<VType>::quantity[idx - 1].x, 0.0);
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-VEC<VType>::quantity[idx + 1].y, VEC<VType>::quantity[idx + 1].x, 0.0) - get_sign(pbc_x) * VType(-VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].y, VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].x, 0.0);
			}
			else {

				sum = VType(-VEC<VType>::quantity[idx + 1].y, VEC<VType>::quantity[idx + 1].x, 0.0);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = get_sign(pbc_x) * VType(-VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].y, VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].x, 0.0) - VType(-VEC<VType>::quantity[idx - 1].y, VEC<VType>::quantity[idx - 1].x, 0.0);
			}
			else {

				sum = VType(VEC<VType>::quantity[idx - 1].y, -VEC<VType>::quantity[idx - 1].x, 0.0);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += VType(0.0, VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z, -VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y) 
			- VType(0.0, VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z, -VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(0.0, VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z, -VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y) 
					- get_sign(pbc_z) * VType(0.0, VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].z, -VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].y);
			}
			else {

				sum += VType(0.0, VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z, -VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += get_sign(pbc_z) * VType(0.0, VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].z, -VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].y) 
					- VType(0.0, VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z, -VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y);
			}
			else {

				sum += VType(0.0, -VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z, VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y);
			}
		}
	}

	return sum;
}

//same as zanisotropic_ngbr_sum but sum normalised values only; for scalar values this is a trivial operation, but for vectors it's not.
template <typename VType>
VType VEC_VC<VType>::yanisotropic_ngbr_dirsum(int idx) const
{
	VType sum = VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		sum = VType(-VEC<VType>::quantity[idx + 1].y, VEC<VType>::quantity[idx + 1].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx + 1]) - VType(-VEC<VType>::quantity[idx - 1].y, VEC<VType>::quantity[idx - 1].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx - 1]);
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-VEC<VType>::quantity[idx + 1].y, VEC<VType>::quantity[idx + 1].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx + 1]) 
					- get_sign(pbc_x) * VType(-VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].y, VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x - 1]);
			}
			else {

				sum = VType(-VEC<VType>::quantity[idx + 1].y, VEC<VType>::quantity[idx + 1].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx + 1]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = get_sign(pbc_x) * VType(-VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].y, VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)]) 
					- VType(-VEC<VType>::quantity[idx - 1].y, VEC<VType>::quantity[idx - 1].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx - 1]);
			}
			else {

				sum = VType(VEC<VType>::quantity[idx - 1].y, -VEC<VType>::quantity[idx - 1].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx - 1]);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += VType(0.0, VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z, -VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y]) 
			- VType(0.0, VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z, -VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(0.0, VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z, -VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y]) 
					- get_sign(pbc_z) * VType(0.0, VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].z, -VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].y) / GetMagnitude(VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y]);
			}
			else {

				sum += VType(0.0, VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z, -VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += get_sign(pbc_z) * VType(0.0, VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].z, -VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].y) / GetMagnitude(VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y]) 
					- VType(0.0, VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z, -VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
			}
			else {

				sum += VType(0.0, -VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z, VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
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
VType VEC_VC<VType>::xanisotropic_ngbr_sum(int idx) const
{
	VType sum = VType();

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x].y, -VEC<VType>::quantity[idx + VEC<VType>::n.x].x, 0.0) - VType(VEC<VType>::quantity[idx - VEC<VType>::n.x].y, -VEC<VType>::quantity[idx - VEC<VType>::n.x].x, 0.0);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x].y, -VEC<VType>::quantity[idx + VEC<VType>::n.x].x, 0.0) - get_sign(pbc_y) * VType(VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y, -VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].x, 0.0);
			}
			else {

				sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x].y, -VEC<VType>::quantity[idx + VEC<VType>::n.x].x, 0.0);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += get_sign(pbc_y) * VType(VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y, -VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].x, 0.0) - VType(VEC<VType>::quantity[idx - VEC<VType>::n.x].y, -VEC<VType>::quantity[idx - VEC<VType>::n.x].x, 0.0);
			}
			else {

				sum += VType(-VEC<VType>::quantity[idx - VEC<VType>::n.x].y, VEC<VType>::quantity[idx - VEC<VType>::n.x].x, 0.0);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x) 
			- VType(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, -VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x) 
					- get_sign(pbc_z) * VType(VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, -VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].x);
			}
			else {

				sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += get_sign(pbc_z) * VType(VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, -VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].x) 
					- VType(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, -VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x);
			}
			else {

				sum += VType(-VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x);
			}
		}
	}

	return sum;
}

//same as zanisotropic_ngbr_sum but sum normalised values only; for scalar values this is a trivial operation, but for vectors it's not.
template <typename VType>
VType VEC_VC<VType>::xanisotropic_ngbr_dirsum(int idx) const
{
	VType sum = VType();

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x].y, -VEC<VType>::quantity[idx + VEC<VType>::n.x].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x])
			- VType(VEC<VType>::quantity[idx - VEC<VType>::n.x].y, -VEC<VType>::quantity[idx - VEC<VType>::n.x].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x]);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x].y, -VEC<VType>::quantity[idx + VEC<VType>::n.x].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x]) 
					- get_sign(pbc_y) * VType(VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y, -VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x]);
			}
			else {

				sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x].y, -VEC<VType>::quantity[idx + VEC<VType>::n.x].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += get_sign(pbc_y) * VType(VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y, -VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x])
					- VType(VEC<VType>::quantity[idx - VEC<VType>::n.x].y, -VEC<VType>::quantity[idx - VEC<VType>::n.x].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x]);
			}
			else {

				sum += VType(-VEC<VType>::quantity[idx - VEC<VType>::n.x].y, VEC<VType>::quantity[idx - VEC<VType>::n.x].x, 0.0) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x]);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y])
			- VType(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, -VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y])
					- get_sign(pbc_z) * VType(VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, -VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].x) / GetMagnitude(VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y]);
			}
			else {

				sum += VType(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, -VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x) / GetMagnitude(VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += get_sign(pbc_z) * VType(VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, -VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y].x) / GetMagnitude(VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1)*VEC<VType>::n.x*VEC<VType>::n.y])
					- VType(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, -VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
			}
			else {

				sum += VType(-VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z, 0.0, VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / GetMagnitude(VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
			}
		}
	}

	return sum;
}