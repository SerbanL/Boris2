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
		sum = (quantity[idx + 1] + quantity[idx - 1]);
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = (quantity[idx + 1] + quantity[idx + n.x - 1]);
			}
			else {

				sum = quantity[idx + 1];
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = (quantity[idx - (n.x - 1)] + quantity[idx - 1]);
			}
			else {

				sum = quantity[idx - 1];
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += (quantity[idx + n.x] + quantity[idx - n.x]);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (quantity[idx + n.x] + quantity[idx + (n.y - 1)*n.x]);
			}
			else {

				sum += quantity[idx + n.x];
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (quantity[idx - (n.y - 1)*n.x] + quantity[idx - n.x]);
			}
			else {

				sum += quantity[idx - n.x];
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += (quantity[idx + n.x*n.y] + quantity[idx - n.x*n.y]);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (quantity[idx + n.x*n.y] + quantity[idx + (n.z - 1) * n.x*n.y]);
			}
			else {

				sum += quantity[idx + n.x*n.y];
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (quantity[idx - (n.z - 1) * n.x*n.y] + quantity[idx - n.x*n.y]);
			}
			else {

				sum += quantity[idx - n.x*n.y];
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
		sum = (quantity[idx + 1] / cu_GetMagnitude(quantity[idx + 1]) + quantity[idx - 1] / cu_GetMagnitude(quantity[idx - 1]));
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = (quantity[idx + 1] / cu_GetMagnitude(quantity[idx + 1]) + quantity[idx + n.x - 1] / cu_GetMagnitude(quantity[idx + n.x - 1]));
			}
			else {

				sum = quantity[idx + 1] / cu_GetMagnitude(quantity[idx + 1]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = (quantity[idx - (n.x - 1)] / cu_GetMagnitude(quantity[idx - (n.x - 1)]) + quantity[idx - 1] / cu_GetMagnitude(quantity[idx - 1]));
			}
			else {

				sum = quantity[idx - 1] / cu_GetMagnitude(quantity[idx - 1]);
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += (quantity[idx + n.x] / cu_GetMagnitude(quantity[idx + n.x]) + quantity[idx - n.x] / cu_GetMagnitude(quantity[idx - n.x]));
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (quantity[idx + n.x] / cu_GetMagnitude(quantity[idx + n.x]) + quantity[idx + (n.y - 1)*n.x] / cu_GetMagnitude(quantity[idx + (n.y - 1)*n.x]));
			}
			else {

				sum += quantity[idx + n.x] / cu_GetMagnitude(quantity[idx + n.x]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (quantity[idx - (n.y - 1)*n.x] / cu_GetMagnitude(quantity[idx - (n.y - 1)*n.x]) + quantity[idx - n.x] / cu_GetMagnitude(quantity[idx - n.x]));
			}
			else {

				sum += quantity[idx - n.x] / cu_GetMagnitude(quantity[idx - n.x]);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += (quantity[idx + n.x*n.y] / cu_GetMagnitude(quantity[idx + n.x*n.y]) + quantity[idx - n.x*n.y] / cu_GetMagnitude(quantity[idx - n.x*n.y]));
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (quantity[idx + n.x*n.y] / cu_GetMagnitude(quantity[idx + n.x*n.y]) + quantity[idx + (n.z - 1) * n.x*n.y] / cu_GetMagnitude(quantity[idx + (n.z - 1) * n.x*n.y]));
			}
			else {

				sum += quantity[idx + n.x*n.y] / cu_GetMagnitude(quantity[idx + n.x*n.y]);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (quantity[idx - (n.z - 1) * n.x*n.y] / cu_GetMagnitude(quantity[idx - (n.z - 1) * n.x*n.y]) + quantity[idx - n.x*n.y] / cu_GetMagnitude(quantity[idx - n.x*n.y]));
			}
			else {

				sum += quantity[idx - n.x*n.y] / cu_GetMagnitude(quantity[idx - n.x*n.y]);
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
		sum = VType(0.0, -quantity[idx + 1].z, quantity[idx + 1].y) - VType(0.0, -quantity[idx - 1].z, quantity[idx - 1].y);
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(0.0, -quantity[idx + 1].z, quantity[idx + 1].y) - VType(0.0, -quantity[idx + n.x - 1].z, quantity[idx + n.x - 1].y);
			}
			else {

				sum = VType(0.0, -quantity[idx + 1].z, quantity[idx + 1].y);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(0.0, -quantity[idx - (n.x - 1)].z, quantity[idx - (n.x - 1)].y) - VType(0.0, -quantity[idx - 1].z, quantity[idx - 1].y);
			}
			else {

				sum = VType(0.0, quantity[idx - 1].z, -quantity[idx - 1].y);
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(quantity[idx + n.x].z, 0.0, -quantity[idx + n.x].x) - VType(quantity[idx - n.x].z, 0.0, -quantity[idx - n.x].x);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(quantity[idx + n.x].z, 0.0, -quantity[idx + n.x].x) - VType(quantity[idx + (n.y - 1)*n.x].z, 0.0, -quantity[idx + (n.y - 1)*n.x].x);
			}
			else {

				sum += VType(quantity[idx + n.x].z, 0.0, -quantity[idx + n.x].x);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(quantity[idx - (n.y - 1)*n.x].z, 0.0, -quantity[idx - (n.y - 1)*n.x].x) - VType(quantity[idx - n.x].z, 0.0, -quantity[idx - n.x].x);
			}
			else {

				sum += VType(-quantity[idx - n.x].z, 0.0, +quantity[idx - n.x].x);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += VType(-quantity[idx + n.x*n.y].y, quantity[idx + n.x*n.y].x, 0.0) - VType(-quantity[idx - n.x*n.y].y, quantity[idx - n.x*n.y].x, 0.0);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(-quantity[idx + n.x*n.y].y, quantity[idx + n.x*n.y].x, 0.0) - VType(-quantity[idx + (n.z - 1) * n.x*n.y].y, quantity[idx + (n.z - 1) * n.x*n.y].x, 0.0);
			}
			else {

				sum += VType(-quantity[idx + n.x*n.y].y, quantity[idx + n.x*n.y].x, 0.0);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(-quantity[idx - (n.z - 1) * n.x*n.y].y, quantity[idx - (n.z - 1) * n.x*n.y].x, 0.0) - VType(-quantity[idx - n.x*n.y].y, quantity[idx - n.x*n.y].x, 0.0);
			}
			else {

				sum += VType(quantity[idx - n.x*n.y].y, -quantity[idx - n.x*n.y].x, 0.0);
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
		sum = VType(0.0, -quantity[idx + 1].z, quantity[idx + 1].y) / cu_GetMagnitude(quantity[idx + 1]) - VType(0.0, -quantity[idx - 1].z, quantity[idx - 1].y) / cu_GetMagnitude(quantity[idx - 1]);
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(0.0, -quantity[idx + 1].z, quantity[idx + 1].y) / cu_GetMagnitude(quantity[idx + 1]) - VType(0.0, -quantity[idx + n.x - 1].z, quantity[idx + n.x - 1].y) / cu_GetMagnitude(quantity[idx + n.x - 1]);
			}
			else {

				sum = VType(0.0, -quantity[idx + 1].z, quantity[idx + 1].y) / cu_GetMagnitude(quantity[idx + 1]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(0.0, -quantity[idx - (n.x - 1)].z, quantity[idx - (n.x - 1)].y) / cu_GetMagnitude(quantity[idx - (n.x - 1)]) - VType(0.0, -quantity[idx - 1].z, quantity[idx - 1].y) / cu_GetMagnitude(quantity[idx - 1]);
			}
			else {

				sum = VType(0.0, quantity[idx - 1].z, -quantity[idx - 1].y) / cu_GetMagnitude(quantity[idx - 1]);
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(quantity[idx + n.x].z, 0.0, -quantity[idx + n.x].x) / cu_GetMagnitude(quantity[idx + n.x]) - VType(quantity[idx - n.x].z, 0.0, -quantity[idx - n.x].x) / cu_GetMagnitude(quantity[idx - n.x]);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(quantity[idx + n.x].z, 0.0, -quantity[idx + n.x].x) / cu_GetMagnitude(quantity[idx + n.x]) - VType(quantity[idx + (n.y - 1)*n.x].z, 0.0, -quantity[idx + (n.y - 1)*n.x].x) / cu_GetMagnitude(quantity[idx + (n.y - 1)*n.x]);
			}
			else {

				sum += VType(quantity[idx + n.x].z, 0.0, -quantity[idx + n.x].x) / cu_GetMagnitude(quantity[idx + n.x]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(quantity[idx - (n.y - 1)*n.x].z, 0.0, -quantity[idx - (n.y - 1)*n.x].x) / cu_GetMagnitude(quantity[idx - (n.y - 1)*n.x]) - VType(quantity[idx - n.x].z, 0.0, -quantity[idx - n.x].x) / cu_GetMagnitude(quantity[idx - n.x]);
			}
			else {

				sum += VType(-quantity[idx - n.x].z, 0.0, +quantity[idx - n.x].x) / cu_GetMagnitude(quantity[idx - n.x]);
			}
		}
	}

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		sum += VType(-quantity[idx + n.x*n.y].y, quantity[idx + n.x*n.y].x, 0.0) / cu_GetMagnitude(quantity[idx + n.x*n.y]) - VType(-quantity[idx - n.x*n.y].y, quantity[idx - n.x*n.y].x, 0.0) / cu_GetMagnitude(quantity[idx - n.x*n.y]);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(-quantity[idx + n.x*n.y].y, quantity[idx + n.x*n.y].x, 0.0) / cu_GetMagnitude(quantity[idx + n.x*n.y]) - VType(-quantity[idx + (n.z - 1) * n.x*n.y].y, quantity[idx + (n.z - 1) * n.x*n.y].x, 0.0) / cu_GetMagnitude(quantity[idx + (n.z - 1) * n.x*n.y]);
			}
			else {

				sum += VType(-quantity[idx + n.x*n.y].y, quantity[idx + n.x*n.y].x, 0.0) / cu_GetMagnitude(quantity[idx + n.x*n.y]);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(-quantity[idx - (n.z - 1) * n.x*n.y].y, quantity[idx - (n.z - 1) * n.x*n.y].x, 0.0) / cu_GetMagnitude(quantity[idx - (n.z - 1) * n.x*n.y]) - VType(-quantity[idx - n.x*n.y].y, quantity[idx - n.x*n.y].x, 0.0) / cu_GetMagnitude(quantity[idx - n.x*n.y]);
			}
			else {

				sum += VType(quantity[idx - n.x*n.y].y, -quantity[idx - n.x*n.y].x, 0.0) / cu_GetMagnitude(quantity[idx - n.x*n.y]);
			}
		}
	}

	return sum;
}

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
		sum = VType(-quantity[idx + 1].z, 0.0, quantity[idx + 1].x) - VType(-quantity[idx - 1].z, 0.0, quantity[idx - 1].x);
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-quantity[idx + 1].z, 0.0, quantity[idx + 1].x) - VType(-quantity[idx + n.x - 1].z, 0.0, quantity[idx + n.x - 1].x);
			}
			else {

				sum = VType(-quantity[idx + 1].z, 0.0, quantity[idx + 1].x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-quantity[idx - (n.x - 1)].z, 0.0, quantity[idx - (n.x - 1)].x) - VType(-quantity[idx - 1].z, 0.0, quantity[idx - 1].x);
			}
			else {

				sum = VType(quantity[idx - 1].z, 0.0, -quantity[idx - 1].x);
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(0.0, -quantity[idx + n.x].z, quantity[idx + n.x].y) - VType(0.0, -quantity[idx - n.x].z, quantity[idx - n.x].y);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(0.0, -quantity[idx + n.x].z, quantity[idx + n.x].y) - VType(0.0, -quantity[idx + (n.y - 1)*n.x].z, quantity[idx + (n.y - 1)*n.x].y);
			}
			else {

				sum += VType(0.0, -quantity[idx + n.x].z, quantity[idx + n.x].y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(0.0, -quantity[idx - (n.y - 1)*n.x].z, quantity[idx - (n.y - 1)*n.x].y) - VType(0.0, -quantity[idx - n.x].z, quantity[idx - n.x].y);
			}
			else {

				sum += VType(0.0, quantity[idx - n.x].z, -quantity[idx - n.x].y);
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
		sum = VType(-quantity[idx + 1].z, 0.0, quantity[idx + 1].x) / cu_GetMagnitude(quantity[idx + 1]) - VType(-quantity[idx - 1].z, 0.0, quantity[idx - 1].x) / cu_GetMagnitude(quantity[idx - 1]);
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-quantity[idx + 1].z, 0.0, quantity[idx + 1].x) / cu_GetMagnitude(quantity[idx + 1]) - VType(-quantity[idx + n.x - 1].z, 0.0, quantity[idx + n.x - 1].x) / cu_GetMagnitude(quantity[idx + n.x - 1]);
			}
			else {

				sum = VType(-quantity[idx + 1].z, 0.0, quantity[idx + 1].x) / cu_GetMagnitude(quantity[idx + 1]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-quantity[idx - (n.x - 1)].z, 0.0, quantity[idx - (n.x - 1)].x) / cu_GetMagnitude(quantity[idx - (n.x - 1)]) - VType(-quantity[idx - 1].z, 0.0, quantity[idx - 1].x) / cu_GetMagnitude(quantity[idx - 1]);
			}
			else {

				sum = VType(quantity[idx - 1].z, 0.0, -quantity[idx - 1].x) / cu_GetMagnitude(quantity[idx - 1]);
			}
		}
	}

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		sum += VType(0.0, -quantity[idx + n.x].z, quantity[idx + n.x].y) / cu_GetMagnitude(quantity[idx + n.x]) - VType(0.0, -quantity[idx - n.x].z, quantity[idx - n.x].y) / cu_GetMagnitude(quantity[idx - n.x]);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(0.0, -quantity[idx + n.x].z, quantity[idx + n.x].y) / cu_GetMagnitude(quantity[idx + n.x]) - VType(0.0, -quantity[idx + (n.y - 1)*n.x].z, quantity[idx + (n.y - 1)*n.x].y) / cu_GetMagnitude(quantity[idx + (n.y - 1)*n.x]);
			}
			else {

				sum += VType(0.0, -quantity[idx + n.x].z, quantity[idx + n.x].y) / cu_GetMagnitude(quantity[idx + n.x]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(0.0, -quantity[idx - (n.y - 1)*n.x].z, quantity[idx - (n.y - 1)*n.x].y) / cu_GetMagnitude(quantity[idx - (n.y - 1)*n.x]) - VType(0.0, -quantity[idx - n.x].z, quantity[idx - n.x].y) / cu_GetMagnitude(quantity[idx - n.x]);
			}
			else {

				sum += VType(0.0, quantity[idx - n.x].z, -quantity[idx - n.x].y) / cu_GetMagnitude(quantity[idx - n.x]);
			}
		}
	}

	return sum;
}
