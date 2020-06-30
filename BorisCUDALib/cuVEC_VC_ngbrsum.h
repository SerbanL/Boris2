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
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1]);
			}
			else {

				sum = cuVEC<VType>::quantity[idx + 1];
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)] + cuVEC<VType>::quantity[idx - 1]);
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
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]);
			}
			else {

				sum += cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x];
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
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
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
			else {

				sum += cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y];
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
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
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = (cuVEC<VType>::quantity[idx + 1] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]) + cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1]));
			}
			else {

				sum = cuVEC<VType>::quantity[idx + 1] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)]) + cuVEC<VType>::quantity[idx - 1] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]));
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
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]) + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]));
			}
			else {

				sum += cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]) + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]));
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
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]) + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y]));
			}
			else {

				sum += cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y]) + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]));
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
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(0.0, -cuVEC<VType>::quantity[idx + 1].z, cuVEC<VType>::quantity[idx + 1].y) - VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].y);
			}
			else {

				sum = VType(0.0, -cuVEC<VType>::quantity[idx + 1].z, cuVEC<VType>::quantity[idx + 1].y);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(0.0, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].z, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].y) - VType(0.0, -cuVEC<VType>::quantity[idx - 1].z, cuVEC<VType>::quantity[idx - 1].y);
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
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x) - VType(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].x);
			}
			else {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].x) - VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x);
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
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) - VType(-cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0);
			}
			else {

				sum += VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(-cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) - VType(-cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0);
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
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(0.0, -cuVEC<VType>::quantity[idx + 1].z, cuVEC<VType>::quantity[idx + 1].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]) - VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1]);
			}
			else {

				sum = VType(0.0, -cuVEC<VType>::quantity[idx + 1].z, cuVEC<VType>::quantity[idx + 1].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(0.0, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].z, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)]) - VType(0.0, -cuVEC<VType>::quantity[idx - 1].z, cuVEC<VType>::quantity[idx - 1].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]);
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
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]) - VType(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]);
			}
			else {

				sum += VType(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]) - VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, 0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
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
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]) - VType(-cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
			else {

				sum += VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCZ) {

				sum += VType(-cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y]) - VType(-cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
			else {

				sum += VType(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].y, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].x, 0.0) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
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
		sum = VType(-cuVEC<VType>::quantity[idx + 1].z, 0.0, cuVEC<VType>::quantity[idx + 1].x) - VType(-cuVEC<VType>::quantity[idx - 1].z, 0.0, cuVEC<VType>::quantity[idx - 1].x);
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-cuVEC<VType>::quantity[idx + 1].z, 0.0, cuVEC<VType>::quantity[idx + 1].x) - VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].z, 0.0, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].x);
			}
			else {

				sum = VType(-cuVEC<VType>::quantity[idx + 1].z, 0.0, cuVEC<VType>::quantity[idx + 1].x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].z, 0.0, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].x) - VType(-cuVEC<VType>::quantity[idx - 1].z, 0.0, cuVEC<VType>::quantity[idx - 1].x);
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
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y) - VType(0.0, -cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y);
			}
			else {

				sum += VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(0.0, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y) - VType(0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y);
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
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-cuVEC<VType>::quantity[idx + 1].z, 0.0, cuVEC<VType>::quantity[idx + 1].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]) - VType(-cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].z, 0.0, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1]);
			}
			else {

				sum = VType(-cuVEC<VType>::quantity[idx + 1].z, 0.0, cuVEC<VType>::quantity[idx + 1].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + 1]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCX) {

				sum = VType(-cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].z, 0.0, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)]) - VType(-cuVEC<VType>::quantity[idx - 1].z, 0.0, cuVEC<VType>::quantity[idx - 1].x) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - 1]);
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
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]) - VType(0.0, -cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]);
			}
			else {

				sum += VType(0.0, -cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise just one contribution.
			if (ngbrFlags[idx] & NF_PBCY) {

				sum += VType(0.0, -cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]) - VType(0.0, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			}
			else {

				sum += VType(0.0, cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].z, -cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cu_GetMagnitude(cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			}
		}
	}

	return sum;
}
