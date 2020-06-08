#pragma once

#include "VEC_VC.h"

//-------------------------------- NEIGHBOR SUM

//calculate 6-point neighbor sum at given index
//missing neighbors not added, including at boundaries, but taking into account pbc
template <typename VType>
VType VEC_VC<VType>::ngbr_sum(int idx) const
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