#pragma once

#include "cuVEC_VC.h"

//---- SECOND ORDER DIFFERENTIALS

//HOMOGENEOUS NEUMANN

//homogeneous second order.
//Use Neumann boundary conditions.
//Returns zero at composite media boundary cells
template <typename VType>
__device__ VType cuVEC_VC<VType>::dxx_neu(int idx) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.x*cuVEC<VType>::h.x;

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		diff = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)] + cuVEC<VType>::quantity[idx - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx - 1] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
__device__ VType cuVEC_VC<VType>::dyy_neu(int idx) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.y*cuVEC<VType>::h.y;

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
__device__ VType cuVEC_VC<VType>::dzz_neu(int idx) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.z*cuVEC<VType>::h.z;

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}

//NON-HOMOGENEOUS NEUMANN

//Use non-homogeneous Neumann boundary conditions.
//Returns zero at composite media boundary cells
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
template <typename VType>
template <typename Class_BDiff>
__device__ VType cuVEC_VC<VType>::dxx_nneu(int idx, const Class_BDiff& bdiff_class) const
{
	VType diff;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.x*cuVEC<VType>::h.x;

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		diff = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx]) - bdiff_class.bdiff(idx).x * cuVEC<VType>::h.x) / hsq;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)] + cuVEC<VType>::quantity[idx - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx - 1] - cuVEC<VType>::quantity[idx]) + bdiff_class.bdiff(idx).x * cuVEC<VType>::h.x) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
template <typename Class_BDiff>
__device__ VType cuVEC_VC<VType>::dyy_nneu(int idx, const Class_BDiff& bdiff_class) const
{
	VType diff;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.y*cuVEC<VType>::h.y;

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) - bdiff_class.bdiff(idx).y * cuVEC<VType>::h.y) / hsq;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) + bdiff_class.bdiff(idx).y * cuVEC<VType>::h.y) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
template <typename Class_BDiff>
__device__ VType cuVEC_VC<VType>::dzz_nneu(int idx, const Class_BDiff& bdiff_class) const
{
	VType diff;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.z*cuVEC<VType>::h.z;

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) - bdiff_class.bdiff(idx).z * cuVEC<VType>::h.z) / hsq;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) + bdiff_class.bdiff(idx).z * cuVEC<VType>::h.z) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
__device__ VType cuVEC_VC<VType>::dxx_nneu(int idx, const cuVAL3<VType>& bdiff) const
{
	VType diff;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.x*cuVEC<VType>::h.x;

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		diff = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx]) - bdiff.x * cuVEC<VType>::h.x) / hsq;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)] + cuVEC<VType>::quantity[idx - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx - 1] - cuVEC<VType>::quantity[idx]) + bdiff.x * cuVEC<VType>::h.x) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
__device__ VType cuVEC_VC<VType>::dyy_nneu(int idx, const cuVAL3<VType>& bdiff) const
{
	VType diff;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.y*cuVEC<VType>::h.y;

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) - bdiff.y * cuVEC<VType>::h.y) / hsq;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) + bdiff.y * cuVEC<VType>::h.y) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
__device__ VType cuVEC_VC<VType>::dzz_nneu(int idx, const cuVAL3<VType>& bdiff) const
{
	VType diff;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.z*cuVEC<VType>::h.z;

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) - bdiff.z * cuVEC<VType>::h.z) / hsq;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) + bdiff.z * cuVEC<VType>::h.z) / hsq;
			}
		}
	}

	return diff;
}

//DIRICHLET, HOMOGENEOUS NEUMANN

//Use Dirichlet boundary conditions, else Neumann boundary conditions (homogeneous).
//Returns zero at composite media boundary cells.
template <typename VType>
__device__ VType cuVEC_VC<VType>::dxx_diri(int idx) const
{
	VType diff;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.x*cuVEC<VType>::h.x;

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		diff = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		//only one neighbor available. Does it use a fixed boundary value (Dirichlet)?
		if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPX) diff = (2 * cuVEC<VType>::quantity[idx + 1] + 4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
			else								   diff = (2 * cuVEC<VType>::quantity[idx - 1] + 4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
		}
		//stencil at surface for homogeneous Neumann boundary condition - correct, do not multiply by 2! 
		//See conservation law approach derivation for heat equation (integrate then use mid-point rule approximations) for cell-centered grids.
		else if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)] + cuVEC<VType>::quantity[idx - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx - 1] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
__device__ VType cuVEC_VC<VType>::dyy_diri(int idx) const
{
	VType diff;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.y*cuVEC<VType>::h.y;

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPY) diff = (2 * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + 4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
			else								   diff = (2 * cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] + 4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
		}
		else if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
__device__ VType cuVEC_VC<VType>::dzz_diri(int idx) const
{
	VType diff;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.z*cuVEC<VType>::h.z;

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) diff = (2 * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + 4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
			else								 diff = (2 * cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] + 4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
		}
		else if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}

//DIRICHLET, NON-HOMOGENEOUS NEUMANN

//Use Dirichlet boundary conditions, else non-homogeneous Neumann boundary conditions.
//Returns zero at composite media boundary cells.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
template <typename VType>
template <typename Class_BDiff>
__device__ VType cuVEC_VC<VType>::dxx_diri_nneu(int idx, const Class_BDiff& bdiff_class) const
{
	VType diff;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.x*cuVEC<VType>::h.x;

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		diff = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		//only one neighbor available. Does it use a fixed boundary value (Dirichlet)?
		if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPX) diff = (2 * cuVEC<VType>::quantity[idx + 1] + 4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
			else								 diff = (2 * cuVEC<VType>::quantity[idx - 1] + 4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
		}
		else if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx]) - bdiff_class.bdiff(idx).x * cuVEC<VType>::h.x) / hsq;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)] + cuVEC<VType>::quantity[idx - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx - 1] - cuVEC<VType>::quantity[idx]) + bdiff_class.bdiff(idx).x * cuVEC<VType>::h.x) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
template <typename Class_BDiff>
__device__ VType cuVEC_VC<VType>::dyy_diri_nneu(int idx, const Class_BDiff& bdiff_class) const
{
	VType diff;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.y*cuVEC<VType>::h.y;

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPY) diff = (2 * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + 4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
			else								 diff = (2 * cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] + 4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
		}
		else if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) - bdiff_class.bdiff(idx).y * cuVEC<VType>::h.y) / hsq;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) + bdiff_class.bdiff(idx).y * cuVEC<VType>::h.y) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
template <typename Class_BDiff>
__device__ VType cuVEC_VC<VType>::dzz_diri_nneu(int idx, const Class_BDiff& bdiff_class) const
{
	VType diff;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.z*cuVEC<VType>::h.z;

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) diff = (2 * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + 4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
			else								 diff = (2 * cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] + 4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
		}
		else if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) - bdiff_class.bdiff(idx).z * cuVEC<VType>::h.z) / hsq;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) + bdiff_class.bdiff(idx).z * cuVEC<VType>::h.z) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
__device__ VType cuVEC_VC<VType>::dxx_diri_nneu(int idx, const cuVAL3<VType>& bdiff) const
{
	VType diff;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.x*cuVEC<VType>::h.x;

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		diff = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		//only one neighbor available. Does it use a fixed boundary value (Dirichlet)?
		if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPX) diff = (2 * cuVEC<VType>::quantity[idx + 1] + 4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
			else								 diff = (2 * cuVEC<VType>::quantity[idx - 1] + 4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
		}
		else if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx]) - bdiff.x * cuVEC<VType>::h.x) / hsq;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)] + cuVEC<VType>::quantity[idx - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx - 1] - cuVEC<VType>::quantity[idx]) + bdiff.x * cuVEC<VType>::h.x) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
__device__ VType cuVEC_VC<VType>::dyy_diri_nneu(int idx, const cuVAL3<VType>& bdiff) const
{
	VType diff;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.y*cuVEC<VType>::h.y;

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPY) diff = (2 * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + 4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
			else								 diff = (2 * cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] + 4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
		}
		else if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) - bdiff.y * cuVEC<VType>::h.y) / hsq;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) + bdiff.y * cuVEC<VType>::h.y) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
__device__ VType cuVEC_VC<VType>::dzz_diri_nneu(int idx, const cuVAL3<VType>& bdiff) const
{
	VType diff;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.z*cuVEC<VType>::h.z;

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) diff = (2 * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + 4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
			else								 diff = (2 * cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] + 4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) - 6 * cuVEC<VType>::quantity[idx]) / hsq;
		}
		else if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) - bdiff.z * cuVEC<VType>::h.z) / hsq;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) + bdiff.z * cuVEC<VType>::h.z) / hsq;
			}
		}
	}

	return diff;
}

//ROBIN BOUNDARY CONDITIONS

//Use Robin boundary conditions (defaulting to Neumann if not set).
//Returns zero at composite media boundary cells.
//The K constant is used in Robin boundary condition calculations, where -K*diff_norm(T) = alpha*(Tboundary - Tambient) is the flux normal to the boundary - K is the thermal conductivity in the heat equation
template <typename VType>
__device__ VType cuVEC_VC<VType>::dxx_robin(int idx, cuBReal K) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.x*cuVEC<VType>::h.x;

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {
		//both neighbors available so must be an inner point along this direction
		diff = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		//only one neighbor available. Does it use a Robin condition?
		if (using_extended_flags && (ngbrFlags2[idx] & NF2_ROBINX)) {

			//Robin cell on +x side of boundary
			if (ngbrFlags2[idx] & NF2_ROBINPX) {

				//use void Robin values
				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff = ((1 + robin_v.i * cuVEC<VType>::h.x / (2 * K)) * cuVEC<VType>::quantity[idx + 1] + robin_v.i * cuVEC<VType>::h.x * robin_v.j / K - (1 + 3 * robin_v.i * cuVEC<VType>::h.x / (2 * K)) * cuVEC<VType>::quantity[idx]) / hsq;
				}
				//don't use void Robin values
				else {

					diff = ((1 + robin_nx.i * cuVEC<VType>::h.x / (2 * K)) * cuVEC<VType>::quantity[idx + 1] + robin_nx.i * cuVEC<VType>::h.x * robin_nx.j / K - (1 + 3 * robin_nx.i * cuVEC<VType>::h.x / (2 * K)) * cuVEC<VType>::quantity[idx]) / hsq;
				}
			}
			//Robin cell on -x side of boundary
			else {

				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff = ((1 + robin_v.i * cuVEC<VType>::h.x / (2 * K)) * cuVEC<VType>::quantity[idx - 1] + robin_v.i * cuVEC<VType>::h.x * robin_v.j / K - (1 + 3 * robin_v.i * cuVEC<VType>::h.x / (2 * K)) * cuVEC<VType>::quantity[idx]) / hsq;
				}
				else {

					diff = ((1 + robin_px.i * cuVEC<VType>::h.x / (2 * K)) * cuVEC<VType>::quantity[idx - 1] + robin_px.i * cuVEC<VType>::h.x * robin_px.j / K - (1 + 3 * robin_px.i * cuVEC<VType>::h.x / (2 * K)) * cuVEC<VType>::quantity[idx]) / hsq;
				}
			}
		}
		//stencil at surface for homogeneous Neumann boundary condition - correct, do not multiply by 2! 
		//See conservation law approach derivation for heat equation (integrate then use mid-point rule approximations) for cell-centered grids.
		else if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)] + cuVEC<VType>::quantity[idx - 1] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx - 1] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
__device__ VType cuVEC_VC<VType>::dyy_robin(int idx, cuBReal K) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.y*cuVEC<VType>::h.y;

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		//only one neighbor available. Does it use a Robin condition?
		if (using_extended_flags && (ngbrFlags2[idx] & NF2_ROBINY)) {

			//Robin cell on +x side of boundary
			if (ngbrFlags2[idx] & NF2_ROBINPY) {

				//use void Robin values
				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff = ((1 + robin_v.i * cuVEC<VType>::h.y / (2 * K)) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + robin_v.i * cuVEC<VType>::h.y * robin_v.j / K - (1 + 3 * robin_v.i * cuVEC<VType>::h.y / (2 * K)) * cuVEC<VType>::quantity[idx]) / hsq;
				}
				//don't use void Robin values
				else {

					diff = ((1 + robin_ny.i * cuVEC<VType>::h.y / (2 * K)) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + robin_ny.i * cuVEC<VType>::h.y * robin_ny.j / K - (1 + 3 * robin_ny.i * cuVEC<VType>::h.y / (2 * K)) * cuVEC<VType>::quantity[idx]) / hsq;
				}
			}
			//Robin cell on -x side of boundary
			else {

				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff = ((1 + robin_v.i * cuVEC<VType>::h.y / (2 * K)) * cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] + robin_v.i * cuVEC<VType>::h.y * robin_v.j / K - (1 + 3 * robin_v.i * cuVEC<VType>::h.y / (2 * K)) * cuVEC<VType>::quantity[idx]) / hsq;
				}
				else {

					diff = ((1 + robin_py.i * cuVEC<VType>::h.y / (2 * K)) * cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] + robin_py.i * cuVEC<VType>::h.y * robin_py.j / K - (1 + 3 * robin_py.i * cuVEC<VType>::h.y / (2 * K)) * cuVEC<VType>::quantity[idx]) / hsq;
				}
			}
		}
		else if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
__device__ VType cuVEC_VC<VType>::dzz_robin(int idx, cuBReal K) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	cuBReal hsq = cuVEC<VType>::h.z*cuVEC<VType>::h.z;

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		//only one neighbor available. Does it use a Robin condition?
		if (using_extended_flags && (ngbrFlags2[idx] & NF2_ROBINZ)) {

			//Robin cell on +x side of boundary
			if (ngbrFlags2[idx] & NF2_ROBINPZ) {

				//use void Robin values
				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff = ((1 + robin_v.i * cuVEC<VType>::h.z / (2 * K)) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + robin_v.i * cuVEC<VType>::h.z * robin_v.j / K - (1 + 3 * robin_v.i * cuVEC<VType>::h.z / (2 * K)) * cuVEC<VType>::quantity[idx]) / hsq;
				}
				//don't use void Robin values
				else {

					diff = ((1 + robin_nz.i * cuVEC<VType>::h.z / (2 * K)) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + robin_nz.i * cuVEC<VType>::h.z * robin_nz.j / K - (1 + 3 * robin_nz.i * cuVEC<VType>::h.z / (2 * K)) * cuVEC<VType>::quantity[idx]) / hsq;
				}
			}
			//Robin cell on -x side of boundary
			else {

				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff = ((1 + robin_v.i * cuVEC<VType>::h.z / (2 * K)) * cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] + robin_v.i * cuVEC<VType>::h.z * robin_v.j / K - (1 + 3 * robin_v.i * cuVEC<VType>::h.z / (2 * K)) * cuVEC<VType>::quantity[idx]) / hsq;
				}
				else {

					diff = ((1 + robin_pz.i * cuVEC<VType>::h.z / (2 * K)) * cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] + robin_pz.i * cuVEC<VType>::h.z * robin_pz.j / K - (1 + 3 * robin_pz.i * cuVEC<VType>::h.z / (2 * K)) * cuVEC<VType>::quantity[idx]) / hsq;
				}
			}
		}
		else if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - 2 * cuVEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}

//-------------------------------- SECOND ORDER MIXED

//Use Neumann boundary conditions(homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
//dxy same as dyx
template <typename VType>
__device__ VType cuVEC_VC<VType>::dxy_neu(int idx) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	if ((ngbrFlags[idx] & NF_XY_FULL) == NF_XY_FULL) {

		//full off-axis XY stencil available
		diff = (cuVEC<VType>::quantity[idx + 1 + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx - 1 - cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx + 1 - cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx - 1 + cuVEC<VType>::n.x]);
	}
	else if (ngbrFlags[idx] & NF_XY_OASTENCIL) {

		//not full, but XY off-axis stencil is available

		//try to get right side first
		//right side (+ - order)
		if ((ngbrFlags[idx] & NF_XY_PXPY) && (ngbrFlags[idx] & NF_XY_PXNY)) {

			//o - o
			diff = (cuVEC<VType>::quantity[idx + 1 + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx + 1 - cuVEC<VType>::n.x]);
		}
		else if ((ngbrFlags[idx] & NF_XY_PXPY) && (ngbrFlags[idx] & NF_NPX)) {

			//o a e
			diff = (cuVEC<VType>::quantity[idx + 1 + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx + 1]);
		}
		else if ((ngbrFlags[idx] & NF_NPX) && (ngbrFlags[idx] & NF_XY_PXNY)) {

			//e a o
			diff = (cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx + 1 - cuVEC<VType>::n.x]);
		}
		else {

			//right side not available, must have center and left.

			//center column in (+ -) order
			if ((ngbrFlags[idx] & NF_NPY) && (ngbrFlags[idx] & NF_NNY)) {

				//a c a
				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			}
			else if (ngbrFlags[idx] & NF_NPY) {

				//a c e
				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]);
			}
			else {

				//e c a (only other possibility)
				diff = (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
			}

			//now get left side (- +) order

			//yup, a goto statement! I could use a lambda closure but I think this solution is better (faster, or certainly not slower, depending on compiler).
			goto dxy_neu_left_side;
		}

		//right side was available. Now get center or left (- + order) : one of them must be available
		if ((ngbrFlags[idx] & NF_NPY) && (ngbrFlags[idx] & NF_NNY)) {

			//a c a
			diff += (cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]);
			return diff / (4 * cuVEC<VType>::h.x*cuVEC<VType>::h.y);	//finished stencil
		}
		else if (ngbrFlags[idx] & NF_NPY) {

			//a c e
			diff += (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]);
			return diff / (4 * cuVEC<VType>::h.x*cuVEC<VType>::h.y);	//finished stencil
		}
		else if (ngbrFlags[idx] & NF_NNY) {

			//e c a
			diff += (cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]);
			return diff / (4 * cuVEC<VType>::h.x*cuVEC<VType>::h.y);	//finished stencil
		}

		//center column not available (otherwise would have returned): must have left side (- +) order.

	dxy_neu_left_side:

		//left side (- + order)
		if ((ngbrFlags[idx] & NF_XY_NXPY) && (ngbrFlags[idx] & NF_XY_NXNY)) {

			//o - o
			diff += (cuVEC<VType>::quantity[idx - 1 - cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx - 1 + cuVEC<VType>::n.x]);
		}
		else if ((ngbrFlags[idx] & NF_NNX) && (ngbrFlags[idx] & NF_XY_NXPY)) {

			//o o e
			diff += (cuVEC<VType>::quantity[idx - 1] - cuVEC<VType>::quantity[idx - 1 + cuVEC<VType>::n.x]);
		}
		else {

			//e o o (only other possibility left)
			diff += (cuVEC<VType>::quantity[idx - 1 - cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx - 1]);
		}
	}

	return diff / (4 * cuVEC<VType>::h.x*cuVEC<VType>::h.y);
}

//dxz same as dzx
template <typename VType>
__device__ VType cuVEC_VC<VType>::dxz_neu(int idx) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	if ((ngbrFlags[idx] & NF_XZ_FULL) == NF_XZ_FULL) {

		//full off-axis XZ stencil available
		diff = (cuVEC<VType>::quantity[idx + 1 + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - 1 - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx + 1 - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx - 1 + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
	}
	else if (ngbrFlags[idx] & NF_XZ_OASTENCIL) {

		//not full, but XZ off-axis stencil is available

		//try to get right side first
		//right side (+ - order)
		if ((ngbrFlags[idx] & NF_XZ_PXPZ) && (ngbrFlags[idx] & NF_XZ_PXNZ)) {

			//o - o
			diff = (cuVEC<VType>::quantity[idx + 1 + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx + 1 - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
		}
		else if ((ngbrFlags[idx] & NF_XZ_PXPZ) && (ngbrFlags[idx] & NF_NPX)) {

			//o a e
			diff = (cuVEC<VType>::quantity[idx + 1 + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx + 1]);
		}
		else if ((ngbrFlags[idx] & NF_NPX) && (ngbrFlags[idx] & NF_XZ_PXNZ)) {

			//e a o
			diff = (cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx + 1 - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
		}
		else {

			//right side not available, must have center and left.

			//center column in (+ -) order
			if ((ngbrFlags[idx] & NF_NPZ) && (ngbrFlags[idx] & NF_NNZ)) {

				//a c a
				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
			else if (ngbrFlags[idx] & NF_NPZ) {

				//a c e
				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]);
			}
			else {

				//e c a (only other possibility)
				diff = (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}

			//now get left side (- +) order

			//yup, a goto statement! I could use a lambda closure but I think this solution is better (faster, or certainly not slower, depending on compiler).
			goto dxz_neu_left_side;
		}

		//right side was available. Now get center or left (- + order) : one of them must be available
		if ((ngbrFlags[idx] & NF_NPZ) && (ngbrFlags[idx] & NF_NNZ)) {

			//a c a
			diff += (cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			return diff / (4 * cuVEC<VType>::h.x*cuVEC<VType>::h.z);	//finished stencil
		}
		else if (ngbrFlags[idx] & NF_NPZ) {

			//a c e
			diff += (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			return diff / (4 * cuVEC<VType>::h.x*cuVEC<VType>::h.z);	//finished stencil
		}
		else if (ngbrFlags[idx] & NF_NNZ) {

			//e c a
			diff += (cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]);
			return diff / (4 * cuVEC<VType>::h.x*cuVEC<VType>::h.z);	//finished stencil
		}

		//center column not available (otherwise would have returned): must have left side (- +) order.

	dxz_neu_left_side:

		//left side (- + order)
		if ((ngbrFlags[idx] & NF_XZ_NXPZ) && (ngbrFlags[idx] & NF_XZ_NXNZ)) {

			//o - o
			diff += (cuVEC<VType>::quantity[idx - 1 - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx - 1 + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
		}
		else if ((ngbrFlags[idx] & NF_NNX) && (ngbrFlags[idx] & NF_XZ_NXPZ)) {

			//o o e
			diff += (cuVEC<VType>::quantity[idx - 1] - cuVEC<VType>::quantity[idx - 1 + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
		}
		else {

			//e o o (only other possibility left)
			diff += (cuVEC<VType>::quantity[idx - 1 - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx - 1]);
		}
	}

	return diff / (4 * cuVEC<VType>::h.x*cuVEC<VType>::h.z);
}

//dyz same as dzy
template <typename VType>
__device__ VType cuVEC_VC<VType>::dyz_neu(int idx) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	if ((ngbrFlags[idx] & NF_YZ_FULL) == NF_YZ_FULL) {

		//full off-axis YZ stencil available
		diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
	}
	else if (ngbrFlags[idx] & NF_YZ_OASTENCIL) {

		//not full, but YZ off-axis stencil is available

		//try to get right side first
		//right side (+ - order)
		if ((ngbrFlags[idx] & NF_YZ_PYPZ) && (ngbrFlags[idx] & NF_YZ_PYNZ)) {

			//o - o
			diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
		}
		else if ((ngbrFlags[idx] & NF_YZ_PYPZ) && (ngbrFlags[idx] & NF_NPY)) {

			//o a e
			diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x]);
		}
		else if ((ngbrFlags[idx] & NF_NPY) && (ngbrFlags[idx] & NF_YZ_PYNZ)) {

			//e a o
			diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
		}
		else {

			//right side not available, must have center and left.

			//center column in (+ -) order
			if ((ngbrFlags[idx] & NF_NPZ) && (ngbrFlags[idx] & NF_NNZ)) {

				//a c a
				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}
			else if (ngbrFlags[idx] & NF_NPZ) {

				//a c e
				diff = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]);
			}
			else {

				//e c a (only other possibility)
				diff = (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			}

			//now get left side (- +) order

			//yup, a goto statement! I could use a lambda closure but I think this solution is better (faster, or certainly not slower, depending on compiler).
			goto dyz_neu_left_side;
		}

		//right side was available. Now get center or left (- + order) : one of them must be available
		if ((ngbrFlags[idx] & NF_NPZ) && (ngbrFlags[idx] & NF_NNZ)) {

			//a c a
			diff += (cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			return diff / (4 * cuVEC<VType>::h.y*cuVEC<VType>::h.z);	//finished stencil
		}
		else if (ngbrFlags[idx] & NF_NPZ) {

			//a c e
			diff += (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
			return diff / (4 * cuVEC<VType>::h.y*cuVEC<VType>::h.z);	//finished stencil
		}
		else if (ngbrFlags[idx] & NF_NNZ) {

			//e c a
			diff += (cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]);
			return diff / (4 * cuVEC<VType>::h.y*cuVEC<VType>::h.z);	//finished stencil
		}

		//center column not available (otherwise would have returned): must have left side (- +) order.

	dyz_neu_left_side:

		//left side (- + order)
		if ((ngbrFlags[idx] & NF_YZ_NYPZ) && (ngbrFlags[idx] & NF_YZ_NYNZ)) {

			//o - o
			diff += (cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
		}
		else if ((ngbrFlags[idx] & NF_NNY) && (ngbrFlags[idx] & NF_YZ_NYPZ)) {

			//o o e
			diff += (cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x + cuVEC<VType>::n.x*cuVEC<VType>::n.y]);
		}
		else {

			//e o o (only other possibility left)
			diff += (cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x - cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]);
		}
	}

	return diff / (4 * cuVEC<VType>::h.y*cuVEC<VType>::h.z);
}