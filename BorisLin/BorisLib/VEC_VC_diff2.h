#pragma once

#include "VEC_VC.h"

//-------------------------------- SECOND ORDER DIFFERENTIALS

//HOMOGENEOUS NEUMANN

//Use Neumann boundary conditions (homogeneous).
//Returns zero at composite media boundary cells.
template <typename VType>
VType VEC_VC<VType>::dxx_neu(int idx) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	double hsq = VEC<VType>::h.x*VEC<VType>::h.x;

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		diff = (VEC<VType>::quantity[idx + 1] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (VEC<VType>::quantity[idx + 1] + get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx - 1] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}

//Use Neumann boundary conditions (homogeneous).
//Returns zero at composite media boundary cells.
template <typename VType>
VType VEC_VC<VType>::dyy_neu(int idx) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	double hsq = VEC<VType>::h.y*VEC<VType>::h.y;

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx - VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}

//Use Neumann boundary conditions (homogeneous).
//Returns zero at composite media boundary cells.
template <typename VType>
VType VEC_VC<VType>::dzz_neu(int idx) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	double hsq = VEC<VType>::h.z*VEC<VType>::h.z;

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / hsq;
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
VType VEC_VC<VType>::dxx_nneu(int idx, const VAL3<VType>& bdiff) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	double hsq = VEC<VType>::h.x*VEC<VType>::h.x;

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		diff = (VEC<VType>::quantity[idx + 1] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (VEC<VType>::quantity[idx + 1] + get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) - bdiff.x * VEC<VType>::h.x) / hsq;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((VEC<VType>::quantity[idx - 1] - VEC<VType>::quantity[idx]) + bdiff.x * VEC<VType>::h.x) / hsq;
			}
		}
	}

	return diff;
}

//Use non-homogeneous Neumann boundary conditions.
//Returns zero at composite media boundary cells
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
template <typename VType>
VType VEC_VC<VType>::dyy_nneu(int idx, const VAL3<VType>& bdiff) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	double hsq = VEC<VType>::h.y*VEC<VType>::h.y;

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) - bdiff.y * VEC<VType>::h.y) / hsq;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((VEC<VType>::quantity[idx - VEC<VType>::n.x] - VEC<VType>::quantity[idx]) + bdiff.y * VEC<VType>::h.y) / hsq;
			}
		}
	}

	return diff;
}

//Use non-homogeneous Neumann boundary conditions.
//Returns zero at composite media boundary cells
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
template <typename VType>
VType VEC_VC<VType>::dzz_nneu(int idx, const VAL3<VType>& bdiff) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();
	
	double hsq = VEC<VType>::h.z*VEC<VType>::h.z;

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) - bdiff.z * VEC<VType>::h.z) / hsq;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) + bdiff.z * VEC<VType>::h.z) / hsq;
			}
		}
	}

	return diff;
}

//DIRICHLET, HOMOGENEOUS NEUMANN

//Use Dirichlet boundary conditions, else Neumann boundary conditions (homogeneous).
//Returns zero at composite media boundary cells.
template <typename VType>
VType VEC_VC<VType>::dxx_diri(int idx) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();
	
	double hsq = VEC<VType>::h.x*VEC<VType>::h.x;

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		diff = (VEC<VType>::quantity[idx + 1] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		//only one neighbor available. Does it use a fixed boundary value (Dirichlet)?
		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPX) diff = (2 * VEC<VType>::quantity[idx + 1] + 4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) - 6 * VEC<VType>::quantity[idx]) / hsq;
			else								   diff = (2 * VEC<VType>::quantity[idx - 1] + 4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) - 6 * VEC<VType>::quantity[idx]) / hsq;
		}
		//stencil at surface for homogeneous Neumann boundary condition - correct, do not multiply by 2! 
		//See conservation law approach derivation for heat equation (integrate then use mid-point rule approximations) for cell-centered grids.
		else if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (VEC<VType>::quantity[idx + 1] + get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx - 1] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}

//Use Dirichlet boundary conditions, else Neumann boundary conditions (homogeneous).
//Returns zero at composite media boundary cells.
template <typename VType>
VType VEC_VC<VType>::dyy_diri(int idx) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	double hsq = VEC<VType>::h.y*VEC<VType>::h.y;

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPY) diff = (2 * VEC<VType>::quantity[idx + VEC<VType>::n.x] + 4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) - 6 * VEC<VType>::quantity[idx]) / hsq;
			else								   diff = (2 * VEC<VType>::quantity[idx - VEC<VType>::n.x] + 4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) - 6 * VEC<VType>::quantity[idx]) / hsq;
		}
		else if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx - VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}

//Use Dirichlet boundary conditions, else Neumann boundary conditions (homogeneous).
//Returns zero at composite media boundary cells.
template <typename VType>
VType VEC_VC<VType>::dzz_diri(int idx) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	double hsq = VEC<VType>::h.z*VEC<VType>::h.z;

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) diff = (2 * VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + 4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) - 6 * VEC<VType>::quantity[idx]) / hsq;
			else								 diff = (2 * VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] + 4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) - 6 * VEC<VType>::quantity[idx]) / hsq;
		}
		else if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / hsq;
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
VType VEC_VC<VType>::dxx_diri_nneu(int idx, const VAL3<VType>& bdiff) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	double hsq = VEC<VType>::h.x*VEC<VType>::h.x;

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		diff = (VEC<VType>::quantity[idx + 1] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		//only one neighbor available. Does it use a fixed boundary value (Dirichlet)?
		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPX) diff = (2 * VEC<VType>::quantity[idx + 1] + 4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) - 6 * VEC<VType>::quantity[idx]) / hsq;
			else								   diff = (2 * VEC<VType>::quantity[idx - 1] + 4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) - 6 * VEC<VType>::quantity[idx]) / hsq;
		}
		else if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (VEC<VType>::quantity[idx + 1] + get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) - bdiff.x * VEC<VType>::h.x) / hsq;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((VEC<VType>::quantity[idx - 1] - VEC<VType>::quantity[idx]) + bdiff.x * VEC<VType>::h.x) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
VType VEC_VC<VType>::dyy_diri_nneu(int idx, const VAL3<VType>& bdiff) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();
	
	double hsq = VEC<VType>::h.y*VEC<VType>::h.y;

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPY) diff = (2 * VEC<VType>::quantity[idx + VEC<VType>::n.x] + 4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) - 6 * VEC<VType>::quantity[idx]) / hsq;
			else								 diff = (2 * VEC<VType>::quantity[idx - VEC<VType>::n.x] + 4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) - 6 * VEC<VType>::quantity[idx]) / hsq;
		}
		else if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) - bdiff.y * VEC<VType>::h.y) / hsq;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((VEC<VType>::quantity[idx - VEC<VType>::n.x] - VEC<VType>::quantity[idx]) + bdiff.y * VEC<VType>::h.y) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
VType VEC_VC<VType>::dzz_diri_nneu(int idx, const VAL3<VType>& bdiff) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	double hsq = VEC<VType>::h.z*VEC<VType>::h.z;

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) diff = (2 * VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + 4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) - 6 * VEC<VType>::quantity[idx]) / hsq;
			else								 diff = (2 * VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] + 4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) - 6 * VEC<VType>::quantity[idx]) / hsq;
		}
		else if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) - bdiff.z * VEC<VType>::h.z) / hsq;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = ((VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) + bdiff.z * VEC<VType>::h.z) / hsq;
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
VType VEC_VC<VType>::dxx_robin(int idx, double K) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	double hsq = VEC<VType>::h.x*VEC<VType>::h.x;

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {
		//both neighbors available so must be an inner point along this direction
		diff = (VEC<VType>::quantity[idx + 1] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		//only one neighbor available. Does it use a Robin condition?
		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_ROBINX)) {

			//Robin cell on +x side of boundary
			if (ngbrFlags2[idx] & NF2_ROBINPX) {

				//use void Robin values
				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff = ((1 + robin_v.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx + 1] + robin_v.i * VEC<VType>::h.x * robin_v.j / K - (1 + 3 * robin_v.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx]) / hsq;
				}
				//don't use void Robin values
				else {

					diff = ((1 + robin_nx.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx + 1] + robin_nx.i * VEC<VType>::h.x * robin_nx.j / K - (1 + 3 * robin_nx.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx]) / hsq;
				}
			}
			//Robin cell on -x side of boundary
			else {

				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff = ((1 + robin_v.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx - 1] + robin_v.i * VEC<VType>::h.x * robin_v.j / K - (1 + 3 * robin_v.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx]) / hsq;
				}
				else {

					diff = ((1 + robin_px.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx - 1] + robin_px.i * VEC<VType>::h.x * robin_px.j / K - (1 + 3 * robin_px.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx]) / hsq;
				}
			}
		}
		//stencil at surface for homogeneous Neumann boundary condition - correct, do not multiply by 2! 
		//See conservation law approach derivation for heat equation (integrate then use mid-point rule approximations) for cell-centered grids.
		else if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (VEC<VType>::quantity[idx + 1] + get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx - 1] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
VType VEC_VC<VType>::dyy_robin(int idx, double K) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	double hsq = VEC<VType>::h.y*VEC<VType>::h.y;

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		//only one neighbor available. Does it use a Robin condition?
		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_ROBINY)) {

			//Robin cell on +x side of boundary
			if (ngbrFlags2[idx] & NF2_ROBINPY) {

				//use void Robin values
				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff = ((1 + robin_v.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx + VEC<VType>::n.x] + robin_v.i * VEC<VType>::h.y * robin_v.j / K - (1 + 3 * robin_v.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx]) / hsq;
				}
				//don't use void Robin values
				else {

					diff = ((1 + robin_ny.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx + VEC<VType>::n.x] + robin_ny.i * VEC<VType>::h.y * robin_ny.j / K - (1 + 3 * robin_ny.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx]) / hsq;
				}
			}
			//Robin cell on -x side of boundary
			else {

				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff = ((1 + robin_v.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx - VEC<VType>::n.x] + robin_v.i * VEC<VType>::h.y * robin_v.j / K - (1 + 3 * robin_v.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx]) / hsq;
				}
				else {

					diff = ((1 + robin_py.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx - VEC<VType>::n.x] + robin_py.i * VEC<VType>::h.y * robin_py.j / K - (1 + 3 * robin_py.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx]) / hsq;
				}
			}
		}
		else if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx - VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}

template <typename VType>
VType VEC_VC<VType>::dzz_robin(int idx, double K) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	double hsq = VEC<VType>::h.z*VEC<VType>::h.z;

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]) / hsq;
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		//only one neighbor available. Does it use a Robin condition?
		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_ROBINZ)) {

			//Robin cell on +x side of boundary
			if (ngbrFlags2[idx] & NF2_ROBINPZ) {

				//use void Robin values
				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff = ((1 + robin_v.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + robin_v.i * VEC<VType>::h.z * robin_v.j / K - (1 + 3 * robin_v.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx]) / hsq;
				}
				//don't use void Robin values
				else {

					diff = ((1 + robin_nz.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + robin_nz.i * VEC<VType>::h.z * robin_nz.j / K - (1 + 3 * robin_nz.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx]) / hsq;
				}
			}
			//Robin cell on -x side of boundary
			else {

				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff = ((1 + robin_v.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] + robin_v.i * VEC<VType>::h.z * robin_v.j / K - (1 + 3 * robin_v.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx]) / hsq;
				}
				else {

					diff = ((1 + robin_pz.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] + robin_pz.i * VEC<VType>::h.z * robin_pz.j / K - (1 + 3 * robin_pz.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx]) / hsq;
				}
			}
		}
		else if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]) / hsq;
			}
			else {

				diff = (VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / hsq;
			}
		}
	}

	return diff;
}


//-------------------------------- SECOND ORDER MIXED

//Use Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
template <typename VType>
VType VEC_VC<VType>::dxy_neu(int idx) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	if ((ngbrFlags[idx] & NF_XY_FULL) == NF_XY_FULL) {

		//full off-axis XY stencil available
		diff = (VEC<VType>::quantity[idx + 1 + VEC<VType>::n.x] + VEC<VType>::quantity[idx - 1 - VEC<VType>::n.x] - VEC<VType>::quantity[idx + 1 - VEC<VType>::n.x] - VEC<VType>::quantity[idx - 1 + VEC<VType>::n.x]);
	}
	else if (ngbrFlags[idx] & NF_XY_OASTENCIL) {

		//not full, but XY off-axis stencil is available

		//try to get right side first
		//right side (+ - order)
		if ((ngbrFlags[idx] & NF_XY_PXPY) && (ngbrFlags[idx] & NF_XY_PXNY)) {

			//o - o
			diff = (VEC<VType>::quantity[idx + 1 + VEC<VType>::n.x] - VEC<VType>::quantity[idx + 1 - VEC<VType>::n.x]);
		}
		else if ((ngbrFlags[idx] & NF_XY_PXPY) && (ngbrFlags[idx] & NF_NPX)) {

			//o a e
			diff = (VEC<VType>::quantity[idx + 1 + VEC<VType>::n.x] - VEC<VType>::quantity[idx + 1]);
		}
		else if ((ngbrFlags[idx] & NF_NPX) && (ngbrFlags[idx] & NF_XY_PXNY)) {

			//e a o
			diff = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx + 1 - VEC<VType>::n.x]);
		}
		else {

			//right side not available, must have center and left.

			//center column in (+ -) order
			if ((ngbrFlags[idx] & NF_NPY) && (ngbrFlags[idx] & NF_NNY)) {

				//a c a
				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]);
			}
			else if (ngbrFlags[idx] & NF_NPY) {

				//a c e
				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]);
			}
			else {

				//e c a (only other possibility)
				diff = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]);
			}

			//now get left side (- +) order

			//yup, a goto statement! I could use a lambda closure but I think this solution is better (faster, or certainly not slower, depending on compiler).
			goto dxy_neu_left_side;
		}

		//right side was available. Now get center or left (- + order) : one of them must be available
		if ((ngbrFlags[idx] & NF_NPY) && (ngbrFlags[idx] & NF_NNY)) {

			//a c a
			diff += (VEC<VType>::quantity[idx - VEC<VType>::n.x] - VEC<VType>::quantity[idx + VEC<VType>::n.x]);
			return diff / (4 * VEC<VType>::h.x*VEC<VType>::h.y);	//finished stencil
		}
		else if (ngbrFlags[idx] & NF_NPY) {

			//a c e
			diff += (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx + VEC<VType>::n.x]);
			return diff / (4 * VEC<VType>::h.x*VEC<VType>::h.y);	//finished stencil
		}
		else if (ngbrFlags[idx] & NF_NNY) {

			//e c a
			diff += (VEC<VType>::quantity[idx - VEC<VType>::n.x] - VEC<VType>::quantity[idx]);
			return diff / (4 * VEC<VType>::h.x*VEC<VType>::h.y);	//finished stencil
		}
		
		//center column not available (otherwise would have returned): must have left side (- +) order.

dxy_neu_left_side:

		//left side (- + order)
		if ((ngbrFlags[idx] & NF_XY_NXPY) && (ngbrFlags[idx] & NF_XY_NXNY)) {

			//o - o
			diff += (VEC<VType>::quantity[idx - 1 - VEC<VType>::n.x] - VEC<VType>::quantity[idx - 1 + VEC<VType>::n.x]);
		}
		else if ((ngbrFlags[idx] & NF_NNX) && (ngbrFlags[idx] & NF_XY_NXPY)) {

			//o o e
			diff += (VEC<VType>::quantity[idx - 1] - VEC<VType>::quantity[idx - 1 + VEC<VType>::n.x]);
		}
		else {

			//e o o (only other possibility left)
			diff += (VEC<VType>::quantity[idx - 1 - VEC<VType>::n.x] - VEC<VType>::quantity[idx - 1]);
		}
	}

	return diff / (4 * VEC<VType>::h.x*VEC<VType>::h.y);
}

//Use Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
template <typename VType>
VType VEC_VC<VType>::dxz_neu(int idx) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	if ((ngbrFlags[idx] & NF_XZ_FULL) == NF_XZ_FULL) {

		//full off-axis XZ stencil available
		diff = (VEC<VType>::quantity[idx + 1 + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - 1 - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx + 1 - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - 1 + VEC<VType>::n.x*VEC<VType>::n.y]);
	}
	else if (ngbrFlags[idx] & NF_XZ_OASTENCIL) {

		//not full, but XZ off-axis stencil is available

		//try to get right side first
		//right side (+ - order)
		if ((ngbrFlags[idx] & NF_XZ_PXPZ) && (ngbrFlags[idx] & NF_XZ_PXNZ)) {

			//o - o
			diff = (VEC<VType>::quantity[idx + 1 + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx + 1 - VEC<VType>::n.x*VEC<VType>::n.y]);
		}
		else if ((ngbrFlags[idx] & NF_XZ_PXPZ) && (ngbrFlags[idx] & NF_NPX)) {

			//o a e
			diff = (VEC<VType>::quantity[idx + 1 + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx + 1]);
		}
		else if ((ngbrFlags[idx] & NF_NPX) && (ngbrFlags[idx] & NF_XZ_PXNZ)) {

			//e a o
			diff = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx + 1 - VEC<VType>::n.x*VEC<VType>::n.y]);
		}
		else {

			//right side not available, must have center and left.

			//center column in (+ -) order
			if ((ngbrFlags[idx] & NF_NPZ) && (ngbrFlags[idx] & NF_NNZ)) {

				//a c a
				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
			}
			else if (ngbrFlags[idx] & NF_NPZ) {

				//a c e
				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]);
			}
			else {

				//e c a (only other possibility)
				diff = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
			}

			//now get left side (- +) order

			//yup, a goto statement! I could use a lambda closure but I think this solution is better (faster, or certainly not slower, depending on compiler).
			goto dxz_neu_left_side;
		}

		//right side was available. Now get center or left (- + order) : one of them must be available
		if ((ngbrFlags[idx] & NF_NPZ) && (ngbrFlags[idx] & NF_NNZ)) {

			//a c a
			diff += (VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y]);
			return diff / (4 * VEC<VType>::h.x*VEC<VType>::h.z);	//finished stencil
		}
		else if (ngbrFlags[idx] & NF_NPZ) {

			//a c e
			diff += (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y]);
			return diff / (4 * VEC<VType>::h.x*VEC<VType>::h.z);	//finished stencil
		}
		else if (ngbrFlags[idx] & NF_NNZ) {

			//e c a
			diff += (VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]);
			return diff / (4 * VEC<VType>::h.x*VEC<VType>::h.z);	//finished stencil
		}

		//center column not available (otherwise would have returned): must have left side (- +) order.

	dxz_neu_left_side:

		//left side (- + order)
		if ((ngbrFlags[idx] & NF_XZ_NXPZ) && (ngbrFlags[idx] & NF_XZ_NXNZ)) {

			//o - o
			diff += (VEC<VType>::quantity[idx - 1 - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - 1 + VEC<VType>::n.x*VEC<VType>::n.y]);
		}
		else if ((ngbrFlags[idx] & NF_NNX) && (ngbrFlags[idx] & NF_XZ_NXPZ)) {

			//o o e
			diff += (VEC<VType>::quantity[idx - 1] - VEC<VType>::quantity[idx - 1 + VEC<VType>::n.x*VEC<VType>::n.y]);
		}
		else {

			//e o o (only other possibility left)
			diff += (VEC<VType>::quantity[idx - 1 - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - 1]);
		}
	}

	return diff / (4 * VEC<VType>::h.x*VEC<VType>::h.z);
}

//Use Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
template <typename VType>
VType VEC_VC<VType>::dyz_neu(int idx) const
{
	VType diff = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	if ((ngbrFlags[idx] & NF_YZ_FULL) == NF_YZ_FULL) {

		//full off-axis YZ stencil available
		diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx + VEC<VType>::n.x - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x + VEC<VType>::n.x*VEC<VType>::n.y]);
	}
	else if (ngbrFlags[idx] & NF_YZ_OASTENCIL) {

		//not full, but YZ off-axis stencil is available

		//try to get right side first
		//right side (+ - order)
		if ((ngbrFlags[idx] & NF_YZ_PYPZ) && (ngbrFlags[idx] & NF_YZ_PYNZ)) {

			//o - o
			diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx + VEC<VType>::n.x - VEC<VType>::n.x*VEC<VType>::n.y]);
		}
		else if ((ngbrFlags[idx] & NF_YZ_PYPZ) && (ngbrFlags[idx] & NF_NPY)) {

			//o a e
			diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx + VEC<VType>::n.x]);
		}
		else if ((ngbrFlags[idx] & NF_NPY) && (ngbrFlags[idx] & NF_YZ_PYNZ)) {

			//e a o
			diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx + VEC<VType>::n.x - VEC<VType>::n.x*VEC<VType>::n.y]);
		}
		else {

			//right side not available, must have center and left.

			//center column in (+ -) order
			if ((ngbrFlags[idx] & NF_NPZ) && (ngbrFlags[idx] & NF_NNZ)) {

				//a c a
				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
			}
			else if (ngbrFlags[idx] & NF_NPZ) {

				//a c e
				diff = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]);
			}
			else {

				//e c a (only other possibility)
				diff = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]);
			}

			//now get left side (- +) order

			//yup, a goto statement! I could use a lambda closure but I think this solution is better (faster, or certainly not slower, depending on compiler).
			goto dyz_neu_left_side;
		}

		//right side was available. Now get center or left (- + order) : one of them must be available
		if ((ngbrFlags[idx] & NF_NPZ) && (ngbrFlags[idx] & NF_NNZ)) {

			//a c a
			diff += (VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y]);
			return diff / (4 * VEC<VType>::h.y*VEC<VType>::h.z);	//finished stencil
		}
		else if (ngbrFlags[idx] & NF_NPZ) {

			//a c e
			diff += (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y]);
			return diff / (4 * VEC<VType>::h.y*VEC<VType>::h.z);	//finished stencil
		}
		else if (ngbrFlags[idx] & NF_NNZ) {

			//e c a
			diff += (VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]);
			return diff / (4 * VEC<VType>::h.y*VEC<VType>::h.z);	//finished stencil
		}

		//center column not available (otherwise would have returned): must have left side (- +) order.

	dyz_neu_left_side:

		//left side (- + order)
		if ((ngbrFlags[idx] & NF_YZ_NYPZ) && (ngbrFlags[idx] & NF_YZ_NYNZ)) {

			//o - o
			diff += (VEC<VType>::quantity[idx - VEC<VType>::n.x - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x + VEC<VType>::n.x*VEC<VType>::n.y]);
		}
		else if ((ngbrFlags[idx] & NF_NNY) && (ngbrFlags[idx] & NF_YZ_NYPZ)) {

			//o o e
			diff += (VEC<VType>::quantity[idx - VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x + VEC<VType>::n.x*VEC<VType>::n.y]);
		}
		else {

			//e o o (only other possibility left)
			diff += (VEC<VType>::quantity[idx - VEC<VType>::n.x - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x]);
		}
	}

	return diff / (4 * VEC<VType>::h.y*VEC<VType>::h.z);
}