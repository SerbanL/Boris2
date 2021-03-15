#pragma once

#include "VEC_VC.h"

//-------------------------------- LAPLACE OPERATOR

//calculate Laplace operator at cell with given index. Use Neumann boundary conditions (homogeneous).
//Returns zero at composite media boundary cells.
template <typename VType>
VType VEC_VC<VType>::delsq_neu(int idx) const
{
	VType diff_x = VType(), diff_y = VType(), diff_z = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		diff_x = (VEC<VType>::quantity[idx + 1] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff_x = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff_x = (VEC<VType>::quantity[idx + 1] + get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_x = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff_x = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_x = (VEC<VType>::quantity[idx - 1] - VEC<VType>::quantity[idx]);
			}
		}
	}

	diff_x /= (VEC<VType>::h.x*VEC<VType>::h.x);

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff_y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff_y = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff_y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff_y = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_y = (VEC<VType>::quantity[idx - VEC<VType>::n.x] - VEC<VType>::quantity[idx]);
			}
		}
	}

	diff_y /= (VEC<VType>::h.y*VEC<VType>::h.y);

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff_z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff_z = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff_z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff_z = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_z = (VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]);
			}
		}
	}

	diff_z /= (VEC<VType>::h.z*VEC<VType>::h.z);

	return (diff_x + diff_y + diff_z);
}

//calculate Laplace operator at cell with given index. Use non-homogeneous Neumann boundary conditions with the specified boundary differential.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
//Returns zero at composite media boundary cells.
template <typename VType>
VType VEC_VC<VType>::delsq_nneu(int idx, const VAL3<VType>& bdiff) const
{
	VType diff_x = VType(), diff_y = VType(), diff_z = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		diff_x = (VEC<VType>::quantity[idx + 1] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff_x = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff_x = (VEC<VType>::quantity[idx + 1] + get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_x = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) - bdiff.x * VEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff_x = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_x = (VEC<VType>::quantity[idx - 1] - VEC<VType>::quantity[idx]) + bdiff.x * VEC<VType>::h.x;
			}
		}
	}

	diff_x /= (VEC<VType>::h.x*VEC<VType>::h.x);

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff_y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff_y = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff_y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) - bdiff.y * VEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff_y = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_y = (VEC<VType>::quantity[idx - VEC<VType>::n.x] - VEC<VType>::quantity[idx]) + bdiff.y * VEC<VType>::h.y;
			}
		}
	}

	diff_y /= (VEC<VType>::h.y*VEC<VType>::h.y);

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff_z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff_z = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff_z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) - bdiff.z * VEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff_z = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_z = (VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) + bdiff.z * VEC<VType>::h.z;
			}
		}
	}

	diff_z /= (VEC<VType>::h.z*VEC<VType>::h.z);

	return (diff_x + diff_y + diff_z);
}

//calculate Laplace operator at cell with given index. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
//Returns zero at composite media boundary cells.
template <typename VType>
VType VEC_VC<VType>::delsq_diri(int idx) const
{
	VType diff_x = VType(), diff_y = VType(), diff_z = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		diff_x = (VEC<VType>::quantity[idx + 1] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff_x = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		//only one neighbor available. Does it use a fixed boundary value (Dirichlet)?
		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPX) diff_x = (2 * VEC<VType>::quantity[idx + 1] + 4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) - 6 * VEC<VType>::quantity[idx]);
			else								   diff_x = (2 * VEC<VType>::quantity[idx - 1] + 4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) - 6 * VEC<VType>::quantity[idx]);
		}
		//stencil at surface for homogeneous Neumann boundary condition - correct, do not multiply by 2! 
		//See conservation law approach derivation for heat equation (integrate then use mid-point rule approximations) for cell-centered grids.
		else if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff_x = (VEC<VType>::quantity[idx + 1] + get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_x = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff_x = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_x = (VEC<VType>::quantity[idx - 1] - VEC<VType>::quantity[idx]);
			}
		}
	}

	diff_x /= (VEC<VType>::h.x*VEC<VType>::h.x);

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff_y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff_y = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPY) diff_y = (2 * VEC<VType>::quantity[idx + VEC<VType>::n.x] + 4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) - 6 * VEC<VType>::quantity[idx]);
			else								 diff_y = (2 * VEC<VType>::quantity[idx - VEC<VType>::n.x] + 4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) - 6 * VEC<VType>::quantity[idx]);
		}
		else if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff_y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff_y = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_y = (VEC<VType>::quantity[idx - VEC<VType>::n.x] - VEC<VType>::quantity[idx]);
			}
		}
	}

	diff_y /= (VEC<VType>::h.y*VEC<VType>::h.y);

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff_z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff_z = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) diff_z = (2 * VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + 4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) - 6 * VEC<VType>::quantity[idx]);
			else								 diff_z = (2 * VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] + 4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) - 6 * VEC<VType>::quantity[idx]);
		}
		else if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff_z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff_z = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_z = (VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]);
			}
		}
	}

	diff_z /= (VEC<VType>::h.z*VEC<VType>::h.z);

	return (diff_x + diff_y + diff_z);
}

//calculate Laplace operator at cell with given index. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
//Returns zero at composite media boundary cells.
template <typename VType>
VType VEC_VC<VType>::delsq_diri_nneu(int idx, const VAL3<VType>& bdiff) const
{
	VType diff_x = VType(), diff_y = VType(), diff_z = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		diff_x = (VEC<VType>::quantity[idx + 1] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff_x = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		//only one neighbor available. Does it use a fixed boundary value (Dirichlet)?
		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPX) diff_x = (2 * VEC<VType>::quantity[idx + 1] + 4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) - 6 * VEC<VType>::quantity[idx]);
			else								   diff_x = (2 * VEC<VType>::quantity[idx - 1] + 4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) - 6 * VEC<VType>::quantity[idx]);
		}
		else if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff_x = (VEC<VType>::quantity[idx + 1] + get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_x = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) - bdiff.x * VEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff_x = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_x = (VEC<VType>::quantity[idx - 1] - VEC<VType>::quantity[idx]) + bdiff.x * VEC<VType>::h.x;
			}
		}
	}


	diff_x /= (VEC<VType>::h.x*VEC<VType>::h.x);

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff_y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff_y = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPY) diff_y = (2 * VEC<VType>::quantity[idx + VEC<VType>::n.x] + 4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) - 6 * VEC<VType>::quantity[idx]);
			else								 diff_y = (2 * VEC<VType>::quantity[idx - VEC<VType>::n.x] + 4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) - 6 * VEC<VType>::quantity[idx]);
		}
		else if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff_y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) - bdiff.y * VEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff_y = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_y = (VEC<VType>::quantity[idx - VEC<VType>::n.x] - VEC<VType>::quantity[idx]) + bdiff.y * VEC<VType>::h.y;
			}
		}
	}

	diff_y /= (VEC<VType>::h.y*VEC<VType>::h.y);

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff_z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff_z = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

			if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) diff_z = (2 * VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + 4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) - 6 * VEC<VType>::quantity[idx]);
			else								 diff_z = (2 * VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] + 4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) - 6 * VEC<VType>::quantity[idx]);
		}
		else if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff_z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) - bdiff.z * VEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff_z = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_z = (VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) + bdiff.z * VEC<VType>::h.z;
			}
		}
	}

	diff_z /= (VEC<VType>::h.z*VEC<VType>::h.z);

	return (diff_x + diff_y + diff_z);
}

//calculate Laplace operator at cell with given index. Use Robin boundary conditions, defaulting to Neumann if not set.
//Returns zero at composite media boundary cells.
template <typename VType>
VType VEC_VC<VType>::delsq_robin(int idx, double K) const
{
	VType diff_x = VType(), diff_y = VType(), diff_z = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x axis
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {
		//both neighbors available so must be an inner point along this direction
		diff_x = (VEC<VType>::quantity[idx + 1] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff_x = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		//only one neighbor available. Does it use a Robin condition?
		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_ROBINX)) {

			//Robin cell on +x side of boundary
			if (ngbrFlags2[idx] & NF2_ROBINPX) {

				//use void Robin values
				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff_x = ((1 + robin_v.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx + 1] + robin_v.i * VEC<VType>::h.x * robin_v.j / K - (1 + 3 * robin_v.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx]);
				}
				//don't use void Robin values
				else {

					diff_x = ((1 + robin_nx.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx + 1] + robin_nx.i * VEC<VType>::h.x * robin_nx.j / K - (1 + 3 * robin_nx.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx]);
				}
			}
			//Robin cell on -x side of boundary
			else {

				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff_x = ((1 + robin_v.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx - 1] + robin_v.i * VEC<VType>::h.x * robin_v.j / K - (1 + 3 * robin_v.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx]);
				}
				else {

					diff_x = ((1 + robin_px.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx - 1] + robin_px.i * VEC<VType>::h.x * robin_px.j / K - (1 + 3 * robin_px.i * VEC<VType>::h.x / (2 * K)) * VEC<VType>::quantity[idx]);
				}
			}
		}
		//stencil at surface for homogeneous Neumann boundary condition - correct, do not multiply by 2! 
		//See conservation law approach derivation for heat equation (integrate then use mid-point rule approximations) for cell-centered grids.
		else if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff_x = (VEC<VType>::quantity[idx + 1] + get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_x = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff_x = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] + VEC<VType>::quantity[idx - 1] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_x = (VEC<VType>::quantity[idx - 1] - VEC<VType>::quantity[idx]);
			}
		}
	}


	diff_x /= (VEC<VType>::h.x*VEC<VType>::h.x);

	//y axis
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff_y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff_y = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		//only one neighbor available. Does it use a Robin condition?
		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_ROBINY)) {

			//Robin cell on +x side of boundary
			if (ngbrFlags2[idx] & NF2_ROBINPY) {

				//use void Robin values
				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff_y = ((1 + robin_v.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx + VEC<VType>::n.x] + robin_v.i * VEC<VType>::h.y * robin_v.j / K - (1 + 3 * robin_v.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx]);
				}
				//don't use void Robin values
				else {

					diff_y = ((1 + robin_ny.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx + VEC<VType>::n.x] + robin_ny.i * VEC<VType>::h.y * robin_ny.j / K - (1 + 3 * robin_ny.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx]);
				}
			}
			//Robin cell on -x side of boundary
			else {

				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff_y = ((1 + robin_v.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx - VEC<VType>::n.x] + robin_v.i * VEC<VType>::h.y * robin_v.j / K - (1 + 3 * robin_v.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx]);
				}
				else {

					diff_y = ((1 + robin_py.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx - VEC<VType>::n.x] + robin_py.i * VEC<VType>::h.y * robin_py.j / K - (1 + 3 * robin_py.i * VEC<VType>::h.y / (2 * K)) * VEC<VType>::quantity[idx]);
				}
			}
		}
		else if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff_y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff_y = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x] + VEC<VType>::quantity[idx - VEC<VType>::n.x] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_y = (VEC<VType>::quantity[idx - VEC<VType>::n.x] - VEC<VType>::quantity[idx]);
			}
		}
	}

	diff_y /= (VEC<VType>::h.y*VEC<VType>::h.y);

	//z axis
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff_z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff_z = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		//only one neighbor available. Does it use a Robin condition?
		if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_ROBINZ)) {

			//Robin cell on +x side of boundary
			if (ngbrFlags2[idx] & NF2_ROBINPZ) {

				//use void Robin values
				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff_z = ((1 + robin_v.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + robin_v.i * VEC<VType>::h.z * robin_v.j / K - (1 + 3 * robin_v.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx]);
				}
				//don't use void Robin values
				else {

					diff_z = ((1 + robin_nz.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + robin_nz.i * VEC<VType>::h.z * robin_nz.j / K - (1 + 3 * robin_nz.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx]);
				}
			}
			//Robin cell on -x side of boundary
			else {

				if (ngbrFlags2[idx] & NF2_ROBINV) {

					diff_z = ((1 + robin_v.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] + robin_v.i * VEC<VType>::h.z * robin_v.j / K - (1 + 3 * robin_v.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx]);
				}
				else {

					diff_z = ((1 + robin_pz.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] + robin_pz.i * VEC<VType>::h.z * robin_pz.j / K - (1 + 3 * robin_pz.i * VEC<VType>::h.z / (2 * K)) * VEC<VType>::quantity[idx]);
				}
			}
		}
		else if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff_z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff_z = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - 2 * VEC<VType>::quantity[idx]);
			}
			else {

				diff_z = (VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]);
			}
		}
	}

	diff_z /= (VEC<VType>::h.z*VEC<VType>::h.z);

	return (diff_x + diff_y + diff_z);
}