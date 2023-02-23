#pragma once

#include "VEC_VC.h"

//-------------------------------- GRADIENT OPERATOR

//gradient operator. Use Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
template <typename VType>
VAL3<VType> VEC_VC<VType>::grad_neu(int idx) const
{
	VAL3<VType> diff = VAL3<VType>();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return diff;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		diff.x = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) diff.x = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / VEC<VType>::h.x;
		else if (ngbrFlags[idx] & NF_NNX) diff.x = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / VEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (VEC<VType>::quantity[idx + 1] - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1]) / (2 * VEC<VType>::h.x);
			}
			else {

				diff.x = 0.5 * (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / VEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
			}
			else {

				diff.x = 0.5 * (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / VEC<VType>::h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / VEC<VType>::h.y;
		else if (ngbrFlags[idx] & NF_NNY) diff.y = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / VEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1) * VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diff.y = 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / VEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1) * VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diff.y = 0.5 * (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / VEC<VType>::h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / VEC<VType>::h.z;
		else if (ngbrFlags[idx] & NF_NNZ) diff.z = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / VEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diff.z = 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / VEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diff.z = 0.5 * (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / VEC<VType>::h.z;
			}
		}
	}

	return diff;
}

//gradient operator. Use non-homogeneous Neumann boundary conditions.
//Can be used at composite media boundaries where sided differentials will be used instead.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
template <typename VType>
VAL3<VType> VEC_VC<VType>::grad_nneu(int idx, const VAL3<VType>& bdiff) const
{
	VAL3<VType> diff = VAL3<VType>();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return diff;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		diff.x = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
	}
	//Is it a CMBND boundary? - if not then use non-homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) diff.x = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / VEC<VType>::h.x;
		else if (ngbrFlags[idx] & NF_NNX) diff.x = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / VEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (VEC<VType>::quantity[idx + 1] - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1]) / (2 * VEC<VType>::h.x);
			}
			else {

				diff.x = 0.5 * ((VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / VEC<VType>::h.x + bdiff.x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
			}
			else {

				diff.x = 0.5 * ((VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / VEC<VType>::h.x + bdiff.x);
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / VEC<VType>::h.y;
		else if (ngbrFlags[idx] & NF_NNY) diff.y = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / VEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1) * VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diff.y = 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / VEC<VType>::h.y + bdiff.y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1) * VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diff.y = 0.5 * ((VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / VEC<VType>::h.y + bdiff.y);
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / VEC<VType>::h.z;
		else if (ngbrFlags[idx] & NF_NNZ) diff.z = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / VEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diff.z = 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / VEC<VType>::h.z + bdiff.z);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diff.z = 0.5 * ((VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / VEC<VType>::h.z + bdiff.z);
			}
		}
	}

	return diff;
}

//gradient operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
template <typename VType>
VAL3<VType> VEC_VC<VType>::grad_diri(int idx) const
{
	VAL3<VType> diff = VAL3<VType>();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return diff;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		diff.x = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPX) diff.x = (VEC<VType>::quantity[idx + 1] + VEC<VType>::quantity[idx] - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx)) / (2 * VEC<VType>::h.x);
		else								   diff.x = (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx) - VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
	}
	//Not Dirichlet, is it a CMBND boundary? - if not this either then use homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) diff.x = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / VEC<VType>::h.x;
		else if (ngbrFlags[idx] & NF_NNX) diff.x = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / VEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (VEC<VType>::quantity[idx + 1] - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1]) / (2 * VEC<VType>::h.x);
			}
			else {

				diff.x = 0.5 * (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / VEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
			}
			else {

				diff.x = 0.5 * (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / VEC<VType>::h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPY) diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + VEC<VType>::quantity[idx] - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx)) / (2 * VEC<VType>::h.y);
		else								 diff.y = (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx) - VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / VEC<VType>::h.y;
		else if (ngbrFlags[idx] & NF_NNY) diff.y = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / VEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1) * VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diff.y = 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / VEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1) * VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diff.y = 0.5 * (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / VEC<VType>::h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx] - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx)) / (2 * VEC<VType>::h.z);
		else								 diff.z = (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) - VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / VEC<VType>::h.z;
		else if (ngbrFlags[idx] & NF_NNZ) diff.z = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / VEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diff.z = 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / VEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diff.z = 0.5 * (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / VEC<VType>::h.z;
			}
		}
	}

	return diff;
}

//gradient operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
//Can be used at composite media boundaries where sided differentials will be used instead.
template <typename VType>
VAL3<VType> VEC_VC<VType>::grad_diri_nneu(int idx, const VAL3<VType>& bdiff) const
{
	VAL3<VType> diff = VAL3<VType>();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return diff;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		diff.x = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPX) diff.x = (VEC<VType>::quantity[idx + 1] + VEC<VType>::quantity[idx] - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx)) / (2 * VEC<VType>::h.x);
		else								   diff.x = (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx) - VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
	}
	//Not Dirichlet, is it a CMBND boundary? - if not this either then use homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) diff.x = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / VEC<VType>::h.x;
		else if (ngbrFlags[idx] & NF_NNX) diff.x = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / VEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (VEC<VType>::quantity[idx + 1] - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1]) / (2 * VEC<VType>::h.x);
			}
			else {

				diff.x = 0.5 * ((VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / VEC<VType>::h.x + bdiff.x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
			}
			else {

				diff.x = 0.5 * ((VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / VEC<VType>::h.x + bdiff.x);
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPY) diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + VEC<VType>::quantity[idx] - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx)) / (2 * VEC<VType>::h.y);
		else								 diff.y = (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx) - VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / VEC<VType>::h.y;
		else if (ngbrFlags[idx] & NF_NNY) diff.y = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / VEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1) * VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diff.y = 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / VEC<VType>::h.y + bdiff.y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1) * VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diff.y = 0.5 * ((VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / VEC<VType>::h.y + bdiff.y);
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx] - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx)) / (2 * VEC<VType>::h.z);
		else								 diff.z = (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) - VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / VEC<VType>::h.z;
		else if (ngbrFlags[idx] & NF_NNZ) diff.z = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / VEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diff.z = 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / VEC<VType>::h.z + bdiff.z);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diff.z = 0.5 * ((VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / VEC<VType>::h.z + bdiff.z);
			}
		}
	}

	return diff;
}

//gradient operator. Use sided differentials (also at composite media boundaries)
template <typename VType>
VAL3<VType> VEC_VC<VType>::grad_sided(int idx) const
{
	VAL3<VType> diff = VAL3<VType>();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return diff;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		diff.x = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
	}
	//use sided differentials if one of the neighbors is present
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (VEC<VType>::quantity[idx + 1] - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1]) / (2 * VEC<VType>::h.x);
			}
			else {


				diff.x = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / VEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
			}
			else {

				diff.x = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / VEC<VType>::h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1) * VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diff.y = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / VEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1) * VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diff.y = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / VEC<VType>::h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diff.z = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / VEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diff.z = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / VEC<VType>::h.z;
			}
		}
	}

	return diff;
}