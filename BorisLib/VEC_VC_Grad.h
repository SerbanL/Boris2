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
		diff.x = (quantity[idx + 1] - quantity[idx - 1]) / (2 * h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) diff.x = (quantity[idx + 1] - quantity[idx]) / h.x;
		if (ngbrFlags[idx] & NF_NNX) diff.x = (quantity[idx] - quantity[idx - 1]) / h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (quantity[idx + 1] - quantity[idx + n.x - 1]) / (2 * h.x);
			}
			else {

				diff.x = 0.5 * (quantity[idx + 1] - quantity[idx]) / h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (quantity[idx - (n.x - 1)] - quantity[idx - 1]) / (2 * h.x);
			}
			else {

				diff.x = 0.5 * (quantity[idx] - quantity[idx - 1]) / h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff.y = (quantity[idx + n.x] - quantity[idx - n.x]) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) diff.y = (quantity[idx + n.x] - quantity[idx]) / h.y;
		if (ngbrFlags[idx] & NF_NNY) diff.y = (quantity[idx] - quantity[idx - n.x]) / h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (quantity[idx + n.x] - quantity[idx + (n.y - 1) * n.x]) / (2 * h.y);
			}
			else {

				diff.y = 0.5 * (quantity[idx + n.x] - quantity[idx]) / h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (quantity[idx - (n.y - 1) * n.x] - quantity[idx - n.x]) / (2 * h.y);
			}
			else {

				diff.y = 0.5 * (quantity[idx] - quantity[idx - n.x]) / h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff.z = (quantity[idx + n.x*n.y] - quantity[idx - n.x*n.y]) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) diff.z = (quantity[idx + n.x*n.y] - quantity[idx]) / h.z;
		if (ngbrFlags[idx] & NF_NNZ) diff.z = (quantity[idx] - quantity[idx - n.x*n.y]) / h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (quantity[idx + n.x*n.y] - quantity[idx + (n.z - 1) * n.x*n.y]) / (2 * h.z);
			}
			else {

				diff.z = 0.5 * (quantity[idx + n.x*n.y] - quantity[idx]) / h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (quantity[idx - (n.z - 1) * n.x*n.y] - quantity[idx - n.x*n.y]) / (2 * h.z);
			}
			else {

				diff.z = 0.5 * (quantity[idx] - quantity[idx - n.x*n.y]) / h.z;
			}
		}
	}

	return diff;
}

//gradient operator. Use non-homogeneous Neumann boundary conditions.
//Can be used at composite media boundaries where sided differentials will be used instead.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
template <typename VType>
VAL3<VType> VEC_VC<VType>::grad_nneu(int idx, VAL3<VType>& bdiff) const
{
	VAL3<VType> diff = VAL3<VType>();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return diff;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		diff.x = (quantity[idx + 1] - quantity[idx - 1]) / (2 * h.x);
	}
	//Is it a CMBND boundary? - if not then use non-homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) diff.x = (quantity[idx + 1] - quantity[idx]) / h.x;
		if (ngbrFlags[idx] & NF_NNX) diff.x = (quantity[idx] - quantity[idx - 1]) / h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (quantity[idx + 1] - quantity[idx + n.x - 1]) / (2 * h.x);
			}
			else {

				diff.x = 0.5 * ((quantity[idx + 1] - quantity[idx]) / h.x + bdiff.x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (quantity[idx - (n.x - 1)] - quantity[idx - 1]) / (2 * h.x);
			}
			else {

				diff.x = 0.5 * ((quantity[idx] - quantity[idx - 1]) / h.x + bdiff.x);
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff.y = (quantity[idx + n.x] - quantity[idx - n.x]) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) diff.y = (quantity[idx + n.x] - quantity[idx]) / h.y;
		if (ngbrFlags[idx] & NF_NNY) diff.y = (quantity[idx] - quantity[idx - n.x]) / h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (quantity[idx + n.x] - quantity[idx + (n.y - 1) * n.x]) / (2 * h.y);
			}
			else {

				diff.y = 0.5 * ((quantity[idx + n.x] - quantity[idx]) / h.y + bdiff.y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (quantity[idx - (n.y - 1) * n.x] - quantity[idx - n.x]) / (2 * h.y);
			}
			else {

				diff.y = 0.5 * ((quantity[idx] - quantity[idx - n.x]) / h.y + bdiff.y);
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff.z = (quantity[idx + n.x*n.y] - quantity[idx - n.x*n.y]) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) diff.z = (quantity[idx + n.x*n.y] - quantity[idx]) / h.z;
		if (ngbrFlags[idx] & NF_NNZ) diff.z = (quantity[idx] - quantity[idx - n.x*n.y]) / h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (quantity[idx + n.x*n.y] - quantity[idx + (n.z - 1) * n.x*n.y]) / (2 * h.z);
			}
			else {

				diff.z = 0.5 * ((quantity[idx + n.x*n.y] - quantity[idx]) / h.z + bdiff.z);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (quantity[idx - (n.z - 1) * n.x*n.y] - quantity[idx - n.x*n.y]) / (2 * h.z);
			}
			else {

				diff.z = 0.5 * ((quantity[idx] - quantity[idx - n.x*n.y]) / h.z + bdiff.z);
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
		diff.x = (quantity[idx + 1] - quantity[idx - 1]) / (2 * h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (ngbrFlags[idx] & NF_DIRICHLETX) {

		if (ngbrFlags[idx] & NF_DIRICHLETPX) diff.x = (quantity[idx + 1] + quantity[idx] - 2 * get_dirichlet_value(NF_DIRICHLETPX, idx)) / (2 * h.x);
		else								 diff.x = (2 * get_dirichlet_value(NF_DIRICHLETNX, idx) - quantity[idx] - quantity[idx - 1]) / (2 * h.x);
	}
	//Not Dirichlet, is it a CMBND boundary? - if not this either then use homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) diff.x = (quantity[idx + 1] - quantity[idx]) / h.x;
		if (ngbrFlags[idx] & NF_NNX) diff.x = (quantity[idx] - quantity[idx - 1]) / h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (quantity[idx + 1] - quantity[idx + n.x - 1]) / (2 * h.x);
			}
			else {

				diff.x = 0.5 * (quantity[idx + 1] - quantity[idx]) / h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (quantity[idx - (n.x - 1)] - quantity[idx - 1]) / (2 * h.x);
			}
			else {

				diff.x = 0.5 * (quantity[idx] - quantity[idx - 1]) / h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff.y = (quantity[idx + n.x] - quantity[idx - n.x]) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETY) {

		if (ngbrFlags[idx] & NF_DIRICHLETPY) diff.y = (quantity[idx + n.x] + quantity[idx] - 2 * get_dirichlet_value(NF_DIRICHLETPY, idx)) / (2 * h.y);
		else								 diff.y = (2 * get_dirichlet_value(NF_DIRICHLETNY, idx) - quantity[idx] - quantity[idx - n.x]) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) diff.y = (quantity[idx + n.x] - quantity[idx]) / h.y;
		if (ngbrFlags[idx] & NF_NNY) diff.y = (quantity[idx] - quantity[idx - n.x]) / h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (quantity[idx + n.x] - quantity[idx + (n.y - 1) * n.x]) / (2 * h.y);
			}
			else {

				diff.y = 0.5 * (quantity[idx + n.x] - quantity[idx]) / h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (quantity[idx - (n.y - 1) * n.x] - quantity[idx - n.x]) / (2 * h.y);
			}
			else {

				diff.y = 0.5 * (quantity[idx] - quantity[idx - n.x]) / h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff.z = (quantity[idx + n.x*n.y] - quantity[idx - n.x*n.y]) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

		if (ngbrFlags[idx] & NF_DIRICHLETPZ) diff.z = (quantity[idx + n.x*n.y] + quantity[idx] - 2 * get_dirichlet_value(NF_DIRICHLETPZ, idx)) / (2 * h.z);
		else								 diff.z = (2 * get_dirichlet_value(NF_DIRICHLETNZ, idx) - quantity[idx] - quantity[idx - n.x*n.y]) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) diff.z = (quantity[idx + n.x*n.y] - quantity[idx]) / h.z;
		if (ngbrFlags[idx] & NF_NNZ) diff.z = (quantity[idx] - quantity[idx - n.x*n.y]) / h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (quantity[idx + n.x*n.y] - quantity[idx + (n.z - 1) * n.x*n.y]) / (2 * h.z);
			}
			else {

				diff.z = 0.5 * (quantity[idx + n.x*n.y] - quantity[idx]) / h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (quantity[idx - (n.z - 1) * n.x*n.y] - quantity[idx - n.x*n.y]) / (2 * h.z);
			}
			else {

				diff.z = 0.5 * (quantity[idx] - quantity[idx - n.x*n.y]) / h.z;
			}
		}
	}

	return diff;
}

//gradient operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
//Can be used at composite media boundaries where sided differentials will be used instead.
template <typename VType>
VAL3<VType> VEC_VC<VType>::grad_diri_nneu(int idx, VAL3<VType>& bdiff) const
{
	VAL3<VType> diff = VAL3<VType>();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return diff;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		diff.x = (quantity[idx + 1] - quantity[idx - 1]) / (2 * h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (ngbrFlags[idx] & NF_DIRICHLETX) {

		if (ngbrFlags[idx] & NF_DIRICHLETPX) diff.x = (quantity[idx + 1] + quantity[idx] - 2 * get_dirichlet_value(NF_DIRICHLETPX, idx)) / (2 * h.x);
		else								 diff.x = (2 * get_dirichlet_value(NF_DIRICHLETNX, idx) - quantity[idx] - quantity[idx - 1]) / (2 * h.x);
	}
	//Not Dirichlet, is it a CMBND boundary? - if not this either then use homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) diff.x = (quantity[idx + 1] - quantity[idx]) / h.x;
		if (ngbrFlags[idx] & NF_NNX) diff.x = (quantity[idx] - quantity[idx - 1]) / h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (quantity[idx + 1] - quantity[idx + n.x - 1]) / (2 * h.x);
			}
			else {

				diff.x = 0.5 * ((quantity[idx + 1] - quantity[idx]) / h.x + bdiff.x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (quantity[idx - (n.x - 1)] - quantity[idx - 1]) / (2 * h.x);
			}
			else {

				diff.x = 0.5 * ((quantity[idx] - quantity[idx - 1]) / h.x + bdiff.x);
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff.y = (quantity[idx + n.x] - quantity[idx - n.x]) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETY) {

		if (ngbrFlags[idx] & NF_DIRICHLETPY) diff.y = (quantity[idx + n.x] + quantity[idx] - 2 * get_dirichlet_value(NF_DIRICHLETPY, idx)) / (2 * h.y);
		else								 diff.y = (2 * get_dirichlet_value(NF_DIRICHLETNY, idx) - quantity[idx] - quantity[idx - n.x]) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) diff.y = (quantity[idx + n.x] - quantity[idx]) / h.y;
		if (ngbrFlags[idx] & NF_NNY) diff.y = (quantity[idx] - quantity[idx - n.x]) / h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (quantity[idx + n.x] - quantity[idx + (n.y - 1) * n.x]) / (2 * h.y);
			}
			else {

				diff.y = 0.5 * ((quantity[idx + n.x] - quantity[idx]) / h.y + bdiff.y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (quantity[idx - (n.y - 1) * n.x] - quantity[idx - n.x]) / (2 * h.y);
			}
			else {

				diff.y = 0.5 * ((quantity[idx] - quantity[idx - n.x]) / h.y + bdiff.y);
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff.z = (quantity[idx + n.x*n.y] - quantity[idx - n.x*n.y]) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

		if (ngbrFlags[idx] & NF_DIRICHLETPZ) diff.z = (quantity[idx + n.x*n.y] + quantity[idx] - 2 * get_dirichlet_value(NF_DIRICHLETPZ, idx)) / (2 * h.z);
		else								 diff.z = (2 * get_dirichlet_value(NF_DIRICHLETNZ, idx) - quantity[idx] - quantity[idx - n.x*n.y]) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) diff.z = (quantity[idx + n.x*n.y] - quantity[idx]) / h.z;
		if (ngbrFlags[idx] & NF_NNZ) diff.z = (quantity[idx] - quantity[idx - n.x*n.y]) / h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (quantity[idx + n.x*n.y] - quantity[idx + (n.z - 1) * n.x*n.y]) / (2 * h.z);
			}
			else {

				diff.z = 0.5 * ((quantity[idx + n.x*n.y] - quantity[idx]) / h.z + bdiff.z);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (quantity[idx - (n.z - 1) * n.x*n.y] - quantity[idx - n.x*n.y]) / (2 * h.z);
			}
			else {

				diff.z = 0.5 * ((quantity[idx] - quantity[idx - n.x*n.y]) / h.z + bdiff.z);
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
		diff.x = (quantity[idx + 1] - quantity[idx - 1]) / (2 * h.x);
	}
	//use sided differentials if one of the neighbors is present
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (quantity[idx + 1] - quantity[idx + n.x - 1]) / (2 * h.x);
			}
			else {


				diff.x = (quantity[idx + 1] - quantity[idx]) / h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diff.x = (quantity[idx - (n.x - 1)] - quantity[idx - 1]) / (2 * h.x);
			}
			else {

				diff.x = (quantity[idx] - quantity[idx - 1]) / h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff.y = (quantity[idx + n.x] - quantity[idx - n.x]) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (quantity[idx + n.x] - quantity[idx + (n.y - 1) * n.x]) / (2 * h.y);
			}
			else {

				diff.y = (quantity[idx + n.x] - quantity[idx]) / h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diff.y = (quantity[idx - (n.y - 1) * n.x] - quantity[idx - n.x]) / (2 * h.y);
			}
			else {

				diff.y = (quantity[idx] - quantity[idx - n.x]) / h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff.z = (quantity[idx + n.x*n.y] - quantity[idx - n.x*n.y]) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (quantity[idx + n.x*n.y] - quantity[idx + (n.z - 1) * n.x*n.y]) / (2 * h.z);
			}
			else {

				diff.z = (quantity[idx + n.x*n.y] - quantity[idx]) / h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diff.z = (quantity[idx - (n.z - 1) * n.x*n.y] - quantity[idx - n.x*n.y]) / (2 * h.z);
			}
			else {

				diff.z = (quantity[idx] - quantity[idx - n.x*n.y]) / h.z;
			}
		}
	}

	return diff;
}