#pragma once

#include "VEC_VC.h"

//-------------------------------- DIVERGENCE OPERATOR

//divergence operator. Use Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
//div operator can be applied if VType is a VAL3<Type>, returning Type
template double VEC_VC<DBL3>::div_neu(int idx) const;

template <typename VType>
double VEC_VC<VType>::div_neu(int idx) const
{
	double div = 0.0;
	
	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if (ngbrFlags[idx] & NF_BOTHX) {

		//inner point along this direction
		div += (quantity[idx + 1].x - quantity[idx - 1].x) / (2 * h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) div += (quantity[idx + 1].x - quantity[idx].x) / h.x;
		if (ngbrFlags[idx] & NF_NNX) div += (quantity[idx].x - quantity[idx - 1].x) / h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (quantity[idx + 1].x - quantity[idx + n.x - 1].x) / (2 * h.x);
			}
			else {

				div += 0.5 * (quantity[idx + 1].x - quantity[idx].x) / h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (quantity[idx - (n.x - 1)].x - quantity[idx - 1].x) / (2 * h.x);
			}
			else {

				div += 0.5 * (quantity[idx].x - quantity[idx - 1].x) / h.x;
			}
		}
	}

	//y direction
	if (ngbrFlags[idx] & NF_BOTHY) {

		div += (quantity[idx + n.x].y - quantity[idx - n.x].y) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) div += (quantity[idx + n.x].y - quantity[idx].y) / h.y;
		if (ngbrFlags[idx] & NF_NNY) div += (quantity[idx].y - quantity[idx - n.x].y) / h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (quantity[idx + n.x].y - quantity[idx + (n.y - 1)*n.x].y) / (2 * h.y);
			}
			else {

				div += 0.5 * (quantity[idx + n.x].y - quantity[idx].y) / h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (quantity[idx - (n.y - 1)*n.x].y - quantity[idx - n.x].y) / (2 * h.y);
			}
			else {

				div += 0.5 * (quantity[idx].y - quantity[idx - n.x].y) / h.y;
			}
		}
	}

	//z direction
	if (ngbrFlags[idx] & NF_BOTHZ) {

		div += (quantity[idx + n.x*n.y].z - quantity[idx - n.x*n.y].z) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) div += (quantity[idx + n.x*n.y].z - quantity[idx].z) / h.z;
		if (ngbrFlags[idx] & NF_NNZ) div += (quantity[idx].z - quantity[idx - n.x*n.y].z) / h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			div += 0.5 * (quantity[idx + n.x*n.y].z - quantity[idx].z) / h.z;
		}
		else {

			div += 0.5 * (quantity[idx].z - quantity[idx - n.x*n.y].z) / h.z;
		}
	}
	
	return div;
}

//divergence operator. Use non-homogeneous Neumann boundary conditions.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
//Can be used at composite media boundaries where sided differentials will be used instead.
//div operator can be applied if VType is a VAL3<Type>, returning Type
template double VEC_VC<DBL3>::div_nneu(int idx, VAL3<DBL3>& bdiff) const;

template <typename VType>
double VEC_VC<VType>::div_nneu(int idx, VAL3<VType>& bdiff) const
{
	double div = 0.0;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if (ngbrFlags[idx] & NF_BOTHX) {

		//inner point along this direction
		div += (quantity[idx + 1].x - quantity[idx - 1].x) / (2 * h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) div += (quantity[idx + 1].x - quantity[idx].x) / h.x;
		if (ngbrFlags[idx] & NF_NNX) div += (quantity[idx].x - quantity[idx - 1].x) / h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (quantity[idx + 1].x - quantity[idx + n.x - 1].x) / (2 * h.x);
			}
			else {

				div += 0.5 * ((quantity[idx + 1].x - quantity[idx].x) / h.x + bdiff.x.x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (quantity[idx - (n.x - 1)].x - quantity[idx - 1].x) / (2 * h.x);
			}
			else {

				div += 0.5 * ((quantity[idx].x - quantity[idx - 1].x) / h.x + bdiff.x.x);
			}
		}
	}

	//y direction
	if (ngbrFlags[idx] & NF_BOTHY) {

		div += (quantity[idx + n.x].y - quantity[idx - n.x].y) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) div += (quantity[idx + n.x].y - quantity[idx].y) / h.y;
		if (ngbrFlags[idx] & NF_NNY) div += (quantity[idx].y - quantity[idx - n.x].y) / h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (quantity[idx + n.x].y - quantity[idx + (n.y - 1)*n.x].y) / (2 * h.y);
			}
			else {

				div += 0.5 * ((quantity[idx + n.x].y - quantity[idx].y) / h.y + bdiff.y.y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (quantity[idx - (n.y - 1)*n.x].y - quantity[idx - n.x].y) / (2 * h.y);
			}
			else {

				div += 0.5 * ((quantity[idx].y - quantity[idx - n.x].y) / h.y + bdiff.y.y);
			}
		}
	}

	//z direction
	if (ngbrFlags[idx] & NF_BOTHZ) {

		div += (quantity[idx + n.x*n.y].z - quantity[idx - n.x*n.y].z) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) div += (quantity[idx + n.x*n.y].z - quantity[idx].z) / h.z;
		if (ngbrFlags[idx] & NF_NNZ) div += (quantity[idx].z - quantity[idx - n.x*n.y].z) / h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			div += 0.5 * ((quantity[idx + n.x*n.y].z - quantity[idx].z) / h.z + bdiff.z.z);
		}
		else {

			div += 0.5 * ((quantity[idx].z - quantity[idx - n.x*n.y].z) / h.z + bdiff.z.z);
		}
	}

	return div;
}

//divergence operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
//div operator can be applied if VType is a VAL3<Type>, returning Type
template double VEC_VC<DBL3>::div_diri(int idx) const;

template <typename VType>
double VEC_VC<VType>::div_diri(int idx) const
{
	double div = 0.0;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if (ngbrFlags[idx] & NF_BOTHX) {

		//inner point along this direction
		div += (quantity[idx + 1].x - quantity[idx - 1].x) / (2 * h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (ngbrFlags[idx] & NF_DIRICHLETX) {

		if (ngbrFlags[idx] & NF_DIRICHLETPX) div += (quantity[idx + 1].x + quantity[idx].x - 2 * get_dirichlet_value(NF_DIRICHLETPX, idx).x) / (2 * h.x);
		else								 div += (2 * get_dirichlet_value(NF_DIRICHLETNX, idx).x - quantity[idx].x - quantity[idx - 1].x) / (2 * h.x);
	}
	//Not Dirichlet, is it a CMBND boundary? - if not this either then use homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) div += (quantity[idx + 1].x - quantity[idx].x) / h.x;
		if (ngbrFlags[idx] & NF_NNX) div += (quantity[idx].x - quantity[idx - 1].x) / h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (quantity[idx + 1].x - quantity[idx + n.x - 1].x) / (2 * h.x);
			}
			else {

				div += 0.5 * (quantity[idx + 1].x - quantity[idx].x) / h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (quantity[idx - (n.x - 1)].x - quantity[idx - 1].x) / (2 * h.x);
			}
			else {

				div += 0.5 * (quantity[idx].x - quantity[idx - 1].x) / h.x;
			}
		}
	}

	//y direction
	if (ngbrFlags[idx] & NF_BOTHY) {

		div += (quantity[idx + n.x].y - quantity[idx - n.x].y) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETY) {

		if (ngbrFlags[idx] & NF_DIRICHLETPY) div += (quantity[idx + n.x].y + quantity[idx].y - 2 * get_dirichlet_value(NF_DIRICHLETPY, idx).y) / (2 * h.y);
		else								 div += (2 * get_dirichlet_value(NF_DIRICHLETNY, idx).y - quantity[idx].y - quantity[idx - n.x].y) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) div += (quantity[idx + n.x].y - quantity[idx].y) / h.y;
		if (ngbrFlags[idx] & NF_NNY) div += (quantity[idx].y - quantity[idx - n.x].y) / h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (quantity[idx + n.x].y - quantity[idx + (n.y - 1)*n.x].y) / (2 * h.y);
			}
			else {

				div += 0.5 * (quantity[idx + n.x].y - quantity[idx].y) / h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (quantity[idx - (n.y - 1)*n.x].y - quantity[idx - n.x].y) / (2 * h.y);
			}
			else {

				div += 0.5 * (quantity[idx].y - quantity[idx - n.x].y) / h.y;
			}
		}
	}

	//z direction
	if (ngbrFlags[idx] & NF_BOTHZ) {

		div += (quantity[idx + n.x*n.y].z - quantity[idx - n.x*n.y].z) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

		if (ngbrFlags[idx] & NF_DIRICHLETPZ) div += (quantity[idx + n.x*n.y].z + quantity[idx].z - 2 * get_dirichlet_value(NF_DIRICHLETPZ, idx).z) / (2 * h.z);
		else								 div += (2 * get_dirichlet_value(NF_DIRICHLETNZ, idx).z - quantity[idx].z - quantity[idx - n.x*n.y].z) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) div += (quantity[idx + n.x*n.y].z - quantity[idx].z) / h.z;
		if (ngbrFlags[idx] & NF_NNZ) div += (quantity[idx].z - quantity[idx - n.x*n.y].z) / h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			div += 0.5 * (quantity[idx + n.x*n.y].z - quantity[idx].z) / h.z;
		}
		else {

			div += 0.5 * (quantity[idx].z - quantity[idx - n.x*n.y].z) / h.z;
		}
	}

	return div;
}

//divergence operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
//Can be used at composite media boundaries where sided differentials will be used instead.
//div operator can be applied if VType is a VAL3<Type>, returning Type
template double VEC_VC<DBL3>::div_diri_nneu(int idx, VAL3<DBL3>& bdiff) const;

template <typename VType>
double VEC_VC<VType>::div_diri_nneu(int idx, VAL3<VType>& bdiff) const
{
	double div = 0.0;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if (ngbrFlags[idx] & NF_BOTHX) {

		//inner point along this direction
		div += (quantity[idx + 1].x - quantity[idx - 1].x) / (2 * h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (ngbrFlags[idx] & NF_DIRICHLETX) {

		if (ngbrFlags[idx] & NF_DIRICHLETPX) div += (quantity[idx + 1].x + quantity[idx].x - 2 * get_dirichlet_value(NF_DIRICHLETPX, idx).x) / (2 * h.x);
		else								 div += (2 * get_dirichlet_value(NF_DIRICHLETNX, idx).x - quantity[idx].x - quantity[idx - 1].x) / (2 * h.x);
	}
	//Not Dirichlet, is it a CMBND boundary? - if not this either then use homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) div += (quantity[idx + 1].x - quantity[idx].x) / h.x;
		if (ngbrFlags[idx] & NF_NNX) div += (quantity[idx].x - quantity[idx - 1].x) / h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (quantity[idx + 1].x - quantity[idx + n.x - 1].x) / (2 * h.x);
			}
			else {

				div += 0.5 * ((quantity[idx + 1].x - quantity[idx].x) / h.x + bdiff.x.x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (quantity[idx - (n.x - 1)].x - quantity[idx - 1].x) / (2 * h.x);
			}
			else {

				div += 0.5 * ((quantity[idx].x - quantity[idx - 1].x) / h.x + bdiff.x.x);
			}
		}
	}

	//y direction
	if (ngbrFlags[idx] & NF_BOTHY) {

		div += (quantity[idx + n.x].y - quantity[idx - n.x].y) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETY) {

		if (ngbrFlags[idx] & NF_DIRICHLETPY) div += (quantity[idx + n.x].y + quantity[idx].y - 2 * get_dirichlet_value(NF_DIRICHLETPY, idx).y) / (2 * h.y);
		else								 div += (2 * get_dirichlet_value(NF_DIRICHLETNY, idx).y - quantity[idx].y - quantity[idx - n.x].y) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) div += (quantity[idx + n.x].y - quantity[idx].y) / h.y;
		if (ngbrFlags[idx] & NF_NNY) div += (quantity[idx].y - quantity[idx - n.x].y) / h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (quantity[idx + n.x].y - quantity[idx + (n.y - 1)*n.x].y) / (2 * h.y);
			}
			else {

				div += 0.5 * ((quantity[idx + n.x].y - quantity[idx].y) / h.y + bdiff.y.y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (quantity[idx - (n.y - 1)*n.x].y - quantity[idx - n.x].y) / (2 * h.y);
			}
			else {

				div += 0.5 * ((quantity[idx].y - quantity[idx - n.x].y) / h.y + bdiff.y.y);
			}
		}
	}

	//z direction
	if (ngbrFlags[idx] & NF_BOTHZ) {

		div += (quantity[idx + n.x*n.y].z - quantity[idx - n.x*n.y].z) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

		if (ngbrFlags[idx] & NF_DIRICHLETPZ) div += (quantity[idx + n.x*n.y].z + quantity[idx].z - 2 * get_dirichlet_value(NF_DIRICHLETPZ, idx).z) / (2 * h.z);
		else								 div += (2 * get_dirichlet_value(NF_DIRICHLETNZ, idx).z - quantity[idx].z - quantity[idx - n.x*n.y].z) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) div += (quantity[idx + n.x*n.y].z - quantity[idx].z) / h.z;
		if (ngbrFlags[idx] & NF_NNZ) div += (quantity[idx].z - quantity[idx - n.x*n.y].z) / h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			div += 0.5 * ((quantity[idx + n.x*n.y].z - quantity[idx].z) / h.z + bdiff.z.z);
		}
		else {

			div += 0.5 * ((quantity[idx].z - quantity[idx - n.x*n.y].z) / h.z + bdiff.z.z);
		}
	}

	return div;
}

//divergence operator. Use sided differentials (also at composite media boundaries)
//div operator can be applied if VType is a VAL3<Type>, returning Type
template double VEC_VC<DBL3>::div_sided(int idx) const;

template <typename VType>
double VEC_VC<VType>::div_sided(int idx) const
{
	double div = 0.0;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if (ngbrFlags[idx] & NF_BOTHX) {

		//inner point along this direction
		div += (quantity[idx + 1].x - quantity[idx - 1].x) / (2 * h.x);
	}
	//use sided differentials if one of the neighbors is present
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (quantity[idx + 1].x - quantity[idx + n.x - 1].x) / (2 * h.x);
			}
			else {

				div += (quantity[idx + 1].x - quantity[idx].x) / h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (quantity[idx - (n.x - 1)].x - quantity[idx - 1].x) / (2 * h.x);
			}
			else {

				div += (quantity[idx].x - quantity[idx - 1].x) / h.x;
			}
		}
	}

	//y direction
	if (ngbrFlags[idx] & NF_BOTHY) {

		div += (quantity[idx + n.x].y - quantity[idx - n.x].y) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (quantity[idx + n.x].y - quantity[idx + (n.y - 1)*n.x].y) / (2 * h.y);
			}
			else {

				div += (quantity[idx + n.x].y - quantity[idx].y) / h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (quantity[idx - (n.y - 1)*n.x].y - quantity[idx - n.x].y) / (2 * h.y);
			}
			else {

				div += (quantity[idx].y - quantity[idx - n.x].y) / h.y;
			}
		}
	}

	//z direction
	if (ngbrFlags[idx] & NF_BOTHZ) {

		div += (quantity[idx + n.x*n.y].z - quantity[idx - n.x*n.y].z) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			div += (quantity[idx + n.x*n.y].z - quantity[idx].z) / h.z;
		}
		else {

			div += (quantity[idx].z - quantity[idx - n.x*n.y].z) / h.z;
		}
	}

	return div;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

//-------------------------------- DIVERGENCE OPERATOR applied after multiplying with unit antisymmetric tensor (epsilon3)

//divergence operator of epsilon3(quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
//Use Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
template DBL3 VEC_VC<DBL3>::diveps3_neu(int idx) const;
template FLT3 VEC_VC<FLT3>::diveps3_neu(int idx) const;

template <typename VType>
VType VEC_VC<VType>::diveps3_neu(int idx) const
{
	//differentials along x, y, z directions
	VType diffx, diffy, diffz;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x direction
	if (ngbrFlags[idx] & NF_BOTHX) {

		//inner point along this direction
		diffx = (quantity[idx + 1] - quantity[idx - 1]) / (2 * h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) diffx = (quantity[idx + 1] - quantity[idx]) / h.x;
		if (ngbrFlags[idx] & NF_NNX) diffx = (quantity[idx] - quantity[idx - 1]) / h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx = (quantity[idx + 1] - quantity[idx + n.x - 1]) / (2 * h.x);
			}
			else {

				diffx = 0.5 * (quantity[idx + 1] - quantity[idx]) / h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx += (quantity[idx - (n.x - 1)] - quantity[idx - 1]) / (2 * h.x);
			}
			else {

				diffx = 0.5 * (quantity[idx] - quantity[idx - 1]) / h.x;
			}
		}
	}

	//y direction
	if (ngbrFlags[idx] & NF_BOTHY) {

		diffy = (quantity[idx + n.x] - quantity[idx - n.x]) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) diffy = (quantity[idx + n.x] - quantity[idx]) / h.y;
		if (ngbrFlags[idx] & NF_NNY) diffy = (quantity[idx] - quantity[idx - n.x]) / h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy += (quantity[idx + n.x] - quantity[idx + (n.y - 1)*n.x]) / (2 * h.y);
			}
			else {

				diffy = 0.5 * (quantity[idx + n.x] - quantity[idx]) / h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy += (quantity[idx - (n.y - 1)*n.x] - quantity[idx - n.x]) / (2 * h.y);
			}
			else {

				diffy = 0.5 * (quantity[idx] - quantity[idx - n.x]) / h.y;
			}
		}
	}

	//z direction
	if (ngbrFlags[idx] & NF_BOTHZ) {

		diffz = (quantity[idx + n.x*n.y] - quantity[idx - n.x*n.y]) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) diffz = (quantity[idx + n.x*n.y] - quantity[idx]) / h.z;
		if (ngbrFlags[idx] & NF_NNZ) diffz = (quantity[idx] - quantity[idx - n.x*n.y]) / h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			diffz = 0.5 * (quantity[idx + n.x*n.y] - quantity[idx]) / h.z;
		}
		else {

			diffz = 0.5 * (quantity[idx] - quantity[idx - n.x*n.y]) / h.z;
		}
	}

	return VType(
		-diffy.z + diffz.y,
		-diffz.x + diffx.z,
		-diffx.y + diffy.x
		);
}

//divergence operator of epsilon3(quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
//Use Dirichlet conditions if set, else Neumann boundary conditions(homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
template DBL3 VEC_VC<DBL3>::diveps3_diri(int idx) const;
template FLT3 VEC_VC<FLT3>::diveps3_diri(int idx) const;

template <typename VType>
VType VEC_VC<VType>::diveps3_diri(int idx) const
{
	//differentials along x, y, z directions
	VType diffx, diffy, diffz;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x direction
	if (ngbrFlags[idx] & NF_BOTHX) {

		//inner point along this direction
		diffx = (quantity[idx + 1] - quantity[idx - 1]) / (2 * h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (ngbrFlags[idx] & NF_DIRICHLETX) {

		if (ngbrFlags[idx] & NF_DIRICHLETPX) diffx = (quantity[idx + 1] + quantity[idx] - 2 * get_dirichlet_value(NF_DIRICHLETPX, idx)) / (2 * h.x);
		else								 diffx = (2 * get_dirichlet_value(NF_DIRICHLETNX, idx) - quantity[idx] - quantity[idx - 1]) / (2 * h.x);
	}
	//Not Dirichlet, is it a CMBND boundary? - if not this either then use homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) diffx = (quantity[idx + 1] - quantity[idx]) / h.x;
		if (ngbrFlags[idx] & NF_NNX) diffx = (quantity[idx] - quantity[idx - 1]) / h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx += (quantity[idx + 1] - quantity[idx + n.x - 1]) / (2 * h.x);
			}
			else {

				diffx = 0.5 * (quantity[idx + 1] - quantity[idx]) / h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx += (quantity[idx - (n.x - 1)] - quantity[idx - 1]) / (2 * h.x);
			}
			else {

				diffx = 0.5 * (quantity[idx] - quantity[idx - 1]) / h.x;
			}
		}
	}

	//y direction
	if (ngbrFlags[idx] & NF_BOTHY) {

		diffy = (quantity[idx + n.x] - quantity[idx - n.x]) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETY) {

		if (ngbrFlags[idx] & NF_DIRICHLETPY) diffy = (quantity[idx + n.x] + quantity[idx] - 2 * get_dirichlet_value(NF_DIRICHLETPY, idx)) / (2 * h.y);
		else								 diffy = (2 * get_dirichlet_value(NF_DIRICHLETNY, idx) - quantity[idx] - quantity[idx - n.x]) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) diffy = (quantity[idx + n.x] - quantity[idx]) / h.y;
		if (ngbrFlags[idx] & NF_NNY) diffy = (quantity[idx] - quantity[idx - n.x]) / h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy += (quantity[idx + n.x] - quantity[idx + (n.y - 1)*n.x]) / (2 * h.y);
			}
			else {

				diffy = 0.5 * (quantity[idx + n.x] - quantity[idx]) / h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy += (quantity[idx - (n.y - 1)*n.x] - quantity[idx - n.x]) / (2 * h.y);
			}
			else {

				diffy = 0.5 * (quantity[idx] - quantity[idx - n.x]) / h.y;
			}
		}
	}

	//z direction
	if (ngbrFlags[idx] & NF_BOTHZ) {

		diffz = (quantity[idx + n.x*n.y] - quantity[idx - n.x*n.y]) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

		if (ngbrFlags[idx] & NF_DIRICHLETPZ) diffz = (quantity[idx + n.x*n.y] + quantity[idx] - 2 * get_dirichlet_value(NF_DIRICHLETPZ, idx)) / (2 * h.z);
		else								 diffz = (2 * get_dirichlet_value(NF_DIRICHLETNZ, idx) - quantity[idx] - quantity[idx - n.x*n.y]) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) diffz = (quantity[idx + n.x*n.y] - quantity[idx]) / h.z;
		if (ngbrFlags[idx] & NF_NNZ) diffz = (quantity[idx] - quantity[idx - n.x*n.y]) / h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			diffz = 0.5 * (quantity[idx + n.x*n.y] - quantity[idx]) / h.z;
		}
		else {

			diffz = 0.5 * (quantity[idx] - quantity[idx - n.x*n.y]) / h.z;
		}
	}

	return VType(
		-diffy.z + diffz.y,
		-diffz.x + diffx.z,
		-diffx.y + diffy.x
	);
}

//divergence operator of epsilon3(quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
//Use sided differentials (also at composite media boundaries)
template DBL3 VEC_VC<DBL3>::diveps3_sided(int idx) const;
template FLT3 VEC_VC<FLT3>::diveps3_sided(int idx) const;

template <typename VType>
VType VEC_VC<VType>::diveps3_sided(int idx) const
{
	//differentials along x, y, z directions
	VType diffx, diffy, diffz;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x direction
	if (ngbrFlags[idx] & NF_BOTHX) {

		//inner point along this direction
		diffx = (quantity[idx + 1] - quantity[idx - 1]) / (2 * h.x);
	}
	//use sided differentials if one of the neighbors is present
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx += (quantity[idx + 1] - quantity[idx + n.x - 1]) / (2 * h.x);
			}
			else {

				diffx = (quantity[idx + 1] - quantity[idx]) / h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx += (quantity[idx - (n.x - 1)] - quantity[idx - 1]) / (2 * h.x);
			}
			else {

				diffx = (quantity[idx] - quantity[idx - 1]) / h.x;
			}
		}
	}

	//y direction
	if (ngbrFlags[idx] & NF_BOTHY) {

		diffy = (quantity[idx + n.x] - quantity[idx - n.x]) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy += (quantity[idx + n.x] - quantity[idx + (n.y - 1)*n.x]) / (2 * h.y);
			}
			else {

				diffy = (quantity[idx + n.x] - quantity[idx]) / h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy += (quantity[idx - (n.y - 1)*n.x] - quantity[idx - n.x]) / (2 * h.y);
			}
			else {

				diffy = (quantity[idx] - quantity[idx - n.x]) / h.y;
			}
		}
	}

	//z direction
	if (ngbrFlags[idx] & NF_BOTHZ) {

		diffz = (quantity[idx + n.x*n.y] - quantity[idx - n.x*n.y]) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			diffz = (quantity[idx + n.x*n.y] - quantity[idx]) / h.z;
		}
		else {

			diffz = (quantity[idx] - quantity[idx - n.x*n.y]) / h.z;
		}
	}

	return VType(
		-diffy.z + diffz.y,
		-diffz.x + diffx.z,
		-diffx.y + diffy.x
	);
}