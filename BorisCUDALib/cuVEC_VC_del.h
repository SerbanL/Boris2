#pragma once

#include "cuVEC_VC.h"

//-------------------------------- LAPLACE OPERATOR

//calculate Laplace operator at cell with given index. Use Neumann boundary conditions (homogeneous).
//Returns zero at composite media boundary cells.
template <typename VType>
__device__ VType cuVEC_VC<VType>::delsq_neu(int idx) const
{
	VType diff_x = VType(), diff_y = VType(), diff_z = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x axis
	if (ngbrFlags[idx] & NF_BOTHX) {
		//both neighbors available so must be an inner point along this direction
		diff_x = (quantity[idx + 1] + quantity[idx - 1] - 2 * quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff_x = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) diff_x = (quantity[idx + 1] - quantity[idx]);
		else						 diff_x = (quantity[idx - 1] - quantity[idx]);
	}

	diff_x /= (h.x*h.x);

	//y axis
	if (ngbrFlags[idx] & NF_BOTHY) {

		diff_y = (quantity[idx + n.x] + quantity[idx - n.x] - 2 * quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff_y = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) diff_y = (quantity[idx + n.x] - quantity[idx]);
		else						 diff_y = (quantity[idx - n.x] - quantity[idx]);
	}

	diff_y /= (h.y*h.y);

	//z axis
	if (ngbrFlags[idx] & NF_BOTHZ) {

		diff_z = (quantity[idx + n.x*n.y] + quantity[idx - n.x*n.y] - 2 * quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff_z = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) diff_z = (quantity[idx + n.x*n.y] - quantity[idx]);
		else						 diff_z = (quantity[idx - n.x*n.y] - quantity[idx]);
	}

	diff_z /= (h.z*h.z);

	return (diff_x + diff_y + diff_z);
}

//calculate Laplace operator at cell with given index. Use non-homogeneous Neumann boundary conditions with the specified boundary differential.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
//Returns zero at composite media boundary cells.
template <typename VType>
template <typename Class_BDiff>
__device__ VType cuVEC_VC<VType>::delsq_nneu(int idx, Class_BDiff& bdiff_class) const
{
	VType diff_x = VType(), diff_y = VType(), diff_z = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x axis
	if (ngbrFlags[idx] & NF_BOTHX) {
		//both neighbors available so must be an inner point along this direction
		diff_x = (quantity[idx + 1] + quantity[idx - 1] - 2 * quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff_x = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPX) diff_x = (quantity[idx + 1] - quantity[idx]) - bdiff_val.x * h.x;
		else						 diff_x = (quantity[idx - 1] - quantity[idx]) + bdiff_val.x * h.x;
	}

	diff_x /= (h.x*h.x);

	//y axis
	if (ngbrFlags[idx] & NF_BOTHY) {

		diff_y = (quantity[idx + n.x] + quantity[idx - n.x] - 2 * quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff_y = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPY) diff_y = (quantity[idx + n.x] - quantity[idx]) - bdiff_val.y * h.y;
		else						 diff_y = (quantity[idx - n.x] - quantity[idx]) + bdiff_val.y * h.y;
	}

	diff_y /= (h.y*h.y);

	//z axis
	if (ngbrFlags[idx] & NF_BOTHZ) {

		diff_z = (quantity[idx + n.x*n.y] + quantity[idx - n.x*n.y] - 2 * quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff_z = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPZ) diff_z = (quantity[idx + n.x*n.y] - quantity[idx]) - bdiff_val.z * h.z;
		else						 diff_z = (quantity[idx - n.x*n.y] - quantity[idx]) + bdiff_val.z * h.z;
	}

	diff_z /= (h.z*h.z);

	return (diff_x + diff_y + diff_z);
}

//calculate Laplace operator at cell with given index. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
//Returns zero at composite media boundary cells.
template <typename VType>
__device__ VType cuVEC_VC<VType>::delsq_diri(int idx) const
{
	VType diff_x = VType(), diff_y = VType(), diff_z = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x axis
	if (ngbrFlags[idx] & NF_BOTHX) {
		//both neighbors available so must be an inner point along this direction
		diff_x = (quantity[idx + 1] + quantity[idx - 1] - 2 * quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff_x = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		//only one neighbor available. Does it use a fixed boundary value (Dirichlet)?
		if (ngbrFlags[idx] & NF_DIRICHLETX) {

			if (ngbrFlags[idx] & NF_DIRICHLETPX) diff_x = (2 * quantity[idx + 1] + 4 * get_dirichlet_value(NF_DIRICHLETPX, idx) - 6 * quantity[idx]);
			else								 diff_x = (2 * quantity[idx - 1] + 4 * get_dirichlet_value(NF_DIRICHLETNX, idx) - 6 * quantity[idx]);
		}
		//stencil at surface for homogeneous Neumann boundary condition - correct, do not multiply by 2! 
		//See conservation law approach derivation for heat equation (integrate then use mid-point rule approximations) for cell-centered grids.
		else if (ngbrFlags[idx] & NF_NPX) diff_x = (quantity[idx + 1] - quantity[idx]);
		else							  diff_x = (quantity[idx - 1] - quantity[idx]);
	}


	diff_x /= (h.x*h.x);

	//y axis
	if (ngbrFlags[idx] & NF_BOTHY) {

		diff_y = (quantity[idx + n.x] + quantity[idx - n.x] - 2 * quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff_y = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_DIRICHLETY) {

			if (ngbrFlags[idx] & NF_DIRICHLETPY) diff_y = (2 * quantity[idx + n.x] + 4 * get_dirichlet_value(NF_DIRICHLETPY, idx) - 6 * quantity[idx]);
			else								 diff_y = (2 * quantity[idx - n.x] + 4 * get_dirichlet_value(NF_DIRICHLETNY, idx) - 6 * quantity[idx]);
		}
		else if (ngbrFlags[idx] & NF_NPY) diff_y = (quantity[idx + n.x] - quantity[idx]);
		else							  diff_y = (quantity[idx - n.x] - quantity[idx]);
	}

	diff_y /= (h.y*h.y);

	//z axis
	if (ngbrFlags[idx] & NF_BOTHZ) {

		diff_z = (quantity[idx + n.x*n.y] + quantity[idx - n.x*n.y] - 2 * quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff_z = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_DIRICHLETZ) {

			if (ngbrFlags[idx] & NF_DIRICHLETPZ) diff_z = (2 * quantity[idx + n.x*n.y] + 4 * get_dirichlet_value(NF_DIRICHLETPZ, idx) - 6 * quantity[idx]);
			else								 diff_z = (2 * quantity[idx - n.x*n.y] + 4 * get_dirichlet_value(NF_DIRICHLETNZ, idx) - 6 * quantity[idx]);
		}
		else if (ngbrFlags[idx] & NF_NPZ) diff_z = (quantity[idx + n.x*n.y] - quantity[idx]);
		else							  diff_z = (quantity[idx - n.x*n.y] - quantity[idx]);
	}

	diff_z /= (h.z*h.z);

	return (diff_x + diff_y + diff_z);
}

//calculate Laplace operator at cell with given index. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
//Returns zero at composite media boundary cells.
template <typename VType>
template <typename Class_BDiff>
__device__ VType cuVEC_VC<VType>::delsq_diri_nneu(int idx, Class_BDiff& bdiff_class) const
{
	VType diff_x = VType(), diff_y = VType(), diff_z = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x axis
	if (ngbrFlags[idx] & NF_BOTHX) {
		//both neighbors available so must be an inner point along this direction
		diff_x = (quantity[idx + 1] + quantity[idx - 1] - 2 * quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff_x = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		//only one neighbor available. Does it use a fixed boundary value (Dirichlet)?
		if (ngbrFlags[idx] & NF_DIRICHLETX) {

			if (ngbrFlags[idx] & NF_DIRICHLETPX) diff_x = (2 * quantity[idx + 1] + 4 * get_dirichlet_value(NF_DIRICHLETPX, idx) - 6 * quantity[idx]);
			else								 diff_x = (2 * quantity[idx - 1] + 4 * get_dirichlet_value(NF_DIRICHLETNX, idx) - 6 * quantity[idx]);
		}
		else {

			cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

			if (ngbrFlags[idx] & NF_NPX) diff_x = (quantity[idx + 1] - quantity[idx]) - bdiff_val.x * h.x;
			else						 diff_x = (quantity[idx - 1] - quantity[idx]) + bdiff_val.x * h.x;
		}
	}


	diff_x /= (h.x*h.x);

	//y axis
	if (ngbrFlags[idx] & NF_BOTHY) {

		diff_y = (quantity[idx + n.x] + quantity[idx - n.x] - 2 * quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff_y = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_DIRICHLETY) {

			if (ngbrFlags[idx] & NF_DIRICHLETPY) diff_y = (2 * quantity[idx + n.x] + 4 * get_dirichlet_value(NF_DIRICHLETPY, idx) - 6 * quantity[idx]);
			else								 diff_y = (2 * quantity[idx - n.x] + 4 * get_dirichlet_value(NF_DIRICHLETNY, idx) - 6 * quantity[idx]);
		}
		else {

			cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

			if (ngbrFlags[idx] & NF_NPY) diff_y = (quantity[idx + n.x] - quantity[idx]) - bdiff_val.y * h.y;
			else						 diff_y = (quantity[idx - n.x] - quantity[idx]) + bdiff_val.y * h.y;
		}
	}

	diff_y /= (h.y*h.y);

	//z axis
	if (ngbrFlags[idx] & NF_BOTHZ) {

		diff_z = (quantity[idx + n.x*n.y] + quantity[idx - n.x*n.y] - 2 * quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff_z = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_DIRICHLETZ) {

			if (ngbrFlags[idx] & NF_DIRICHLETPZ) diff_z = (2 * quantity[idx + n.x*n.y] + 4 * get_dirichlet_value(NF_DIRICHLETPZ, idx) - 6 * quantity[idx]);
			else								 diff_z = (2 * quantity[idx - n.x*n.y] + 4 * get_dirichlet_value(NF_DIRICHLETNZ, idx) - 6 * quantity[idx]);
		}
		else {

			cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

			if (ngbrFlags[idx] & NF_NPZ) diff_z = (quantity[idx + n.x*n.y] - quantity[idx]) - bdiff_val.z * h.z;
			else						 diff_z = (quantity[idx - n.x*n.y] - quantity[idx]) + bdiff_val.z * h.z;
		}
	}

	diff_z /= (h.z*h.z);

	return (diff_x + diff_y + diff_z);
}

//calculate Laplace operator at cell with given index. Use Robin boundary conditions, defaulting to Neumann if not set.
//Returns zero at composite media boundary cells.
template <typename VType>
__device__ VType cuVEC_VC<VType>::delsq_robin(int idx, cuReal K) const
{
	VType diff_x = VType(), diff_y = VType(), diff_z = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x axis
	if (ngbrFlags[idx] & NF_BOTHX) {
		//both neighbors available so must be an inner point along this direction
		diff_x = (quantity[idx + 1] + quantity[idx - 1] - 2 * quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		diff_x = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		//only one neighbor available. Does it use a Robin condition?
		if (ngbrFlags[idx] & NF_ROBINX) {

			//Robin cell on +x side of boundary
			if (ngbrFlags[idx] & NF_ROBINPX) {

				//use void Robin values
				if (ngbrFlags[idx] & NF_ROBINV) {

					diff_x = ((1 + robin_v.i * h.x / (2 * K)) * quantity[idx + 1] + robin_v.i * h.x * robin_v.j / K - (1 + 3 * robin_v.i * h.x / (2 * K)) * quantity[idx]);
				}
				//don't use void Robin values
				else {

					diff_x = ((1 + robin_nx.i * h.x / (2 * K)) * quantity[idx + 1] + robin_nx.i * h.x * robin_nx.j / K - (1 + 3 * robin_nx.i * h.x / (2 * K)) * quantity[idx]);
				}
			}
			//Robin cell on -x side of boundary
			else {

				if (ngbrFlags[idx] & NF_ROBINV) {

					diff_x = ((1 + robin_v.i * h.x / (2 * K)) * quantity[idx - 1] + robin_v.i * h.x * robin_v.j / K - (1 + 3 * robin_v.i * h.x / (2 * K)) * quantity[idx]);
				}
				else {

					diff_x = ((1 + robin_px.i * h.x / (2 * K)) * quantity[idx - 1] + robin_px.i * h.x * robin_px.j / K - (1 + 3 * robin_px.i * h.x / (2 * K)) * quantity[idx]);
				}
			}
		}
		//stencil at surface for homogeneous Neumann boundary condition - correct, do not multiply by 2! 
		//See conservation law approach derivation for heat equation (integrate then use mid-point rule approximations) for cell-centered grids.
		else if (ngbrFlags[idx] & NF_NPX) diff_x = (quantity[idx + 1] - quantity[idx]);
		else							  diff_x = (quantity[idx - 1] - quantity[idx]);
	}


	diff_x /= (h.x*h.x);

	//y axis
	if (ngbrFlags[idx] & NF_BOTHY) {

		diff_y = (quantity[idx + n.x] + quantity[idx - n.x] - 2 * quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		diff_y = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		//only one neighbor available. Does it use a Robin condition?
		if (ngbrFlags[idx] & NF_ROBINY) {

			//Robin cell on +x side of boundary
			if (ngbrFlags[idx] & NF_ROBINPY) {

				//use void Robin values
				if (ngbrFlags[idx] & NF_ROBINV) {

					diff_y = ((1 + robin_v.i * h.y / (2 * K)) * quantity[idx + n.x] + robin_v.i * h.y * robin_v.j / K - (1 + 3 * robin_v.i * h.y / (2 * K)) * quantity[idx]);
				}
				//don't use void Robin values
				else {

					diff_y = ((1 + robin_ny.i * h.y / (2 * K)) * quantity[idx + n.x] + robin_ny.i * h.y * robin_ny.j / K - (1 + 3 * robin_ny.i * h.y / (2 * K)) * quantity[idx]);
				}
			}
			//Robin cell on -x side of boundary
			else {

				if (ngbrFlags[idx] & NF_ROBINV) {

					diff_y = ((1 + robin_v.i * h.y / (2 * K)) * quantity[idx - n.x] + robin_v.i * h.y * robin_v.j / K - (1 + 3 * robin_v.i * h.y / (2 * K)) * quantity[idx]);
				}
				else {

					diff_y = ((1 + robin_py.i * h.y / (2 * K)) * quantity[idx - n.x] + robin_py.i * h.y * robin_py.j / K - (1 + 3 * robin_py.i * h.y / (2 * K)) * quantity[idx]);
				}
			}
		}
		else if (ngbrFlags[idx] & NF_NPY) diff_y = (quantity[idx + n.x] - quantity[idx]);
		else							  diff_y = (quantity[idx - n.x] - quantity[idx]);
	}

	diff_y /= (h.y*h.y);

	//z axis
	if (ngbrFlags[idx] & NF_BOTHZ) {

		diff_z = (quantity[idx + n.x*n.y] + quantity[idx - n.x*n.y] - 2 * quantity[idx]);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		diff_z = 0;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		//only one neighbor available. Does it use a Robin condition?
		if (ngbrFlags[idx] & NF_ROBINZ) {

			//Robin cell on +x side of boundary
			if (ngbrFlags[idx] & NF_ROBINPZ) {

				//use void Robin values
				if (ngbrFlags[idx] & NF_ROBINV) {

					diff_z = ((1 + robin_v.i * h.z / (2 * K)) * quantity[idx + n.x*n.y] + robin_v.i * h.z * robin_v.j / K - (1 + 3 * robin_v.i * h.z / (2 * K)) * quantity[idx]);
				}
				//don't use void Robin values
				else {

					diff_z = ((1 + robin_nz.i * h.z / (2 * K)) * quantity[idx + n.x*n.y] + robin_nz.i * h.z * robin_nz.j / K - (1 + 3 * robin_nz.i * h.z / (2 * K)) * quantity[idx]);
				}
			}
			//Robin cell on -x side of boundary
			else {

				if (ngbrFlags[idx] & NF_ROBINV) {

					diff_z = ((1 + robin_v.i * h.z / (2 * K)) * quantity[idx - n.x*n.y] + robin_v.i * h.z * robin_v.j / K - (1 + 3 * robin_v.i * h.z / (2 * K)) * quantity[idx]);
				}
				else {

					diff_z = ((1 + robin_pz.i * h.z / (2 * K)) * quantity[idx - n.x*n.y] + robin_pz.i * h.z * robin_pz.j / K - (1 + 3 * robin_pz.i * h.z / (2 * K)) * quantity[idx]);
				}
			}
		}
		else if (ngbrFlags[idx] & NF_NPZ) diff_z = (quantity[idx + n.x*n.y] - quantity[idx]);
		else							  diff_z = (quantity[idx - n.x*n.y] - quantity[idx]);
	}

	diff_z /= (h.z*h.z);

	return (diff_x + diff_y + diff_z);
}

