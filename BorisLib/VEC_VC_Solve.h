#pragma once

#include "VEC_VC.h"

//-------------------------------- LAPLACE EQUATION

template <typename VType>
DBL2 VEC_VC<VType>::IterateLaplace_SOR(double relaxation_param)
{
	//get maximum cell side
	double h_max = maximum(h.x, h.y, h.z);

	//get weights
	double w_x = (h_max / h.x) * (h_max / h.x);
	double w_y = (h_max / h.y) * (h_max / h.y);
	double w_z = (h_max / h.z) * (h_max / h.z);

	magnitude_reduction.new_minmax_reduction();
	magnitude_reduction2.new_minmax_reduction();

	//red-black : two passes will be taken
	int rb = 0;
	while (rb < 2) {

		#pragma omp parallel for
		for (int idx_jk = 0; idx_jk < n.y * n.z; idx_jk++) {

			int j = idx_jk % n.y;
			int k = (idx_jk / n.y) % n.z;

			//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
			bool red_nudge = (((j % 2) == 1 && (k % 2) == 0) || (((j % 2) == 0 && (k % 2) == 1)));

			//For red pass (first) i starts from red_nudge. For black pass (second) i starts from !red_nudge.
			for (int i = (1 - rb) * red_nudge + rb * (!red_nudge); i < n.x; i += 2) {

				int idx = i + j * n.x + k * n.x*n.y;

				//calculate new value only in non-empty; also skip if indicated as a composite media boundary condition cell
				if ((ngbrFlags[idx] & NF_CMBND) || !(ngbrFlags[idx] & NF_NOTEMPTY)) continue;

				VType weighted_sum = VType();
				double total_weight = 0;

				//x direction
				if (ngbrFlags[idx] & NF_BOTHX) {

					total_weight += 2 * w_x;
					weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
				}
				else if (ngbrFlags[idx] & NF_DIRICHLETX) {

					total_weight += 6 * w_x;

					if (ngbrFlags[idx] & NF_DIRICHLETPX) weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
					else								 weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
				}
				else if (ngbrFlags[idx] & NF_NGBRX) {

					total_weight += w_x;

					if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * quantity[idx + 1];
					else						 weighted_sum += w_x * quantity[idx - 1];
				}

				//y direction
				if (ngbrFlags[idx] & NF_BOTHY) {

					total_weight += 2 * w_y;
					weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
				}
				else if (ngbrFlags[idx] & NF_DIRICHLETY) {

					total_weight += 6 * w_y;

					if (ngbrFlags[idx] & NF_DIRICHLETPY) weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
					else								 weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
				}
				else if (ngbrFlags[idx] & NF_NGBRY) {

					total_weight += w_y;

					if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * quantity[idx + n.x];
					else						 weighted_sum += w_y * quantity[idx - n.x];
				}

				//z direction
				if (ngbrFlags[idx] & NF_BOTHZ) {

					total_weight += 2 * w_z;
					weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
				}
				else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

					total_weight += 6 * w_z;

					if (ngbrFlags[idx] & NF_DIRICHLETPZ) weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
					else								 weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
				}
				else if (ngbrFlags[idx] & NF_NGBRZ) {

					total_weight += w_z;

					if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * quantity[idx + n.x*n.y];
					else						 weighted_sum += w_z * quantity[idx - n.x*n.y];
				}

				//advance using SOR formula
				VType old_value = quantity[idx];
				quantity[idx] = quantity[idx] * (1 - relaxation_param) + relaxation_param * (weighted_sum / total_weight);

				magnitude_reduction.reduce_max(GetMagnitude(old_value - quantity[idx]));
				magnitude_reduction2.reduce_max(GetMagnitude(quantity[idx]));
			}
		}

		rb++;
	}

	return DBL2(magnitude_reduction.maximum(), magnitude_reduction2.maximum());
}

//-------------------------------- POISSON EQUATION

template <typename VType>
template <typename Owner>
DBL2 VEC_VC<VType>::IteratePoisson_SOR(std::function<VType(const Owner&, int)> Poisson_RHS, Owner& instance, double relaxation_param)
{
	//get maximum cell side
	double h_max_sq = maximum(h.x, h.y, h.z);
	h_max_sq *= h_max_sq;

	//get weights
	double w_x = h_max_sq / (h.x*h.x);
	double w_y = h_max_sq / (h.y*h.y);
	double w_z = h_max_sq / (h.z*h.z);

	magnitude_reduction.new_minmax_reduction();
	magnitude_reduction2.new_minmax_reduction();

	//red-black : two passes will be taken
	int rb = 0;
	while (rb < 2) {

		#pragma omp parallel for
		for(int idx_jk = 0; idx_jk < n.y * n.z; idx_jk++) {

			int j = idx_jk % n.y;
			int k = (idx_jk / n.y) % n.z;

			//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
			bool red_nudge = (((j % 2) == 1 && (k % 2) == 0) || (((j % 2) == 0 && (k % 2) == 1)));

			//For red pass (first) i starts from red_nudge. For black pass (second) i starts from !red_nudge.
			for (int i = (1 - rb) * red_nudge + rb * (!red_nudge); i < n.x; i += 2) {

				int idx = i + j * n.x + k * n.x*n.y;

				//calculate new value only in non-empty cells with non-fixed values; also skip if indicated as a composite media boundary condition cell
				if ((ngbrFlags[idx] & NF_CMBND) || !(ngbrFlags[idx] & NF_NOTEMPTY)) continue;

				VType weighted_sum = VType(0);
				double total_weight = 0;

				//x direction
				if (ngbrFlags[idx] & NF_BOTHX) {

					total_weight += 2 * w_x;
					weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
				}
				else if (ngbrFlags[idx] & NF_DIRICHLETX) {

					total_weight += 6 * w_x;

					if (ngbrFlags[idx] & NF_DIRICHLETPX) weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
					else								 weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
				}
				else if (ngbrFlags[idx] & NF_NGBRX) {

					total_weight += w_x;

					if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * quantity[idx + 1];
					else						 weighted_sum += w_x * quantity[idx - 1];
				}

				//y direction
				if (ngbrFlags[idx] & NF_BOTHY) {

					total_weight += 2 * w_y;
					weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
				}
				else if (ngbrFlags[idx] & NF_DIRICHLETY) {

					total_weight += 6 * w_y;

					if (ngbrFlags[idx] & NF_DIRICHLETPY) weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
					else								 weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
				}
				else if (ngbrFlags[idx] & NF_NGBRY) {

					total_weight += w_y;

					if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * quantity[idx + n.x];
					else						 weighted_sum += w_y * quantity[idx - n.x];
				}

				//z direction
				if (ngbrFlags[idx] & NF_BOTHZ) {

					total_weight += 2 * w_z;
					weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
				}
				else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

					total_weight += 6 * w_z;

					if (ngbrFlags[idx] & NF_DIRICHLETPZ) weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
					else								 weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
				}
				else if (ngbrFlags[idx] & NF_NGBRZ) {

					total_weight += w_z;

					if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * quantity[idx + n.x*n.y];
					else						 weighted_sum += w_z * quantity[idx - n.x*n.y];
				}

				//advance using SOR formula
				VType old_value = quantity[idx];
				quantity[idx] = quantity[idx] * (1 - relaxation_param) + relaxation_param * ((weighted_sum - h_max_sq*Poisson_RHS(instance, idx)) / total_weight);

				magnitude_reduction.reduce_max(GetMagnitude(old_value - quantity[idx]));
				magnitude_reduction2.reduce_max(GetMagnitude(quantity[idx]));
			}
		}

		rb++;
	}

	return DBL2(magnitude_reduction.maximum(), magnitude_reduction2.maximum());
}

//This solves delsq V = F + M * V : For M use Tensor_RHS (For VType double M returns type double, For VType DBL3 M returns DBL33)
//For Poisson equation we need a function to specify the RHS of the equation delsq V = F : use Poisson_RHS
//F must be a member const method of Owner taking an index value (the index ranges over this VEC) and returning a double value : F(index) evaluated at the index-th cell.
//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
//Dirichlet boundary conditions used, defaulting to Neumann boundary conditions where not set, and composite media boundary cells skipped (use boundary conditions to set cmbnd cells after calling this)
template <typename VType>
template <typename Owner, typename MType>
DBL2 VEC_VC<VType>::IteratePoisson_SOR(std::function<VType(const Owner&, int)> Poisson_RHS, std::function<MType(const Owner&, int)> Tensor_RHS, Owner& instance, double relaxation_param = 1.9)
{
	//get maximum cell side
	double h_max_sq = maximum(h.x, h.y, h.z);
	h_max_sq *= h_max_sq;

	//get weights
	double w_x = h_max_sq / (h.x*h.x);
	double w_y = h_max_sq / (h.y*h.y);
	double w_z = h_max_sq / (h.z*h.z);

	magnitude_reduction.new_minmax_reduction();
	magnitude_reduction2.new_minmax_reduction();

	//red-black : two passes will be taken
	int rb = 0;
	while (rb < 2) {

#pragma omp parallel for
		for (int idx_jk = 0; idx_jk < n.y * n.z; idx_jk++) {

			int j = idx_jk % n.y;
			int k = (idx_jk / n.y) % n.z;

			//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
			bool red_nudge = (((j % 2) == 1 && (k % 2) == 0) || (((j % 2) == 0 && (k % 2) == 1)));

			//For red pass (first) i starts from red_nudge. For black pass (second) i starts from !red_nudge.
			for (int i = (1 - rb) * red_nudge + rb * (!red_nudge); i < n.x; i += 2) {

				int idx = i + j * n.x + k * n.x*n.y;

				//calculate new value only in non-empty cells with non-fixed values; also skip if indicated as a composite media boundary condition cell
				if ((ngbrFlags[idx] & NF_CMBND) || !(ngbrFlags[idx] & NF_NOTEMPTY)) continue;

				VType weighted_sum = VType(0);
				double total_weight = 0;

				//x direction
				if (ngbrFlags[idx] & NF_BOTHX) {

					total_weight += 2 * w_x;
					weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
				}
				else if (ngbrFlags[idx] & NF_DIRICHLETX) {

					total_weight += 6 * w_x;

					if (ngbrFlags[idx] & NF_DIRICHLETPX) weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
					else								 weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
				}
				else if (ngbrFlags[idx] & NF_NGBRX) {

					total_weight += w_x;

					if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * quantity[idx + 1];
					else						 weighted_sum += w_x * quantity[idx - 1];
				}

				//y direction
				if (ngbrFlags[idx] & NF_BOTHY) {

					total_weight += 2 * w_y;
					weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
				}
				else if (ngbrFlags[idx] & NF_DIRICHLETY) {

					total_weight += 6 * w_y;

					if (ngbrFlags[idx] & NF_DIRICHLETPY) weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
					else								 weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
				}
				else if (ngbrFlags[idx] & NF_NGBRY) {

					total_weight += w_y;

					if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * quantity[idx + n.x];
					else						 weighted_sum += w_y * quantity[idx - n.x];
				}

				//z direction
				if (ngbrFlags[idx] & NF_BOTHZ) {

					total_weight += 2 * w_z;
					weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
				}
				else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

					total_weight += 6 * w_z;

					if (ngbrFlags[idx] & NF_DIRICHLETPZ) weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
					else								 weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
				}
				else if (ngbrFlags[idx] & NF_NGBRZ) {

					total_weight += w_z;

					if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * quantity[idx + n.x*n.y];
					else						 weighted_sum += w_z * quantity[idx - n.x*n.y];
				}

				MType tensor = total_weight * ident<MType>() + h_max_sq * Tensor_RHS(instance, idx);

				//advance using SOR formula
				VType old_value = quantity[idx];
				quantity[idx] = quantity[idx] * (1 - relaxation_param) + relaxation_param * inverse<MType>(tensor) * (weighted_sum - h_max_sq * Poisson_RHS(instance, idx));

				magnitude_reduction.reduce_max(GetMagnitude(old_value - quantity[idx]));
				magnitude_reduction2.reduce_max(GetMagnitude(quantity[idx]));
			}
		}

		rb++;
	}

	return DBL2(magnitude_reduction.maximum(), magnitude_reduction2.maximum());
}

template <typename VType>
template <typename Owner>
DBL2 VEC_VC<VType>::IteratePoisson_aSOR(std::function<VType(const Owner&, int)> Poisson_RHS, Owner& instance, bool start_iters, double err_limit)
{
	DBL2 error_parts = IteratePoisson_SOR(Poisson_RHS, instance, aSOR_damping);

	double error = (error_parts.j > 0 ? error_parts.i / error_parts.j : error_parts.i);

	//prepare start of a sequence of iterations but don't adjust damping at the start
	if (start_iters) {

		aSOR_lasterror = error;
		aSOR_lastgrad = 0;
		return error_parts;
	}

	//adjust damping
	double grad_lnerror = (log(error) - log(aSOR_lasterror)) / aSOR_damping;
	adjust_aSOR_damping(grad_lnerror, error, err_limit);

	//save parameters from this iteration
	aSOR_lasterror = error;
	aSOR_lastgrad = grad_lnerror;

	//return the current error
	return error_parts;
}

//Poisson equation solved using SOR, but using non-homogeneous Neumann boundary condition - this is evaluated using the bdiff call-back method.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
template <typename VType>
template <typename Owner>
DBL2 VEC_VC<VType>::IteratePoisson_SOR(std::function<VType(const Owner&, int)> Poisson_RHS, std::function<VAL3<VType>(const Owner&, int)> bdiff, Owner& instance, double relaxation_param)
{
	//get maximum cell side
	double h_max_sq = maximum(h.x, h.y, h.z);
	h_max_sq *= h_max_sq;

	//get weights
	double w_x = h_max_sq / (h.x*h.x);
	double w_y = h_max_sq / (h.y*h.y);
	double w_z = h_max_sq / (h.z*h.z);

	magnitude_reduction.new_minmax_reduction();
	magnitude_reduction2.new_minmax_reduction();

	//red-black : two passes will be taken
	int rb = 0;
	while (rb < 2) {

#pragma omp parallel for
		for (int idx_jk = 0; idx_jk < n.y * n.z; idx_jk++) {

			int j = idx_jk % n.y;
			int k = (idx_jk / n.y) % n.z;

			//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
			bool red_nudge = (((j % 2) == 1 && (k % 2) == 0) || (((j % 2) == 0 && (k % 2) == 1)));

			//For red pass (first) i starts from red_nudge. For black pass (second) i starts from !red_nudge.
			for (int i = (1 - rb) * red_nudge + rb * (!red_nudge); i < n.x; i += 2) {

				int idx = i + j * n.x + k * n.x*n.y;

				//calculate new value only in non-empty cells with non-fixed values; also skip if indicated as a composite media boundary condition cell
				if ((ngbrFlags[idx] & NF_CMBND) || !(ngbrFlags[idx] & NF_NOTEMPTY)) continue;

				VType weighted_sum = VType(0);
				double total_weight = 0;

				//x direction
				if (ngbrFlags[idx] & NF_BOTHX) {

					total_weight += 2 * w_x;
					weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
				}
				else if (ngbrFlags[idx] & NF_DIRICHLETX) {

					total_weight += 6 * w_x;

					if (ngbrFlags[idx] & NF_DIRICHLETPX) weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
					else								 weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
				}
				else if (ngbrFlags[idx] & NF_NGBRX) {

					total_weight += w_x;

					if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * (quantity[idx + 1] - bdiff(instance, idx).x * h.x);
					else						 weighted_sum += w_x * (quantity[idx - 1] + bdiff(instance, idx).x * h.x);
				}

				//y direction
				if (ngbrFlags[idx] & NF_BOTHY) {

					total_weight += 2 * w_y;
					weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
				}
				else if (ngbrFlags[idx] & NF_DIRICHLETY) {

					total_weight += 6 * w_y;

					if (ngbrFlags[idx] & NF_DIRICHLETPY) weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
					else								 weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
				}
				else if (ngbrFlags[idx] & NF_NGBRY) {

					total_weight += w_y;

					if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * (quantity[idx + n.x] - bdiff(instance, idx).y * h.y);
					else						 weighted_sum += w_y * (quantity[idx - n.x] + bdiff(instance, idx).y * h.y);
				}
				
				//z direction
				if (ngbrFlags[idx] & NF_BOTHZ) {

					total_weight += 2 * w_z;
					weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
				}
				else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

					total_weight += 6 * w_z;

					if (ngbrFlags[idx] & NF_DIRICHLETPZ) weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
					else								 weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
				}
				else if (ngbrFlags[idx] & NF_NGBRZ) {

					total_weight += w_z;

					if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * (quantity[idx + n.x*n.y] - bdiff(instance, idx).z * h.z);
					else						 weighted_sum += w_z * (quantity[idx - n.x*n.y] + bdiff(instance, idx).z * h.z);
				}

				//advance using SOR formula
				VType old_value = quantity[idx];
				quantity[idx] = quantity[idx] * (1 - relaxation_param) + relaxation_param * ((weighted_sum - h_max_sq * Poisson_RHS(instance, idx)) / total_weight);

				magnitude_reduction.reduce_max(GetMagnitude(old_value - quantity[idx]));
				magnitude_reduction2.reduce_max(GetMagnitude(quantity[idx]));
			}
		}

		rb++;
	}

	return DBL2(magnitude_reduction.maximum(), magnitude_reduction2.maximum());
}

//This solves delsq V = F + M * V : For M use Tensor_RHS (For VType double M returns type double, For VType DBL3 M returns DBL33)
//Poisson equation solved using SOR, but using non-homogeneous Neumann boundary condition - this is evaluated using the bdiff call-back method.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
template <typename VType>
template <typename Owner, typename MType>
DBL2 VEC_VC<VType>::IteratePoisson_SOR(std::function<VType(const Owner&, int)> Poisson_RHS, std::function<MType(const Owner&, int)> Tensor_RHS, std::function<VAL3<VType>(const Owner&, int)> bdiff, Owner& instance, double relaxation_param = 1.9)
{
	//get maximum cell side
	double h_max_sq = maximum(h.x, h.y, h.z);
	h_max_sq *= h_max_sq;

	//get weights
	double w_x = h_max_sq / (h.x*h.x);
	double w_y = h_max_sq / (h.y*h.y);
	double w_z = h_max_sq / (h.z*h.z);

	magnitude_reduction.new_minmax_reduction();
	magnitude_reduction2.new_minmax_reduction();

	//red-black : two passes will be taken
	int rb = 0;
	while (rb < 2) {

#pragma omp parallel for
		for (int idx_jk = 0; idx_jk < n.y * n.z; idx_jk++) {

			int j = idx_jk % n.y;
			int k = (idx_jk / n.y) % n.z;

			//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
			bool red_nudge = (((j % 2) == 1 && (k % 2) == 0) || (((j % 2) == 0 && (k % 2) == 1)));

			//For red pass (first) i starts from red_nudge. For black pass (second) i starts from !red_nudge.
			for (int i = (1 - rb) * red_nudge + rb * (!red_nudge); i < n.x; i += 2) {

				int idx = i + j * n.x + k * n.x*n.y;

				//calculate new value only in non-empty cells with non-fixed values; also skip if indicated as a composite media boundary condition cell
				if ((ngbrFlags[idx] & NF_CMBND) || !(ngbrFlags[idx] & NF_NOTEMPTY)) continue;

				VType weighted_sum = VType(0);
				double total_weight = 0;

				//x direction
				if (ngbrFlags[idx] & NF_BOTHX) {

					total_weight += 2 * w_x;
					weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
				}
				else if (ngbrFlags[idx] & NF_DIRICHLETX) {

					total_weight += 6 * w_x;

					if (ngbrFlags[idx] & NF_DIRICHLETPX) weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
					else								 weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
				}
				else if (ngbrFlags[idx] & NF_NGBRX) {

					total_weight += w_x;

					if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * (quantity[idx + 1] - bdiff(instance, idx).x * h.x);
					else						 weighted_sum += w_x * (quantity[idx - 1] + bdiff(instance, idx).x * h.x);
				}

				//y direction
				if (ngbrFlags[idx] & NF_BOTHY) {

					total_weight += 2 * w_y;
					weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
				}
				else if (ngbrFlags[idx] & NF_DIRICHLETY) {

					total_weight += 6 * w_y;

					if (ngbrFlags[idx] & NF_DIRICHLETPY) weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
					else								 weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
				}
				else if (ngbrFlags[idx] & NF_NGBRY) {

					total_weight += w_y;

					if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * (quantity[idx + n.x] - bdiff(instance, idx).y * h.y);
					else						 weighted_sum += w_y * (quantity[idx - n.x] + bdiff(instance, idx).y * h.y);
				}

				//z direction
				if (ngbrFlags[idx] & NF_BOTHZ) {

					total_weight += 2 * w_z;
					weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
				}
				else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

					total_weight += 6 * w_z;

					if (ngbrFlags[idx] & NF_DIRICHLETPZ) weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
					else								 weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
				}
				else if (ngbrFlags[idx] & NF_NGBRZ) {

					total_weight += w_z;

					if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * (quantity[idx + n.x*n.y] - bdiff(instance, idx).z * h.z);
					else						 weighted_sum += w_z * (quantity[idx - n.x*n.y] + bdiff(instance, idx).z * h.z);
				}

				MType tensor = total_weight * ident<MType>() + h_max_sq * Tensor_RHS(instance, idx);

				//advance using SOR formula
				VType old_value = quantity[idx];
				quantity[idx] = quantity[idx] * (1 - relaxation_param) + relaxation_param * inverse<MType>(tensor) * (weighted_sum - h_max_sq * Poisson_RHS(instance, idx));

				magnitude_reduction.reduce_max(GetMagnitude(old_value - quantity[idx]));
				magnitude_reduction2.reduce_max(GetMagnitude(quantity[idx]));
			}
		}

		rb++;
	}

	return DBL2(magnitude_reduction.maximum(), magnitude_reduction2.maximum());
}

//Poisson equation solved using adaptive SOR algorithm, using non-homogeneous Neumann boundary condition - this is evaluated using the bdiff call-back method.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
template <typename VType>
template <typename Owner>
DBL2 VEC_VC<VType>::IteratePoisson_aSOR(std::function<VType(const Owner&, int)> Poisson_RHS, std::function<VAL3<VType>(const Owner&, int)> bdiff, Owner& instance, bool start_iters, double err_limit)
{
	DBL2 error_parts = IteratePoisson_SOR(Poisson_RHS, bdiff, instance, aSOR_damping);

	double error = (error_parts.j > 0 ? error_parts.i / error_parts.j : error_parts.i);

	//prepare start of a sequence of iterations but don't adjust damping at the start
	if (start_iters) {

		aSOR_lasterror = error;
		aSOR_lastgrad = 0;
		return error_parts;
	}

	//adjust damping
	double grad_lnerror = (log(error) - log(aSOR_lasterror)) / aSOR_damping;
	adjust_aSOR_damping(grad_lnerror, error, err_limit);

	//save parameters from this iteration
	aSOR_lasterror = error;
	aSOR_lastgrad = grad_lnerror;

	//return the current error
	return error_parts;
}

//adjust aSOR damping based on current gradient of ln(error), if error > err_limit.
template <typename VType>
void VEC_VC<VType>::adjust_aSOR_damping(double grad_lnerror, double error, double err_limit)
{

#define ASOR_SPIKEVAL	0.4
#define ASOR_EXPONENT	2.1
#define ASOR_BIAS	0.02
#define ASOR_NUDGE	1.018
#define	ASOR_MINDAMPING	0.2
#define	ASOR_MAXDAMPING	2.0

	//apply full adjustment mechanism only if error is above threshold : below this cannot apply the normal mechanism due to "numerical noise"
	if (error > err_limit) {

		//positive gradient - should decrease damping
		if (grad_lnerror >= 0) {

			//avoid spikes : do nothing if simple spike detected
			if (aSOR_lastgrad <= 0 && grad_lnerror > ASOR_SPIKEVAL) return;

			//decrease damping using formula : larger g results in bigger decrease
			aSOR_damping *= exp(-grad_lnerror * ASOR_EXPONENT - ASOR_BIAS);
		}
		//negative g - might be able to do better by increasing damping
		else {

			aSOR_damping *= ASOR_NUDGE;
		}
	}
	else {

		//error is below threshold, but don't want to be stuck with a low damping value : give it a gentle increase
		aSOR_damping *= ASOR_NUDGE;
	}

	//make sure damping is within bounds
	if (aSOR_damping < ASOR_MINDAMPING) aSOR_damping = ASOR_MINDAMPING;
	if (aSOR_damping > ASOR_MAXDAMPING) aSOR_damping = ASOR_MAXDAMPING;
}
