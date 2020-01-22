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

	//need to check for DIRICHLET flags which are held in the extended ngbrFlags (may be empty if not set)
	bool using_extended_flags = ngbrFlags2.size();

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
				if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

					total_weight += 2 * w_x;
					weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
				}
				else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

					total_weight += 6 * w_x;

					if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

						weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
					}
					else {

						weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
					}
				}
				else if (ngbrFlags[idx] & NF_NGBRX) {

					total_weight += w_x;

					if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * quantity[idx + 1];
					else						 weighted_sum += w_x * quantity[idx - 1];
				}

				//y direction
				if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

					total_weight += 2 * w_y;
					weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
				}
				else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

					total_weight += 6 * w_y;

					if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

						weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
					}
					else {

						weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
					}
				}
				else if (ngbrFlags[idx] & NF_NGBRY) {

					total_weight += w_y;

					if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * quantity[idx + n.x];
					else						 weighted_sum += w_y * quantity[idx - n.x];
				}

				//z direction
				if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

					total_weight += 2 * w_z;
					weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
				}
				else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

					total_weight += 6 * w_z;

					if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

						weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
					}
					else {

						weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
					}
				}
				else if (ngbrFlags[idx] & NF_NGBRZ) {

					total_weight += w_z;

					if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * quantity[idx + n.x*n.y];
					else						 weighted_sum += w_z * quantity[idx - n.x*n.y];
				}

				//advance using SOR equation
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

	//need to check for DIRICHLET flags which are held in the extended ngbrFlags (may be empty if not set)
	bool using_extended_flags = ngbrFlags2.size();

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
				if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

					total_weight += 2 * w_x;
					weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
				}
				else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

					total_weight += 6 * w_x;

					if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

						weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
					}
					else {

						weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
					}
				}
				else if (ngbrFlags[idx] & NF_NGBRX) {

					total_weight += w_x;

					if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * quantity[idx + 1];
					else						 weighted_sum += w_x * quantity[idx - 1];
				}

				//y direction
				if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

					total_weight += 2 * w_y;
					weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
				}
				else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

					total_weight += 6 * w_y;

					if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

						weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
					}
					else {

						weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
					}
				}
				else if (ngbrFlags[idx] & NF_NGBRY) {

					total_weight += w_y;

					if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * quantity[idx + n.x];
					else						 weighted_sum += w_y * quantity[idx - n.x];
				}

				//z direction
				if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

					total_weight += 2 * w_z;
					weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
				}
				else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

					total_weight += 6 * w_z;

					if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

						weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
					}
					else {

						weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
					}
				}
				else if (ngbrFlags[idx] & NF_NGBRZ) {

					total_weight += w_z;

					if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * quantity[idx + n.x*n.y];
					else						 weighted_sum += w_z * quantity[idx - n.x*n.y];
				}

				//advance using SOR equation
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

	//need to check for DIRICHLET flags which are held in the extended ngbrFlags (may be empty if not set)
	bool using_extended_flags = ngbrFlags2.size();

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
				if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

					total_weight += 2 * w_x;
					weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
				}
				else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

					total_weight += 6 * w_x;

					if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

						weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
					}
					else {

						weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
					}
				}
				else if (ngbrFlags[idx] & NF_NGBRX) {

					total_weight += w_x;

					if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * quantity[idx + 1];
					else						 weighted_sum += w_x * quantity[idx - 1];
				}

				//y direction
				if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

					total_weight += 2 * w_y;
					weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
				}
				else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

					total_weight += 6 * w_y;

					if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

						weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
					}
					else {

						weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
					}
				}
				else if (ngbrFlags[idx] & NF_NGBRY) {

					total_weight += w_y;

					if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * quantity[idx + n.x];
					else						 weighted_sum += w_y * quantity[idx - n.x];
				}

				//z direction
				if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

					total_weight += 2 * w_z;
					weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
				}
				else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

					total_weight += 6 * w_z;

					if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

						weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
					}
					else {

						weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
					}
				}
				else if (ngbrFlags[idx] & NF_NGBRZ) {

					total_weight += w_z;

					if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * quantity[idx + n.x*n.y];
					else						 weighted_sum += w_z * quantity[idx - n.x*n.y];
				}

				MType tensor = total_weight * ident<MType>() + h_max_sq * Tensor_RHS(instance, idx);

				//advance using SOR equation
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

	//need to check for DIRICHLET flags which are held in the extended ngbrFlags (may be empty if not set)
	bool using_extended_flags = ngbrFlags2.size();

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
				if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

					total_weight += 2 * w_x;
					weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
				}
				else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

					total_weight += 6 * w_x;

					if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

						weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
					}
					else {

						weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
					}
				}
				else if (ngbrFlags[idx] & NF_NGBRX) {

					total_weight += w_x;

					if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * (quantity[idx + 1] - bdiff(instance, idx).x * h.x);
					else						 weighted_sum += w_x * (quantity[idx - 1] + bdiff(instance, idx).x * h.x);
				}

				//y direction
				if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

					total_weight += 2 * w_y;
					weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
				}
				else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

					total_weight += 6 * w_y;

					if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

						weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
					}
					else {

						weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
					}
				}
				else if (ngbrFlags[idx] & NF_NGBRY) {

					total_weight += w_y;

					if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * (quantity[idx + n.x] - bdiff(instance, idx).y * h.y);
					else						 weighted_sum += w_y * (quantity[idx - n.x] + bdiff(instance, idx).y * h.y);
				}
				
				//z direction
				if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

					total_weight += 2 * w_z;
					weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
				}
				else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

					total_weight += 6 * w_z;

					if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

						weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
					}
					else {

						weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
					}
				}
				else if (ngbrFlags[idx] & NF_NGBRZ) {

					total_weight += w_z;

					if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * (quantity[idx + n.x*n.y] - bdiff(instance, idx).z * h.z);
					else						 weighted_sum += w_z * (quantity[idx - n.x*n.y] + bdiff(instance, idx).z * h.z);
				}

				//advance using SOR equation
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

	//need to check for DIRICHLET flags which are held in the extended ngbrFlags (may be empty if not set)
	bool using_extended_flags = ngbrFlags2.size();

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
				if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

					total_weight += 2 * w_x;
					weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
				}
				else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

					total_weight += 6 * w_x;

					if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

						weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
					}
					else {

						weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
					}
				}
				else if (ngbrFlags[idx] & NF_NGBRX) {

					total_weight += w_x;

					if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * (quantity[idx + 1] - bdiff(instance, idx).x * h.x);
					else						 weighted_sum += w_x * (quantity[idx - 1] + bdiff(instance, idx).x * h.x);
				}

				//y direction
				if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

					total_weight += 2 * w_y;
					weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
				}
				else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

					total_weight += 6 * w_y;

					if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

						weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
					}
					else {

						weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
					}
				}
				else if (ngbrFlags[idx] & NF_NGBRY) {

					total_weight += w_y;

					if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * (quantity[idx + n.x] - bdiff(instance, idx).y * h.y);
					else						 weighted_sum += w_y * (quantity[idx - n.x] + bdiff(instance, idx).y * h.y);
				}

				//z direction
				if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

					total_weight += 2 * w_z;
					weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
				}
				else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

					total_weight += 6 * w_z;

					if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

						weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
					}
					else {

						weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
					}
				}
				else if (ngbrFlags[idx] & NF_NGBRZ) {

					total_weight += w_z;

					if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * (quantity[idx + n.x*n.y] - bdiff(instance, idx).z * h.z);
					else						 weighted_sum += w_z * (quantity[idx - n.x*n.y] + bdiff(instance, idx).z * h.z);
				}

				MType tensor = total_weight * ident<MType>() + h_max_sq * Tensor_RHS(instance, idx);

				//advance using SOR equation
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
