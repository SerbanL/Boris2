#pragma once

#include "VEC_VC.h"
#include "VEC_VC_cmbnd.h"

//Testing CG Solver - not ready to be used properly, just solving Laplace equation in a single mesh.
//Cannot use with multiple meshes with CMBND conditions since this doesn't form a system where for Ax = b, A is a symmetric and positive definite matrix (the CMBND conditions makes A non-symmetric).
//Might be worth having a go at bi-conjugate gradient solver since this doesn't require A to be symmetric.

//-------------------------------- Conjugate Gradient Solver

//CG Solver flow:

//The problem to solve is represented by a linear system as A x = b
//x is the data in the VEC_VC, A is the operator matrix to apply to it, which must be positive definite (this is the case for Laplace and Poisson equations)
//b contains all data that is not multiplied by an x (the constants).
//For example with Laplace equation, if there are no Dirichlet boundary conditions then A is just the Laplace operator and b is zero - trivial problem of course.
//Non-trivial problem: If we have Dirichlet boundary conditions then A is obtained from the Laplace operator with Dirichlet boundary conditions, but all Dirichlet values (together with their multipliers) are moved to the RHS, i.e. into b, as they are just constants.

//NOTE : b - Ax, is simply the negative of the delsq operator (with any required boundary conditions) applied to x at a certain point in space.

//Then we have:

//1. Prepare CG Solver so we have starting values for p and r : p0, r0.
//
// Method 1: 
// p0 = r0 = b - A x0, where x0 is the current data in the VEC_VC

//2. Calculate coefficient alpha as:

// alpha = (r_T * r) / (pT * A p)
// Here r_T is the transpose of r, thus r_T * r is just the sum of squares for data in r.

//3. Update x values as:

// x_new = x_current + alpha * p

//4a. Update r values as (using updated x values):

// r_new = b - A x = r_current - alpha * A p

//4b. Calculate beta values using r_current and r_new as:

// beta = (r_new_T * r_new) / (r_T * r)

//5. Update p values as (using updated r values):

//p_new = r + beta * p_current

template <typename VType>
class CGSolve {

	friend VEC_VC<VType>;

private:

	//The VEC_VC holding this CG Solver
	VEC_VC<VType>* pV;

	//auxiliary data used by CG solver - see description above
	std::vector<VType> p, r;

	//r_T * r and p_T * Ap values calculated at the end of each iteration so it can be used for the next
	double rT_r, pT_Ap;

	//is the cg solver primed and ready to iterate?
	bool primed = false;

private:

	//Apply A operator to p for a Laplace equation and:
	//Use Dirichlet conditions if set, else Neumann boundary conditions(homogeneous).
	//Returns zero at composite media boundary cells.
	//All boundary info stored in *pV, p is just used to store values in a simple std::vector
	VType A_p_Laplace_diri(int idx);

	VType b_minus_A(int idx);

	//Prepare solver to iterate - called if primed flag is false
	void PrimeSolver(void);

public:

	CGSolve(VEC_VC<VType>* pV_) :
		pV(pV_)
	{
		p.resize(pV->n.dim());
		r.resize(pV->n.dim());
	}

	//----LAPLACE EQUATION with Dirichlet boundary conditions where set, and homogeneous Neumann boundary conditions.

	DBL2 IterateSolver(void);
};


//-------------------------------------------------------------------------------------------------

//Apply A operator to p for a Laplace equation and:
//Use Dirichlet conditions if set, else Neumann boundary conditions(homogeneous).
//Returns zero at composite media boundary cells.
//All boundary info stored in *pV, p is just used to store values in a simple std::vector
template <typename VType>
VType CGSolve<VType>::A_p_Laplace_diri(int idx)
{
	DBL3& h = pV->h;
	SZ3& n = pV->n;

	VType diff_x = VType(), diff_y = VType(), diff_z = VType();

	if (!(pV->ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	bool using_extended_flags = pV->ngbrFlags2.size();

	//x axis
	if ((pV->ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//both neighbors available so must be an inner point along this direction
		diff_x = (p[idx + 1] + p[idx - 1] - 2 * p[idx]) / (h.x*h.x);
	}
	else if (pV->ngbrFlags[idx] & NF_CMBNDX) {

		diff_x = 0;
	}
	else if (pV->ngbrFlags[idx] & NF_NGBRX) {

		//only one neighbor available. Does it use a fixed boundary value (Dirichlet)?
		if (using_extended_flags && (pV->ngbrFlags2[idx] & NF2_DIRICHLETX)) {

			if (pV->ngbrFlags2[idx] & NF2_DIRICHLETPX) diff_x = (2 * p[idx + 1] - 6 * p[idx]) / (h.x*h.x);
			else									 diff_x = (2 * p[idx - 1] - 6 * p[idx]) / (h.x*h.x);
		}
		//stencil at surface for homogeneous Neumann boundary condition - correct, do not multiply by 2! 
		//See conservation law approach derivation for heat equation (integrate then use mid-point rule approximations) for cell-centered grids.
		else if (pV->ngbrFlags[idx] & NF_NPX) diff_x = (p[idx + 1] - p[idx]) / (h.x*h.x);
		else								  diff_x = (p[idx - 1] - p[idx]) / (h.x*h.x);
	}

	//y axis
	if ((pV->ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff_y = (p[idx + n.x] + p[idx - n.x] - 2 * p[idx]) / (h.y*h.y);
	}
	else if (pV->ngbrFlags[idx] & NF_CMBNDY) {

		diff_y = 0;
	}
	else if (pV->ngbrFlags[idx] & NF_NGBRY) {

		if (using_extended_flags && (pV->ngbrFlags2[idx] & NF2_DIRICHLETY)) {

			if (pV->ngbrFlags2[idx] & NF2_DIRICHLETPY) diff_y = (2 * p[idx + n.x] - 6 * p[idx]) / (h.y*h.y);
			else									 diff_y = (2 * p[idx - n.x] - 6 * p[idx]) / (h.y*h.y);
		}
		else if (pV->ngbrFlags[idx] & NF_NPY) diff_y = (p[idx + n.x] - p[idx]) / (h.y*h.y);
		else								  diff_y = (p[idx - n.x] - p[idx]) / (h.y*h.y);
	}

	//z axis
	if ((pV->ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff_z = (p[idx + n.x*n.y] + p[idx - n.x*n.y] - 2 * p[idx]) / (h.z*h.z);
	}
	else if (pV->ngbrFlags[idx] & NF_CMBNDZ) {

		diff_z = 0;
	}
	else if (pV->ngbrFlags[idx] & NF_NGBRZ) {

		if (using_extended_flags && (pV->ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

			if (pV->ngbrFlags2[idx] & NF2_DIRICHLETPZ) diff_z = (2 * p[idx + n.x*n.y] - 6 * p[idx]) / (h.z*h.z);
			else									 diff_z = (2 * p[idx - n.x*n.y] - 6 * p[idx]) / (h.z*h.z);
		}
		else if (pV->ngbrFlags[idx] & NF_NPZ) diff_z = (p[idx + n.x*n.y] - p[idx]) / (h.z*h.z);
		else								  diff_z = (p[idx - n.x*n.y] - p[idx]) / (h.z*h.z);
	}

	return (diff_x + diff_y + diff_z);
}

template <typename VType>
VType CGSolve<VType>::b_minus_A(int idx)
{
	DBL3& h = pV->h;
	SZ3& n = pV->n;

	VType diff_x = VType(), diff_y = VType(), diff_z = VType();

	if (!(pV->ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	bool using_extended_flags = pV->ngbrFlags2.size();

	//x axis
	if ((pV->ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {
		//both neighbors available so must be an inner point along this direction
		diff_x = (pV->quantity[idx + 1] + pV->quantity[idx - 1] - 2 * pV->quantity[idx]);
	}
	else if (pV->ngbrFlags[idx] & NF_CMBNDX) {

		diff_x = 0;
	}
	else if (pV->ngbrFlags[idx] & NF_NGBRX) {

		//only one neighbor available. Does it use a fixed boundary value (Dirichlet)?
		if (using_extended_flags && (pV->ngbrFlags2[idx] & NF2_DIRICHLETX)) {

			if (pV->ngbrFlags2[idx] & NF2_DIRICHLETPX) diff_x = (2 * pV->quantity[idx + 1] + 4 * pV->get_dirichlet_value(NF2_DIRICHLETPX, idx) - 6 * pV->quantity[idx]);
			else									 diff_x = (2 * pV->quantity[idx - 1] + 4 * pV->get_dirichlet_value(NF2_DIRICHLETNX, idx) - 6 * pV->quantity[idx]);
		}
		//stencil at surface for homogeneous Neumann boundary condition - correct, do not multiply by 2! 
		//See conservation law approach derivation for heat equation (integrate then use mid-point rule approximations) for cell-centered grids.
		else if (pV->ngbrFlags[idx] & NF_NPX) diff_x = (pV->quantity[idx + 1] - pV->quantity[idx]);
		else								  diff_x = (pV->quantity[idx - 1] - pV->quantity[idx]);
	}

	diff_x /= (h.x*h.x);

	//y axis
	if ((pV->ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diff_y = (pV->quantity[idx + n.x] + pV->quantity[idx - n.x] - 2 * pV->quantity[idx]);
	}
	else if (pV->ngbrFlags[idx] & NF_CMBNDY) {

		diff_y = 0;
	}
	else if (pV->ngbrFlags[idx] & NF_NGBRY) {

		if (using_extended_flags && (pV->ngbrFlags2[idx] & NF2_DIRICHLETY)) {

			if (pV->ngbrFlags2[idx] & NF2_DIRICHLETPY) diff_y = (2 * pV->quantity[idx + n.x] + 4 * pV->get_dirichlet_value(NF2_DIRICHLETPY, idx) - 6 * pV->quantity[idx]);
			else									 diff_y = (2 * pV->quantity[idx - n.x] + 4 * pV->get_dirichlet_value(NF2_DIRICHLETNY, idx) - 6 * pV->quantity[idx]);
		}
		else if (pV->ngbrFlags[idx] & NF_NPY) diff_y = (pV->quantity[idx + n.x] - pV->quantity[idx]);
		else								  diff_y = (pV->quantity[idx - n.x] - pV->quantity[idx]);
	}

	diff_y /= (h.y*h.y);

	//z axis
	if ((pV->ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diff_z = (pV->quantity[idx + n.x*n.y] + pV->quantity[idx - n.x*n.y] - 2 * pV->quantity[idx]);
	}
	else if (pV->ngbrFlags[idx] & NF_CMBNDZ) {

		diff_z = 0;
	}
	else if (pV->ngbrFlags[idx] & NF_NGBRZ) {

		if (using_extended_flags && (pV->ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

			if (pV->ngbrFlags2[idx] & NF2_DIRICHLETPZ) diff_z = (2 * pV->quantity[idx + n.x*n.y] + 4 * pV->get_dirichlet_value(NF2_DIRICHLETPZ, idx) - 6 * pV->quantity[idx]);
			else									 diff_z = (2 * pV->quantity[idx - n.x*n.y] + 4 * pV->get_dirichlet_value(NF2_DIRICHLETNZ, idx) - 6 * pV->quantity[idx]);
		}
		else if (pV->ngbrFlags[idx] & NF_NPZ) diff_z = (pV->quantity[idx + n.x*n.y] - pV->quantity[idx]);
		else								  diff_z = (pV->quantity[idx - n.x*n.y] - pV->quantity[idx]);
	}

	diff_z /= (h.z*h.z);

	return -(diff_x + diff_y + diff_z);
}

template <typename VType>
void CGSolve<VType>::PrimeSolver(void)
{
	// Prepare CG Solver so we have starting values for p and r : p0, r0.
	//
	// p0 = r0 = b - A x0, where x0 is the current data in the VEC_VC

	//if called first time then p and r will not have the right size: allocate memory
	if (p.size() != pV->n.dim()) {

		p.resize(pV->n.dim());
		r.resize(pV->n.dim());
	}

#pragma omp parallel for
	for (int idx = 0; idx < pV->n.dim(); idx++) {

		//this is b - Ax
		//VType value = -1 * pV->delsq_diri(idx);
		VType value = b_minus_A(idx);

		p[idx] = value;
		r[idx] = value;
	}

	//find pT_Ap and rT_r:

	double rT_r_new = 0.0;
	double pT_Ap_new = 0.0;

#pragma omp parallel for reduction(+:rT_r_new, pT_Ap_new)
	for (int idx = 0; idx < pV->n.dim(); idx++) {

		rT_r_new += r[idx] * r[idx];
		pT_Ap_new += p[idx] * A_p_Laplace_diri(idx);
	}

	rT_r = rT_r_new;
	pT_Ap = pT_Ap_new;

	primed = true;
}

template <typename VType>
DBL2 CGSolve<VType>::IterateSolver(void)
{
	if (!primed || p.size() != pV->n.dim()) PrimeSolver();

	pV->magnitude_reduction.new_minmax_reduction();
	pV->magnitude_reduction2.new_minmax_reduction();

	//rT_r and pT_Ap must be available from previous iteration
	double alpha = rT_r / pT_Ap;

	//3. Update x values
#pragma omp parallel for
	for (int idx = 0; idx < pV->n.dim(); idx++) {

		if (!(pV->ngbrFlags[idx] & NF_NOTEMPTY)) continue;

		//3. Update x values
		VType old_value = pV->quantity[idx];

		pV->quantity[idx] += alpha * p[idx];

		pV->magnitude_reduction.reduce_max(GetMagnitude(old_value - pV->quantity[idx]));
		pV->magnitude_reduction2.reduce_max(GetMagnitude(pV->quantity[idx]));
	}

	//Update r values (using updated x values):
#pragma omp parallel for
	for (int idx = 0; idx < pV->n.dim(); idx++) {

		//Update r values (using updated x values):
		r[idx] = b_minus_A(idx);
	}

	//Calculate beta value:

	//calculate new rT_r value
	double rT_r_new = 0.0;

#pragma omp parallel for reduction(+:rT_r_new)
	for (int idx = 0; idx < pV->n.dim(); idx++) {

		rT_r_new += r[idx] * r[idx];
	}

	double beta = rT_r_new / rT_r;

	rT_r = rT_r_new;

	//Update p values (using updated r values):

#pragma omp parallel for
	for (int idx = 0; idx < pV->n.dim(); idx++) {

		p[idx] = r[idx] + beta * p[idx];
	}

	//update pT_Ap

	double pT_Ap_new = 0.0;

#pragma omp parallel for reduction(+:pT_Ap_new)
	for (int idx = 0; idx < pV->n.dim(); idx++) {

		pT_Ap_new += p[idx] * A_p_Laplace_diri(idx);
	}

	pT_Ap = pT_Ap_new;

	return DBL2(pV->magnitude_reduction.maximum(), pV->magnitude_reduction2.maximum());
}