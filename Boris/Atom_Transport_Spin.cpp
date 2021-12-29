#include "stdafx.h"
#include "Atom_Transport.h"

#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "SuperMesh.h"
#include "Atom_MeshParamsControl.h"

#if COMPILECUDA == 1
#include "Atom_TransportCUDA.h"
#endif

//-------------------Calculation Methods used by Spin Current Solver only

//---------------Charge part

DBL2 Atom_Transport::IterateSpinSolver_Charge_SOR(double damping)
{
	//TO DO
	//return paMesh->V.IteratePoisson_SOR<Atom_Transport>(&Atom_Transport::Evaluate_SpinSolver_delsqV_RHS, *this, damping);
	return DBL2();
}

//before iterating the spin solver (charge part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
void Atom_Transport::PrimeSpinSolver_Charge(void)
{
	//TO DO
	/*
	//Update dM_dt values if needed
	if (Need_dM_dt_Calculation()) {

#pragma omp parallel for
		for (int idx = 0; idx < dM_dt.linear_size(); idx++) {

			if (paMesh->M.is_not_empty(idx)) {

				dM_dt[idx] = paMesh->dMdt(idx);
			}
		}
	}

	//the rest are terms to calculate in delsq_V_fixed
	if (!Need_delsq_V_fixed_Precalculation()) return;

	bool cppgmr_enabled = IsNZ(paMesh->betaD.get0());
	bool cpump_enabled = IsNZ(paMesh->cpump_eff.get0());

#pragma omp parallel for
	for (int idx = 0; idx < delsq_V_fixed.linear_size(); idx++) {

		delsq_V_fixed[idx] = 0.0;

		if (paMesh->V.is_not_empty(idx)) {

			if (cppgmr_enabled || cpump_enabled) {

				double Ms = paMesh->Ms;
				paMesh->update_parameters_ecoarse(idx, paMesh->Ms, Ms);

				int idx_M = paMesh->M.position_to_cellidx(paMesh->V.cellidx_to_position(idx));
				DBL33 grad_m = paMesh->M.grad_neu(idx_M) / Ms;
				DBL3 m = paMesh->M[idx_M] / Ms;

				//CPP-GMR contribution
				if (cppgmr_enabled) {

					double De = paMesh->De;
					double betaD = paMesh->betaD;
					paMesh->update_parameters_ecoarse(idx, paMesh->De, De, paMesh->betaD, betaD);

					DBL33 grad_S = paMesh->S.grad_neu(idx);
					DBL3 delsq_S = paMesh->S.delsq_neu(idx);
					double div_grad_S_m = (grad_S.i * grad_m.i) + (grad_S.j * grad_m.j) + (grad_S.k * grad_m.k) + (m * delsq_S);

					delsq_V_fixed[idx] += div_grad_S_m * betaD * De / (MUB_E * paMesh->elC[idx]);
				}

				//Charge pumping pre-calculation
				if (cpump_enabled) {

					DBL33 grad_dm_dt = dM_dt.grad_neu(idx_M) / Ms;
					DBL3 dm_dt = dM_dt[idx_M] / Ms;

					double P = paMesh->P;
					paMesh->update_parameters_ecoarse(idx, paMesh->P, P);

					DBL3 dx_m = grad_m.x;
					DBL3 dy_m = grad_m.y;
					DBL3 dxx_m = paMesh->M.dxx_neu(idx_M) / Ms;
					DBL3 dyy_m = paMesh->M.dyy_neu(idx_M) / Ms;

					delsq_V_fixed[idx] += (paMesh->cpump_eff.get0() * P * HBAR_E / 2) * ((grad_dm_dt.x ^ dx_m) + (grad_dm_dt.y ^ dy_m) + (dm_dt ^ (dxx_m + dyy_m))) * m;
				}
			}
		}
	}
	*/
}

//call-back method for Poisson equation for V
//VERIFIED - CORRECT
double Atom_Transport::Evaluate_SpinSolver_delsqV_RHS(int idx) const
{
	//The Poisson solver calls this method to evaluate the RHS of this equation
	double value = 0.0;

	//TO DO
	/*
	//homogeneous Neumann boundary condition applies to V in magnetic meshes
	DBL3 grad_V = paMesh->V.grad_diri(idx);

	//1. principal term : always present
	value = -(grad_V * paMesh->elC.grad_sided(idx)) / paMesh->elC[idx];

	//2. topological Hall effect contribution
	if (IsNZ(paMesh->the_eff.get0())) {

		double Ms = paMesh->Ms;
		double P = paMesh->P;
		double n_density = paMesh->n_density;
		paMesh->update_parameters_ecoarse(idx, paMesh->Ms, Ms, paMesh->P, P, paMesh->n_density, n_density);

		int idx_M = paMesh->M.position_to_cellidx(paMesh->V.cellidx_to_position(idx));
		DBL3 m = paMesh->M[idx_M] / Ms;

		DBL33 grad_m = paMesh->M.grad_neu(idx_M) / Ms;
		DBL3 dx_m = grad_m.x;
		DBL3 dy_m = grad_m.y;
		DBL3 dxy_m = paMesh->M.dxy_neu(idx_M) / Ms;
		DBL3 dxx_m = paMesh->M.dxx_neu(idx_M) / Ms;
		DBL3 dyy_m = paMesh->M.dyy_neu(idx_M) / Ms;

		DBL3 B_the = DBL3(
			((dxy_m ^ dy_m) + (dx_m ^ dyy_m)) * m,
			-1.0 * ((dxx_m ^ dy_m) + (dx_m ^ dxy_m)) * m,
			0.0);

		value -= (paMesh->the_eff.get0() * P * paMesh->elC[idx] * HBAR_E / (ECHARGE * n_density)) * (grad_V * B_the);
	}

	//additional fixed contributions if needed (e.g. CPP-GMR and charge pumping)
	if (delsq_V_fixed.linear_size()) value += delsq_V_fixed[idx];
	*/
	return value;
}

//---------------Spin part

//solve for spin accumulation using Poisson equation for delsq_S, solved using SOR algorithm
DBL2 Atom_Transport::IterateSpinSolver_Spin_SOR(double damping)
{
	//TO DO
	//no SHE contribution. Note, SHE is not included in magnetic meshes.
	//return paMesh->S.IteratePoisson_SOR<Atom_Transport>(&Atom_Transport::Evaluate_SpinSolver_delsqS_RHS, *this, damping);
	return DBL2();
}

//before iterating the spin solver (spin part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
void Atom_Transport::PrimeSpinSolver_Spin(void)
{
	//TO DO
	/*
	//the rest are terms to calculate in delsq_S_fixed
	if (!Need_delsq_S_fixed_Precalculation()) return;

	bool cpump_enabled = IsNZ(paMesh->cpump_eff.get0());
	bool the_enabled = IsNZ(paMesh->the_eff.get0());
	bool she_enabled = IsNZ(paMesh->SHA.get0());

#pragma omp parallel for
	for (int idx = 0; idx < delsq_S_fixed.linear_size(); idx++) {

		delsq_S_fixed[idx] = DBL3();

		if (paMesh->V.is_not_empty(idx)) {

			if (stsolve == STSOLVE_FERROMAGNETIC) {

				//magnetic mesh

				double Ms = paMesh->Ms;
				double P = paMesh->P;
				double De = paMesh->De;
				paMesh->update_parameters_ecoarse(idx, paMesh->Ms, Ms, paMesh->P, P, paMesh->De, De);

				//term due to drift (non-uniformity of M term, and delsq V contribution - non-uniformity of E term)

				//find grad M and M at the M cell in which the current S cell center is
				int idx_M = paMesh->M.position_to_cellidx(paMesh->V.cellidx_to_position(idx));

				DBL3 m = paMesh->M[idx_M] / Ms;
				DBL33 grad_m = paMesh->M.grad_neu(idx_M) / Ms;
				DBL3 E_dot_del_m = grad_m | paMesh->E[idx];

				//E_dot_del_m term is very important, but Evaluate_SpinSolver_delsqV_RHS term could be neglected in most cases especially if E is uniform.
				delsq_S_fixed[idx] += (P * MUB_E * paMesh->elC[idx] / De) * (Evaluate_SpinSolver_delsqV_RHS(idx) * m - E_dot_del_m);

				//charge pumping and topological Hall effect
				if (cpump_enabled || the_enabled) {

					DBL3 dx_m = grad_m.x;
					DBL3 dy_m = grad_m.y;
					DBL3 dxy_m = paMesh->M.dxy_neu(idx_M) / Ms;
					DBL3 dxx_m = paMesh->M.dxx_neu(idx_M) / Ms;
					DBL3 dyy_m = paMesh->M.dyy_neu(idx_M) / Ms;

					if (cpump_enabled) {

						DBL3 dmdt = dM_dt[idx_M] / Ms;
						DBL33 grad_dm_dt = dM_dt.grad_neu(idx_M) / Ms;

						delsq_S_fixed[idx] += paMesh->cpump_eff.get0() * (paMesh->elC[idx] * HBAR_E * MUB_E / (2 * De)) * ((grad_dm_dt.x ^ dx_m) + (grad_dm_dt.y ^ dy_m) + (dmdt ^ (dxx_m + dyy_m)));
					}

					if (the_enabled) {

						double n_density = paMesh->n_density;
						paMesh->update_parameters_ecoarse(idx, paMesh->n_density, n_density);

						delsq_S_fixed[idx] += paMesh->the_eff.get0() * (HBAR_E * MUB_E * paMesh->elC[idx] * paMesh->elC[idx] / (ECHARGE * n_density * De)) * (paMesh->E[idx].x * ((dxy_m ^ dy_m) + (dx_m ^ dyy_m)) - paMesh->E[idx].y * ((dxx_m ^ dy_m) + (dx_m ^ dxy_m)));
					}
				}
			}

			//terms occuring only in non-magnetic meshes
			else {

				//1. SHA term (this is negligible in most cases, even if E is non-uniform, but might as well include it) 
				if (she_enabled) {

					double SHA = paMesh->SHA;
					double De = paMesh->De;
					paMesh->update_parameters_ecoarse(idx, paMesh->SHA, SHA, paMesh->De, De);

					//Check boundary conditions for this term : should be Dirichlet with 0 Jc value normal to the boundary except for electrodes.
					delsq_S_fixed[idx] += (SHA * paMesh->elC[idx] * MUB_E / De) * paMesh->E.diveps3_sided(idx);
				}
			}
		}
	}
	*/
}

//call-back method for Poisson equation for S
//VERIFIED - CORRECT
DBL3 Atom_Transport::Evaluate_SpinSolver_delsqS_RHS(int idx) const
{
	DBL3 delsq_S_RHS;
	//TO DO

	/*
	double l_sf = paMesh->l_sf;
	paMesh->update_parameters_ecoarse(idx, paMesh->l_sf, l_sf);

	//Contributions which apply equally in ferromagnetic and non-ferromagnetic meshes
	//longitudinal S decay term
	delsq_S_RHS = (paMesh->S[idx] / (l_sf * l_sf));

	//Terms occuring only in magnetic meshes
	if (stsolve == STSOLVE_FERROMAGNETIC) {

		double Ms = paMesh->Ms;
		double l_ex = paMesh->l_ex;
		double l_ph = paMesh->l_ph;
		paMesh->update_parameters_ecoarse(idx, paMesh->Ms, Ms, paMesh->l_ex, l_ex, paMesh->l_ph, l_ph);

		//find grad M and M at the M cell in which the current S cell center is
		int idx_M = paMesh->M.position_to_cellidx(paMesh->V.cellidx_to_position(idx));
		DBL3 m = paMesh->M[idx_M] / Ms;

		//transverse S decay terms
		delsq_S_RHS += ((paMesh->S[idx] ^ m) / (l_ex * l_ex) + (m ^ (paMesh->S[idx] ^ m)) / (l_ph * l_ph));
	}

	//additional fixed contributions if needed
	if (delsq_S_fixed.linear_size()) delsq_S_RHS += delsq_S_fixed[idx];
	*/
	return delsq_S_RHS;
}

//-------------------

//Non-homogeneous Neumann boundary condition for V' - call-back method for Poisson equation for V
//VERIFIED - CORRECT
DBL3 Atom_Transport::NHNeumann_Vdiff(int idx) const
{
	return DBL3();
}

//Non-homogeneous Neumann boundary condition for S' - call-back method for Poisson equation for S
//VERIFIED - CORRECT
DBL33 Atom_Transport::NHNeumann_Sdiff(int idx) const
{
	return DBL33();
}

//-------------------CMBND computation methods

//Spin transport : V and S
//CMBND for V
//flux = a + b V' at the interface, b = -sigma, a = betaD * (De*e/muB) * (grad S)m + (SHA*De*e/muB) * curl S + charge pumping + topological Hall effect
//Note, the topological Hall effect term includes E, thus gradients in V, but we can include these in the a term for 2 reasons:
//1. these CMBND functions especially with the topological Hall effect enabled is used for interfaces along z direction normally, and here Ez is zero 
//(for such interfaces both charge pumping and topological Hall effect have zero contribution to z direction charge current)
//2. even if the interface is along x or y we can still use the previously calculated E field, and the solution will converge to the same value (but might take more iterations).
//VERIFIED - CORRECT
double Atom_Transport::afunc_st_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	//TO DO

	double a = 0.0;

	/*
	DBL3 u = shift.normalized() * -1;

	bool cppgmr_enabled = IsNZ(paMesh->betaD.get0());
	bool cpump_enabled = IsNZ(paMesh->cpump_eff.get0()) && IsZ(shift.z);
	bool the_enabled = IsNZ(paMesh->the_eff.get0()) && IsZ(shift.z);

	if (stsolve == STSOLVE_FERROMAGNETIC && (cppgmr_enabled || cpump_enabled || the_enabled)) {

		//magnetic mesh

		double Ms = paMesh->Ms;
		paMesh->update_parameters_atposition(relpos_m1, paMesh->Ms, Ms);

		DBL3 m1 = paMesh->M.weighted_average(relpos_m1, stencil) / Ms;
		DBL3 m2 = paMesh->M.weighted_average(relpos_m1 + shift, stencil) / Ms;

		//1. CPP-GMR contribution
		if (cppgmr_enabled) {

			double betaD = paMesh->betaD;
			double De = paMesh->De;
			paMesh->update_parameters_atposition(relpos_m1, paMesh->betaD, betaD, paMesh->De, De);

			int idx_S1 = paMesh->S.position_to_cellidx(relpos_m1);
			int idx_S2 = paMesh->S.position_to_cellidx(relpos_m1 + shift);

			//value a1
			DBL33 grad_S1 = paMesh->S.grad_neu(idx_S1);

			double a1 = ((grad_S1 * m1) * betaD * De / MUB_E) * u;

			//value a2
			DBL33 grad_S2 = paMesh->S.grad_neu(idx_S2);

			double a2 = ((grad_S2 * m2) * betaD * De / MUB_E) * u;

			//final interpolated a value
			a += (1.5 * a1 - 0.5 * a2);
		}

		//2. Charge pumping
		//3. Topological Hall effect
		if (cpump_enabled || the_enabled) {

			double P = paMesh->P;
			double n_density = paMesh->n_density;
			paMesh->update_parameters_atposition(relpos_m1, paMesh->P, P, paMesh->n_density, n_density);

			int idx_M1 = paMesh->M.position_to_cellidx(relpos_m1);
			int idx_M2 = paMesh->M.position_to_cellidx(relpos_m1 + shift);

			DBL33 grad_m1 = paMesh->M.grad_neu(idx_M1) / Ms;
			DBL33 grad_m2 = paMesh->M.grad_neu(idx_M2) / Ms;

			int idx_V1 = paMesh->V.position_to_cellidx(relpos_m1);
			//do not read off the E field directly as it's only calculated after the spin solver (charge part) has relaxed
			DBL3 E1 = -1.0 * paMesh->V.grad_diri(idx_V1);

			int idx_V2 = paMesh->V.position_to_cellidx(relpos_m1 + shift);
			DBL3 E2 = -1.0 * paMesh->V.grad_diri(idx_V2);

			double sigma_1 = paMesh->elC.weighted_average(relpos_m1, stencil);
			double sigma_2 = paMesh->elC.weighted_average(relpos_m1 + shift, stencil);

			//topological Hall effect contribution
			if (the_enabled) {

				//value a1
				double Bz_the_1 = (grad_m1.x ^ grad_m1.y) * m1;
				double a1 = paMesh->the_eff.get0() * (-P * sigma_1 * sigma_1 * HBAR_E / (ECHARGE * n_density)) * DBL3(E1.y * Bz_the_1, -E1.x * Bz_the_1, 0.0) * u;

				//value a2
				double Bz_the_2 = (grad_m2.x ^ grad_m2.y) * m2;
				double a2 = paMesh->the_eff.get0() * (-P * sigma_2 * sigma_2 * HBAR_E / (ECHARGE * n_density)) * DBL3(E1.y * Bz_the_1, -E1.x * Bz_the_1, 0.0) * u;

				//final interpolated a value
				a += (1.5 * a1 - 0.5 * a2);
			}

			//charge pumping contribution
			if (cpump_enabled) {

				//value a1
				DBL3 dm_dt_1 = dM_dt.weighted_average(relpos_m1, stencil) / Ms;
				double a1 = paMesh->cpump_eff.get0() * (P * sigma_1 * HBAR_E / 2) * DBL3((dm_dt_1 ^ grad_m1.x) * m1, (dm_dt_1 ^ grad_m1.y) * m1, 0.0) * u;

				//value a2
				DBL3 dm_dt_2 = dM_dt.weighted_average(relpos_m1 + shift, stencil) / Ms;
				double a2 = paMesh->cpump_eff.get0() * (P * sigma_2 * HBAR_E / 2) * DBL3((dm_dt_2 ^ grad_m2.x) * m2, (dm_dt_2 ^ grad_m2.y) * m2, 0.0) * u;

				//final interpolated a value
				a += (1.5 * a1 - 0.5 * a2);
			}
		}
	}
	else {

		//non-magnetic mesh

		//1. ISHE contribution
		if (IsNZ(paMesh->iSHA.get0()) && stsolve != STSOLVE_NONE) {

			double iSHA = paMesh->iSHA;
			double SHA = paMesh->SHA;
			double De = paMesh->De;
			paMesh->update_parameters_atposition(relpos_m1, paMesh->iSHA, iSHA, paMesh->SHA, SHA, paMesh->De, De);

			int idx_S1 = paMesh->V.position_to_cellidx(relpos_m1);
			int idx_S2 = paMesh->V.position_to_cellidx(relpos_m1 + shift);

			//need to find value at boundary so use interpolation

			//value a1
			double sigma_1 = paMesh->elC.weighted_average(relpos_m1, stencil);
			//do not read off the E field directly as it's only calculated after the spin solver (charge part) has relaxed
			DBL3 E1 = -1.0 * paMesh->V.grad_diri_nneu(idx_S1, (iSHA * De / (MUB_E * sigma_1)) * paMesh->S.curl_neu(idx_S1));
			double a1 = (iSHA * De / MUB_E) * paMesh->S.curl_nneu(idx_S1, epsilon3(E1) * (SHA * sigma_1 * MUB_E / De)) * u;

			//value a2
			double sigma_2 = paMesh->elC.weighted_average(relpos_m1 + shift, stencil);
			DBL3 E2 = -1.0 * paMesh->V.grad_diri_nneu(idx_S2, (iSHA * De / (MUB_E * sigma_2)) * paMesh->S.curl_neu(idx_S2));
			double a2 = (iSHA * De / MUB_E) * paMesh->S.curl_nneu(idx_S2, epsilon3(E2) * (SHA * sigma_2 * MUB_E / De)) * u;

			//final interpolated a value
			a += (1.5 * a1 - 0.5 * a2);
		}
	}
	*/
	return a;
}

//flux = a + b V' at the interface, b = -sigma, a = betaD * (De*e/muB) * (grad S)m + (SHA*De*e/muB) * curl S + charge pumping + topological Hall effect
//Note, the topological Hall effect term includes E, thus gradients in V, but we can include these in the a term for 2 reasons:
//1. these CMBND functions especially with the topological Hall effect enabled is used for interfaces along z direction normally, and here Ez is zero 
//(for such interfaces both charge pumping and topological Hall effect have zero contribution to z direction charge current)
//2. even if the interface is along x or y we can still use the previously calculated E field, and the solution will converge to the same value (but might take more iterations).
//VERIFIED - CORRECT
double Atom_Transport::afunc_st_V_pri(int cell1_idx, int cell2_idx, DBL3 shift) const
{
	double a = 0.0;

	//TO DO

	/*
	DBL3 u = shift.normalized() * -1;

	bool cppgmr_enabled = IsNZ(paMesh->betaD.get0());
	bool cpump_enabled = IsNZ(paMesh->cpump_eff.get0()) && IsZ(shift.z);
	bool the_enabled = IsNZ(paMesh->the_eff.get0()) && IsZ(shift.z);

	if (stsolve == STSOLVE_FERROMAGNETIC && (cppgmr_enabled || cpump_enabled || the_enabled)) {

		//magnetic mesh

		double Ms = paMesh->Ms;
		paMesh->update_parameters_ecoarse(cell1_idx, paMesh->Ms, Ms);

		int idx_M1 = paMesh->M.position_to_cellidx(paMesh->V.cellidx_to_position(cell1_idx));
		int idx_M2 = paMesh->M.position_to_cellidx(paMesh->V.cellidx_to_position(cell2_idx));

		DBL3 m1 = paMesh->M[idx_M1] / Ms;
		DBL3 m2 = paMesh->M[idx_M2] / Ms;

		//1. CPP-GMR contribution
		if (cppgmr_enabled) {

			double betaD = paMesh->betaD;
			double De = paMesh->De;
			paMesh->update_parameters_ecoarse(cell1_idx, paMesh->betaD, betaD, paMesh->De, De);

			//value a1
			DBL33 grad_S1 = paMesh->S.grad_neu(cell1_idx);

			double a1 = ((grad_S1 * m1) * betaD * De / MUB_E) * u;

			//value a2
			DBL33 grad_S2 = paMesh->S.grad_neu(cell2_idx);

			double a2 = ((grad_S2 * m2) * betaD * De / MUB_E) * u;

			//final interpolated a value
			a += (1.5 * a1 - 0.5 * a2);
		}

		//2. Charge pumping
		//3. Topological Hall effect
		if (cpump_enabled || the_enabled) {

			double P = paMesh->P;
			double n_density = paMesh->n_density;
			paMesh->update_parameters_ecoarse(cell1_idx, paMesh->P, P, paMesh->n_density, n_density);

			DBL33 grad_m1 = paMesh->M.grad_neu(idx_M1) / Ms;
			DBL33 grad_m2 = paMesh->M.grad_neu(idx_M2) / Ms;

			//do not read off the E field directly as it's only calculated after the spin solver (charge part) has relaxed
			DBL3 E1 = -1.0 * paMesh->V.grad_diri(cell1_idx);
			DBL3 E2 = -1.0 * paMesh->V.grad_diri(cell2_idx);

			double sigma_1 = paMesh->elC[cell1_idx];
			double sigma_2 = paMesh->elC[cell2_idx];

			//topological Hall effect contribution
			if (the_enabled) {

				//value a1
				double Bz_the_1 = (grad_m1.x ^ grad_m1.y) * m1;
				double a1 = paMesh->the_eff.get0() * (-P * sigma_1 * sigma_1 * HBAR_E / (ECHARGE * n_density)) * DBL3(E1.y * Bz_the_1, -E1.x * Bz_the_1, 0.0) * u;

				//value a2
				double Bz_the_2 = (grad_m2.x ^ grad_m2.y) * m2;
				double a2 = paMesh->the_eff.get0() * (-P * sigma_2 * sigma_2 * HBAR_E / (ECHARGE * n_density)) * DBL3(E1.y * Bz_the_1, -E1.x * Bz_the_1, 0.0) * u;

				//final interpolated a value
				a += (1.5 * a1 - 0.5 * a2);
			}

			//charge pumping contribution
			if (cpump_enabled) {

				//value a1
				DBL3 dm_dt_1 = dM_dt[idx_M1] / Ms;
				double a1 = paMesh->cpump_eff.get0() * (P * sigma_1 * HBAR_E / 2) * DBL3((dm_dt_1 ^ grad_m1.x) * m1, (dm_dt_1 ^ grad_m1.y) * m1, 0.0) * u;

				//value a2
				DBL3 dm_dt_2 = dM_dt[idx_M2] / Ms;
				double a2 = paMesh->cpump_eff.get0() * (P * sigma_2 * HBAR_E / 2) * DBL3((dm_dt_2 ^ grad_m2.x) * m2, (dm_dt_2 ^ grad_m2.y) * m2, 0.0) * u;

				//final interpolated a value
				a += (1.5 * a1 - 0.5 * a2);
			}
		}
	}
	else {

		//non-magnetic mesh

		//1. ISHE contribution
		if (IsNZ(paMesh->iSHA.get0()) && stsolve != STSOLVE_NONE) {

			double iSHA = paMesh->iSHA;
			double SHA = paMesh->SHA;
			double De = paMesh->De;
			paMesh->update_parameters_atposition(cell1_idx, paMesh->iSHA, iSHA, paMesh->SHA, SHA, paMesh->De, De);

			//need to find value at boundary so use interpolation

			//value a1

			//do not read off the E field directly as it's only calculated after the spin solver (charge part) has relaxed
			DBL3 E1 = -1.0 * paMesh->V.grad_diri_nneu(cell1_idx, (iSHA * De / (MUB_E * paMesh->elC[cell1_idx])) * paMesh->S.curl_neu(cell1_idx));
			double a1 = (iSHA * De / MUB_E) * paMesh->S.curl_nneu(cell1_idx, epsilon3(E1) * (SHA * paMesh->elC[cell1_idx] * MUB_E / De)) * u;

			//value a2
			DBL3 E2 = -1.0 * paMesh->V.grad_diri_nneu(cell2_idx, (iSHA * De / (MUB_E * paMesh->elC[cell2_idx])) * paMesh->S.curl_neu(cell2_idx));
			double a2 = (iSHA * De / MUB_E) * paMesh->S.curl_nneu(cell2_idx, epsilon3(E2) * (SHA * paMesh->elC[cell2_idx] * MUB_E / De)) * u;

			//final interpolated a value
			a += (1.5 * a1 - 0.5 * a2);
		}
	}
	*/
	return a;
}

//VERIFIED - CORRECT
double Atom_Transport::bfunc_st_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	return -(1.5 * paMesh->elC.weighted_average(relpos_m1, stencil) - 0.5 * paMesh->elC.weighted_average(relpos_m1 + shift, stencil));
}

//VERIFIED - CORRECT
double Atom_Transport::bfunc_st_V_pri(int cell1_idx, int cell2_idx) const
{
	return -(1.5 * paMesh->elC[cell1_idx] - 0.5 * paMesh->elC[cell2_idx]);
}

//second order differential of V along the shift axis
//this is simply Evaluate_SpinSolver_delsqV_RHS from which we subtract second order differentials orthogonal to the shift axis
//VERIFIED - CORRECT
double Atom_Transport::diff2_st_V_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const
{
	int cellm1_idx = paMesh->V.position_to_cellidx(relpos_m1);

	return Evaluate_SpinSolver_delsqV_RHS(cellm1_idx);
}

//second order differential of V along the shift axis
//this is simply Evaluate_SpinSolver_delsqV_RHS from which we subtract second order differentials orthogonal to the shift axis
//VERIFIED - CORRECT
double Atom_Transport::diff2_st_V_pri(int cell1_idx, DBL3 shift) const
{
	return Evaluate_SpinSolver_delsqV_RHS(cell1_idx);
}

//CMBND for S
//flux = a + b S' at the interface, b = -De, a = -(muB*P*sigma/e) * E ox m + (SHA*sigma*muB/e)*epsE + charge pumping + topological Hall effect
//VERIFIED - CORRECT
DBL3 Atom_Transport::afunc_st_S_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	//need to find value at boundary so use interpolation
	//TO DO
	
	/*
	//unit vector perpendicular to interface (pointing from secondary to primary mesh)
	DBL3 u = shift.normalized() * -1;

	//values on secondary side
	if (stsolve == STSOLVE_FERROMAGNETIC) {

		double Ms = paMesh->Ms;
		double De = paMesh->De;
		double P = paMesh->P;
		paMesh->update_parameters_atposition(relpos_m1, paMesh->Ms, Ms, paMesh->De, De, paMesh->P, P);

		//a1 value
		DBL3 m1 = paMesh->M.weighted_average(relpos_m1, stencil) / Ms;
		DBL3 E1 = paMesh->E.weighted_average(relpos_m1, stencil);
		double sigma_1 = paMesh->elC.weighted_average(relpos_m1, stencil);

		int idx_S1 = paMesh->S.position_to_cellidx(relpos_m1);
		DBL33 grad_S1 = paMesh->S.grad_neu(idx_S1);

		DBL3 a1 = -MUB_E * ((E1 | m1) | u) * (P * sigma_1);

		//a2 value
		DBL3 m2 = paMesh->M.weighted_average(relpos_m1 + shift, stencil) / Ms;
		DBL3 E2 = paMesh->E.weighted_average(relpos_m1 + shift, stencil);
		double sigma_2 = paMesh->elC.weighted_average(relpos_m1 + shift, stencil);

		int idx_S2 = paMesh->S.position_to_cellidx(relpos_m1 + shift);
		DBL33 grad_S2 = paMesh->S.grad_neu(idx_S2);

		DBL3 a2 = -MUB_E * ((E2 | m2) | u) * (P * sigma_2);

		bool cpump_enabled = IsNZ(paMesh->cpump_eff.get0());
		bool the_enabled = IsNZ(paMesh->the_eff.get0());

		if (IsZ(shift.z) && (cpump_enabled || the_enabled)) {

			int idx_M1 = paMesh->M.position_to_cellidx(relpos_m1);
			int idx_M2 = paMesh->M.position_to_cellidx(relpos_m1 + shift);

			DBL33 grad_m1 = paMesh->M.grad_neu(idx_M1) / Ms;
			DBL33 grad_m2 = paMesh->M.grad_neu(idx_M2) / Ms;

			//topological Hall effect contribution
			if (the_enabled) {

				double n_density = paMesh->n_density;
				paMesh->update_parameters_atposition(relpos_m1, paMesh->n_density, n_density);

				DBL3 B1 = (grad_m1.x ^ grad_m1.y);
				a1 += paMesh->the_eff.get0() * (HBAR_E * MUB_E * sigma_1 * sigma_1 / (ECHARGE * n_density)) * ((-E1.y * B1 * u.x) + (E1.x * B1 * u.y));

				DBL3 B2 = (grad_m2.x ^ grad_m2.y);
				a2 += paMesh->the_eff.get0() * (HBAR_E * MUB_E * sigma_2 * sigma_2 / (ECHARGE * n_density)) * ((-E2.y * B2 * u.x) + (E2.x * B2 * u.y));
			}

			//charge pumping contribution
			if (cpump_enabled) {

				//value a1
				DBL3 dm_dt_1 = dM_dt.weighted_average(relpos_m1, stencil) / Ms;
				a1 += paMesh->cpump_eff.get0() * (HBAR_E * MUB_E * sigma_1 / 2) * (((dm_dt_1 ^ grad_m1.x) * u.x) + ((dm_dt_1 ^ grad_m1.y) * u.y));

				DBL3 dm_dt_2 = dM_dt.weighted_average(relpos_m1 + shift, stencil) / Ms;
				a2 += paMesh->cpump_eff.get0() * (HBAR_E * MUB_E * sigma_2 / 2) * (((dm_dt_2 ^ grad_m2.x) * u.x) + ((dm_dt_2 ^ grad_m2.y) * u.y));
			}
		}

		return (1.5 * a1 - 0.5 * a2);
	}
	*/
	return DBL3();
}

//flux = a + b S' at the interface, b = -De, a = -(muB*P*sigma/e) * E ox m + (SHA*sigma*muB/e)*epsE + charge pumping + topological Hall effect
//VERIFIED - CORRECT
DBL3 Atom_Transport::afunc_st_S_pri(int cell1_idx, int cell2_idx, DBL3 shift) const
{
	//need to find value at boundary so use interpolation

	//TO DO
	
	/*
	//unit vector perpendicular to interface (pointing from secondary to primary mesh)
	DBL3 u = shift.normalized() * -1;

	//values on secondary side
	if (stsolve == STSOLVE_FERROMAGNETIC) {

		double Ms = paMesh->Ms;
		double De = paMesh->De;
		double P = paMesh->P;
		paMesh->update_parameters_ecoarse(cell1_idx, paMesh->Ms, Ms, paMesh->De, De, paMesh->P, P);

		int idx_M1 = paMesh->M.position_to_cellidx(paMesh->S.cellidx_to_position(cell1_idx));
		int idx_M2 = paMesh->M.position_to_cellidx(paMesh->S.cellidx_to_position(cell2_idx));

		//a1 value
		DBL3 m1 = paMesh->M[idx_M1] / Ms;
		DBL3 E1 = paMesh->E[cell1_idx];
		double sigma_1 = paMesh->elC[cell1_idx];

		DBL33 grad_S1 = paMesh->S.grad_neu(cell1_idx);

		DBL3 a1 = -MUB_E * ((E1 | m1) | u) * (P * sigma_1);

		//a2 value
		DBL3 m2 = paMesh->M[idx_M2] / Ms;
		DBL3 E2 = paMesh->E[cell2_idx];
		double sigma_2 = paMesh->elC[cell2_idx];

		DBL33 grad_S2 = paMesh->S.grad_neu(cell2_idx);

		DBL3 a2 = -MUB_E * ((E2 | m2) | u) * (P * sigma_2);

		bool cpump_enabled = IsNZ(paMesh->cpump_eff.get0());
		bool the_enabled = IsNZ(paMesh->the_eff.get0());

		if (IsZ(shift.z) && (cpump_enabled || the_enabled)) {

			DBL33 grad_m1 = paMesh->M.grad_neu(idx_M1) / Ms;
			DBL33 grad_m2 = paMesh->M.grad_neu(idx_M2) / Ms;

			//topological Hall effect contribution
			if (the_enabled) {

				double n_density = paMesh->n_density;
				paMesh->update_parameters_ecoarse(cell1_idx, paMesh->n_density, n_density);

				DBL3 B1 = (grad_m1.x ^ grad_m1.y);
				a1 += paMesh->the_eff.get0() * (HBAR_E * MUB_E * sigma_1 * sigma_1 / (ECHARGE * n_density)) * ((-E1.y * B1 * u.x) + (E1.x * B1 * u.y));

				DBL3 B2 = (grad_m2.x ^ grad_m2.y);
				a2 += paMesh->the_eff.get0() * (HBAR_E * MUB_E * sigma_2 * sigma_2 / (ECHARGE * n_density)) * ((-E2.y * B2 * u.x) + (E2.x * B2 * u.y));
			}

			//charge pumping contribution
			if (cpump_enabled) {

				//value a1
				DBL3 dm_dt_1 = dM_dt[idx_M1] / Ms;
				a1 += paMesh->cpump_eff.get0() * (HBAR_E * MUB_E * sigma_1 / 2) * (((dm_dt_1 ^ grad_m1.x) * u.x) + ((dm_dt_1 ^ grad_m1.y) * u.y));

				DBL3 dm_dt_2 = dM_dt[idx_M2] / Ms;
				a2 += paMesh->cpump_eff.get0() * (HBAR_E * MUB_E * sigma_2 / 2) * (((dm_dt_2 ^ grad_m2.x) * u.x) + ((dm_dt_2 ^ grad_m2.y) * u.y));
			}
		}

		return (1.5 * a1 - 0.5 * a2);
	}
	*/
	return DBL3();
}

//flux = a + b S' at the interface, b = -De
//VERIFIED - CORRECT
double Atom_Transport::bfunc_st_S_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	//TO DO
	/*
	double De = paMesh->De;
	paMesh->update_parameters_atposition(relpos_m1, paMesh->De, De);

	return -1.0 * De;
	*/
	return 0.0;
}

//flux = a + b S' at the interface, b = -De
//VERIFIED - CORRECT
double Atom_Transport::bfunc_st_S_pri(int cell1_idx, int cell2_idx) const
{
	//TO DO
	/*
	double De = paMesh->De;
	paMesh->update_parameters_ecoarse(cell1_idx, paMesh->De, De);

	return -1.0 * De;
	*/
	return 0.0;
}

//second order differential of S along the shift axis
//this is simply Evaluate_SpinSolver_delsqS_RHS from which we subtract second order differentials orthogonal to the shift axis
//VERIFIED - CORRECT
DBL3 Atom_Transport::diff2_st_S_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const
{
	//TO DO
	/*
	int cellm1_idx = paMesh->S.position_to_cellidx(relpos_m1);

	return Evaluate_SpinSolver_delsqS_RHS(cellm1_idx);
	*/
	return 0.0;
}

//second order differential of S along the shift axis
//this is simply Evaluate_SpinSolver_delsqS_RHS from which we subtract second order differentials orthogonal to the shift axis
//VERIFIED - CORRECT
DBL3 Atom_Transport::diff2_st_S_pri(int cell1_idx, DBL3 shift) const
{
	return Evaluate_SpinSolver_delsqS_RHS(cell1_idx);
}

//-------------------

//VERIFIED - CORRECT
double Atom_Transport::cfunc_sec(DBL3 relpos, DBL3 stencil) const
{
	//TO DO
	/*
	double De = paMesh->De;
	paMesh->update_parameters_atposition(relpos, paMesh->De, De);

	return De / (paMesh->elC.weighted_average(relpos, stencil) * MUB_E);
	*/
	return 0.0;
}

//VERIFIED - CORRECT
double Atom_Transport::cfunc_pri(int cell_idx) const
{
	//TO DO
	/*
	double De = paMesh->De;
	paMesh->update_parameters_ecoarse(cell_idx, paMesh->De, De);

	return De / (paMesh->elC[cell_idx] * MUB_E);
	*/
	return 0.0;
}

//-------------------

//calculate the field resulting from a spin accumulation (spin current solver enabled) so a spin accumulation torque results when using the LLG or LLB equations
//VERIFIED - CORRECT
void Atom_Transport::CalculateSAField(void)
{
	//TO DO
	/*
	if (stsolve != STSOLVE_FERROMAGNETIC) return;

#pragma omp parallel for
	for (int idx = 0; idx < paMesh->M.linear_size(); idx++) {

		if (paMesh->M.is_not_empty(idx)) {

			double Ms = paMesh->Ms;
			double grel = paMesh->grel;
			double De = paMesh->De;
			double ts_eff = paMesh->ts_eff;
			double l_ex = paMesh->l_ex;
			double l_ph = paMesh->l_ph;
			paMesh->update_parameters_mcoarse(idx, paMesh->Ms, Ms, paMesh->grel, grel, paMesh->De, De, paMesh->ts_eff, ts_eff, paMesh->l_ex, l_ex, paMesh->l_ph, l_ph);

			if (IsNZ((double)grel)) {

				DBL3 Sa = paMesh->S.weighted_average(paMesh->M.cellidx_to_position(idx), paMesh->h);

				paMesh->Heff[idx] += (De * ts_eff / (GAMMA * grel * Ms)) * (Sa / (l_ex * l_ex) + (paMesh->M[idx] ^ Sa) / (l_ph * l_ph * Ms));
			}
		}
	}
	*/
}

//Calculate interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set)
//VERIFIED - CORRECT
void Atom_Transport::CalculateSAInterfaceField(Atom_Transport* ptrans_sec, CMBNDInfo& contact)
{
	//TO DO
	/*
	//the top contacting mesh sets G values
	bool isGInterface_Enabled = ((contact.IsPrimaryTop() && paMesh->GInterface_Enabled()) || (!contact.IsPrimaryTop() && ptrans_sec->GInterface_Enabled()));

	if (isGInterface_Enabled && stsolve == STSOLVE_FERROMAGNETIC && ptrans_sec->Get_STSolveType() == STSOLVE_NORMALMETAL) {

		//interface conductance method with F being the primary mesh (N-F) contact : calculate and set spin torque

		//convert the cells box from S mesh to M mesh
		INT3 mbox_start = paMesh->M.cellidx_from_position(paMesh->S.cellidx_to_position(contact.cells_box.s) + paMesh->meshRect.s);
		INT3 mbox_end = paMesh->M.cellidx_from_position(paMesh->S.cellidx_to_position(contact.cells_box.e - INT3(1)) + paMesh->meshRect.s);

		if ((mbox_end.i - mbox_start.i) == 0) mbox_end.i = mbox_start.i + 1;
		if ((mbox_end.j - mbox_start.j) == 0) mbox_end.j = mbox_start.j + 1;
		if ((mbox_end.k - mbox_start.k) == 0) mbox_end.k = mbox_start.k + 1;

		INT3 box_sizes = mbox_end - mbox_start;

		//the cellsize perpendicular to the contact (in the M mesh)
		double dh = (DBL3(contact.cell_shift) & paMesh->h).norm();

		Mesh* pMesh_sec = ptrans_sec->pMesh;

		//primary cells in this contact
#pragma omp parallel for
		for (int box_idx = 0; box_idx < box_sizes.dim(); box_idx++) {

			int i = (box_idx % box_sizes.x) + mbox_start.i;
			int j = ((box_idx / box_sizes.x) % box_sizes.y) + mbox_start.j;
			int k = (box_idx / (box_sizes.x * box_sizes.y)) + mbox_start.k;

			//index of magnetic cell 1
			int mcell1_idx = i + j * paMesh->n.x + k * paMesh->n.x * paMesh->n.y;

			if (paMesh->M.is_empty(mcell1_idx)) continue;

			double Ms = paMesh->Ms;
			double grel = paMesh->grel;
			double tsi_eff = paMesh->tsi_eff;
			paMesh->update_parameters_mcoarse(mcell1_idx, paMesh->grel, grel, paMesh->Ms, Ms, paMesh->tsi_eff, tsi_eff);

			if (IsNZ((double)grel)) {

				//position at interface relative to primary mesh
				DBL3 mhshift_primary = contact.hshift_primary.normalized() & paMesh->h;
				DBL3 relpos_interf = ((DBL3(i, j, k) + DBL3(0.5)) & paMesh->h) + mhshift_primary / 2;

				DBL3 relpos_1 = relpos_interf - contact.hshift_primary / 2;

				DBL3 relpos_m1 = paMesh->meshRect.s - pMesh_sec->meshRect.s + relpos_interf + contact.hshift_secondary / 2;

				DBL3 stencil = paMesh->h - mod(mhshift_primary) + mod(contact.hshift_secondary);

				//S values
				DBL3 S_1 = paMesh->S.weighted_average(relpos_1, stencil);
				DBL3 S_2 = paMesh->S.weighted_average(relpos_1 - contact.hshift_primary, stencil);
				DBL3 S_m1 = pMesh_sec->S.weighted_average(relpos_m1, stencil);
				DBL3 S_m2 = pMesh_sec->S.weighted_average(relpos_m1 + contact.hshift_secondary, stencil);

				//c values
				double c_m1 = ptrans_sec->cfunc_sec(relpos_m1, stencil);
				double c_m2 = ptrans_sec->cfunc_sec(relpos_m1 + contact.hshift_secondary, stencil);
				double c_1 = cfunc_sec(relpos_1, stencil);
				double c_2 = cfunc_sec(relpos_1 - contact.hshift_primary, stencil);

				//Calculate S drop at the interface
				DBL3 Vs_F = 1.5 * c_1 * S_1 - 0.5 * c_2 * S_2;
				DBL3 Vs_N = 1.5 * c_m1 * S_m1 - 0.5 * c_m2 * S_m2;
				DBL3 dVs = Vs_F - Vs_N;

				//Get G values from top contacting mesh
				DBL2 Gmix;
				if (contact.IsPrimaryTop()) {

					Gmix = paMesh->Gmix;
					paMesh->update_parameters_mcoarse(mcell1_idx, paMesh->Gmix, Gmix);
				}
				else {

					Gmix = pMesh_sec->Gmix;
					pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gmix, Gmix);
				}

				double gI = (2.0 * GMUB_2E / dh) * Gmix.j / (-GAMMA * grel * Ms);
				double gR = (2.0 * GMUB_2E / dh) * Gmix.i / (-GAMMA * grel * Ms);

				paMesh->Heff[mcell1_idx] += tsi_eff * (gI * dVs + gR * (paMesh->M[mcell1_idx] ^ dVs) / Ms);
			}
		}
	}
	*/
}

#endif