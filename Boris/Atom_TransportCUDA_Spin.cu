#include "Atom_TransportCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "Atom_MeshCUDA.h"
#include "SuperMeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

//-------------------Calculation Methods : Iterate Spin-Charge Solver

void Atom_TransportCUDA::IterateSpinSolver_Charge_SOR(cu_obj<cuBReal>& damping, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value, bool use_NNeu)
{
	//use_NNeu not needed here, but implementing TransportBaseCUDA interface
	paMeshCUDA->V()->IteratePoisson_SOR(paMeshCUDA->n_e.dim(), (Atom_TransportCUDA_Spin_V_Funcs&)poisson_Spin_V, damping, max_error, max_value);
}

//------------------- PRIME SPIN-CHARGE SOLVER

__global__ void Atom_Get_dM_dt_Kernel(Atom_TransportCUDA_Spin_S_Funcs& poisson_Spin_S)
{
	cuVEC_VC<cuReal3>& dM_dt = *poisson_Spin_S.pdM_dt;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < dM_dt.linear_size()) {

		if (dM_dt.is_not_empty(idx)) {

			//TO DO
			//dM_dt[idx] = poisson_Spin_S.pcuDiffEq->dMdt(idx);
		}
	}
}

__global__ void Atom_PrimeSpinSolver_Charge_Kernel(ManagedAtom_MeshCUDA& cuaMesh, Atom_TransportCUDA_Spin_V_Funcs& poisson_Spin_V)
{
	cuVEC_VC<cuBReal>& V = *cuaMesh.pV;
	cuVEC_VC<cuBReal>& elC = *cuaMesh.pelC;
	cuVEC_VC<cuReal3>& S = *cuaMesh.pS;
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;

	cuVEC<cuBReal>& delsq_V_fixed = *poisson_Spin_V.pdelsq_V_fixed;
	cuVEC_VC<cuReal3>& dM_dt = *poisson_Spin_V.pdM_dt;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < delsq_V_fixed.linear_size()) {

		delsq_V_fixed[idx] = 0.0;

		//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
		if (V.is_not_empty(idx)) {

			bool cppgmr_enabled = cuIsNZ(cuaMesh.pbetaD->get0());
			bool cpump_enabled = cuIsNZ(cuaMesh.pcpump_eff->get0());

			//TO DO
			/*
			if (cppgmr_enabled || cpump_enabled) {

				cuBReal Ms = *cuaMesh.pMs;
				cuaMesh.update_parameters_ecoarse(idx, *cuaMesh.pMs, Ms);

				int idx_M = M.position_to_cellidx(V.cellidx_to_position(idx));
				cuReal33 grad_m = M.grad_neu(idx_M) / Ms;
				cuReal3 m = M[idx_M] / Ms;

				//CPP-GMR contribution
				if (cppgmr_enabled) {

					cuBReal De = *cuaMesh.pDe;
					cuBReal betaD = *cuaMesh.pbetaD;
					cuaMesh.update_parameters_ecoarse(idx, *cuaMesh.pDe, De, *cuaMesh.pbetaD, betaD);

					cuReal33 grad_S = S.grad_neu(idx);
					cuReal3 delsq_S = S.delsq_neu(idx);
					cuBReal div_grad_S_m = (grad_S.i * grad_m.i) + (grad_S.j * grad_m.j) + (grad_S.k * grad_m.k) + (m * delsq_S);

					delsq_V_fixed[idx] += div_grad_S_m * betaD * De / ((cuBReal)MUB_E * elC[idx]);
				}

				//Charge pumping pre-calculation
				if (cpump_enabled) {

					cuReal33 grad_dm_dt = dM_dt.grad_neu(idx_M) / Ms;
					cuReal3 dm_dt = dM_dt[idx_M] / Ms;

					cuBReal P = *cuaMesh.pP;
					cuaMesh.update_parameters_ecoarse(idx, *cuaMesh.pP, P);

					cuReal3 dx_m = grad_m.x;
					cuReal3 dy_m = grad_m.y;
					cuReal3 dxx_m = M.dxx_neu(idx_M) / Ms;
					cuReal3 dyy_m = M.dyy_neu(idx_M) / Ms;

					delsq_V_fixed[idx] += (cuaMesh.pcpump_eff->get0() * P * (cuBReal)HBAR_E / 2) * ((grad_dm_dt.x ^ dx_m) + (grad_dm_dt.y ^ dy_m) + (dm_dt ^ (dxx_m + dyy_m))) * m;
				}
			}
			*/
		}
	}
}

//before iterating the spin solver (charge part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
void Atom_TransportCUDA::PrimeSpinSolver_Charge(void)
{
	//Update dM_dt values if needed
	if (Need_dM_dt_Calculation()) {

		Atom_Get_dM_dt_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (poisson_Spin_S);
	}

	//the rest are terms to calculate in delsq_V_fixed
	if (Need_delsq_V_fixed_Precalculation()) {

		Atom_PrimeSpinSolver_Charge_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, poisson_Spin_V);
	}
}

//-------------------Calculation Methods : Iterate Spin-Spin Solver

//solve for spin accumulation using Poisson equation for delsq_S, solved using SOR algorithm
void Atom_TransportCUDA::IterateSpinSolver_Spin_SOR(cu_obj<cuBReal>& damping, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value, bool use_NNeu)
{
	//use_NNeu not needed here, but implementing TransportBaseCUDA interface
	paMeshCUDA->S()->IteratePoisson_SOR(paMeshCUDA->n_e.dim(), (Atom_TransportCUDA_Spin_S_Funcs&)poisson_Spin_S, damping, max_error, max_value);
}

//------------------- PRIME SPIN-SPIN SOLVER

__global__ void Atom_PrimeSpinSolver_Spin_Kernel(ManagedAtom_MeshCUDA& cuaMesh, Atom_TransportCUDA_Spin_S_Funcs& poisson_Spin_S)
{
	cuVEC_VC<cuBReal>& V = *cuaMesh.pV;
	cuVEC_VC<cuBReal>& elC = *cuaMesh.pelC;
	cuVEC_VC<cuReal3>& E = *cuaMesh.pE;
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;

	cuVEC<cuReal3>& delsq_S_fixed = *poisson_Spin_S.pdelsq_S_fixed;
	cuVEC_VC<cuReal3>& dM_dt = *poisson_Spin_S.pdM_dt;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < delsq_S_fixed.linear_size()) {

		delsq_S_fixed[idx] = 0.0;

		//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
		if (V.is_not_empty(idx)) {

			bool cpump_enabled = cuIsNZ(cuaMesh.pcpump_eff->get0());
			bool the_enabled = cuIsNZ(cuaMesh.pthe_eff->get0());
			bool she_enabled = cuIsNZ(cuaMesh.pSHA->get0());

			//TO DO
			/*
			if (poisson_Spin_S.stsolve == STSOLVE_FERROMAGNETIC) {

				//magnetic mesh

				cuBReal Ms = *cuaMesh.pMs;
				cuBReal P = *cuaMesh.pP;
				cuBReal De = *cuaMesh.pDe;
				cuaMesh.update_parameters_ecoarse(idx, *cuaMesh.pMs, Ms, *cuaMesh.pP, P, *cuaMesh.pDe, De);

				//term due to drift (non-uniformity of M term, and delsq V contribution - non-uniformity of E term)

				//find grad M and M at the M cell in which the current S cell center is
				int idx_M = M.position_to_cellidx(V.cellidx_to_position(idx));

				cuReal3 m = M[idx_M] / Ms;
				cuReal33 grad_m = M.grad_neu(idx_M) / Ms;
				cuReal3 E_dot_del_m = grad_m | E[idx];

				//E_dot_del_m term is very important, but Evaluate_SpinSolver_delsqV_RHS term could be neglected in most cases especially if E is uniform.
				delsq_S_fixed[idx] += (P * (cuBReal)MUB_E * elC[idx] / De) * (poisson_Spin_S.pPoisson_Spin_V->Poisson_RHS(idx) * m - E_dot_del_m);

				//charge pumping and topological Hall effect
				if (cpump_enabled || the_enabled) {

					cuReal3 dx_m = grad_m.x;
					cuReal3 dy_m = grad_m.y;
					cuReal3 dxy_m = M.dxy_neu(idx_M) / Ms;
					cuReal3 dxx_m = M.dxx_neu(idx_M) / Ms;
					cuReal3 dyy_m = M.dyy_neu(idx_M) / Ms;

					if (cpump_enabled) {

						cuReal3 dmdt = dM_dt[idx_M] / Ms;
						cuReal33 grad_dm_dt = dM_dt.grad_neu(idx_M) / Ms;

						delsq_S_fixed[idx] += cuaMesh.pcpump_eff->get0() * (elC[idx] * (cuBReal)HBAR_E * (cuBReal)MUB_E / (2 * De)) * ((grad_dm_dt.x ^ dx_m) + (grad_dm_dt.y ^ dy_m) + (dmdt ^ (dxx_m + dyy_m)));
					}

					if (the_enabled) {

						cuBReal n_density = *cuaMesh.pn_density;
						cuaMesh.update_parameters_ecoarse(idx, *cuaMesh.pn_density, n_density);

						delsq_S_fixed[idx] += cuaMesh.pthe_eff->get0() * ((cuBReal)HBAR_E * (cuBReal)MUB_E * elC[idx] * elC[idx] / ((cuBReal)ECHARGE * n_density * De)) * (E[idx].x * ((dxy_m ^ dy_m) + (dx_m ^ dyy_m)) - E[idx].y * ((dxx_m ^ dy_m) + (dx_m ^ dxy_m)));
					}
				}
			}
			*/
		}
	}
}

//before iterating the spin solver (spin part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
void Atom_TransportCUDA::PrimeSpinSolver_Spin(void)
{
	if (Need_delsq_S_fixed_Precalculation()) {

		Atom_PrimeSpinSolver_Spin_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, poisson_Spin_S);
	}
}

//--------------------------------------------------------------- Effective field from spin accumulation

__global__ void Atom_CalculateSAField_Kernel(ManagedAtom_MeshCUDA& cuaMesh)
{
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC_VC<cuReal3>& S = *cuaMesh.pS;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx)) {

			//TO DO
			/*
			cuBReal De = *cuaMesh.pDe;
			cuBReal ts_eff = *cuaMesh.pts_eff;
			cuBReal grel = *cuaMesh.pgrel;
			cuBReal Ms = *cuaMesh.pMs;
			cuBReal l_ex = *cuaMesh.pl_ex;
			cuBReal l_ph = *cuaMesh.pl_ph;
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pgrel, grel, *cuaMesh.pMs, Ms, *cuaMesh.pDe, De, *cuaMesh.pts_eff, ts_eff, *cuaMesh.pl_ex, l_ex, *cuaMesh.pl_ph, l_ph);

			if (cuIsNZ((cuBReal)grel)) {

				cuReal3 Sa = S.weighted_average(M.cellidx_to_position(idx), M.h);

				Heff[idx] += (De * ts_eff / ((cuBReal)GAMMA * grel * Ms)) * (Sa / (l_ex * l_ex) + (M[idx] ^ Sa) / (l_ph * l_ph * Ms));
			}
			*/
		}
	}
}

//Spin accumulation field
void Atom_TransportCUDA::CalculateSAField(void)
{
	if (stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

		Atom_CalculateSAField_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh);
	}
}

//--------------------------------------------------------------- Effective field from interface spin accumulation drop

__global__ void Atom_CalculateSAInterfaceField_Kernel(CMBNDInfoCUDA& contact, Atom_TransportCUDA_Spin_S_Funcs& cmbndFuncs_sec, Atom_TransportCUDA_Spin_S_Funcs& cmbndFuncs_pri)
{
	cuVEC<cuReal3>& Heff1 = *cmbndFuncs_pri.pcuaMesh->pHeff1;
	cuVEC_VC<cuReal3>& M1 = *cmbndFuncs_pri.pcuaMesh->pM1;
	cuVEC_VC<cuReal3>& S_pri = *cmbndFuncs_pri.pcuaMesh->pS;
	cuVEC_VC<cuReal3>& S_sec = *cmbndFuncs_sec.pcuaMesh->pS;

	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	//interface conductance method with F being the primary mesh : calculate and set spin torque

	//convert the cells box from S mesh to M mesh
	cuINT3 mbox_start = M1.cellidx_from_position(S_pri.cellidx_to_position(contact.cells_box.s) + M1.rect.s);
	cuINT3 mbox_end = M1.cellidx_from_position(S_pri.cellidx_to_position(contact.cells_box.e - cuINT3(1)) + M1.rect.s);

	if ((mbox_end.i - mbox_start.i) == 0) mbox_end.i = mbox_start.i + 1;
	if ((mbox_end.j - mbox_start.j) == 0) mbox_end.j = mbox_start.j + 1;
	if ((mbox_end.k - mbox_start.k) == 0) mbox_end.k = mbox_start.k + 1;

	cuINT3 box_sizes = mbox_end - mbox_start;
	
	//TO DO
	/*
	if (box_idx < box_sizes.dim()) {

		//the cellsize perpendicular to the contact (in the M mesh)
		cuBReal dh = (cuReal3(contact.cell_shift) & M.h).norm();

		int i = (box_idx % box_sizes.x) + mbox_start.i;
		int j = ((box_idx / box_sizes.x) % box_sizes.y) + mbox_start.j;
		int k = (box_idx / (box_sizes.x * box_sizes.y)) + mbox_start.k;

		//index of magnetic cell 1
		int mcell1_idx = i + j * M.n.x + k * M.n.x*M.n.y;

		if (M.is_empty(mcell1_idx)) return;

		cuBReal grel = *cmbndFuncs_pri.pcuaMesh->pgrel;
		cuBReal Ms = *cmbndFuncs_pri.pcuaMesh->pMs;
		cuBReal tsi_eff = *cmbndFuncs_pri.pcuaMesh->ptsi_eff;
		cmbndFuncs_pri.pcuaMesh->update_parameters_mcoarse(mcell1_idx, *cmbndFuncs_pri.pcuaMesh->pgrel, grel, *cmbndFuncs_pri.pcuaMesh->pMs, Ms, *cmbndFuncs_pri.pcuaMesh->ptsi_eff, tsi_eff);

		if (cuIsNZ((cuBReal)grel)) {

			//position at interface relative to primary mesh
			cuReal3 mhshift_primary = contact.hshift_primary.normalized() & M.h;
			cuReal3 relpos_interf = ((cuReal3(i, j, k) + cuReal3(0.5)) & M.h) + mhshift_primary / 2;

			cuReal3 relpos_1 = relpos_interf - contact.hshift_primary / 2;

			cuReal3 relpos_m1 = S_pri.rect.s - S_sec.rect.s + relpos_interf + contact.hshift_secondary / 2;

			cuReal3 stencil = M.h - cu_mod(mhshift_primary) + cu_mod(contact.hshift_secondary);

			//S values
			cuReal3 S_1 = S_pri.weighted_average(relpos_1, stencil);
			cuReal3 S_2 = S_pri.weighted_average(relpos_1 - contact.hshift_primary, stencil);
			cuReal3 S_m1 = S_sec.weighted_average(relpos_m1, stencil);
			cuReal3 S_m2 = S_sec.weighted_average(relpos_m1 + contact.hshift_secondary, stencil);

			//c values
			cuBReal c_m1 = cmbndFuncs_sec.c_func_sec(relpos_m1, stencil);
			cuBReal c_m2 = cmbndFuncs_sec.c_func_sec(relpos_m1 + contact.hshift_secondary, stencil);
			cuBReal c_1 = cmbndFuncs_pri.c_func_sec(relpos_1, stencil);
			cuBReal c_2 = cmbndFuncs_pri.c_func_sec(relpos_1 - contact.hshift_primary, stencil);

			//Calculate S drop at the interface
			cuReal3 Vs_F = 1.5 * c_1 * S_1 - 0.5 * c_2 * S_2;
			cuReal3 Vs_N = 1.5 * c_m1 * S_m1 - 0.5 * c_m2 * S_m2;
			cuReal3 dVs = Vs_F - Vs_N;

			//Get G values from top contacting mesh
			cuReal2 Gmix;
			if (contact.IsPrimaryTop()) {

				Gmix = *cmbndFuncs_pri.pcuaMesh->pGmix;
				cmbndFuncs_pri.pcuaMesh->update_parameters_mcoarse(mcell1_idx, *cmbndFuncs_pri.pcuaMesh->pGmix, Gmix);
			}
			else {

				Gmix = *cmbndFuncs_sec.pcuaMesh->pGmix;
				cmbndFuncs_sec.pcuaMesh->update_parameters_atposition(relpos_m1, *cmbndFuncs_sec.pcuaMesh->pGmix, Gmix);
			}

			cuBReal gI = (2.0 * (cuBReal)GMUB_2E / dh) * cuReal2(Gmix).j / (-(cuBReal)GAMMA * grel * Ms);
			cuBReal gR = (2.0 * (cuBReal)GMUB_2E / dh) * cuReal2(Gmix).i / (-(cuBReal)GAMMA * grel * Ms);

			Heff[mcell1_idx] += tsi_eff * (gI * dVs + gR * (M[mcell1_idx] ^ dVs) / Ms);
		}
	}
	*/
}

//Calculate the field resulting from interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set)
void Atom_TransportCUDA::CalculateSAInterfaceField(Atom_TransportCUDA* ptrans_sec, CMBNDInfoCUDA& contactCUDA, bool primary_top)
{
	//the top contacting mesh sets G values
	bool isGInterface_Enabled = ((primary_top && paMeshCUDA->GInterface_Enabled()) || (!primary_top && ptrans_sec->paMeshCUDA->GInterface_Enabled()));
	//TO DO
	/*
	if (stsolve == STSOLVE_FERROMAGNETIC && ptrans_sec->Get_STSolveType() == STSOLVE_NORMALMETAL && isGInterface_Enabled) {

		Atom_CalculateSAInterfaceField_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (contactCUDA, ptrans_sec->poisson_Spin_S, poisson_Spin_S);
	}
	*/
}

#endif

#endif