#include "TransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_TRANSPORT

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"
#include "SuperMeshCUDA.h"
#include "MeshParamsControlCUDA.h"

//-------------------Calculation Methods

void TransportCUDA::IterateSpinSolver_Charge_aSOR(bool start_iters, cuBReal conv_error, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value, bool use_NNeu)
{
	if (use_NNeu) {

		//iSHE enabled, must use non-homogeneous Neumann boundary condition for grad V -> Note homogeneous Neumann boundary conditions apply when calculating S differentials here (due to Jc.n = 0 at boundaries)
		pMeshCUDA->V()->IteratePoisson_NNeu_aSOR(pMeshCUDA->n_e.dim(), (TransportCUDA_Spin_V_Funcs&)poisson_Spin_V, start_iters, conv_error, max_error, max_value);
	}
	else {

		//no iSHE contribution. Note, iSHE is not included in magnetic meshes.
		pMeshCUDA->V()->IteratePoisson_aSOR(pMeshCUDA->n_e.dim(), (TransportCUDA_Spin_V_Funcs&)poisson_Spin_V, start_iters, conv_error, max_error, max_value);
	}
}

void TransportCUDA::IterateSpinSolver_Charge_SOR(cu_obj<cuBReal>& damping, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value, bool use_NNeu)
{
	if (use_NNeu) {

		//iSHE enabled, must use non-homogeneous Neumann boundary condition for grad V -> Note homogeneous Neumann boundary conditions apply when calculating S differentials here (due to Jc.n = 0 at boundaries)
		pMeshCUDA->V()->IteratePoisson_NNeu_SOR(pMeshCUDA->n_e.dim(), (TransportCUDA_Spin_V_Funcs&)poisson_Spin_V, damping, max_error, max_value);
	}
	else {

		//no iSHE contribution. Note, iSHE is not included in magnetic meshes.
		pMeshCUDA->V()->IteratePoisson_SOR(pMeshCUDA->n_e.dim(), (TransportCUDA_Spin_V_Funcs&)poisson_Spin_V, damping, max_error, max_value);
	}
}

//solve for spin accumulation using Poisson equation for delsq_S, solved using SOR algorithm
void TransportCUDA::IterateSpinSolver_Spin_aSOR(bool start_iters, cuBReal conv_error, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value, bool use_NNeu)
{
	if (use_NNeu) {

		//SHE enabled, must use non-homogeneous Neumann boundary condition for grad S
		pMeshCUDA->S()->IteratePoisson_NNeu_aSOR(pMeshCUDA->n_e.dim(), (TransportCUDA_Spin_S_Funcs&)poisson_Spin_S, start_iters, conv_error, max_error, max_value);
	}
	else {

		//no SHE contribution. Note, SHE is not included in magnetic meshes.
		pMeshCUDA->S()->IteratePoisson_aSOR(pMeshCUDA->n_e.dim(), (TransportCUDA_Spin_S_Funcs&)poisson_Spin_S, start_iters, conv_error, max_error, max_value);
	}
}

//solve for spin accumulation using Poisson equation for delsq_S, solved using SOR algorithm
void TransportCUDA::IterateSpinSolver_Spin_SOR(cu_obj<cuBReal>& damping, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value, bool use_NNeu)
{
	if (use_NNeu) {

		//SHE enabled, must use non-homogeneous Neumann boundary condition for grad S
		pMeshCUDA->S()->IteratePoisson_NNeu_SOR(pMeshCUDA->n_e.dim(), (TransportCUDA_Spin_S_Funcs&)poisson_Spin_S, damping, max_error, max_value);
	}
	else {

		//no SHE contribution. Note, SHE is not included in magnetic meshes.
		pMeshCUDA->S()->IteratePoisson_SOR(pMeshCUDA->n_e.dim(), (TransportCUDA_Spin_S_Funcs&)poisson_Spin_S, damping, max_error, max_value);
	}
}

//--------------------------------------------------------------- Effective field from spin accumulation

__global__ void CalculateSAField_Kernel(ManagedMeshCUDA& cuMesh)
{
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& S = *cuMesh.pS;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			cuBReal De = *cuMesh.pDe;
			cuBReal ts_eff = *cuMesh.pts_eff;
			cuBReal grel = *cuMesh.pgrel;
			cuBReal Ms = *cuMesh.pMs;
			cuBReal l_ex = *cuMesh.pl_ex;
			cuBReal l_ph = *cuMesh.pl_ph;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pgrel, grel, *cuMesh.pMs, Ms, *cuMesh.pDe, De, *cuMesh.pts_eff, ts_eff, *cuMesh.pl_ex, l_ex, *cuMesh.pl_ph, l_ph);

			if (cuIsNZ((cuBReal)grel)) {

				cuReal3 Sa = S.weighted_average(M.cellidx_to_position(idx), M.h);

				Heff[idx] += (De * ts_eff / ((cuBReal)GAMMA * grel * Ms)) * (Sa / (l_ex * l_ex) + (M[idx] ^ Sa) / (l_ph * l_ph * Ms));
			}
		}
	}
}

//Spin accumulation field
void TransportCUDA::CalculateSAField(void)
{
	CalculateSAField_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh);
}

//--------------------------------------------------------------- Effective field from interface spin accumulation drop

__global__ void CalculateSAInterfaceField_Kernel(CMBNDInfoCUDA& contact, TransportCUDA_Spin_S_Funcs& cmbndFuncs_sec, TransportCUDA_Spin_S_Funcs& cmbndFuncs_pri)
{
	cuVEC<cuReal3>& Heff = *cmbndFuncs_pri.pcuMesh->pHeff;
	cuVEC_VC<cuReal3>& M = *cmbndFuncs_pri.pcuMesh->pM;
	cuVEC_VC<cuReal3>& S_pri = *cmbndFuncs_pri.pcuMesh->pS;
	cuVEC_VC<cuReal3>& S_sec = *cmbndFuncs_sec.pcuMesh->pS;

	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	//interface conductance method with F being the primary mesh : calculate and set spin torque

	//convert the cells box from S mesh to M mesh
	cuINT3 mbox_start = M.cellidx_from_position(S_pri.cellidx_to_position(contact.cells_box.s) + M.rect.s);
	cuINT3 mbox_end = M.cellidx_from_position(S_pri.cellidx_to_position(contact.cells_box.e - cuINT3(1)) + M.rect.s);

	if ((mbox_end.i - mbox_start.i) == 0) mbox_end.i = mbox_start.i + 1;
	if ((mbox_end.j - mbox_start.j) == 0) mbox_end.j = mbox_start.j + 1;
	if ((mbox_end.k - mbox_start.k) == 0) mbox_end.k = mbox_start.k + 1;

	cuINT3 box_sizes = mbox_end - mbox_start;

	if (box_idx < box_sizes.dim()) {

		//the cellsize perpendicular to the contact (in the M mesh)
		cuBReal dh = (cuReal3(contact.cell_shift) & M.h).norm();

		int i = (box_idx % box_sizes.x) + mbox_start.i;
		int j = ((box_idx / box_sizes.x) % box_sizes.y) + mbox_start.j;
		int k = (box_idx / (box_sizes.x * box_sizes.y)) + mbox_start.k;

		//index of magnetic cell 1
		int mcell1_idx = i + j * M.n.x + k * M.n.x*M.n.y;

		if (M.is_empty(mcell1_idx)) return;

		cuBReal grel = *cmbndFuncs_pri.pcuMesh->pgrel;
		cuBReal Ms = *cmbndFuncs_pri.pcuMesh->pMs;
		cuBReal tsi_eff = *cmbndFuncs_pri.pcuMesh->ptsi_eff;
		cmbndFuncs_pri.pcuMesh->update_parameters_mcoarse(mcell1_idx, *cmbndFuncs_pri.pcuMesh->pgrel, grel, *cmbndFuncs_pri.pcuMesh->pMs, Ms, *cmbndFuncs_pri.pcuMesh->ptsi_eff, tsi_eff);

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

				Gmix = *cmbndFuncs_pri.pcuMesh->pGmix;
				cmbndFuncs_pri.pcuMesh->update_parameters_mcoarse(mcell1_idx, *cmbndFuncs_pri.pcuMesh->pGmix, Gmix);
			}
			else {

				Gmix = *cmbndFuncs_sec.pcuMesh->pGmix;
				cmbndFuncs_sec.pcuMesh->update_parameters_atposition(relpos_m1, *cmbndFuncs_sec.pcuMesh->pGmix, Gmix);
			}

			cuBReal gI = (2.0 * (cuBReal)GMUB_2E / dh) * cuReal2(Gmix).j / (-(cuBReal)GAMMA * grel * Ms);
			cuBReal gR = (2.0 * (cuBReal)GMUB_2E / dh) * cuReal2(Gmix).i / (-(cuBReal)GAMMA * grel * Ms);

			Heff[mcell1_idx] += tsi_eff * (gI * dVs + gR * (M[mcell1_idx] ^ dVs) / Ms);
		}
	}
}

//Calculate the field resulting from interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set)
void TransportCUDA::CalculateSAInterfaceField(TransportCUDA* ptrans_sec, CMBNDInfoCUDA& contactCUDA, bool primary_top)
{
	//the top contacting mesh sets G values
	bool GInterface_Enabled = ((primary_top && pMeshCUDA->GInterface_Enabled()) || (!primary_top && ptrans_sec->pMeshCUDA->GInterface_Enabled()));
	
	if (pMeshCUDA->MComputation_Enabled() && !ptrans_sec->pMeshCUDA->Magnetisation_Enabled() && GInterface_Enabled) {

		CalculateSAInterfaceField_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (contactCUDA, ptrans_sec->poisson_Spin_S, poisson_Spin_S);
	}
}

//-------------------Display Calculation Methods

//SPIN CURRENT

__global__ void GetSpinCurrent_Kernel(int component, cuVEC<cuReal3>& displayVEC, ManagedMeshCUDA& cuMesh, TransportCUDA_Spin_S_Funcs& poisson_Spin_S)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& S = *cuMesh.pS;
	cuVEC_VC<cuReal3>& Jc = *cuMesh.pJc;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < S.linear_size()) {

		cuReal33 Js = cuReal33();

		if (S.is_not_empty(idx)) {

			cuBReal De = *cuMesh.pDe;
			cuBReal P = *cuMesh.pP;
			cuBReal Ms = *cuMesh.pMs;
			cuBReal SHA = *cuMesh.pSHA;
			cuBReal betaD = *cuMesh.pbetaD;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pP, P, *cuMesh.pDe, De, *cuMesh.pbetaD, betaD, *cuMesh.pSHA, SHA);

			if (M.linear_size()) {

				//magnetic mesh terms

				//1. drift
				int idx_M = M.position_to_cellidx(S.cellidx_to_position(idx));

				cuReal3 Mval = M[idx_M];
				cuReal33 grad_S = S.grad_neu(idx);

				Js = (Jc[idx] | Mval) * (P / Ms) * (-(cuBReal)MUB_E);

				//2. diffusion with homogeneous Neumann boundary condition
				Js -= grad_S * De;

				//3. CPP-GMR term
				if (cuIsNZ((cuBReal)betaD)) {

					cuReal3 delsq_S = S.delsq_neu(idx);

					Js += ((grad_S * Mval) | Mval) * P * betaD * De / (Ms * Ms);
				}
			}
			else {

				//non-magnetic mesh terms

				//1. SHE contribution
				Js = cu_epsilon3(Jc[idx]) * SHA * (cuBReal)MUB_E;

				//2. diffusion with non-homogeneous Neumann boundary condition
				Js -= S.grad_nneu(idx, poisson_Spin_S) * De;
			}
		}

		switch (component) {

		case 0:
			displayVEC[idx] = Js.x;
			break;
		case 1:
			displayVEC[idx] = Js.y;
			break;
		case 2:
			displayVEC[idx] = Js.z;
			break;
		}
	}
}

//SPIN TORQUE

__global__ void GetSpinTorque_Kernel(cuVEC<cuReal3>& displayVEC, ManagedMeshCUDA& cuMesh)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& S = *cuMesh.pS;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < M.linear_size()) {

		if (M.is_empty(idx)) {

			displayVEC[idx] = cuReal3();
			return;
		}

		cuBReal De = *cuMesh.pDe;
		cuBReal ts_eff = *cuMesh.pts_eff;
		cuBReal Ms = *cuMesh.pMs;
		cuBReal l_ex = *cuMesh.pl_ex;
		cuBReal l_ph = *cuMesh.pl_ph;
		cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pDe, De, *cuMesh.pts_eff, ts_eff, *cuMesh.pl_ex, l_ex, *cuMesh.pl_ph, l_ph);

		cuReal3 Sav = S.weighted_average(M.cellidx_to_position(idx), M.h);

		displayVEC[idx] = ts_eff * ((Sav ^ M[idx]) * De / (Ms * l_ex * l_ex) + (M[idx] ^ (Sav ^ M[idx])) * De / (Ms * Ms * l_ph * l_ph));
	}
}

//SPIN INTERFACE TORQUE

__global__ void CalculateDisplaySAInterfaceTorque_Kernel(CMBNDInfoCUDA& contact, TransportCUDA_Spin_S_Funcs& cmbndFuncs_sec, TransportCUDA_Spin_S_Funcs& cmbndFuncs_pri, cuVEC<cuReal3>& displayVEC)
{
	cuVEC_VC<cuReal3>& M = *cmbndFuncs_pri.pcuMesh->pM;
	cuVEC_VC<cuReal3>& S_pri = *cmbndFuncs_pri.pcuMesh->pS;
	cuVEC_VC<cuReal3>& S_sec = *cmbndFuncs_sec.pcuMesh->pS;

	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	//interface conductance method with F being the primary mesh : calculate and set spin torque

	//convert the cells box from S mesh to M mesh
	cuINT3 mbox_start = M.cellidx_from_position(S_pri.cellidx_to_position(contact.cells_box.s) + M.rect.s);
	cuINT3 mbox_end = M.cellidx_from_position(S_pri.cellidx_to_position(contact.cells_box.e - cuINT3(1)) + M.rect.s);

	if ((mbox_end.i - mbox_start.i) == 0) mbox_end.i = mbox_start.i + 1;
	if ((mbox_end.j - mbox_start.j) == 0) mbox_end.j = mbox_start.j + 1;
	if ((mbox_end.k - mbox_start.k) == 0) mbox_end.k = mbox_start.k + 1;

	cuINT3 box_sizes = mbox_end - mbox_start;

	if (box_idx < box_sizes.dim()) {

		//the cellsize perpendicular to the contact (in the M mesh)
		cuBReal dh = (cuReal3(contact.cell_shift) & M.h).norm();

		int i = (box_idx % box_sizes.x) + mbox_start.i;
		int j = ((box_idx / box_sizes.x) % box_sizes.y) + mbox_start.j;
		int k = (box_idx / (box_sizes.x * box_sizes.y)) + mbox_start.k;

		//index of magnetic cell 1
		int mcell1_idx = i + j * M.n.x + k * M.n.x*M.n.y;

		if (M.is_empty(mcell1_idx)) return;

		cuBReal Ms = *cmbndFuncs_pri.pcuMesh->pMs;
		cuBReal tsi_eff = *cmbndFuncs_pri.pcuMesh->ptsi_eff;
		cmbndFuncs_pri.pcuMesh->update_parameters_mcoarse(mcell1_idx, *cmbndFuncs_pri.pcuMesh->pMs, Ms, *cmbndFuncs_pri.pcuMesh->ptsi_eff, tsi_eff);

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

			Gmix = *cmbndFuncs_pri.pcuMesh->pGmix;
			cmbndFuncs_pri.pcuMesh->update_parameters_mcoarse(mcell1_idx, *cmbndFuncs_pri.pcuMesh->pGmix, Gmix);
		}
		else {

			Gmix = *cmbndFuncs_sec.pcuMesh->pGmix;
			cmbndFuncs_sec.pcuMesh->update_parameters_atposition(relpos_m1, *cmbndFuncs_sec.pcuMesh->pGmix, Gmix);
		}

		cuBReal gI = (2.0 * (cuBReal)GMUB_2E / dh) * Gmix.j / Ms;
		cuBReal gR = (2.0 * (cuBReal)GMUB_2E / dh) * Gmix.i / Ms;

		displayVEC[mcell1_idx] += tsi_eff * (gI * (M[mcell1_idx] ^ dVs) + gR * (M[mcell1_idx] ^ (M[mcell1_idx] ^ dVs)) / Ms);
	}
}

//Launchers

//prepare displayVEC ready for calculation of display quantity
bool TransportCUDA::PrepareDisplayVEC(DBL3 cellsize)
{
	if (pSMeshCUDA->SolveSpinCurrent() && pMeshCUDA->EComputation_Enabled()) {

		//make sure memory is allocated to the correct size
		displayVEC()->assign(cellsize, pMeshCUDA->meshRect, cuReal3(0.0));

		return true;
	}
	else displayVEC()->clear();

	return false;
}

//return x, y, or z component of spin current (component = 0, 1, or 2)
cu_obj<cuVEC<cuReal3>>& TransportCUDA::GetSpinCurrent(int component)
{
	if (!PrepareDisplayVEC(pMeshCUDA->h_e)) return displayVEC;

	GetSpinCurrent_Kernel <<< (pMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (component, displayVEC, pMeshCUDA->cuMesh, poisson_Spin_S);

	return displayVEC;
}

//return spin torque computed from spin accumulation
cu_obj<cuVEC<cuReal3>>& TransportCUDA::GetSpinTorque(void)
{
	if (!PrepareDisplayVEC(pMeshCUDA->h)) return displayVEC;

	GetSpinTorque_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (displayVEC, pMeshCUDA->cuMesh);

	return displayVEC;
}

//Calculate the interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set), accumulating result in displayVEC
void TransportCUDA::CalculateDisplaySAInterfaceTorque(TransportCUDA* ptrans_sec, CMBNDInfoCUDA& contactCUDA, bool primary_top)
{
	//the top contacting mesh sets G values
	bool GInterface_Enabled = ((primary_top && pMeshCUDA->GInterface_Enabled()) || (!primary_top && ptrans_sec->pMeshCUDA->GInterface_Enabled()));

	if (pMeshCUDA->MComputation_Enabled() && !ptrans_sec->pMeshCUDA->Magnetisation_Enabled() && GInterface_Enabled) {

		CalculateDisplaySAInterfaceTorque_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (contactCUDA, ptrans_sec->poisson_Spin_S, poisson_Spin_S, displayVEC);
	}
}

#endif

#endif