#include "Atom_TransportCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "Atom_MeshCUDA.h"
#include "SuperMeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

//-------------------Display Calculation Methods

//SPIN CURRENT

__global__ void Atom_GetSpinCurrent_Kernel(int component, cuVEC<cuReal3>& displayVEC, ManagedAtom_MeshCUDA& cuaMesh, TransportCUDA_Spin_S_Funcs& poisson_Spin_S)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC_VC<cuReal3>& S = *cuaMesh.pS;
	cuVEC_VC<cuReal3>& E = *cuaMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuaMesh.pelC;

	cuVEC_VC<cuReal3>& dM_dt = *poisson_Spin_S.pdM_dt;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (idx < S.linear_size()) {

		bool cpump_enabled = cuIsNZ(cuaMesh.pcpump_eff->get0());
		bool the_enabled = cuIsNZ(cuaMesh.pthe_eff->get0());

		cuReal33 Js = cuReal33();

		if (S.is_not_empty(idx)) {
			
			if (poisson_Spin_S.stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

				//magnetic mesh terms

				cuBReal mu_s = *cuaMesh.pmu_s;
				cuBReal P = *cuaMesh.pP;
				cuBReal De = *cuaMesh.pDe;
				cuaMesh.update_parameters_ecoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pP, P, *cuaMesh.pDe, De);

				//1. drift
				int idx_M = M1.position_to_cellidx(S.cellidx_to_position(idx));

				cuReal3 Mval = M1[idx_M];
				cuReal33 grad_S = S.grad_neu(idx);

				Js = (E[idx] | Mval) * (P * elC[idx] / mu_s) * (-(cuBReal)MUB_E);

				//2. diffusion with homogeneous Neumann boundary condition
				Js -= grad_S * De;

				//3. charge pumping
				//4. topological Hall effect

				if (component != 2 && (cpump_enabled || the_enabled)) {

					cuReal33 grad_m = M1.grad_neu(idx_M) / mu_s;

					//topological Hall effect contribution
					if (the_enabled) {

						cuBReal n_density = *cuaMesh.pn_density;
						cuaMesh.update_parameters_ecoarse(idx, *cuaMesh.pn_density, n_density);

						cuReal3 B = (grad_m.x ^ grad_m.y);
						Js += cuaMesh.pthe_eff->get0() * ((cuBReal)HBAR_E * (cuBReal)MUB_E * elC[idx] * elC[idx] / ((cuBReal)ECHARGE * n_density)) * cuReal33(-E[idx].y * B, E[idx].x * B, cuReal3());
					}

					//charge pumping contribution
					if (cpump_enabled) {

						//value a1
						cuReal3 dm_dt = dM_dt[idx_M] / mu_s;
						Js += cuaMesh.pcpump_eff->get0() * ((cuBReal)HBAR_E * (cuBReal)MUB_E * elC[idx] / 2) * cuReal33(dm_dt ^ grad_m.x, dm_dt ^ grad_m.y, cuReal3());
					}
				}
			}
			else {

				//non-magnetic mesh terms

				cuBReal De = *cuaMesh.pDe;
				cuBReal SHA = *cuaMesh.pSHA;
				cuaMesh.update_parameters_ecoarse(idx, *cuaMesh.pDe, De, *cuaMesh.pSHA, SHA);

				//1. SHE contribution
				Js = cu_epsilon3(E[idx]) * SHA * elC[idx] * (cuBReal)MUB_E;

				//2. diffusion with non-homogeneous Neumann boundary condition
				Js -= S.grad_nneu(idx, cu_epsilon3(E[idx]) * (SHA * elC[idx] * (cuBReal)MUB_E / De)) * De;
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

__global__ void Atom_GetSpinTorque_Kernel(cuVEC<cuReal3>& displayVEC, ManagedAtom_MeshCUDA& cuaMesh)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC_VC<cuReal3>& S = *cuaMesh.pS;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < M1.linear_size()) {

		if (M1.is_empty(idx)) {

			displayVEC[idx] = cuReal3();
			return;
		}

		cuBReal De = *cuaMesh.pDe;
		cuBReal ts_eff = *cuaMesh.pts_eff;
		cuBReal mu_s = *cuaMesh.pmu_s;
		cuBReal l_ex = *cuaMesh.pl_ex;
		cuBReal l_ph = *cuaMesh.pl_ph;
		cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pDe, De, *cuaMesh.pts_eff, ts_eff, *cuaMesh.pl_ex, l_ex, *cuaMesh.pl_ph, l_ph);

		cuReal3 Sav = S.weighted_average(M1.cellidx_to_position(idx), M1.h);

		displayVEC[idx] = ts_eff * ((Sav ^ M1[idx]) * De / (mu_s * l_ex * l_ex) + (M1[idx] ^ (Sav ^ M1[idx])) * De / (mu_s * mu_s * l_ph * l_ph));
	}
}

//SPIN INTERFACE TORQUE

__global__ void Atom_CalculateDisplaySAInterfaceTorque_Kernel(
	CMBNDInfoCUDA& contact, 
	TransportCUDA_Spin_S_Funcs& cmbndFuncs_sec, TransportCUDA_Spin_S_Funcs& cmbndFuncs_pri, 
	cuVEC<cuReal3>& displayVEC)
{
	//primary mesh is atomistic (cmbndFuncs_pri.pcuaMesh), and secondary is MeshCUDA (cmbndFuncs_sec.pcuMesh) since secondary is non-magnetic

	cuVEC_VC<cuReal3>& M1 = *cmbndFuncs_pri.pcuaMesh->pM1;
	cuVEC_VC<cuReal3>& S_pri = *cmbndFuncs_pri.pcuaMesh->pS;
	cuVEC_VC<cuReal3>& S_sec = *cmbndFuncs_sec.pcuMesh->pS;

	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	//interface conductance method with F being the primary mesh : calculate and set spin torque

	//convert the cells box from S mesh to M mesh
	cuINT3 mbox_start = M1.cellidx_from_position(S_pri.cellidx_to_position(contact.cells_box.s) + M1.rect.s);
	cuINT3 mbox_end = M1.cellidx_from_position(S_pri.cellidx_to_position(contact.cells_box.e - cuINT3(1)) + M1.rect.s) + cuINT3(1);

	if ((mbox_end.i - mbox_start.i) == 0) mbox_end.i = mbox_start.i + 1;
	if ((mbox_end.j - mbox_start.j) == 0) mbox_end.j = mbox_start.j + 1;
	if ((mbox_end.k - mbox_start.k) == 0) mbox_end.k = mbox_start.k + 1;

	cuINT3 box_sizes = mbox_end - mbox_start;
	
	if (box_idx < box_sizes.dim()) {

		//the cellsize perpendicular to the contact (in the M mesh)
		cuBReal dh = (cuReal3(contact.cell_shift) & M1.h).norm();

		int i = (box_idx % box_sizes.x) + mbox_start.i;
		int j = ((box_idx / box_sizes.x) % box_sizes.y) + mbox_start.j;
		int k = (box_idx / (box_sizes.x * box_sizes.y)) + mbox_start.k;

		//index of magnetic cell 1
		int mcell1_idx = i + j * M1.n.x + k * M1.n.x*M1.n.y;

		if (M1.is_empty(mcell1_idx)) return;

		cuBReal mu_s = *cmbndFuncs_pri.pcuaMesh->pmu_s;
		cuBReal tsi_eff = *cmbndFuncs_pri.pcuaMesh->ptsi_eff;
		cmbndFuncs_pri.pcuaMesh->update_parameters_mcoarse(mcell1_idx, *cmbndFuncs_pri.pcuaMesh->pmu_s, mu_s, *cmbndFuncs_pri.pcuaMesh->ptsi_eff, tsi_eff);

		//position at interface relative to primary mesh
		cuReal3 mhshift_primary = contact.hshift_primary.normalized() & M1.h;
		cuReal3 relpos_interf = ((cuReal3(i, j, k) + cuReal3(0.5)) & M1.h) + mhshift_primary / 2;

		cuReal3 relpos_1 = relpos_interf - contact.hshift_primary / 2;

		cuReal3 relpos_m1 = S_pri.rect.s - S_sec.rect.s + relpos_interf + contact.hshift_secondary / 2;

		cuReal3 stencil_sec = M1.h - cu_mod(mhshift_primary) + cu_mod(contact.hshift_secondary);
		cuReal3 stencil_pri = M1.h - cu_mod(mhshift_primary) + cu_mod(contact.hshift_primary);

		//S values
		cuReal3 S_1 = S_pri.weighted_average(relpos_1, stencil_pri);
		cuReal3 S_2 = S_pri.weighted_average(relpos_1 - contact.hshift_primary, stencil_pri);
		cuReal3 S_m1 = S_sec.weighted_average(relpos_m1, stencil_sec);
		cuReal3 S_m2 = S_sec.weighted_average(relpos_m1 + contact.hshift_secondary, stencil_sec);

		//c values
		cuBReal c_1 = cmbndFuncs_pri.c_func_sec(relpos_1, stencil_pri);
		cuBReal c_2 = cmbndFuncs_pri.c_func_sec(relpos_1 - contact.hshift_primary, stencil_pri);
		cuBReal c_m1 = cmbndFuncs_sec.c_func_sec(relpos_m1, stencil_sec);
		cuBReal c_m2 = cmbndFuncs_sec.c_func_sec(relpos_m1 + contact.hshift_secondary, stencil_sec);

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

			Gmix = *cmbndFuncs_sec.pcuMesh->pGmix;
			cmbndFuncs_sec.pcuMesh->update_parameters_atposition(relpos_m1, *cmbndFuncs_sec.pcuMesh->pGmix, Gmix);
		}

		cuBReal gI = (2.0 * (cuBReal)GMUB_2E / dh) * Gmix.j / mu_s;
		cuBReal gR = (2.0 * (cuBReal)GMUB_2E / dh) * Gmix.i / mu_s;

		displayVEC[mcell1_idx] += tsi_eff * (gI * (M1[mcell1_idx] ^ dVs) + gR * (M1[mcell1_idx] ^ (M1[mcell1_idx] ^ dVs)) / mu_s);
	}
}

//return x, y, or z component of spin current (component = 0, 1, or 2)
cu_obj<cuVEC<cuReal3>>& Atom_TransportCUDA::GetSpinCurrent(int component)
{
	if (!PrepareDisplayVEC(paMeshCUDA->h_e)) return displayVEC;

	if (stsolve != STSOLVE_NONE) {

		Atom_GetSpinCurrent_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (component, displayVEC, paMeshCUDA->cuaMesh, poisson_Spin_S);
	}

	return displayVEC;
}

//return spin torque computed from spin accumulation
cu_obj<cuVEC<cuReal3>>& Atom_TransportCUDA::GetSpinTorque(void)
{
	if (!PrepareDisplayVEC(paMeshCUDA->h)) return displayVEC;
	
	if (stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

		Atom_GetSpinTorque_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (displayVEC, paMeshCUDA->cuaMesh);
	}

	return displayVEC;
}

//Calculate the interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set), accumulating result in displayVEC
void Atom_TransportCUDA::CalculateDisplaySAInterfaceTorque(TransportBaseCUDA* ptrans_sec, CMBNDInfoCUDA& contactCUDA, bool primary_top)
{
	//the top contacting mesh sets G values
	bool isGInterface_Enabled = ((primary_top && GInterface_Enabled()) || (!primary_top && ptrans_sec->GInterface_Enabled()));

	if (isGInterface_Enabled && stsolve == STSOLVE_FERROMAGNETIC_ATOM && (ptrans_sec->Get_STSolveType() == STSOLVE_NORMALMETAL || ptrans_sec->Get_STSolveType() == STSOLVE_TUNNELING)) {

		//primary mesh is atomistic (poisson_Spin_S), and secondary is MeshCUDA (ptrans_sec->poisson_Spin_S) since secondary is non-magnetic
		Atom_CalculateDisplaySAInterfaceTorque_Kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> 
			(contactCUDA, ptrans_sec->poisson_Spin_S, poisson_Spin_S, displayVEC);
	}
}

#endif

#endif