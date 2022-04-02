#include "Atom_TransportCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "Atom_MeshCUDA.h"
#include "SuperMeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

//-------------------Display Calculation Methods

//--------------------------------------------------------------- Current Density

//Current density when only charge solver is used
__global__ void Atom_CalculateCurrentDensity_Charge_Kernel(cuVEC_VC<cuReal3>& Jc, cuVEC_VC<cuBReal>& V, cuVEC_VC<cuBReal>& elC)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Jc.linear_size()) {

		//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
		if (V.is_not_empty(idx)) {

			Jc[idx] = -elC[idx] * V.grad_diri(idx);
		}
		else Jc[idx] = cuReal3(0.0);
	}
}

//if transport solver disabled we need to set displayVEC_VC directly from E and elC as Jc = elC * E
__global__ void Atom_CalculateFixedCurrentDensity_Charge_Kernel(cuVEC_VC<cuReal3>& Jc, cuVEC_VC<cuReal3>& E, cuVEC_VC<cuBReal>& elC)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Jc.linear_size()) {

		//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
		if (elC.is_not_empty(idx)) {

			Jc[idx] = elC[idx] * E[idx];
		}
		else Jc[idx] = cuReal3(0.0);
	}
}

__global__ void Atom_CalculateCurrentDensity_Spin_Kernel(cuVEC_VC<cuReal3>& Jc, ManagedAtom_MeshCUDA& cuaMesh, TransportCUDA_Spin_V_Funcs& poisson_Spin_V)
{
	cuVEC<cuReal3>& E = *cuaMesh.pE;
	cuVEC_VC<cuBReal>& V = *cuaMesh.pV;
	cuVEC_VC<cuBReal>& elC = *cuaMesh.pelC;
	cuVEC_VC<cuReal3>& S = *cuaMesh.pS;
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Jc.linear_size()) {

		//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
		if (V.is_not_empty(idx)) {

			bool cppgmr_enabled = cuIsNZ(cuaMesh.pbetaD->get0());
			bool cpump_enabled = cuIsNZ(cuaMesh.pcpump_eff->get0());
			bool the_enabled = cuIsNZ(cuaMesh.pthe_eff->get0());

			//magnetic mesh

			cuReal3 grad_V = V.grad_diri(idx);

			//1. principal term : always present
			Jc[idx] = -elC[idx] * grad_V;

			//additional contributions if enabled
			if (cppgmr_enabled || cpump_enabled || the_enabled) {

				cuBReal mu_s = *cuaMesh.pmu_s;
				cuaMesh.update_parameters_ecoarse(idx, *cuaMesh.pmu_s, mu_s);

				int idx_M = M1.position_to_cellidx(S.cellidx_to_position(idx));

				cuReal3 m = M1[idx_M] / mu_s;
				cuReal33 grad_S = S.grad_neu(idx);		//homogeneous Neumann since SHA = 0 in magnetic meshes

				//2. CPP-GMR contribution
				if (cppgmr_enabled) {

					cuBReal betaD = *cuaMesh.pbetaD;
					cuBReal De = *cuaMesh.pDe;
					cuaMesh.update_parameters_ecoarse(idx, *cuaMesh.pbetaD, betaD, *cuaMesh.pDe, De);

					Jc[idx] += (grad_S * m) * betaD * De / (cuBReal)MUB_E;
				}

				//3. topological Hall effect contribution
				//4. charge pumping contribution
				if (cpump_enabled || the_enabled) {

					cuBReal P = *cuaMesh.pP;
					cuBReal n_density = *cuaMesh.pn_density;
					cuaMesh.update_parameters_ecoarse(idx, *cuaMesh.pP, P, *cuaMesh.pn_density, n_density);

					cuReal33 grad_M = M1.grad_neu(idx_M);
					cuReal3 dx_m = grad_M.x / mu_s;
					cuReal3 dy_m = grad_M.y / mu_s;

					//topological Hall effect contribution
					if (the_enabled) {

						cuBReal Bz_the = (dx_m ^ dy_m) * m;
						Jc[idx] += cuaMesh.pthe_eff->get0() * (P * elC[idx] * (cuBReal)HBAR_E / ((cuBReal)ECHARGE * n_density)) * elC[idx] * cuReal3(grad_V.y * Bz_the, -grad_V.x *Bz_the, 0.0);
					}

					//charge pumping contribution
					if (cpump_enabled) {

						cuReal3 dm_dt = (*poisson_Spin_V.pdM_dt)[idx_M] / mu_s;
						Jc[idx] += cuaMesh.pcpump_eff->get0() * (P * elC[idx] * (cuBReal)HBAR_E / 2) * cuReal3((dm_dt ^ dx_m) * m, (dm_dt ^ dy_m) * m, 0.0);
					}
				}
			}
		}
		else Jc[idx] = cuReal3(0);
	}
}

//-------------------Calculation Methods : Charge Current Density

//calculate charge current density over the mesh
cu_obj<cuVEC_VC<cuReal3>>& Atom_TransportCUDA::GetChargeCurrent(void)
{
	if (!PrepareDisplayVEC_VC(paMeshCUDA->h_e)) return displayVEC_VC;

	if (!pSMeshCUDA->DisabledTransportSolver()) {

		if (stsolve == STSOLVE_NONE) {

			Atom_CalculateCurrentDensity_Charge_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (displayVEC_VC, paMeshCUDA->V, paMeshCUDA->elC);
		}
		else {

			Atom_CalculateCurrentDensity_Spin_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (displayVEC_VC, paMeshCUDA->cuaMesh, poisson_Spin_V);
		}
	}
	else {

		//if transport solver disabled we need to set displayVEC_VC directly from E and elC as Jc = elC * E
		Atom_CalculateFixedCurrentDensity_Charge_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (displayVEC_VC, paMeshCUDA->E, paMeshCUDA->elC);
	}

	return displayVEC_VC;
}

#endif

#endif