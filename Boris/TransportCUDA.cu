#include "TransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "MeshCUDA.h"
#include "SuperMeshCUDA.h"

#include "BorisCUDALib.cuh"

#include "MeshParamsControlCUDA.h"

//--------------------------------------------------------------- Electrical Conductivity with AMR

__global__ void CalculateElectricalConductivity_AMR_Kernel(ManagedMeshCUDA& cuMesh)
{
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;
	cuVEC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < elC.linear_size()) {

		if (elC.is_not_empty(idx)) {

			cuBReal elecCond = *cuMesh.pelecCond;
			cuBReal amrPercentage = *cuMesh.pamrPercentage;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.pelecCond, elecCond, *cuMesh.pamrPercentage, amrPercentage);

			//get current density value at this conductivity cell
			cuReal3 jc_value = cu_normalize(elC[idx] * E[idx]);

			//get M value (M is on n, h mesh so could be different)
			cuReal3 m_value = cu_normalize(M[elC.cellidx_to_position(idx)]);

			cuBReal dotproduct = jc_value * m_value;

			elC[idx] = elecCond / (1 + amrPercentage * dotproduct*dotproduct / 100);
		}
	}
}

//calculate electrical conductivity with AMR present
void TransportCUDA::CalculateElectricalConductivity_AMR(void)
{
	CalculateElectricalConductivity_AMR_Kernel <<< (pMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh);
}

//--------------------------------------------------------------- Electrical Conductivity without AMR

__global__ void CalculateElectricalConductivity_NoAMR_Kernel(ManagedMeshCUDA& cuMesh)
{
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < elC.linear_size()) {

		if (elC.is_not_empty(idx)) {

			cuBReal elecCond = *cuMesh.pelecCond;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.pelecCond, elecCond);

			elC[idx] = elecCond;
		}
	}
}

//calculate electrical conductivity without AMR
void TransportCUDA::CalculateElectricalConductivity_NoAMR(void)
{
	CalculateElectricalConductivity_NoAMR_Kernel <<< (pMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh);
}

//--------------------------------------------------------------- Electric Field

//Current density when only charge solver is used
__global__ void CalculateElectricField_Charge_Kernel(cuVEC<cuReal3>& E, cuVEC_VC<cuBReal>& V)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < V.linear_size()) {

		//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
		if (V.is_not_empty(idx)) {

			E[idx] = -1.0 * V.grad_diri(idx);
		}
		else E[idx] = cuReal3(0.0);
	}
}

__global__ void CalculateElectricField_Spin_withISHE_Kernel(ManagedMeshCUDA& cuMesh, TransportCUDA_Spin_V_Funcs& poisson_Spin_V, TransportCUDA_Spin_S_Funcs& poisson_Spin_S)
{
	cuVEC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& V = *cuMesh.pV;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;
	cuVEC_VC<cuReal3>& S = *cuMesh.pS;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < E.linear_size()) {

		//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
		if (V.is_not_empty(idx)) {

			if (cuIsNZ(cuMesh.piSHA->get0())) {

				//ISHE enabled - use nonhomogeneous Neumann boundary conditions
				cuBReal iSHA = *cuMesh.piSHA;
				cuBReal De = *cuMesh.pDe;
				cuMesh.update_parameters_ecoarse(idx, *cuMesh.piSHA, iSHA, *cuMesh.pDe, De);

				E[idx] = -1.0 * V.grad_diri_nneu(idx, (iSHA * De / ((cuBReal)MUB_E * elC[idx])) * S.curl_neu(idx));
			}
			else {

				//iSHA is zero so don't need to calculate ISHE
				E[idx] = -1.0 * V.grad_diri(idx);
			}
		}
		else E[idx] = cuReal3(0);
	}
}

//-------------------Calculation Methods : Electric Field

//calculate electric field as the negative gradient of V
void TransportCUDA::CalculateElectricField(void)
{
	if (stsolve == STSOLVE_NONE || stsolve == STSOLVE_FERROMAGNETIC) {

		CalculateElectricField_Charge_Kernel << < (pMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->E, pMeshCUDA->V);
	}
	else {

		CalculateElectricField_Spin_withISHE_Kernel << < (pMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, poisson_Spin_V, poisson_Spin_S);
	}
}

#endif

#endif