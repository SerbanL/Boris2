#include "TMRCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TMR

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"
#include "SuperMeshCUDA.h"
#include "MeshParamsControlCUDA.h"

//-------------------Display Calculation Methods

//--------------------------------------------------------------- Current Density

//Current density when only charge solver is used
__global__ void TMR_CalculateCurrentDensity_Charge_Kernel(cuVEC_VC<cuReal3>& Jc, cuVEC_VC<cuBReal>& V, cuVEC_VC<cuBReal>& elC)
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
__global__ void TMR_CalculateFixedCurrentDensity_Charge_Kernel(cuVEC_VC<cuReal3>& Jc, cuVEC_VC<cuReal3>& E, cuVEC_VC<cuBReal>& elC)
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

//-------------------Calculation Methods : Charge Current Density

//calculate charge current density over the mesh
cu_obj<cuVEC_VC<cuReal3>>& TMRCUDA::GetChargeCurrent(void)
{
	if (!PrepareDisplayVEC_VC(pMeshCUDA->h_e)) return displayVEC_VC;

	if (!pSMeshCUDA->DisabledTransportSolver()) {

		TMR_CalculateCurrentDensity_Charge_Kernel <<< (pMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (displayVEC_VC, pMeshCUDA->V, pMeshCUDA->elC);
	}
	else {

		//if transport solver disabled we need to set displayVEC_VC directly from E and elC as Jc = elC * E
		TMR_CalculateFixedCurrentDensity_Charge_Kernel <<< (pMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (displayVEC_VC, pMeshCUDA->E, pMeshCUDA->elC);
	}

	return displayVEC_VC;
}

#endif

#endif