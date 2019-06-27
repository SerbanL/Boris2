#include "HeatCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_HEAT

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"
#include "MeshParamsControlCUDA.h"

//-------------------Calculation Methods

__global__ void IterateHeatEquation_Kernel(ManagedMeshCUDA& cuMesh, cuReal* heatEq_RHS)
{
	cuVEC_VC<cuReal>& Temp = *cuMesh.pTemp; 
	cuVEC_VC<cuReal>& elC = *cuMesh.pelC;
	cuVEC<cuReal3>& Jc = *cuMesh.pJc;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Temp.linear_size()) {

		if (!Temp.is_not_empty(idx) || !Temp.is_not_cmbnd(idx)) return;

		cuReal density = *cuMesh.pdensity;
		cuReal shc = *cuMesh.pshc;
		cuReal thermCond = *cuMesh.pthermCond;
		cuMesh.update_parameters_tcoarse(idx, *cuMesh.pdensity, density, *cuMesh.pshc, shc, *cuMesh.pthermCond, thermCond);

		cuReal cro = density * shc;
		cuReal K = thermCond;

		//heat equation with Robin boundaries (based on Newton's law of cooling)
		heatEq_RHS[idx] = Temp.delsq_robin(idx, K) * K / cro;

		//add Joule heating if set
		if (Jc.linear_size()) {

			cuINT3 n_t = Temp.n;

			int i = idx % n_t.x;
			int j = (idx / n_t.x) % n_t.y;
			int k = idx / (n_t.x * n_t.y);

			cuReal3 Jc_value = Jc.weighted_average(cuINT3(i, j, k), Temp.h);
			cuReal elC_value = elC.weighted_average(cuINT3(i, j, k), Temp.h);

			//add Joule heating source term
			heatEq_RHS[idx] += (Jc_value * Jc_value) / (cro * elC_value);
		}
	}
}

__global__ void TemperatureFTCS_Kernel(cuVEC_VC<cuReal>& Temp, cuReal* heatEq_RHS, cuReal dT)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Temp.linear_size()) {

		Temp[idx] += dT * heatEq_RHS[idx];
	}
}

void HeatCUDA::IterateHeatEquation(cuReal dT)
{
	//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
	IterateHeatEquation_Kernel <<< (pMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, heatEq_RHS);

	//2. Now use forward time to advance by dT
	TemperatureFTCS_Kernel <<< (pMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->Temp, heatEq_RHS, dT);
}

#endif

#endif