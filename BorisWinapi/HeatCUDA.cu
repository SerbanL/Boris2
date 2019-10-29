#include "HeatCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_HEAT

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"
#include "MeshParamsControlCUDA.h"

//-------------------Calculation Methods

__global__ void IterateHeatEquation_Kernel(ManagedMeshCUDA& cuMesh, cuBReal* heatEq_RHS)
{
	cuVEC_VC<cuBReal>& Temp = *cuMesh.pTemp; 
	cuVEC_VC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Temp.linear_size()) {

		if (!Temp.is_not_empty(idx) || !Temp.is_not_cmbnd(idx)) return;

		cuBReal density = *cuMesh.pdensity;
		cuBReal shc = *cuMesh.pshc;
		cuBReal thermCond = *cuMesh.pthermCond;
		cuMesh.update_parameters_tcoarse(idx, *cuMesh.pdensity, density, *cuMesh.pshc, shc, *cuMesh.pthermCond, thermCond);

		cuBReal cro = density * shc;
		cuBReal K = thermCond;

		//heat equation with Robin boundaries (based on Newton's law of cooling)
		heatEq_RHS[idx] = Temp.delsq_robin(idx, K) * K / cro;

		//add Joule heating if set
		if (E.linear_size()) {

			cuReal3 position = Temp.cellidx_to_position(idx);

			cuReal3 E_value = E.weighted_average(position, Temp.h);
			cuBReal elC_value = elC.weighted_average(position, Temp.h);

			//add Joule heating source term
			heatEq_RHS[idx] += (elC_value * E_value * E_value) / cro;
		}

		//add heat source contribution if set
		if (cuIsNZ(cuMesh.pQ->get0())) {

			cuBReal Q = *cuMesh.pQ;
			cuMesh.update_parameters_tcoarse(idx, *cuMesh.pQ, Q);

			heatEq_RHS[idx] += Q / cro;
		}
	}
}

__global__ void TemperatureFTCS_Kernel(cuVEC_VC<cuBReal>& Temp, cuBReal* heatEq_RHS, cuBReal dT)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Temp.linear_size()) {

		Temp[idx] += dT * heatEq_RHS[idx];
	}
}

void HeatCUDA::IterateHeatEquation(cuBReal dT)
{
	//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
	IterateHeatEquation_Kernel <<< (pMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, heatEq_RHS);

	//2. Now use forward time to advance by dT
	TemperatureFTCS_Kernel <<< (pMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->Temp, heatEq_RHS, dT);
}

//-------------------Setters

//non-uniform temperature setting
__global__ void SetBaseTemperature_Nonuniform_Kernel(ManagedMeshCUDA& cuMesh, cuBReal Temperature)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuVEC_VC<cuBReal>& Temp = *cuMesh.pTemp;

	if (idx < Temp.linear_size()) {

		if (Temp.is_not_empty(idx)) {

			cuBReal cT = *cuMesh.pcT;
			cuMesh.update_parameters_tcoarse(idx, *cuMesh.pcT, cT);

			Temp[idx] = cT * Temperature;
		}
	}
}

//set Temp non-uniformly as specified through the cT mesh parameter
void HeatCUDA::SetBaseTemperature_Nonuniform(cuBReal Temperature)
{
	SetBaseTemperature_Nonuniform_Kernel <<< (pMeshCUDA->n_t.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, Temperature);
}

//set Temp uniformly to base temperature
void HeatCUDA::SetBaseTemperature(cuBReal Temperature)
{
	pMeshCUDA->Temp()->setnonempty(Temperature);
}


#endif

#endif