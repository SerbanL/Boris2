#include "Atom_TransportCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"
#include "SuperMeshCUDA.h"

#include "BorisCUDALib.cuh"

#include "Atom_MeshParamsControlCUDA.h"

//--------------------------------------------------------------- Electrical Conductivity with AMR

__global__ void Atom_CalculateElectricalConductivity_AMR_Kernel(ManagedAtom_MeshCUDA& cuMesh)
{
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;
	cuVEC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuReal3>& M1 = *cuMesh.pM1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < elC.linear_size()) {

		if (elC.is_not_empty(idx)) {

			cuBReal elecCond = *cuMesh.pelecCond;
			cuBReal amrPercentage = *cuMesh.pamrPercentage;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.pelecCond, elecCond, *cuMesh.pamrPercentage, amrPercentage);

			//get current density value at this conductivity cell
			cuReal3 Jc_value = elC[idx] * E[idx];

			//get M value (M is on n, h mesh so could be different)
			cuReal3 M_value = M1[elC.cellidx_to_position(idx)];

			cuBReal magnitude = Jc_value.norm() * M_value.norm();
			cuBReal dotproduct = 0.0;

			if (cuIsNZ(magnitude)) dotproduct = (Jc_value * M_value) / magnitude;

			elC[idx] = elecCond / (1 + amrPercentage * dotproduct*dotproduct / 100);
		}
	}
}

//calculate electrical conductivity with AMR present
void Atom_TransportCUDA::CalculateElectricalConductivity_AMR(void)
{
	Atom_CalculateElectricalConductivity_AMR_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh);
}

//--------------------------------------------------------------- Electrical Conductivity without AMR

__global__ void Atom_CalculateElectricalConductivity_NoAMR_Kernel(ManagedAtom_MeshCUDA& cuMesh)
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
void Atom_TransportCUDA::CalculateElectricalConductivity_NoAMR(void)
{
	Atom_CalculateElectricalConductivity_NoAMR_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh);
}

//--------------------------------------------------------------- Electric Field

//Current density when only charge solver is used
__global__ void Atom_CalculateElectricField_Charge_Kernel(cuVEC<cuReal3>& E, cuVEC_VC<cuBReal>& V)
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

//-------------------Calculation Methods : Electric Field

//calculate electric field as the negative gradient of V
void Atom_TransportCUDA::CalculateElectricField(void)
{
	Atom_CalculateElectricField_Charge_Kernel <<< (paMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->E, paMeshCUDA->V);
}

#endif

#endif