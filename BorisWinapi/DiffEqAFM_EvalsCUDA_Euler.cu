#include "DiffEqAFMCUDA.h"


#if COMPILECUDA == 1
#ifdef ODE_EVAL_EULER
#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

#include "MeshParamsControlCUDA.h"

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//----------------------------------------- EVALUATIONS: Euler

__global__ void RunEuler_Kernel_withReductions(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	cuReal3 mxh = cuReal3();
	cuReal3 dmdt = cuReal3();

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];
			(*cuDiffEq.psM1_2)[idx] = (*cuMesh.pM2)[idx];

			if (!cuMesh.pM->is_skipcell(idx)) {

				//obtain average normalized torque term
				cuBReal Mnorm = (*cuMesh.pM)[idx].norm();
				mxh = ((*cuMesh.pM)[idx] ^ (*cuMesh.pHeff)[idx]) / (Mnorm * Mnorm);

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs2;
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, rhs2);

				//Now estimate magnetization for the next time step
				(*cuMesh.pM)[idx] += rhs * dT;
				(*cuMesh.pM2)[idx] += rhs2 * dT;

				if (*cuDiffEq.prenormalize) {

					cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
					cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
					(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
					(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
				}

				//obtain maximum normalized dmdt term
				dmdt = ((*cuMesh.pM)[idx] - (*cuDiffEq.psM1)[idx]) / (dT * (cuBReal)GAMMA * Mnorm * Mnorm);
			}
			else {

				cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
				(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	//only reduce for dmdt (and mxh) if grel is not zero (if it's zero this means magnetisation dynamics is disabled in this mesh)
	if (cuMesh.pgrel->get0()) {

		reduction_avg(0, 1, &mxh, *cuDiffEq.pmxh_av, *cuDiffEq.pavpoints);
		reduction_avg(0, 1, &dmdt, *cuDiffEq.pdmdt_av, *cuDiffEq.pavpoints2);
	}
}

__global__ void RunEuler_Kernel(ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal dT = *cuDiffEq.pdT;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			//Save current magnetization
			(*cuDiffEq.psM1)[idx] = (*cuMesh.pM)[idx];
			(*cuDiffEq.psM1_2)[idx] = (*cuMesh.pM2)[idx];

			if (!cuMesh.pM->is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				cuReal3 rhs2;
				cuReal3 rhs = (cuDiffEq.*(cuDiffEq.pODEFunc))(idx, rhs2);

				//Now estimate magnetization for the next time step
				(*cuMesh.pM)[idx] += rhs * dT;
				(*cuMesh.pM2)[idx] += rhs2 * dT;

				if (*cuDiffEq.prenormalize) {

					cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
					cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
					(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
					(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
				}
			}
			else {

				cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
				cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);
				(*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
				(*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
			}
		}
	}
}

//----------------------------------------- DifferentialEquationCUDA Launchers

//EULER

void DifferentialEquationAFMCUDA::RunEuler(bool calculate_mxh, bool calculate_dmdt)
{
	if (calculate_mxh || calculate_dmdt) {

		RunEuler_Kernel_withReductions << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
	else {

		RunEuler_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuDiffEq, pMeshCUDA->cuMesh);
	}
}

#endif
#endif
#endif