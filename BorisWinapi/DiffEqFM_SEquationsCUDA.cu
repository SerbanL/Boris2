#include "DiffEqFMCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_FERROMAGNETIC

#include "BorisCUDALib.h"
#include "BorisCUDALib.cuh"

//---------------------------------------- OTHER CALCULATION METHODS : GENERATE THERMAL cuVECs

//----------------------------------------

__global__ void GenerateThermalField_Kernel(cuBorisRand& prng, ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal grel = cuMesh.pgrel->get0();

	if (idx < cuDiffEq.pH_Thermal->linear_size() && cuIsNZ(grel)) {

		cuReal3 position = cuDiffEq.pH_Thermal->cellidx_to_position(idx);

		if (cuMesh.pM->is_not_empty(position) && !cuMesh.pM->is_skipcell(position)) {

			cuReal3 h = cuDiffEq.pH_Thermal->h;
			cuBReal dT = *cuDiffEq.pdT;

			cuBReal Temperature;

			if (cuMesh.pTemp->linear_size()) {

				//get temperature at centre of idx M cell
				Temperature = (*cuMesh.pTemp)[position];
			}
			else Temperature = (*cuMesh.pbase_temperature);

			//do not include any damping here - this will be included in the stochastic equations
			cuBReal Hth_const = sqrt(2 * (cuBReal)BOLTZMANN * Temperature / ((cuBReal)GAMMA * grel * h.dim() * (cuBReal)MU0 * cuMesh.pMs->get0() * dT));
			
			(*cuDiffEq.pH_Thermal)[idx] = Hth_const * cuReal3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));
		}
	}
}

//called when using stochastic equations
void DifferentialEquationFMCUDA::GenerateThermalField(void)
{
	GenerateThermalField_Kernel <<< (pMeshCUDA->n_s.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (prng, cuDiffEq, pMeshCUDA->cuMesh);
}

//----------------------------------------

__global__ void GenerateThermalField_and_Torque_Kernel(cuBorisRand& prng, ManagedDiffEqFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal grel = cuMesh.pgrel->get0();

	if (idx < cuDiffEq.pH_Thermal->linear_size() && cuIsNZ(grel)) {

		cuReal3 position = cuDiffEq.pH_Thermal->cellidx_to_position(idx);

		if (cuMesh.pM->is_not_empty(position) && !cuMesh.pM->is_skipcell(position)) {

			cuReal3 h = cuDiffEq.pH_Thermal->h;
			cuBReal dT = *cuDiffEq.pdT;

			cuBReal Temperature;

			if (cuMesh.pTemp->linear_size()) {
				
				//get temperature at centre of idx M cell
				Temperature = (*cuMesh.pTemp)[position];
			}
			else Temperature = (*cuMesh.pbase_temperature);
			
			//1. Thermal Field
			//do not include any damping here - this will be included in the stochastic equations
			cuBReal Hth_const = sqrt(2 * (cuBReal)BOLTZMANN * Temperature / ((cuBReal)GAMMA * grel * h.dim() * (cuBReal)MU0 * cuMesh.pMs->get0() * dT));		

			(*cuDiffEq.pH_Thermal)[idx] = Hth_const * cuReal3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));
			
			//2. Thermal Torque
			//do not include any damping here - this will be included in the stochastic equations
			cuBReal Tth_const = sqrt(2 * (cuBReal)BOLTZMANN * Temperature * (cuBReal)GAMMA * grel * cuMesh.pMs->get0() / ((cuBReal)MU0 * h.dim() * dT));
			
			(*cuDiffEq.pTorque_Thermal)[idx] = Tth_const * cuReal3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));
		}
	}
}

void DifferentialEquationFMCUDA::GenerateThermalField_and_Torque(void)
{
	GenerateThermalField_and_Torque_Kernel <<< (pMeshCUDA->n_s.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (prng, cuDiffEq, pMeshCUDA->cuMesh);
}

#endif
#endif