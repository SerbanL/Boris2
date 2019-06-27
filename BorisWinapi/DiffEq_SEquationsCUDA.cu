#include "DiffEqCUDA.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.h"
#include "BorisCUDALib.cuh"

//---------------------------------------- OTHER CALCULATION METHODS : GENERATE THERMAL cuVECs

//----------------------------------------

__global__ void GenerateThermalField_Kernel(cuBorisRand& prng, ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal grel = cuMesh.pgrel->get0();

	if (idx < cuDiffEq.pH_Thermal->linear_size() && cuIsNZ(grel)) {

		if (cuMesh.pM->is_not_empty(idx) && !cuMesh.pM->is_skipcell(idx)) {

			cuReal3 h = cuMesh.pM->h;
			cuReal dT = *cuDiffEq.pdT;

			cuReal Temperature;

			if (cuMesh.pTemp->linear_size()) {

				//get temperature at centre of idx M cell
				Temperature = (*cuMesh.pTemp)[cuMesh.pM->cellidx_to_position(idx)];
			}
			else Temperature = (*cuMesh.pbase_temperature);

			//do not include any damping here - this will be included in the stochastic equations
			cuReal Mag = prng.rand() * sqrt(2 * (cuReal)BOLTZMANN * Temperature / ((cuReal)GAMMA * grel * h.dim() * (cuReal)MU0 * cuMesh.pMs->get0() * dT));
			cuReal theta = prng.rand() * (cuReal)TWO_PI;
			cuReal phi = prng.rand() * (cuReal)TWO_PI;

			(*cuDiffEq.pH_Thermal)[idx] = Mag * cuReal3(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
		}
	}
}

//called when using stochastic equations
void DifferentialEquationCUDA::GenerateThermalField(void)
{
	GenerateThermalField_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (prng, cuDiffEq, pMeshCUDA->cuMesh);
}

//----------------------------------------

__global__ void GenerateThermalField_and_Torque_Kernel(cuBorisRand& prng, ManagedDiffEqCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal grel = cuMesh.pgrel->get0();

	if (idx < cuDiffEq.pH_Thermal->linear_size() && cuIsNZ(grel)) {

		if (cuMesh.pM->is_not_empty(idx) && !cuMesh.pM->is_skipcell(idx)) {

			cuReal3 h = cuMesh.pM->h;
			cuReal dT = *cuDiffEq.pdT;

			cuReal Temperature;

			if (cuMesh.pTemp->linear_size()) {
				
				//get temperature at centre of idx M cell
				Temperature = (*cuMesh.pTemp)[cuMesh.pM->cellidx_to_position(idx)];
			}
			else Temperature = (*cuMesh.pbase_temperature);

			//1. Thermal Field
			//do not include any damping here - this will be included in the stochastic equations
			cuReal Mag = prng.rand() * sqrt(2 * (cuReal)BOLTZMANN * Temperature / ((cuReal)GAMMA * grel * h.dim() * (cuReal)MU0 * cuMesh.pMs->get0() * dT));
			cuReal theta = prng.rand() * (cuReal)TWO_PI;
			cuReal phi = prng.rand() * (cuReal)TWO_PI;

			(*cuDiffEq.pH_Thermal)[idx] = Mag * cuReal3(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));

			//2. Thermal Torque
			//do not include any damping here - this will be included in the stochastic equations
			Mag = prng.rand() * sqrt(2 * (cuReal)BOLTZMANN * Temperature * (cuReal)GAMMA * grel * cuMesh.pMs->get0() / ((cuReal)MU0 * h.dim() * dT));
			theta = prng.rand() * (cuReal)TWO_PI;
			phi = prng.rand() * (cuReal)TWO_PI;

			(*cuDiffEq.pTorque_Thermal)[idx] = Mag * cuReal3(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
		}
	}
}

void DifferentialEquationCUDA::GenerateThermalField_and_Torque(void)
{
	GenerateThermalField_and_Torque_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (prng, cuDiffEq, pMeshCUDA->cuMesh);
}

#endif