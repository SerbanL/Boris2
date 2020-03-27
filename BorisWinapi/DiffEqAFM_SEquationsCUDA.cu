#include "DiffEqAFMCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

#include "BorisCUDALib.h"
#include "BorisCUDALib.cuh"

//---------------------------------------- OTHER CALCULATION METHODS : GENERATE THERMAL cuVECs

//----------------------------------------

__global__ void GenerateThermalField_Kernel(cuBorisRand& prng, ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal2 grel = cuMesh.pgrel_AFM->get0();

	if (idx < cuDiffEq.pH_Thermal->linear_size() && cuIsNZ(grel.i + grel.j)) {

		cuReal3 h = cuDiffEq.pH_Thermal->h;
		cuBReal dT = *cuDiffEq.pdT;

		cuBReal Temperature;

		if (cuMesh.pTemp->linear_size()) {

			//get temperature at centre of idx M cell
			Temperature = (*cuMesh.pTemp)[cuDiffEq.pH_Thermal->cellidx_to_position(idx)];
		}
		else Temperature = (*cuMesh.pbase_temperature);

		//do not include any damping here - this will be included in the stochastic equations
		cuBReal Hth_const =  sqrt(2 * (cuBReal)BOLTZMANN * Temperature / ((cuBReal)GAMMA * grel.i * h.dim() * (cuBReal)MU0 * cuMesh.pMs_AFM->get0().i * dT));

		(*cuDiffEq.pH_Thermal)[idx] = Hth_const * cuReal3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));

		cuBReal Hth_const_2 = sqrt(2 * (cuBReal)BOLTZMANN * Temperature / ((cuBReal)GAMMA * grel.j * h.dim() * (cuBReal)MU0 * cuMesh.pMs_AFM->get0().j * dT));

		(*cuDiffEq.pH_Thermal_2)[idx] = Hth_const_2 * cuReal3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));
	}
}

//called when using stochastic equations
void DifferentialEquationAFMCUDA::GenerateThermalField(void)
{
	GenerateThermalField_Kernel <<< (pMeshCUDA->n_s.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (prng, cuDiffEq, pMeshCUDA->cuMesh);
}

//----------------------------------------

__global__ void GenerateThermalField_and_Torque_Kernel(cuBorisRand& prng, ManagedDiffEqAFMCUDA& cuDiffEq, ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal2 grel = cuMesh.pgrel_AFM->get0();

	if (idx < cuDiffEq.pH_Thermal->linear_size() && cuIsNZ(grel.i + grel.j)) {

		cuReal3 h = cuDiffEq.pH_Thermal->h;
		cuBReal dT = *cuDiffEq.pdT;

		cuBReal Temperature;

		if (cuMesh.pTemp->linear_size()) {
				
			//get temperature at centre of idx M cell
			Temperature = (*cuMesh.pTemp)[cuDiffEq.pH_Thermal->cellidx_to_position(idx)];
		}
		else Temperature = (*cuMesh.pbase_temperature);

		//A

		//1. Thermal Field
		//do not include any damping here - this will be included in the stochastic equations
		cuBReal Hth_const = sqrt(2 * (cuBReal)BOLTZMANN * Temperature / ((cuBReal)GAMMA * grel.i * h.dim() * (cuBReal)MU0 * cuMesh.pMs_AFM->get0().i * dT));

		(*cuDiffEq.pH_Thermal)[idx] = Hth_const * cuReal3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));

		//2. Thermal Torque
		//do not include any damping here - this will be included in the stochastic equations
		cuBReal Tth_const = sqrt(2 * (cuBReal)BOLTZMANN * Temperature * (cuBReal)GAMMA * grel.i * cuMesh.pMs_AFM->get0().i / ((cuBReal)MU0 * h.dim() * dT));

		(*cuDiffEq.pTorque_Thermal)[idx] = Tth_const * cuReal3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));

		//B

		//1. Thermal Field
		//do not include any damping here - this will be included in the stochastic equations
		cuBReal Hth_const_2 = sqrt(2 * (cuBReal)BOLTZMANN * Temperature / ((cuBReal)GAMMA * grel.j * h.dim() * (cuBReal)MU0 * cuMesh.pMs_AFM->get0().j * dT));

		(*cuDiffEq.pH_Thermal_2)[idx] = Hth_const_2 * cuReal3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));

		//2. Thermal Torque
		//do not include any damping here - this will be included in the stochastic equations
		cuBReal Tth_const_2 = sqrt(2 * (cuBReal)BOLTZMANN * Temperature * (cuBReal)GAMMA * grel.j * cuMesh.pMs_AFM->get0().j / ((cuBReal)MU0 * h.dim() * dT));

		(*cuDiffEq.pTorque_Thermal_2)[idx] = Tth_const_2 * cuReal3(prng.rand_gauss(0, 1), prng.rand_gauss(0, 1), prng.rand_gauss(0, 1));
	}
}

void DifferentialEquationAFMCUDA::GenerateThermalField_and_Torque(void)
{
	GenerateThermalField_and_Torque_Kernel <<< (pMeshCUDA->n_s.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (prng, cuDiffEq, pMeshCUDA->cuMesh);
}

#endif
#endif