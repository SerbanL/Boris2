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
			cuReal3 Jc_value = elC[idx] * E[idx];

			//get M value (M is on n, h mesh so could be different)
			cuReal3 M_value = M[elC.cellidx_to_position(idx)];

			cuBReal magnitude = Jc_value.norm() * M_value.norm();
			cuBReal dotproduct = 0.0;

			if (cuIsNZ(magnitude)) dotproduct = (Jc_value * M_value) / magnitude;

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

//--------------------------------------------------------------- Electrode Current

__global__ void CalculateElectrodeCurrent_nX_Side_Kernel(cuVEC_VC<cuBReal>& elC, cuVEC_VC<cuBReal>& V, cuBReal& current, cuBox electrode_box)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int stride = electrode_box.e.j - electrode_box.s.j;

	//negative x side of box, parse j and k
	cuINT3 ijk = cuINT3(electrode_box.s.i - 1, (idx % stride) + electrode_box.s.j, (idx / stride) + electrode_box.s.k);

	cuReal3 h_e = V.h;

	cuBReal current_ = 0.0;

	if (idx < stride * (electrode_box.e.k - electrode_box.s.k)) {

		current_ = -elC[ijk] * (V[ijk] - V[ijk - cuINT3(1, 0, 0)]) * h_e.y * h_e.z / h_e.x;
	}

	reduction_sum(0, 1, &current_, current);
}

__global__ void CalculateElectrodeCurrent_pX_Side_Kernel(cuVEC_VC<cuBReal>& elC, cuVEC_VC<cuBReal>& V, cuBReal& current, cuBox electrode_box)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int stride = electrode_box.e.j - electrode_box.s.j;

	//positive x side of box, parse j and k
	cuINT3 ijk = cuINT3(electrode_box.e.i, (idx % stride) + electrode_box.s.j, (idx / stride) + electrode_box.s.k);

	cuReal3 h_e = V.h;

	cuBReal current_ = 0.0;

	if (idx < stride * (electrode_box.e.k - electrode_box.s.k)) {

		current_ = elC[ijk] * (V[ijk + cuINT3(1, 0, 0)] - V[ijk]) * h_e.y * h_e.z / h_e.x;
	}

	reduction_sum(0, 1, &current_, current);
}

__global__ void CalculateElectrodeCurrent_nY_Side_Kernel(cuVEC_VC<cuBReal>& elC, cuVEC_VC<cuBReal>& V, cuBReal& current, cuBox electrode_box)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int stride = electrode_box.e.i - electrode_box.s.i;

	//negative y side of box, parse i and k
	cuINT3 ijk = cuINT3((idx % stride) + electrode_box.s.i, electrode_box.s.j - 1, (idx / stride) + electrode_box.s.k);

	cuReal3 h_e = V.h;

	cuBReal current_ = 0.0;

	if (idx < stride * (electrode_box.e.k - electrode_box.s.k)) {

		current_ = -elC[ijk] * (V[ijk] - V[ijk - cuINT3(0, 1, 0)]) * h_e.x * h_e.z / h_e.y;
	}

	reduction_sum(0, 1, &current_, current);
}

__global__ void CalculateElectrodeCurrent_pY_Side_Kernel(cuVEC_VC<cuBReal>& elC, cuVEC_VC<cuBReal>& V, cuBReal& current, cuBox electrode_box)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int stride = electrode_box.e.i - electrode_box.s.i;

	//positive y side of box, parse i and k
	cuINT3 ijk = cuINT3((idx % stride) + electrode_box.s.i, electrode_box.e.j, (idx / stride) + electrode_box.s.k);

	cuReal3 h_e = V.h;

	cuBReal current_ = 0.0;

	if (idx < stride * (electrode_box.e.k - electrode_box.s.k)) {

		current_ = elC[ijk] * (V[ijk + cuINT3(0, 1, 0)] - V[ijk]) * h_e.x * h_e.z / h_e.y;
	}

	reduction_sum(0, 1, &current_, current);
}

__global__ void CalculateElectrodeCurrent_nZ_Side_Kernel(cuVEC_VC<cuBReal>& elC, cuVEC_VC<cuBReal>& V, cuBReal& current, cuBox electrode_box)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int stride = electrode_box.e.i - electrode_box.s.i;

	//negative z side of box, parse i and j
	cuINT3 ijk = cuINT3((idx % stride) + electrode_box.s.i, (idx / stride) + electrode_box.s.j, electrode_box.s.k - 1);

	cuReal3 h_e = V.h;

	cuBReal current_ = 0.0;

	if (idx < stride * (electrode_box.e.j - electrode_box.s.j)) {

		current_ = -elC[ijk] * (V[ijk] - V[ijk - cuINT3(0, 0, 1)]) * h_e.x * h_e.y / h_e.z;
	}

	reduction_sum(0, 1, &current_, current);
}

__global__ void CalculateElectrodeCurrent_pZ_Side_Kernel(cuVEC_VC<cuBReal>& elC, cuVEC_VC<cuBReal>& V, cuBReal& current, cuBox electrode_box)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	int stride = electrode_box.e.i - electrode_box.s.i;

	//positive z side of box, parse i and j
	cuINT3 ijk = cuINT3((idx % stride) + electrode_box.s.i, (idx / stride) + electrode_box.s.j, electrode_box.e.k);

	cuReal3 h_e = V.h;

	cuBReal current_ = 0.0;

	if (idx < stride * (electrode_box.e.j - electrode_box.s.j)) {

		current_ = elC[ijk] * (V[ijk + cuINT3(0, 0, 1)] - V[ijk]) * h_e.x * h_e.y / h_e.z;
	}

	reduction_sum(0, 1, &current_, current);
}

cuBReal TransportCUDA::CalculateElectrodeCurrent(cuBox electrode_box)
{
	//calculate current from current density in cells just next to the box
	//Normally there is only one side of the box we can use so it's easier to separate into multiple kernels - one per side.

	//Obtain the current by reduction in the energy value

	ZeroEnergy();

	//cells on -x side
	if (electrode_box.s.i > 1) {

		size_t ker_size = (electrode_box.e.j - electrode_box.s.j) * (electrode_box.e.k - electrode_box.s.k);

		CalculateElectrodeCurrent_nX_Side_Kernel <<< (ker_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->elC, pMeshCUDA->V, energy, electrode_box);
	}

	//cells on +x side
	if (electrode_box.e.i + 1 < pMeshCUDA->n_e.i) {

		size_t ker_size = (electrode_box.e.j - electrode_box.s.j) * (electrode_box.e.k - electrode_box.s.k);

		CalculateElectrodeCurrent_pX_Side_Kernel <<< (ker_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->elC, pMeshCUDA->V, energy, electrode_box);
	}
	
	//cells on -y side
	if (electrode_box.s.j > 1) {

		size_t ker_size = (electrode_box.e.i - electrode_box.s.i) * (electrode_box.e.k - electrode_box.s.k);

		CalculateElectrodeCurrent_nY_Side_Kernel <<< (ker_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->elC, pMeshCUDA->V, energy, electrode_box);
	}

	//cells on +y side
	if (electrode_box.e.j + 1 < pMeshCUDA->n_e.j) {

		size_t ker_size = (electrode_box.e.i - electrode_box.s.i) * (electrode_box.e.k - electrode_box.s.k);

		CalculateElectrodeCurrent_pY_Side_Kernel <<< (ker_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->elC, pMeshCUDA->V, energy, electrode_box);
	}

	//cells on -z side
	if (electrode_box.s.k > 1) {

		size_t ker_size = (electrode_box.e.i - electrode_box.s.i) * (electrode_box.e.j - electrode_box.s.j);

		CalculateElectrodeCurrent_nZ_Side_Kernel <<< (ker_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->elC, pMeshCUDA->V, energy, electrode_box);
	}

	//cells on +z side
	if (electrode_box.e.k + 1 < pMeshCUDA->n_e.k) {

		size_t ker_size = (electrode_box.e.i - electrode_box.s.i) * (electrode_box.e.j - electrode_box.s.j);

		CalculateElectrodeCurrent_pZ_Side_Kernel <<< (ker_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->elC, pMeshCUDA->V, energy, electrode_box);
	}
	
	//energy has the current value; reset it after as we don't want to count it to the total energy density
	double current = energy.to_cpu();
	ZeroEnergy();

	return current;
}

#endif

#endif