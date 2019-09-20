#include "TransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_TRANSPORT

#include "MeshCUDA.h"
#include "SuperMeshCUDA.h"

#include "BorisCUDALib.cuh"

#include "MeshParamsControlCUDA.h"

//--------------------------------------------------------------- Electrical Conductivity with AMR

__global__ void CalculateElectricalConductivity_AMR_Kernel(ManagedMeshCUDA& cuMesh)
{
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;
	cuVEC<cuReal3>& Jc = *cuMesh.pJc;
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < elC.linear_size()) {

		if (elC.is_not_empty(idx)) {

			cuBReal elecCond = *cuMesh.pelecCond;
			cuBReal amrPercentage = *cuMesh.pamrPercentage;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.pelecCond, elecCond, *cuMesh.pamrPercentage, amrPercentage);

			//get current density value at this conductivity cell
			cuReal3 Jc_value = Jc[idx];

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

//--------------------------------------------------------------- Current Density

//Current density when only charge solver is used
__global__ void CalculateCurrentDensity_Charge_Kernel(cuVEC<cuReal3>& Jc, cuVEC_VC<cuBReal>& V, cuVEC_VC<cuBReal>& elC)
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

__global__ void CalculateCurrentDensity_Spin_Kernel(ManagedMeshCUDA& cuMesh, TransportCUDA_Spin_V_Funcs& poisson_Spin_V, TransportCUDA_Spin_S_Funcs& poisson_Spin_S)
{
	cuVEC<cuReal3>& Jc = *cuMesh.pJc;
	cuVEC_VC<cuBReal>& V = *cuMesh.pV;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;
	cuVEC_VC<cuReal3>& S = *cuMesh.pS;
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Jc.linear_size()) {

		//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
		if (V.is_not_empty(idx)) {

			cuBReal iSHA = *cuMesh.piSHA;
			cuBReal De = *cuMesh.pDe;
			cuBReal betaD = *cuMesh.pbetaD;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.piSHA, iSHA, *cuMesh.pDe, De, *cuMesh.pbetaD, betaD);

			//1. Ohm's law contribution
			if (cuIsZ((cuBReal)iSHA) || M.linear_size()) {

				//no iSHE contribution. Note, iSHE is not included in magnetic meshes.
				Jc[idx] = -elC[idx] * V.grad_diri(idx);
			}
			else {

				//iSHE enabled, must use non-homogeneous Neumann boundary condition for grad V -> Note homogeneous Neumann boundary conditions apply when calculating S differentials here (due to Jc.n = 0 at boundaries)
				Jc[idx] = -elC[idx] * V.grad_diri_nneu(idx, poisson_Spin_V);

				//must also add iSHE contribution -> here we must use non-homogeneous Neumann boundary conditions when calculating S differentials
				Jc[idx] += (iSHA * De / (cuBReal)MUB_E) * S.curl_nneu(idx, poisson_Spin_S);
			}

			//2. CPP-GMR contribution
			if (M.linear_size() && cuIsNZ((cuBReal)betaD)) {

				cuBReal Ms = *cuMesh.pMs;
				cuMesh.update_parameters_ecoarse(idx, *cuMesh.pMs, Ms);

				int idx_M = M.position_to_cellidx(S.cellidx_to_position(idx));

				cuReal3 Mval = M[idx_M];
				cuReal33 grad_S = S.grad_neu(idx);		//homogeneous Neumann since SHA = 0 in magnetic meshes

				Jc[idx] += (grad_S * Mval) * betaD * De / ((cuBReal)MUB_E * Ms);
			}
		}
		else Jc[idx] = cuReal3(0);
	}
}

//-------------------Calculation Methods

//calculate charge current density over the mesh
void TransportCUDA::CalculateCurrentDensity(void)
{
	if (!pSMeshCUDA->SolveSpinCurrent()) {

		CalculateCurrentDensity_Charge_Kernel <<< (pMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->Jc, pMeshCUDA->V, pMeshCUDA->elC);
	}
	else {

		CalculateCurrentDensity_Spin_Kernel <<< (pMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, poisson_Spin_V, poisson_Spin_S);
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