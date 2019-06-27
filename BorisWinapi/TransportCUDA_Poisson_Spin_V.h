#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_TRANSPORT

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "ManagedMeshCUDA.h"
#include "TransportCUDA_Poisson_Spin_S.h"

#include "MeshParamsControlCUDA.h"

class MeshCUDA;

//This is held as a cu_obj managed class in TransportCUDA modules
//It provides methods and access to mesh data for use in cuVEC_VC methods.
//The methods have fixed names, e.g. Poisson_RHS is used by Poisson solvers to evaluate the r.h.s. of the Poisson equation
//The a_func, b_func and diff2_func methods are used to set CMBND conditions based on the continuity of a quantity and a flux.
//If V is the potential, then the flux is the function f(V) = a_func + b_func * V', where the V' differential direction is perpendicular to the interface.
//This particular class is used for charge transport within the spin current solver.
class TransportCUDA_Spin_V_Funcs {

public:

	//managed mesh for access to all required mesh VECs and material parameters
	ManagedMeshCUDA* pcuMesh;

	//need this (held in TransportCUDA) to access non-homogeneous Neumann boundary condition for S
	TransportCUDA_Spin_S_Funcs* pPoisson_Spin_S;

public:

	__host__ void construct_cu_obj(void) {}
	__host__ void destruct_cu_obj(void) {}

	BError set_pointers(MeshCUDA* pMeshCUDA, cu_obj<TransportCUDA_Spin_S_Funcs>& poisson_Spin_S);

	//this evaluates the Poisson RHS when solving the Poisson equation on V (in the context of full spin solver)
	__device__ cuReal Poisson_RHS(int idx)
	{
		cuVEC_VC<cuReal>& V = *pcuMesh->pV;
		cuVEC_VC<cuReal>& elC = *pcuMesh->pelC;
		cuVEC_VC<cuReal3>& S = *pcuMesh->pS;
		cuVEC_VC<cuReal3>& M = *pcuMesh->pM;

		cuReal iSHA = *pcuMesh->piSHA;
		cuReal betaD = *pcuMesh->pbetaD;
		pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->piSHA, iSHA, *pcuMesh->pbetaD, betaD);

		//We are solving the Poisson equation del_sq V = -grad sigma * grad V / sigma + betaD*De*e * del ((grad S) m) / sigma * muB

		//The Poisson solver calls this method to evaluate the RHS of this equation
			
		cuReal value = 0.0;

		if (cuIsZ((cuReal)iSHA) || M.linear_size()) {

			//no iSHE contribution. Note, iSHE is not included in magnetic meshes.
			value = -(V.grad_diri(idx) * elC.grad_sided(idx)) / elC[idx];
		}
		else {

			//iSHE enabled, must use non-homogeneous Neumann boundary condition for grad V -> Note homogeneous Neumann boundary conditions apply when calculating S differentials here (due to Jc.n = 0 at boundaries)
			value = -(V.grad_diri_nneu(idx, *this) * elC.grad_sided(idx)) / elC[idx];
		}

		//CPP-GMR contribution in magnetic meshes
		if (M.linear_size() && cuIsNZ((cuReal)betaD)) {

			cuReal Ms = *pcuMesh->pMs;
			cuReal De = *pcuMesh->pDe;
			pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->pMs, Ms, *pcuMesh->pDe, De);

			int idx_M = M.position_to_cellidx(S.cellidx_to_position(idx));

			cuReal33 grad_M = M.grad_neu(idx_M);
			cuReal33 grad_S = S.grad_neu(idx);
			cuReal3 delsq_S = S.delsq_neu(idx);
			cuReal div_grad_S_M = (grad_S.i * grad_M.i) + (grad_S.j * grad_M.j) + (grad_S.k * grad_M.k) + (M[idx_M] * delsq_S);

			value += div_grad_S_M * betaD * De / ((cuReal)MUB_E * elC[idx] * Ms);
		}

		return value;
	}

	//boundary differential of V for non-homogeneous Neumann boundary conditions
	__device__ cuVAL3<cuReal> bdiff(int idx)
	{
		cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
		cuVEC_VC<cuReal>& elC = *pcuMesh->pelC;
		cuVEC_VC<cuReal3>& S = *pcuMesh->pS;

		if (M.linear_size()) return cuReal3();

		cuReal De = *pcuMesh->pDe;
		cuReal iSHA = *pcuMesh->piSHA;
		pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->piSHA, iSHA, *pcuMesh->pDe, De);

		return (iSHA * De / ((cuReal)MUB_E * elC[idx])) * S.curl_neu(idx);
	}

	//Functions used for calculating CMBND values

	//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface	
	__device__ cuReal a_func_pri(int cell1_idx, int cell2_idx, cuReal3 shift)
	{
		cuVEC_VC<cuReal3>& S = *pcuMesh->pS;
		cuVEC_VC<cuReal3>& Jc = *pcuMesh->pJc;
		cuVEC_VC<cuReal3>& M = *pcuMesh->pM;

		TransportCUDA_Spin_S_Funcs& poisson_Spin_S = *pPoisson_Spin_S;

		cuReal a = 0.0;

		cuReal3 u = shift.normalized() * -1;
		
		cuReal iSHA = *pcuMesh->piSHA;
		cuReal betaD = *pcuMesh->pbetaD;
		pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pbetaD, betaD, *pcuMesh->piSHA, iSHA);

		if (M.linear_size() && cuIsNZ((cuReal)betaD)) {

			cuReal Ms = *pcuMesh->pMs;
			cuReal De = *pcuMesh->pDe;
			pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pMs, Ms, *pcuMesh->pDe, De);

			//need to find value at boundary so use interpolation

			//value a1
			int idx_M1 = M.position_to_cellidx(S.cellidx_to_position(cell1_idx));
			cuReal3 M1 = M[idx_M1];
			cuReal33 grad_S1 = S.grad_neu(cell1_idx);

			cuReal a1 = ((grad_S1 * M1) * betaD * De / ((cuReal)MUB_E * Ms)) * u;

			//value a2
			int idx_M2 = M.position_to_cellidx(S.cellidx_to_position(cell2_idx));
			cuReal3 M2 = M[idx_M2];
			cuReal33 grad_S2 = S.grad_neu(cell2_idx);

			cuReal a2 = ((grad_S2 * M2) * betaD * De / ((cuReal)MUB_E * Ms)) * u;

			//final interpolated a value
			a = (1.5 * a1 - 0.5 * a2);
		}

		//iSHE contribution
		if (cuIsNZ((cuReal)iSHA) && !M.linear_size()) {

			cuReal De = *pcuMesh->pDe;
			pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pDe, De);

			//need to find value at boundary so use interpolation

			//value a1
			cuReal3 Jc1 = Jc[cell1_idx];

			cuReal a1 = (iSHA * De / (cuReal)MUB_E) * S.curl_nneu(cell1_idx, poisson_Spin_S) * u;

			//value a2
			cuReal3 Jc2 = Jc[cell2_idx];

			cuReal a2 = (iSHA * De / (cuReal)MUB_E) * S.curl_nneu(cell2_idx, poisson_Spin_S) * u;

			//final interpolated a value
			a += (1.5 * a1 - 0.5 * a2);
		}

		return a;
	}

	__device__ cuReal a_func_sec(cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil)
	{
		cuVEC_VC<cuReal3>& S = *pcuMesh->pS;
		cuVEC_VC<cuReal3>& Jc = *pcuMesh->pJc;
		cuVEC_VC<cuReal3>& M = *pcuMesh->pM;

		TransportCUDA_Spin_S_Funcs& poisson_Spin_S = *pPoisson_Spin_S;

		cuReal a = 0.0;

		cuReal3 u = shift.normalized() * -1;

		int idx_S1 = S.position_to_cellidx(relpos_m1);
		int idx_S2 = S.position_to_cellidx(relpos_m1 + shift);

		cuReal iSHA = *pcuMesh->piSHA;
		cuReal betaD = *pcuMesh->pbetaD;
		pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pbetaD, betaD, *pcuMesh->piSHA, iSHA);

		//CPP-GMR term
		if (M.linear_size() && cuIsNZ((cuReal)betaD)) {

			cuReal Ms = *pcuMesh->pMs;
			cuReal De = *pcuMesh->pDe;
			pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pMs, Ms, *pcuMesh->pDe, De);

			//need to find value at boundary so use interpolation

			//value a1
			cuReal3 M1 = M.weighted_average(relpos_m1, stencil);
			cuReal33 grad_S1 = S.grad_neu(idx_S1);

			cuReal a1 = ((grad_S1 * M1) * betaD * De / ((cuReal)MUB_E * Ms)) * u;

			//value a2
			cuReal3 M2 = M.weighted_average(relpos_m1 + shift, stencil);
			cuReal33 grad_S2 = S.grad_neu(idx_S2);

			cuReal a2 = ((grad_S2 * M2) * betaD * De / ((cuReal)MUB_E * Ms)) * u;

			//final interpolated a value
			a = (1.5 * a1 - 0.5 * a2);
		}

		//iSHE contribution
		if (cuIsNZ((cuReal)iSHA) && !M.linear_size()) {

			cuReal De = *pcuMesh->pDe;
			pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pDe, De);

			//need to find value at boundary so use interpolation

			//value a1
			cuReal3 Jc1 = Jc.weighted_average(relpos_m1, stencil);

			cuReal a1 = (iSHA * De / (cuReal)MUB_E) * S.curl_nneu(idx_S1, poisson_Spin_S) * u;

			//value a2
			cuReal3 Jc2 = Jc.weighted_average(relpos_m1 + shift, stencil);

			cuReal a2 = (iSHA * De / (cuReal)MUB_E) * S.curl_nneu(idx_S2, poisson_Spin_S) * u;

			//final interpolated a value
			a += (1.5 * a1 - 0.5 * a2);
		}

		return a;
	}

	//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface	
	__device__ cuReal b_func_pri(int cell1_idx, int cell2_idx)
	{
		cuVEC_VC<cuReal>& elC = *pcuMesh->pelC;

		return (-1.5 * elC[cell1_idx] + 0.5 * elC[cell2_idx]);
	}

	__device__ cuReal b_func_sec(cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil)
	{
		cuVEC_VC<cuReal>& elC = *pcuMesh->pelC;

		return (-1.5 * elC.weighted_average(relpos_m1, stencil) + 0.5 * elC.weighted_average(relpos_m1 + shift, stencil));
	}

	//second order differential of V at cells either side of the boundary; delsq V = -grad V * grad elC / elC
	__device__ cuReal diff2_pri(int cell1_idx)
	{
		cuVEC_VC<cuReal>& V = *pcuMesh->pV;
		cuVEC_VC<cuReal>& elC = *pcuMesh->pelC;
		cuVEC_VC<cuReal3>& S = *pcuMesh->pS;
		cuVEC_VC<cuReal3>& M = *pcuMesh->pM;

		cuReal iSHA = *pcuMesh->piSHA;
		cuReal betaD = *pcuMesh->pbetaD;
		pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->piSHA, iSHA, *pcuMesh->pbetaD, betaD);

		cuReal value = 0.0;

		if (cuIsZ((cuReal)iSHA) || M.linear_size()) {

			//no iSHE contribution. Note, iSHE is not included in magnetic meshes.
			value = -(V.grad_diri(cell1_idx) * elC.grad_sided(cell1_idx)) / elC[cell1_idx];
		}
		else {

			//iSHE enabled, must use non-homogeneous Neumann boundary condition for grad V -> Note homogeneous Neumann boundary conditions apply when calculating S differentials here (due to Jc.n = 0 at boundaries)
			value = -(V.grad_diri_nneu(cell1_idx, *this) * elC.grad_sided(cell1_idx)) / elC[cell1_idx];
		}

		if (M.linear_size() && cuIsNZ((cuReal)betaD)) {

			cuReal Ms = *pcuMesh->pMs;
			cuReal De = *pcuMesh->pDe;
			pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pMs, Ms, *pcuMesh->pDe, De);

			int idx_M = M.position_to_cellidx(S.cellidx_to_position(cell1_idx));

			cuReal33 grad_S = S.grad_neu(cell1_idx);
			cuReal33 grad_M = M.grad_neu(idx_M);
			cuReal3 delsq_S = S.delsq_neu(cell1_idx);

			cuReal div_grad_S_M = (grad_S.i * grad_M.i) + (grad_S.j * grad_M.j) + (grad_S.k * grad_M.k) + (M[idx_M] * delsq_S);

			value += div_grad_S_M * betaD * De / ((cuReal)MUB_E * Ms * elC[cell1_idx]);
		}

		return value;
	}

	__device__ cuReal diff2_sec(cuReal3 relpos_m1, cuReal3 stencil)
	{
		cuVEC_VC<cuReal>& V = *pcuMesh->pV;
		cuVEC_VC<cuReal>& elC = *pcuMesh->pelC;
		cuVEC_VC<cuReal3>& S = *pcuMesh->pS;
		cuVEC_VC<cuReal3>& M = *pcuMesh->pM;

		cuReal iSHA = *pcuMesh->piSHA;
		cuReal betaD = *pcuMesh->pbetaD;
		pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->piSHA, iSHA, *pcuMesh->pbetaD, betaD);

		int cellm1_idx = V.position_to_cellidx(relpos_m1);

		cuReal value = 0.0;

		if (cuIsZ((cuReal)iSHA) || M.linear_size()) {

			//no iSHE contribution. Note, iSHE is not included in magnetic meshes.
			value = -(V.grad_diri(cellm1_idx) * elC.grad_sided(cellm1_idx)) / elC[cellm1_idx];
		}
		else {

			//iSHE enabled, must use non-homogeneous Neumann boundary condition for grad V -> Note homogeneous Neumann boundary conditions apply when calculating S differentials here (due to Jc.n = 0 at boundaries)
			value = -(V.grad_diri_nneu(cellm1_idx, *this) * elC.grad_sided(cellm1_idx)) / elC[cellm1_idx];
		}

		if (M.linear_size() && cuIsNZ((cuReal)betaD)) {

			cuReal Ms = *pcuMesh->pMs;
			cuReal De = *pcuMesh->pDe;
			pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pMs, Ms, *pcuMesh->pDe, De);

			int idx_M = M.position_to_cellidx(relpos_m1);

			cuReal33 grad_S = S.grad_neu(cellm1_idx);
			cuReal33 grad_M = M.grad_neu(idx_M);
			cuReal3 delsq_S = S.delsq_neu(cellm1_idx);

			cuReal div_grad_S_M = (grad_S.i * grad_M.i) + (grad_S.j * grad_M.j) + (grad_S.k * grad_M.k) + (M[idx_M] * delsq_S);

			value += div_grad_S_M * betaD * De / ((cuReal)MUB_E * Ms * elC[cellm1_idx]);
		}

		return value;
	}
};

#endif

#endif
