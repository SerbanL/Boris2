#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_TRANSPORT

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "ManagedMeshCUDA.h"

#include "MeshParamsControlCUDA.h"

class MeshCUDA;
class ManagedDiffEqCUDA;
class DifferentialEquationCUDA;

//This is held as a cu_obj managed class in TransportCUDA modules
//It provides methods and access to mesh data for use in cuVEC_VC methods.
//The methods have fixed names, e.g. Poisson_RHS is used by Poisson solvers to evaluate the r.h.s. of the Poisson equation
//The a_func, b_func and diff2_func methods are used to set CMBND conditions based on the continuity of a quantity and a flux.
//If V is the potential, then the flux is the function f(V) = a_func + b_func * V', where the V' differential direction is perpendicular to the interface.
//This particular class is used for spin transport within the spin current solver.
class TransportCUDA_Spin_S_Funcs {

public:

	//managed mesh for access to all required mesh VECs and material parameters
	ManagedMeshCUDA* pcuMesh;

	//need this to obtain dMdt in this mesh when calculating spin pumping - not used here but this object is passed to STransportCUDA_GInterf_S_NF_Funcs to use there.
	//this will be a nullptr if not a magnetic mesh
	ManagedDiffEqCUDA* pcuDiffEq;

public:

	__host__ void construct_cu_obj(void) {}
	__host__ void destruct_cu_obj(void) {}

	BError set_pointers(MeshCUDA* pMeshCUDA, DifferentialEquationCUDA* pdiffEqCUDA);

	//this evaluates the Poisson RHS when solving the Poisson equation on S
	__device__ cuReal3 Poisson_RHS(int idx)
	{
		cuReal3 delsq_S_RHS = cuReal3();

		cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
		cuVEC_VC<cuReal3>& S = *pcuMesh->pS;
		cuVEC_VC<cuReal3>& Jc = *pcuMesh->pJc;

		cuReal De = *pcuMesh->pDe;
		cuReal SHA = *pcuMesh->pSHA;
		cuReal l_sf = *pcuMesh->pl_sf;
		pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->pDe, De, *pcuMesh->pSHA, SHA, *pcuMesh->pl_sf, l_sf);

		//Terms occuring only in ferromagnetic meshes (where SHA = 0)
		if (M.linear_size()) {
			
			cuReal Ms = *pcuMesh->pMs;
			cuReal P = *pcuMesh->pP;
			cuReal betaD = *pcuMesh->pbetaD;
			cuReal l_ex = *pcuMesh->pl_ex;
			cuReal l_ph = *pcuMesh->pl_ph;
			pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->pMs, Ms, *pcuMesh->pP, P, *pcuMesh->pbetaD, betaD, *pcuMesh->pl_ex, l_ex, *pcuMesh->pl_ph, l_ph);

			//1. term due to non-uniformity of M

			//find grad M and M at the M cell in which the current S cell center is
			int idx_M = M.position_to_cellidx(S.cellidx_to_position(idx));

			cuReal3 Mval = M[idx_M];
			cuReal33 grad_M = M.grad_neu(idx_M);
			cuReal3 Jc_dot_del_M = grad_M | Jc[idx];

			delsq_S_RHS -= P * (cuReal)MUB_E * Jc_dot_del_M / (Ms * De);

			//2. transverse S decay terms

			delsq_S_RHS += ((S[idx] ^ Mval) / (Ms * l_ex * l_ex) + (Mval ^ (S[idx] ^ Mval)) / (Ms * Ms * l_ph * l_ph));

			//3. CPP-GMR term
			if (cuIsNZ((cuReal)betaD)) {

				cuReal33 grad_S = S.grad_neu(idx);
				cuReal3 delsq_S = S.delsq_neu(idx);

				cuReal3 grad_S_M_dot_del_M = grad_M | (grad_S * Mval);

				delsq_S_RHS += (grad_S_M_dot_del_M + Mval * ((grad_M.i * grad_S.i) + (grad_M.j * grad_S.j) + (grad_M.k * grad_S.k) + (Mval * delsq_S))) * betaD * P / (Ms * Ms);
			}
		}
		//terms occuring only in non-ferromagnetic meshes (i.e. those with SHA not zero)
		else {

			//1. SHA term. 

			//Check boundary conditions for this term : should be Dirichlet with 0 Jc value normal to the boundary except for electrodes.
			delsq_S_RHS += (SHA * (cuReal)MUB_E / De) * Jc.diveps3_sided(idx);
		}

		//Contributions which apply equally in ferromagnetic and non-ferromagnetic meshes

		//1. longitudinal S decay term

		delsq_S_RHS += (S[idx] / (l_sf * l_sf));

		return delsq_S_RHS;
	}

	//boundary differential of S for non-homogeneous Neumann boundary conditions
	__device__ cuVAL3<cuReal3> bdiff(int idx)
	{
		cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
		cuVEC_VC<cuReal3>& Jc = *pcuMesh->pJc;

		if (M.linear_size()) return cuReal33();

		cuReal De = *pcuMesh->pDe;
		cuReal SHA = *pcuMesh->pSHA;
		pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->pSHA, SHA, *pcuMesh->pDe, De);

		return cu_epsilon3(Jc[idx]) * (SHA * (cuReal)MUB_E / De);
	}

	//Functions used for calculating CMBND values

	//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface	
	__device__ cuReal3 a_func_pri(int cell1_idx, int cell2_idx, cuReal3 shift)
	{
		cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
		cuVEC_VC<cuReal3>& Jc = *pcuMesh->pJc;
		cuVEC_VC<cuReal3>& S = *pcuMesh->pS;

		cuReal3 u = shift.normalized() * -1;

		//values on primary side
		if (M.linear_size()) {

			cuReal Ms = *pcuMesh->pMs;
			cuReal P = *pcuMesh->pP;
			cuReal betaD = *pcuMesh->pbetaD;
			cuReal De = *pcuMesh->pDe;
			pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pMs, Ms, *pcuMesh->pDe, De, *pcuMesh->pP, P, *pcuMesh->pbetaD, betaD);

			//magnetic mesh
			int idx_M1 = M.position_to_cellidx(Jc.cellidx_to_position(cell1_idx));
			cuReal3 M1 = M[idx_M1];
			cuReal3 Jc1 = Jc[cell1_idx];

			cuReal3 a1 = ((Jc1 | M1) | u) * (P / Ms) * (-(cuReal)MUB_E);

			cuReal33 grad_S1 = S.grad_neu(cell1_idx);

			a1 += (((grad_S1 * M1) | M1) * betaD * P * De / (Ms * Ms)) | u;

			//magnetic mesh
			int idx_M2 = M.position_to_cellidx(Jc.cellidx_to_position(cell2_idx));
			cuReal3 M2 = M[idx_M2];
			cuReal3 Jc2 = Jc[cell2_idx];

			cuReal3 a2 = ((Jc2 | M2) | u) * (P / Ms) * (-(cuReal)MUB_E);

			cuReal33 grad_S2 = S.grad_neu(cell2_idx);

			a2 += (((grad_S2 * M2) | M2) * betaD * P * De / (Ms * Ms)) | u;

			return (1.5 * a1 - 0.5 * a2);
		}
		else {

			cuReal SHA = *pcuMesh->pSHA;
			pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pSHA, SHA);

			//non-magnetic mesh
			cuReal3 a1 = (cu_epsilon3(Jc[cell1_idx]) | u) * SHA * (cuReal)MUB_E;
			cuReal3 a2 = (cu_epsilon3(Jc[cell2_idx]) | u) * SHA * (cuReal)MUB_E;

			return (1.5 * a1 - 0.5 * a2);
		}
	}

	__device__ cuReal3 a_func_sec(cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil)
	{
		cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
		cuVEC_VC<cuReal3>& Jc = *pcuMesh->pJc;
		cuVEC_VC<cuReal3>& S = *pcuMesh->pS;

		//unit vector perpendicular to interface (pointing from secondary to primary mesh)
		cuReal3 u = shift.normalized() * -1;

		//values on secondary side
		if (M.linear_size()) {

			cuReal Ms = *pcuMesh->pMs;
			cuReal P = *pcuMesh->pP;
			cuReal betaD = *pcuMesh->pbetaD;
			cuReal De = *pcuMesh->pDe;
			pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pMs, Ms, *pcuMesh->pDe, De, *pcuMesh->pP, P, *pcuMesh->pbetaD, betaD);

			//a1 value
			cuReal3 M1 = M.weighted_average(relpos_m1, stencil);
			cuReal3 Jc1 = Jc.weighted_average(relpos_m1, stencil);
			int idx_S1 = S.position_to_cellidx(relpos_m1);
			cuReal33 grad_S1 = S.grad_neu(idx_S1);

			cuReal3 a1 =
				-(cuReal)MUB_E * ((Jc1 | M1) | u) * (P / Ms)
				+ ((((grad_S1 * M1) | M1) * betaD * P * De / (Ms * Ms)) | u);

			//a2 value
			cuReal3 M2 = M.weighted_average(relpos_m1 + shift, stencil);
			cuReal3 Jc2 = Jc.weighted_average(relpos_m1 + shift, stencil);
			int idx_S2 = S.position_to_cellidx(relpos_m1 + shift);
			cuReal33 grad_S2 = S.grad_neu(idx_S2);

			cuReal3 a2 =
				-(cuReal)MUB_E * ((Jc2 | M2) | u) * (P / Ms)
				+ ((((grad_S2 * M2) | M2) * betaD * P * De / (Ms * Ms)) | u);

			return (1.5 * a1 - 0.5 * a2);
		}
		else {

			cuReal SHA = *pcuMesh->pSHA;
			pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pSHA, SHA);

			//non-magnetic mesh
			cuReal3 a1 = (cu_epsilon3(Jc.weighted_average(relpos_m1, stencil)) | u) * SHA * (cuReal)MUB_E;
			cuReal3 a2 = (cu_epsilon3(Jc.weighted_average(relpos_m1 + shift, stencil)) | u) * SHA * (cuReal)MUB_E;

			return (1.5 * a1 - 0.5 * a2);
		}
	}

	//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface	
	__device__ cuReal b_func_pri(int cell1_idx, int cell2_idx)
	{
		cuReal De = *pcuMesh->pDe;
		pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pDe, De);

		return -1.0 * De;
	}

	__device__ cuReal b_func_sec(cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil)
	{
		cuReal De = *pcuMesh->pDe;
		pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pDe, De);

		return -1.0 * De;
	}

	//second order differential of V at cells either side of the boundary; delsq V = -grad V * grad elC / elC
	__device__ cuReal3 diff2_pri(int cell1_idx)
	{
		cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
		cuVEC_VC<cuReal3>& S = *pcuMesh->pS;
		cuVEC_VC<cuReal3>& Jc = *pcuMesh->pJc;

		if (M.linear_size()) {

			cuReal Ms = *pcuMesh->pMs;
			cuReal De = *pcuMesh->pDe;
			cuReal P = *pcuMesh->pP;
			cuReal betaD = *pcuMesh->pbetaD;
			cuReal l_sf = *pcuMesh->pl_sf;
			cuReal l_ex = *pcuMesh->pl_ex;
			cuReal l_ph = *pcuMesh->pl_ph;
			pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pMs, Ms, *pcuMesh->pDe, De, *pcuMesh->pP, P, *pcuMesh->pbetaD, betaD, *pcuMesh->pl_ex, l_ex, *pcuMesh->pl_ph, l_ph, *pcuMesh->pl_sf, l_sf);

			cuReal3 value = S[cell1_idx] / (l_sf * l_sf);

			//add to value for magnetic mesh

			int idx_M = M.position_to_cellidx(S.cellidx_to_position(cell1_idx));

			cuReal3 Sval = S[cell1_idx];
			cuReal3 Mval = M[idx_M];

			value += (Sval ^ Mval) / (Ms * l_ex * l_ex) + (Mval ^ (Sval ^ Mval)) / (Ms * Ms * l_ph * l_ph);

			//find grad M and M at the M cell in which the current S cell center is
			cuReal33 grad_M = M.grad_neu(idx_M);
			cuReal3 Jc_dot_del_M = grad_M | Jc[cell1_idx];

			value -= P * (cuReal)MUB_E * Jc_dot_del_M / (Ms * De);

			if (cuIsNZ((cuReal)betaD)) {

				cuReal33 grad_S = S.grad_neu(cell1_idx);
				cuReal3 delsq_S = S.delsq_neu(cell1_idx);

				cuReal3 grad_S_M_dot_del_M = grad_M | (grad_S * Mval);

				value += (grad_S_M_dot_del_M + Mval * ((grad_M.i * grad_S.i) + (grad_M.j * grad_S.j) + (grad_M.k * grad_S.k) + (Mval * delsq_S))) * betaD * P / (Ms * Ms);
			}

			return value;
		}
		else {

			cuReal SHA = *pcuMesh->pSHA;
			cuReal De = *pcuMesh->pDe;
			cuReal l_sf = *pcuMesh->pl_sf;
			pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pSHA, SHA, *pcuMesh->pDe, De, *pcuMesh->pl_sf, l_sf);

			cuReal3 value = S[cell1_idx] / (l_sf * l_sf);

			//add to value for non-magnetic mesh
			value += (SHA * (cuReal)MUB_E / De) * Jc.diveps3_sided(cell1_idx);

			return value;
		}
	}

	__device__ cuReal3 diff2_sec(cuReal3 relpos_m1, cuReal3 stencil)
	{
		cuVEC_VC<cuReal3>& M = *pcuMesh->pM;
		cuVEC_VC<cuReal3>& S = *pcuMesh->pS;
		cuVEC_VC<cuReal3>& Jc = *pcuMesh->pJc;

		if (M.linear_size()) {

			cuReal Ms = *pcuMesh->pMs;
			cuReal De = *pcuMesh->pDe;
			cuReal P = *pcuMesh->pP;
			cuReal betaD = *pcuMesh->pbetaD;
			cuReal l_sf = *pcuMesh->pl_sf;
			cuReal l_ex = *pcuMesh->pl_ex;
			cuReal l_ph = *pcuMesh->pl_ph;
			pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pMs, Ms, *pcuMesh->pDe, De, *pcuMesh->pP, P, *pcuMesh->pbetaD, betaD, *pcuMesh->pl_ex, l_ex, *pcuMesh->pl_ph, l_ph, *pcuMesh->pl_sf, l_sf);

			cuReal3 value = S.weighted_average(relpos_m1, stencil) / (l_sf * l_sf);

			//add to value for magnetic mesh

			cuReal3 Sval = S.weighted_average(relpos_m1, stencil);
			cuReal3 Mval = M.weighted_average(relpos_m1, stencil);

			value += (Sval ^ Mval) / (Ms * l_ex * l_ex) + (Mval ^ (Sval ^ Mval)) / (Ms * Ms * l_ph * l_ph);

			//find grad M and M at the M cell in which the current S cell center is
			int idx_M = M.position_to_cellidx(relpos_m1);
			cuReal33 grad_M = M.grad_neu(idx_M);
			cuReal3 Jc_dot_del_M = grad_M | Jc.weighted_average(relpos_m1, stencil);

			value -= P * (cuReal)MUB_E * Jc_dot_del_M / (Ms * De);

			if (cuIsNZ((cuReal)betaD)) {

				int idx_S = S.position_to_cellidx(relpos_m1);
				cuReal33 grad_S = S.grad_neu(idx_S);
				cuReal3 delsq_S = S.delsq_neu(idx_S);

				cuReal3 grad_S_M_dot_del_M = grad_M | (grad_S * Mval);

				value += (grad_S_M_dot_del_M + Mval * ((grad_M.i * grad_S.i) + (grad_M.j * grad_S.j) + (grad_M.k * grad_S.k) + (Mval * delsq_S))) * betaD * P / (Ms * Ms);
			}

			return value;
		}
		else {

			cuReal SHA = *pcuMesh->pSHA;
			cuReal De = *pcuMesh->pDe;
			cuReal l_sf = *pcuMesh->pl_sf;
			pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pSHA, SHA, *pcuMesh->pDe, De, *pcuMesh->pl_sf, l_sf);

			cuReal3 value = S.weighted_average(relpos_m1, stencil) / (l_sf * l_sf);

			//add to value for non-magnetic mesh
			int idx_Jc = Jc.position_to_cellidx(relpos_m1);

			value += (SHA * (cuReal)MUB_E / De) * Jc.diveps3_sided(idx_Jc);

			return value;
		}
	}

	//multiply spin accumulation by these to obtain spin potential, i.e. Vs = (De / elC) * (e/muB) * S, evaluated at the boundary
	__device__ cuReal c_func_pri(int cell_idx)
	{
		cuVEC_VC<cuReal>& elC = *pcuMesh->pelC;
		
		cuReal De = *pcuMesh->pDe;
		pcuMesh->update_parameters_ecoarse(cell_idx, *pcuMesh->pDe, De);

		return De / (elC[cell_idx] * (cuReal)MUB_E);
	}

	__device__ cuReal c_func_sec(cuReal3 relpos, cuReal3 stencil)
	{
		cuVEC_VC<cuReal>& elC = *pcuMesh->pelC;
		
		cuReal De = *pcuMesh->pDe;
		pcuMesh->update_parameters_atposition(relpos, *pcuMesh->pDe, De);

		return De / (elC.weighted_average(relpos, stencil) * (cuReal)MUB_E);
	}
};

#endif

#endif

