#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "ManagedMeshCUDA.h"

#include "MeshParamsControlCUDA.h"

#include "TransportCUDA_Poisson_Spin_V.h"

#include "Transport_Defs.h"

class MeshCUDA;
class ManagedDiffEqFMCUDA;
class DifferentialEquationFMCUDA;
class TransportCUDA;

//This is held as a cu_obj managed class in TransportCUDA modules
//It provides methods and access to mesh data for use in cuVEC_VC methods.
//The methods have fixed names, e.g. Poisson_RHS is used by Poisson solvers to evaluate the r.h.s. of the Poisson equation
//The a_func, b_func and diff2_func methods are used to set CMBND conditions based on the continuity of a quantity and a flux.
//If V is the potential, then the flux is the function f(V) = a_func + b_func * V', where the V' differential direction is perpendicular to the interface.
//This particular class is used for spin transport within the spin current solver.
class TransportCUDA_Spin_S_Funcs {

public:

	//spin transport solver type (see Transport_Defs.h) : copy of stsolve in TransportCUDA, but on the gpu so we can use it in device code
	int stsolve;

	//managed mesh for access to all required mesh VECs and material parameters
	ManagedMeshCUDA* pcuMesh;

	//need this to obtain dMdt in this mesh when calculating spin pumping - not used here but this object is passed to STransportCUDA_GInterf_S_NF_Funcs to use there.
	//this will be a nullptr if not a magnetic mesh
	ManagedDiffEqFMCUDA* pcuDiffEq;

	//need this (held in TransportCUDA) to access delsq V (Poisson_RHS) calculation
	TransportCUDA_Spin_V_Funcs* pPoisson_Spin_V;

	//dM_dt VEC when we need to do vector calculus operations on it
	//points to cuVEC in TransportCUDA
	cuVEC_VC<cuReal3>* pdM_dt;

	//for Poisson equations for S some values are fixed during relaxation, so pre-calculate them and store here to re-use.
	////points to cuVEC in TransportCUDA
	cuVEC<cuReal3>* pdelsq_S_fixed;

public:

	__host__ void construct_cu_obj(void) {}
	__host__ void destruct_cu_obj(void) {}

	BError set_pointers(MeshCUDA* pMeshCUDA, DifferentialEquationFMCUDA* pdiffEqCUDA, TransportCUDA* pTransportCUDA);

	__host__ void set_stsolve(int stsolve_) { set_gpu_value(stsolve, stsolve_); }

	//this evaluates the Poisson RHS when solving the Poisson equation on S
	__device__ cuReal3 Poisson_RHS(int idx)
	{
		cuReal3 delsq_S_RHS;

		cuVEC_VC<cuBReal>& V = *pcuMesh->pV;
		cuVEC_VC<cuReal3>& S = *pcuMesh->pS;
		cuVEC_VC<cuReal3>& M = *pcuMesh->pM;

		cuBReal l_sf = *pcuMesh->pl_sf;
		pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->pl_sf, l_sf);

		//Contributions which apply equally in ferromagnetic and non-ferromagnetic meshes
		//longitudinal S decay term
		delsq_S_RHS = (S[idx] / (l_sf * l_sf));

		//Terms occuring only in magnetic meshes
		if (stsolve == STSOLVE_FERROMAGNETIC) {

			cuBReal Ms = *pcuMesh->pMs;
			cuBReal l_ex = *pcuMesh->pl_ex;
			cuBReal l_ph = *pcuMesh->pl_ph;
			pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->pMs, Ms, *pcuMesh->pl_ex, l_ex, *pcuMesh->pl_ph, l_ph);

			//find grad M and M at the M cell in which the current S cell center is
			int idx_M = M.position_to_cellidx(V.cellidx_to_position(idx));
			cuReal3 m = M[idx_M] / Ms;

			//transverse S decay terms
			delsq_S_RHS += ((S[idx] ^ m) / (l_ex * l_ex) + (m ^ (S[idx] ^ m)) / (l_ph * l_ph));
		}

		//additional fixed contributions if needed
		if (pdelsq_S_fixed->linear_size()) delsq_S_RHS += (*pdelsq_S_fixed)[idx];

		return delsq_S_RHS;
	}

	//boundary differential of S for non-homogeneous Neumann boundary conditions
	__device__ cuVAL3<cuReal3> bdiff(int idx)
	{
		if (stsolve == STSOLVE_FERROMAGNETIC) return cuReal33();

		cuVEC_VC<cuReal3>& E = *pcuMesh->pE;
		cuVEC_VC<cuBReal>& elC = *pcuMesh->pelC;

		cuBReal De = *pcuMesh->pDe;
		cuBReal SHA = *pcuMesh->pSHA;
		pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->pSHA, SHA, *pcuMesh->pDe, De);

		return cu_epsilon3(E[idx]) * (SHA * elC[idx] * (cuBReal)MUB_E / De);
	}

	//Functions used for calculating CMBND values

	//CMBND for S
	//flux = a + b S' at the interface, b = -De, a = -(muB*P*sigma/e) * E ox m + (SHA*sigma*muB/e)*epsE + charge pumping + topological Hall effect	
	__device__ cuReal3 a_func_pri(int cell1_idx, int cell2_idx, cuReal3 shift)
	{
		//need to find value at boundary so use interpolation

		cuVEC_VC<cuReal3>& E = *pcuMesh->pE;
		cuVEC_VC<cuBReal>& elC = *pcuMesh->pelC;
		cuVEC_VC<cuReal3>& S = *pcuMesh->pS;
		cuVEC_VC<cuReal3>& M = *pcuMesh->pM;

		//unit vector perpendicular to interface (pointing from secondary to primary mesh)
		cuReal3 u = shift.normalized() * -1;

		//values on secondary side
		if (stsolve == STSOLVE_FERROMAGNETIC) {

			cuBReal Ms = *pcuMesh->pMs;
			cuBReal De = *pcuMesh->pDe;
			cuBReal P = *pcuMesh->pP;
			pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pMs, Ms, *pcuMesh->pDe, De, *pcuMesh->pP, P);

			int idx_M1 = M.position_to_cellidx(S.cellidx_to_position(cell1_idx));
			int idx_M2 = M.position_to_cellidx(S.cellidx_to_position(cell2_idx));

			//a1 value
			cuReal3 m1 = M[idx_M1] / Ms;
			cuReal3 E1 = E[cell1_idx];
			cuBReal sigma_1 = elC[cell1_idx];

			cuReal33 grad_S1 = S.grad_neu(cell1_idx);

			cuReal3 a1 = -(cuBReal)MUB_E * ((E1 | m1) | u) * (P * sigma_1);

			//a2 value
			cuReal3 m2 = M[idx_M2] / Ms;
			cuReal3 E2 = E[cell2_idx];
			cuBReal sigma_2 = elC[cell2_idx];

			cuReal33 grad_S2 = S.grad_neu(cell2_idx);

			cuReal3 a2 = -(cuBReal)MUB_E * ((E2 | m2) | u) * (P * sigma_2);

			bool cpump_enabled = cuIsNZ(pcuMesh->pcpump_eff->get0());
			bool the_enabled = cuIsNZ(pcuMesh->pthe_eff->get0());

			if (cuIsZ(shift.z) && (cpump_enabled || the_enabled)) {

				cuReal33 grad_m1 = M.grad_neu(idx_M1) / Ms;
				cuReal33 grad_m2 = M.grad_neu(idx_M2) / Ms;

				//topological Hall effect contribution
				if (the_enabled) {

					cuBReal n_density = *pcuMesh->pn_density;
					pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pn_density, n_density);

					cuReal3 B1 = (grad_m1.x ^ grad_m1.y);
					a1 += pcuMesh->pthe_eff->get0() * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_1 * sigma_1 / ((cuBReal)ECHARGE * n_density)) * ((-E1.y * B1 * u.x) + (E1.x * B1 * u.y));

					cuReal3 B2 = (grad_m2.x ^ grad_m2.y);
					a2 += pcuMesh->pthe_eff->get0() * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_2 * sigma_2 / ((cuBReal)ECHARGE * n_density)) * ((-E2.y * B2 * u.x) + (E2.x * B2 * u.y));
				}

				//charge pumping contribution
				if (cpump_enabled) {

					//value a1
					cuReal3 dm_dt_1 = (*pdM_dt)[idx_M1] / Ms;
					a1 += pcuMesh->pcpump_eff->get0() * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_1 / 2) * (((dm_dt_1 ^ grad_m1.x) * u.x) + ((dm_dt_1 ^ grad_m1.y) * u.y));

					cuReal3 dm_dt_2 = (*pdM_dt)[idx_M2] / Ms;
					a2 += pcuMesh->pcpump_eff->get0() * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_2 / 2) * (((dm_dt_2 ^ grad_m2.x) * u.x) + ((dm_dt_2 ^ grad_m2.y) * u.y));
				}
			}

			return (1.5 * a1 - 0.5 * a2);
		}
		else {

			cuBReal SHA = *pcuMesh->pSHA;
			pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pSHA, SHA);

			//non-magnetic mesh
			cuReal3 a1 = (cu_epsilon3(E[cell1_idx]) | u) * SHA * elC[cell1_idx] * (cuBReal)MUB_E;
			cuReal3 a2 = (cu_epsilon3(E[cell2_idx]) | u) * SHA * elC[cell2_idx] * (cuBReal)MUB_E;

			return (1.5 * a1 - 0.5 * a2);
		}
	}

	//flux = a + b S' at the interface, b = -De, a = -(muB*P*sigma/e) * E ox m + (SHA*sigma*muB/e)*epsE + charge pumping + topological Hall effect
	__device__ cuReal3 a_func_sec(cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil)
	{
		//need to find value at boundary so use interpolation

		cuVEC_VC<cuReal3>& E = *pcuMesh->pE;
		cuVEC_VC<cuBReal>& elC = *pcuMesh->pelC;
		cuVEC_VC<cuReal3>& S = *pcuMesh->pS;
		cuVEC_VC<cuReal3>& M = *pcuMesh->pM;

		//unit vector perpendicular to interface (pointing from secondary to primary mesh)
		cuReal3 u = shift.normalized() * -1;

		//values on secondary side
		if (stsolve == STSOLVE_FERROMAGNETIC) {

			cuBReal Ms = *pcuMesh->pMs;
			cuBReal De = *pcuMesh->pDe;
			cuBReal P = *pcuMesh->pP;
			pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pMs, Ms, *pcuMesh->pDe, De, *pcuMesh->pP, P);

			//a1 value
			cuReal3 m1 = M.weighted_average(relpos_m1, stencil) / Ms;
			cuReal3 E1 = E.weighted_average(relpos_m1, stencil);
			cuBReal sigma_1 = elC.weighted_average(relpos_m1, stencil);

			int idx_S1 = S.position_to_cellidx(relpos_m1);
			cuReal33 grad_S1 = S.grad_neu(idx_S1);

			cuReal3 a1 = -(cuBReal)MUB_E * ((E1 | m1) | u) * (P * sigma_1);

			//a2 value
			cuReal3 m2 = M.weighted_average(relpos_m1 + shift, stencil) / Ms;
			cuReal3 E2 = E.weighted_average(relpos_m1 + shift, stencil);
			cuBReal sigma_2 = elC.weighted_average(relpos_m1 + shift, stencil);

			int idx_S2 = S.position_to_cellidx(relpos_m1 + shift);
			cuReal33 grad_S2 = S.grad_neu(idx_S2);

			cuReal3 a2 = -(cuBReal)MUB_E * ((E2 | m2) | u) * (P * sigma_2);

			bool cpump_enabled = cuIsNZ(pcuMesh->pcpump_eff->get0());
			bool the_enabled = cuIsNZ(pcuMesh->pthe_eff->get0());

			if (cuIsZ(shift.z) && (cpump_enabled || the_enabled)) {

				int idx_M1 = M.position_to_cellidx(relpos_m1);
				int idx_M2 = M.position_to_cellidx(relpos_m1 + shift);

				cuReal33 grad_m1 = M.grad_neu(idx_M1) / Ms;
				cuReal33 grad_m2 = M.grad_neu(idx_M2) / Ms;

				//topological Hall effect contribution
				if (the_enabled) {

					cuBReal n_density = *pcuMesh->pn_density;
					pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pn_density, n_density);

					cuReal3 B1 = (grad_m1.x ^ grad_m1.y);
					a1 += pcuMesh->pthe_eff->get0() * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_1 * sigma_1 / ((cuBReal)ECHARGE * n_density)) * ((-E1.y * B1 * u.x) + (E1.x * B1 * u.y));

					cuReal3 B2 = (grad_m2.x ^ grad_m2.y);
					a2 += pcuMesh->pthe_eff->get0() * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_2 * sigma_2 / ((cuBReal)ECHARGE * n_density)) * ((-E2.y * B2 * u.x) + (E2.x * B2 * u.y));
				}

				//charge pumping contribution
				if (cpump_enabled) {

					//value a1
					cuReal3 dm_dt_1 = pdM_dt->weighted_average(relpos_m1, stencil) / Ms;
					a1 += pcuMesh->pcpump_eff->get0() * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_1 / 2) * (((dm_dt_1 ^ grad_m1.x) * u.x) + ((dm_dt_1 ^ grad_m1.y) * u.y));

					cuReal3 dm_dt_2 = pdM_dt->weighted_average(relpos_m1 + shift, stencil) / Ms;
					a2 += pcuMesh->pcpump_eff->get0() * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_2 / 2) * (((dm_dt_2 ^ grad_m2.x) * u.x) + ((dm_dt_2 ^ grad_m2.y) * u.y));
				}
			}

			return (1.5 * a1 - 0.5 * a2);
		}
		else {

			cuBReal SHA = *pcuMesh->pSHA;
			pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pSHA, SHA);

			//non-magnetic mesh
			cuReal3 a1 = (cu_epsilon3(E.weighted_average(relpos_m1, stencil)) | u) * SHA * elC.weighted_average(relpos_m1, stencil) * (cuBReal)MUB_E;
			cuReal3 a2 = (cu_epsilon3(E.weighted_average(relpos_m1 + shift, stencil)) | u) * SHA * elC.weighted_average(relpos_m1 + shift, stencil) * (cuBReal)MUB_E;

			return (1.5 * a1 - 0.5 * a2);
		}
	}

	//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface	
	__device__ cuBReal b_func_pri(int cell1_idx, int cell2_idx)
	{
		cuBReal De = *pcuMesh->pDe;
		pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pDe, De);

		return -1.0 * De;
	}

	__device__ cuBReal b_func_sec(cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil)
	{
		cuBReal De = *pcuMesh->pDe;
		pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pDe, De);

		return -1.0 * De;
	}

	//second order differential of S along the shift axis
	//this is simply Evaluate_SpinSolver_delsqS_RHS from which we subtract second order differentials orthogonal to the shift axis
	__device__ cuReal3 diff2_pri(int cell1_idx, cuReal3 shift)
	{
		return Poisson_RHS(cell1_idx);
	}

	//second order differential of S along the shift axis
	//this is simply Evaluate_SpinSolver_delsqS_RHS from which we subtract second order differentials orthogonal to the shift axis
	__device__ cuReal3 diff2_sec(cuReal3 relpos_m1, cuReal3 stencil, cuReal3 shift)
	{
		cuVEC_VC<cuReal3>& S = *pcuMesh->pS;

		int cellm1_idx = S.position_to_cellidx(relpos_m1);

		return Poisson_RHS(cellm1_idx);
	}

	//multiply spin accumulation by these to obtain spin potential, i.e. Vs = (De / elC) * (e/muB) * S, evaluated at the boundary
	__device__ cuBReal c_func_pri(int cell_idx)
	{
		cuVEC_VC<cuBReal>& elC = *pcuMesh->pelC;
		
		cuBReal De = *pcuMesh->pDe;
		pcuMesh->update_parameters_ecoarse(cell_idx, *pcuMesh->pDe, De);

		return De / (elC[cell_idx] * (cuBReal)MUB_E);
	}

	__device__ cuBReal c_func_sec(cuReal3 relpos, cuReal3 stencil)
	{
		cuVEC_VC<cuBReal>& elC = *pcuMesh->pelC;
		
		cuBReal De = *pcuMesh->pDe;
		pcuMesh->update_parameters_atposition(relpos, *pcuMesh->pDe, De);

		return De / (elC.weighted_average(relpos, stencil) * (cuBReal)MUB_E);
	}
};

#endif

#endif

