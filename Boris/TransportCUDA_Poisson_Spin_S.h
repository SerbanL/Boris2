#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "ManagedMeshCUDA.h"
#include "MeshParamsControlCUDA.h"

#include "ManagedAtom_MeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

#include "TransportCUDA_Poisson_Spin_V.h"

#include "Transport_Defs.h"

class MeshCUDA;
class Atom_MeshCUDA;
class TransportBaseCUDA;

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

	//managed mesh for access to all required mesh VECs and material parameters
	ManagedAtom_MeshCUDA* pcuaMesh;

	//need this (held in TransportCUDA) to access delsq V (Poisson_RHS) calculation
	TransportCUDA_Spin_V_Funcs* pPoisson_Spin_V;

	//dM_dt VEC (charge pumping and spin pumping)
	//points to cuVEC in TransportCUDA
	cuVEC_VC<cuReal3>* pdM_dt;

	//for Poisson equations for S some values are fixed during relaxation, so pre-calculate them and store here to re-use.
	////points to cuVEC in TransportCUDA
	cuVEC<cuReal3>* pdelsq_S_fixed;

public:

	__host__ void construct_cu_obj(void) {}
	__host__ void destruct_cu_obj(void) {}

	//for modules held in micromagnetic meshes
	BError set_pointers_transport(MeshCUDA* pMeshCUDA, TransportBaseCUDA* pTransportCUDA);
	//for modules held in atomistic meshes
	BError set_pointers_atomtransport(Atom_MeshCUDA* paMeshCUDA, TransportBaseCUDA* paTransportCUDA);

	__host__ void set_stsolve(int stsolve_) { set_gpu_value(stsolve, stsolve_); }

	//this evaluates the Poisson RHS when solving the Poisson equation on S
	__device__ cuReal3 Poisson_RHS(int idx)
	{
		cuReal3 delsq_S_RHS;

		cuVEC_VC<cuBReal>& V = (pcuMesh ? *pcuMesh->pV : *pcuaMesh->pV);
		cuVEC_VC<cuReal3>& S = (pcuMesh ? *pcuMesh->pS : *pcuaMesh->pS);
		cuVEC_VC<cuReal3>& M = (pcuMesh ? *pcuMesh->pM : *pcuaMesh->pM1);

		cuBReal l_sf = (pcuMesh ? *pcuMesh->pl_sf : *pcuaMesh->pl_sf);
		if (pcuMesh) pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->pl_sf, l_sf);
		else pcuaMesh->update_parameters_ecoarse(idx, *pcuaMesh->pl_sf, l_sf);

		if (stsolve == STSOLVE_TUNNELING) {

			cuBReal elecCond = *pcuMesh->pelecCond;
			pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->pelecCond, elecCond);

			//Use elecCond to mark out metallic pinholes.
			//tunelling : l_sf tends to infinity
			if (elecCond > 0.0) delsq_S_RHS = (S[idx] / (l_sf * l_sf));
			return delsq_S_RHS;
		}

		//Contributions which apply equally in ferromagnetic and non-ferromagnetic meshes
		//longitudinal S decay term
		delsq_S_RHS = (S[idx] / (l_sf * l_sf));

		//Terms occuring only in magnetic meshes
		if (stsolve == STSOLVE_FERROMAGNETIC || stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

			cuBReal l_ex = (pcuMesh ? *pcuMesh->pl_ex : *pcuaMesh->pl_ex);
			cuBReal l_ph = (pcuMesh ? *pcuMesh->pl_ph : *pcuaMesh->pl_ph);
			if (pcuMesh)pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->pl_ex, l_ex, *pcuMesh->pl_ph, l_ph);
			else pcuaMesh->update_parameters_ecoarse(idx, *pcuaMesh->pl_ex, l_ex, *pcuaMesh->pl_ph, l_ph);

			int idx_M = M.position_to_cellidx(V.cellidx_to_position(idx));
			cuReal3 m = cu_normalize(M[idx_M]);

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
		if (stsolve == STSOLVE_FERROMAGNETIC || stsolve == STSOLVE_FERROMAGNETIC_ATOM || stsolve == STSOLVE_TUNNELING) return cuReal33();

		//only pcuMesh here, but best to check
		if (pcuMesh) {

			cuVEC_VC<cuReal3>& E = *pcuMesh->pE;
			cuVEC_VC<cuBReal>& elC = *pcuMesh->pelC;

			cuBReal De = *pcuMesh->pDe;
			cuBReal SHA = *pcuMesh->pSHA;
			pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->pSHA, SHA, *pcuMesh->pDe, De);

			return cu_epsilon3(E[idx]) * (SHA * elC[idx] * (cuBReal)MUB_E / De);
		}
		else return cuVAL3<cuReal3>();
	}

	//Functions used for calculating CMBND values

	//CMBND for S
	//flux = a + b S' at the interface, b = -De, a = -(muB*P*sigma/e) * E ox m + (SHA*sigma*muB/e)*epsE + charge pumping + topological Hall effect	
	__device__ cuReal3 a_func_pri(int cell1_idx, int cell2_idx, cuReal3 shift)
	{
		if (stsolve == STSOLVE_TUNNELING) return cuReal3();

		//need to find value at boundary so use interpolation

		cuVEC_VC<cuReal3>& E = (pcuMesh ? *pcuMesh->pE : *pcuaMesh->pE);
		cuVEC_VC<cuBReal>& elC = (pcuMesh ? *pcuMesh->pelC : *pcuaMesh->pelC);
		cuVEC_VC<cuReal3>& S = (pcuMesh ? *pcuMesh->pS : *pcuaMesh->pS);
		cuVEC_VC<cuReal3>& M = (pcuMesh ? *pcuMesh->pM : *pcuaMesh->pM1);

		//unit vector perpendicular to interface (pointing from secondary to primary mesh)
		cuReal3 u = shift.normalized() * -1;

		//values on secondary side
		if (stsolve == STSOLVE_FERROMAGNETIC || stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

			cuBReal De = (pcuMesh ? *pcuMesh->pDe : *pcuaMesh->pDe);
			cuBReal P = (pcuMesh ? *pcuMesh->pP : *pcuaMesh->pP);
			if (pcuMesh) pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pDe, De, *pcuMesh->pP, P);
			else pcuaMesh->update_parameters_ecoarse(cell1_idx, *pcuaMesh->pDe, De, *pcuaMesh->pP, P);

			int idx_M1 = M.position_to_cellidx(S.cellidx_to_position(cell1_idx));
			int idx_M2 = M.position_to_cellidx(S.cellidx_to_position(cell2_idx));

			//a1 value
			cuReal3 m1 = cu_normalize(M[idx_M1]);
			cuReal3 E1 = E[cell1_idx];
			cuBReal sigma_1 = elC[cell1_idx];

			cuReal33 grad_S1 = S.grad_neu(cell1_idx);

			cuReal3 a1 = -(cuBReal)MUB_E * ((E1 | m1) | u) * (P * sigma_1);

			//a2 value
			cuReal3 m2 = cu_normalize(M[idx_M2]);
			cuReal3 E2 = E[cell2_idx];
			cuBReal sigma_2 = elC[cell2_idx];

			cuReal33 grad_S2 = S.grad_neu(cell2_idx);

			cuReal3 a2 = -(cuBReal)MUB_E * ((E2 | m2) | u) * (P * sigma_2);

			cuBReal cpump_eff0 = (pcuMesh ? pcuMesh->pcpump_eff->get0() : pcuaMesh->pcpump_eff->get0());
			bool cpump_enabled = cuIsNZ(cpump_eff0);

			cuBReal the_eff0 = (pcuMesh ? pcuMesh->pthe_eff->get0() : pcuaMesh->pthe_eff->get0());
			bool the_enabled = cuIsNZ(the_eff0);

			if (cuIsZ(shift.z) && (cpump_enabled || the_enabled)) {

				cuReal33 grad_m1 = cu_normalize(M.grad_neu(idx_M1), M[idx_M1]);
				cuReal33 grad_m2 = cu_normalize(M.grad_neu(idx_M2), M[idx_M2]);

				//topological Hall effect contribution
				if (the_enabled) {

					cuBReal n_density = (pcuMesh ? *pcuMesh->pn_density : *pcuaMesh->pn_density);
					if (pcuMesh) pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pn_density, n_density);
					else pcuaMesh->update_parameters_ecoarse(cell1_idx, *pcuaMesh->pn_density, n_density);

					cuReal3 B1 = (grad_m1.x ^ grad_m1.y);
					a1 += the_eff0 * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_1 * sigma_1 / ((cuBReal)ECHARGE * n_density)) * ((-E1.y * B1 * u.x) + (E1.x * B1 * u.y));

					cuReal3 B2 = (grad_m2.x ^ grad_m2.y);
					a2 += the_eff0 * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_2 * sigma_2 / ((cuBReal)ECHARGE * n_density)) * ((-E2.y * B2 * u.x) + (E2.x * B2 * u.y));
				}

				//charge pumping contribution
				if (cpump_enabled) {

					//value a1
					cuReal3 dm_dt_1 = cu_normalize((*pdM_dt)[idx_M1], M[idx_M1]);
					a1 += cpump_eff0 * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_1 / 2) * (((dm_dt_1 ^ grad_m1.x) * u.x) + ((dm_dt_1 ^ grad_m1.y) * u.y));

					cuReal3 dm_dt_2 = cu_normalize((*pdM_dt)[idx_M2], M[idx_M2]);
					a2 += cpump_eff0 * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_2 / 2) * (((dm_dt_2 ^ grad_m2.x) * u.x) + ((dm_dt_2 ^ grad_m2.y) * u.y));
				}
			}

			return (1.5 * a1 - 0.5 * a2);
		}
		else {

			//non-magnetic mesh. pcuMesh only, but best to check
			if (pcuMesh) {

				cuBReal SHA = *pcuMesh->pSHA;
				pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pSHA, SHA);

				//non-magnetic mesh
				cuReal3 a1 = (cu_epsilon3(E[cell1_idx]) | u) * SHA * elC[cell1_idx] * (cuBReal)MUB_E;
				cuReal3 a2 = (cu_epsilon3(E[cell2_idx]) | u) * SHA * elC[cell2_idx] * (cuBReal)MUB_E;

				return (1.5 * a1 - 0.5 * a2);
			}
			else return cuReal3();
		}
	}

	//flux = a + b S' at the interface, b = -De, a = -(muB*P*sigma/e) * E ox m + (SHA*sigma*muB/e)*epsE + charge pumping + topological Hall effect
	__device__ cuReal3 a_func_sec(cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil)
	{
		if (stsolve == STSOLVE_TUNNELING) return cuReal3();

		//need to find value at boundary so use interpolation

		cuVEC_VC<cuReal3>& E = (pcuMesh ? *pcuMesh->pE : *pcuaMesh->pE);
		cuVEC_VC<cuBReal>& elC = (pcuMesh ? *pcuMesh->pelC : *pcuaMesh->pelC);
		cuVEC_VC<cuReal3>& S = (pcuMesh ? *pcuMesh->pS : *pcuaMesh->pS);
		cuVEC_VC<cuReal3>& M = (pcuMesh ? *pcuMesh->pM : *pcuaMesh->pM1);

		//unit vector perpendicular to interface (pointing from secondary to primary mesh)
		cuReal3 u = shift.normalized() * -1;

		//values on secondary side
		if (stsolve == STSOLVE_FERROMAGNETIC || stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

			cuBReal De = (pcuMesh ? *pcuMesh->pDe : *pcuaMesh->pDe);
			cuBReal P = (pcuMesh ? *pcuMesh->pP : *pcuaMesh->pP);
			if (pcuMesh) pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pDe, De, *pcuMesh->pP, P);
			else pcuaMesh->update_parameters_atposition(relpos_m1, *pcuaMesh->pDe, De, *pcuaMesh->pP, P);

			//a1 value
			cuReal3 m1 = cu_normalize(M.weighted_average(relpos_m1, stencil));
			cuReal3 E1 = E.weighted_average(relpos_m1, stencil);
			cuBReal sigma_1 = elC.weighted_average(relpos_m1, stencil);

			int idx_S1 = S.position_to_cellidx(relpos_m1);
			cuReal33 grad_S1 = S.grad_neu(idx_S1);

			cuReal3 a1 = -(cuBReal)MUB_E * ((E1 | m1) | u) * (P * sigma_1);

			//a2 value
			cuReal3 m2 = cu_normalize(M.weighted_average(relpos_m1 + shift, stencil));
			cuReal3 E2 = E.weighted_average(relpos_m1 + shift, stencil);
			cuBReal sigma_2 = elC.weighted_average(relpos_m1 + shift, stencil);

			int idx_S2 = S.position_to_cellidx(relpos_m1 + shift);
			cuReal33 grad_S2 = S.grad_neu(idx_S2);

			cuReal3 a2 = -(cuBReal)MUB_E * ((E2 | m2) | u) * (P * sigma_2);

			cuBReal cpump_eff0 = (pcuMesh ? pcuMesh->pcpump_eff->get0() : pcuaMesh->pcpump_eff->get0());
			bool cpump_enabled = cuIsNZ(cpump_eff0);

			cuBReal the_eff0 = (pcuMesh ? pcuMesh->pthe_eff->get0() : pcuaMesh->pthe_eff->get0());
			bool the_enabled = cuIsNZ(the_eff0);

			if (cuIsZ(shift.z) && (cpump_enabled || the_enabled)) {

				int idx_M1 = M.position_to_cellidx(relpos_m1);
				int idx_M2 = M.position_to_cellidx(relpos_m1 + shift);

				cuReal33 grad_m1 = cu_normalize(M.grad_neu(idx_M1), M[idx_M1]);
				cuReal33 grad_m2 = cu_normalize(M.grad_neu(idx_M2), M[idx_M2]);

				//topological Hall effect contribution
				if (the_enabled) {

					cuBReal n_density = (pcuMesh ? *pcuMesh->pn_density : *pcuaMesh->pn_density);
					if (pcuMesh) pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pn_density, n_density);
					else pcuaMesh->update_parameters_atposition(relpos_m1, *pcuaMesh->pn_density, n_density);

					cuReal3 B1 = (grad_m1.x ^ grad_m1.y);
					a1 += the_eff0 * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_1 * sigma_1 / ((cuBReal)ECHARGE * n_density)) * ((-E1.y * B1 * u.x) + (E1.x * B1 * u.y));

					cuReal3 B2 = (grad_m2.x ^ grad_m2.y);
					a2 += the_eff0 * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_2 * sigma_2 / ((cuBReal)ECHARGE * n_density)) * ((-E2.y * B2 * u.x) + (E2.x * B2 * u.y));
				}

				//charge pumping contribution
				if (cpump_enabled) {

					//value a1
					cuReal3 dm_dt_1 = cu_normalize(pdM_dt->weighted_average(relpos_m1, stencil), M[idx_M1]);
					a1 += cpump_eff0 * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_1 / 2) * (((dm_dt_1 ^ grad_m1.x) * u.x) + ((dm_dt_1 ^ grad_m1.y) * u.y));

					cuReal3 dm_dt_2 = cu_normalize(pdM_dt->weighted_average(relpos_m1 + shift, stencil), M[idx_M2]);
					a2 += cpump_eff0 * ((cuBReal)HBAR_E * (cuBReal)MUB_E * sigma_2 / 2) * (((dm_dt_2 ^ grad_m2.x) * u.x) + ((dm_dt_2 ^ grad_m2.y) * u.y));
				}
			}

			return (1.5 * a1 - 0.5 * a2);
		}
		else {

			//non-magnetic mesh. pcuMesh only, but best to check
			if (pcuMesh) {

				cuBReal SHA = *pcuMesh->pSHA;
				pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pSHA, SHA);

				//non-magnetic mesh
				cuReal3 a1 = (cu_epsilon3(E.weighted_average(relpos_m1, stencil)) | u) * SHA * elC.weighted_average(relpos_m1, stencil) * (cuBReal)MUB_E;
				cuReal3 a2 = (cu_epsilon3(E.weighted_average(relpos_m1 + shift, stencil)) | u) * SHA * elC.weighted_average(relpos_m1 + shift, stencil) * (cuBReal)MUB_E;

				return (1.5 * a1 - 0.5 * a2);
			}
			else return cuReal3();
		}
	}

	//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface	
	__device__ cuBReal b_func_pri(int cell1_idx, int cell2_idx)
	{
		cuBReal De = (pcuMesh ? *pcuMesh->pDe : *pcuaMesh->pDe);
		cuBReal elecCond = (pcuMesh ? *pcuMesh->pelecCond : *pcuaMesh->pelecCond);
		if (pcuMesh) pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pDe, De, *pcuMesh->pelecCond, elecCond);
		else pcuaMesh->update_parameters_ecoarse(cell1_idx, *pcuaMesh->pDe, De, *pcuaMesh->pelecCond, elecCond);

		//metallic conduction : use De
		if (elecCond > 0.0) return -De;
		//tunelling : use 1
		else return -1;
	}

	__device__ cuBReal b_func_sec(cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil)
	{
		cuBReal De = (pcuMesh ? *pcuMesh->pDe : *pcuaMesh->pDe);
		cuBReal elecCond = (pcuMesh ? *pcuMesh->pelecCond : *pcuaMesh->pelecCond);
		if (pcuMesh) pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pDe, De, *pcuMesh->pelecCond, elecCond);
		else pcuaMesh->update_parameters_atposition(relpos_m1, *pcuaMesh->pDe, De, *pcuaMesh->pelecCond, elecCond);

		//metallic conduction : use De
		if (elecCond > 0.0) return -De;
		//tunelling : use 1
		else return -1;
	}

	//second order differential of S along the shift axis
	//this is simply Evaluate_SpinSolver_delsqS_RHS from which we subtract second order differentials orthogonal to the shift axis
	__device__ cuReal3 diff2_pri(int cell1_idx, cuReal3 shift)
	{
		cuReal3 delsq_S_RHS;

		cuVEC_VC<cuBReal>& V = (pcuMesh ? *pcuMesh->pV : *pcuaMesh->pV);
		cuVEC_VC<cuReal3>& S = (pcuMesh ? *pcuMesh->pS : *pcuaMesh->pS);
		cuVEC_VC<cuReal3>& M = (pcuMesh ? *pcuMesh->pM : *pcuaMesh->pM1);
		cuVEC_VC<cuReal3>& E = (pcuMesh ? *pcuMesh->pE : *pcuaMesh->pE);
		cuVEC_VC<cuBReal>& elC = (pcuMesh ? *pcuMesh->pelC : *pcuaMesh->pelC);

		cuBReal l_sf = (pcuMesh ? *pcuMesh->pl_sf : *pcuaMesh->pl_sf);
		if (pcuMesh) pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pl_sf, l_sf);
		else pcuaMesh->update_parameters_ecoarse(cell1_idx, *pcuaMesh->pl_sf, l_sf);

		if (stsolve == STSOLVE_TUNNELING) {

			cuBReal elecCond = *pcuMesh->pelecCond;
			pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pelecCond, elecCond);

			//Use elecCond to mark out metallic pinholes.
			//tunelling : l_sf tends to infinity
			if (elecCond > 0.0) delsq_S_RHS = (S[cell1_idx] / (l_sf * l_sf));
			return delsq_S_RHS;
		}

		//Contributions which apply equally in ferromagnetic and non-ferromagnetic meshes
		//longitudinal S decay term
		delsq_S_RHS = (S[cell1_idx] / (l_sf * l_sf));

		//Terms occuring only in magnetic meshes
		if (stsolve == STSOLVE_FERROMAGNETIC || stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

			cuBReal l_ex = (pcuMesh ? *pcuMesh->pl_ex : *pcuaMesh->pl_ex);
			cuBReal l_ph = (pcuMesh ? *pcuMesh->pl_ph : *pcuaMesh->pl_ph);
			if (pcuMesh)pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pl_ex, l_ex, *pcuMesh->pl_ph, l_ph);
			else pcuaMesh->update_parameters_ecoarse(cell1_idx, *pcuaMesh->pl_ex, l_ex, *pcuaMesh->pl_ph, l_ph);

			int idx_M = M.position_to_cellidx(V.cellidx_to_position(cell1_idx));
			cuReal3 m = cu_normalize(M[idx_M]);

			//transverse S decay terms
			delsq_S_RHS += ((S[cell1_idx] ^ m) / (l_ex * l_ex) + (m ^ (S[cell1_idx] ^ m)) / (l_ph * l_ph));
		}

		//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
		if (V.is_not_empty(cell1_idx)) {

			bool cpump_enabled = cuIsNZ(pcuMesh ? pcuMesh->pcpump_eff->get0() : pcuaMesh->pcpump_eff->get0());
			bool the_enabled = cuIsNZ(pcuMesh ? pcuMesh->pthe_eff->get0() : pcuaMesh->pthe_eff->get0());
			bool she_enabled = cuIsNZ(pcuMesh ? pcuMesh->pSHA->get0() : pcuaMesh->pSHA->get0());

			if (stsolve == STSOLVE_FERROMAGNETIC) {

				//magnetic mesh

				cuBReal P = (pcuMesh ? *pcuMesh->pP : *pcuaMesh->pP);
				cuBReal De = (pcuMesh ? *pcuMesh->pDe : *pcuaMesh->pDe);
				if (pcuMesh) pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pP, P, *pcuMesh->pDe, De);
				else pcuaMesh->update_parameters_ecoarse(cell1_idx, *pcuaMesh->pP, P, *pcuaMesh->pDe, De);

				//term due to drift (non-uniformity of M term, and delsq V contribution - non-uniformity of E term)

				//find grad M and M at the M cell in which the current S cell center is
				int idx_M = M.position_to_cellidx(V.cellidx_to_position(cell1_idx));

				cuReal3 m = cu_normalize(M[idx_M]);
				cuReal33 grad_m = cu_normalize(M.grad_neu(idx_M), M[idx_M]);
				cuReal3 E_dot_del_m = grad_m | E[cell1_idx];

				//E_dot_del_m term is very important, but Evaluate_SpinSolver_delsqV_RHS term could be neglected in most cases especially if E is uniform.
				delsq_S_RHS += (P * (cuBReal)MUB_E * elC[cell1_idx] / De) * (pPoisson_Spin_V->diff2_pri(cell1_idx, shift) * m - E_dot_del_m);

				//charge pumping and topological Hall effect
				if (cpump_enabled || the_enabled) {

					cuReal3 dx_m = grad_m.x;
					cuReal3 dy_m = grad_m.y;
					cuReal3 dxy_m = cu_normalize(M.dxy_neu(idx_M), M[idx_M]);
					cuReal3 dxx_m = cu_normalize(M.dxx_neu(idx_M), M[idx_M]);
					cuReal3 dyy_m = cu_normalize(M.dyy_neu(idx_M), M[idx_M]);

					if (cpump_enabled) {

						cuReal3 dmdt = cu_normalize((*pdM_dt)[idx_M], M[idx_M]);
						cuReal33 grad_dm_dt = cu_normalize((*pdM_dt).grad_neu(idx_M), M[idx_M]);

						cuBReal cpump_eff = (pcuMesh ? pcuMesh->pcpump_eff->get0() : pcuaMesh->pcpump_eff->get0());
						delsq_S_RHS += cpump_eff * (elC[cell1_idx] * (cuBReal)HBAR_E * (cuBReal)MUB_E / (2 * De)) * ((grad_dm_dt.x ^ dx_m) + (grad_dm_dt.y ^ dy_m) + (dmdt ^ (dxx_m + dyy_m)));
					}

					if (the_enabled) {

						cuBReal n_density = (pcuMesh ? *pcuMesh->pn_density : *pcuaMesh->pn_density);
						if (pcuMesh) pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pn_density, n_density);
						else pcuaMesh->update_parameters_ecoarse(cell1_idx, *pcuaMesh->pn_density, n_density);

						cuBReal the_eff = (pcuMesh ? pcuMesh->pthe_eff->get0() : pcuaMesh->pthe_eff->get0());
						delsq_S_RHS += the_eff * ((cuBReal)HBAR_E * (cuBReal)MUB_E * elC[cell1_idx] * elC[cell1_idx] / ((cuBReal)ECHARGE * n_density * De)) * (E[cell1_idx].x * ((dxy_m ^ dy_m) + (dx_m ^ dyy_m)) - E[cell1_idx].y * ((dxx_m ^ dy_m) + (dx_m ^ dxy_m)));
					}
				}
			}

			//terms occuring only in non-magnetic meshes
			else {

				//1. SHA term (this is negligible in most cases, even if E is non-uniform, but might as well include it) 
				if (she_enabled) {

					cuBReal SHA = (pcuMesh ? *pcuMesh->pSHA : *pcuaMesh->pSHA);
					cuBReal De = (pcuMesh ? *pcuMesh->pDe : *pcuaMesh->pDe);
					if (pcuMesh) pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pSHA, SHA, *pcuMesh->pDe, De);
					else pcuaMesh->update_parameters_ecoarse(cell1_idx, *pcuaMesh->pSHA, SHA, *pcuaMesh->pDe, De);

					//Check boundary conditions for this term : should be Dirichlet with 0 Jc value normal to the boundary except for electrodes.
					delsq_S_RHS += (SHA * elC[cell1_idx] * (cuBReal)MUB_E / De) * E.diveps3_sided(cell1_idx);
				}
			}
		}

		return delsq_S_RHS;
	}

	//second order differential of S along the shift axis
	//this is simply Evaluate_SpinSolver_delsqS_RHS from which we subtract second order differentials orthogonal to the shift axis
	__device__ cuReal3 diff2_sec(cuReal3 relpos_m1, cuReal3 stencil, cuReal3 shift)
	{
		cuVEC_VC<cuReal3>& S = (pcuMesh ? *pcuMesh->pS : *pcuaMesh->pS);
		int cellm1_idx = S.position_to_cellidx(relpos_m1);
		return diff2_pri(cellm1_idx, shift);
	}

	//multiply spin accumulation by these to obtain spin potential, i.e. Vs = (De / elC) * (e/muB) * S, evaluated at the boundary
	__device__ cuBReal c_func_pri(int cell_idx)
	{
		cuVEC_VC<cuBReal>& elC = (pcuMesh ? *pcuMesh->pelC : *pcuaMesh->pelC);
		
		cuBReal De = (pcuMesh ? *pcuMesh->pDe : *pcuaMesh->pDe);
		cuBReal elecCond = (pcuMesh ? *pcuMesh->pelecCond : *pcuaMesh->pelecCond);
		if (pcuMesh) pcuMesh->update_parameters_ecoarse(cell_idx, *pcuMesh->pDe, De, *pcuMesh->pelecCond, elecCond);
		else pcuaMesh->update_parameters_ecoarse(cell_idx, *pcuaMesh->pDe, De, *pcuaMesh->pelecCond, elecCond);

		//metallic conduction : use De
		if (elecCond > 0.0) return De / (elC[cell_idx] * (cuBReal)MUB_E);
		//tunelling : use 1
		else return 1.0 / (elC[cell_idx] * (cuBReal)MUB_E);
	}

	__device__ cuBReal c_func_sec(cuReal3 relpos, cuReal3 stencil)
	{
		cuVEC_VC<cuBReal>& elC = (pcuMesh ? *pcuMesh->pelC : *pcuaMesh->pelC);
		
		cuBReal De = (pcuMesh ? *pcuMesh->pDe : *pcuaMesh->pDe);
		cuBReal elecCond = (pcuMesh ? *pcuMesh->pelecCond : *pcuaMesh->pelecCond);
		if (pcuMesh) pcuMesh->update_parameters_atposition(relpos, *pcuMesh->pDe, De, *pcuMesh->pelecCond, elecCond);
		else pcuaMesh->update_parameters_atposition(relpos, *pcuaMesh->pDe, De, *pcuaMesh->pelecCond, elecCond);

		//metallic conduction : use De
		if (elecCond > 0.0) return De / (elC.weighted_average(relpos, stencil) * (cuBReal)MUB_E);
		//tunelling : use 1
		else return 1.0 / (elC.weighted_average(relpos, stencil) * (cuBReal)MUB_E);
	}
};

#endif

#endif