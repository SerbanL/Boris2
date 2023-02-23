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

//This is held as a cu_obj managed class in TransportCUDA modules
//It provides methods and access to mesh data for use in cuVEC_VC methods.
//The methods have fixed names, e.g. Poisson_RHS is used by Poisson solvers to evaluate the r.h.s. of the Poisson equation
//The a_func, b_func and diff2_func methods are used to set CMBND conditions based on the continuity of a quantity and a flux.
//If V is the potential, then the flux is the function f(V) = a_func + b_func * V', where the V' differential direction is perpendicular to the interface.
//This particular class is used for charge transport only.
class TransportCUDA_V_Funcs {

public:

	//managed mesh for access to all required mesh VECs and material parameters
	ManagedMeshCUDA* pcuMesh;

	//managed mesh for access to all required mesh VECs and material parameters
	ManagedAtom_MeshCUDA* pcuaMesh;

	//flags set
	bool is_thermoelectric_mesh;
	bool is_open_potential;

public:

	__host__ void construct_cu_obj(void) {}
	__host__ void destruct_cu_obj(void) {}

	//for modules held in micromagnetic meshes
	BError set_pointers_transport(MeshCUDA* pMeshCUDA);
	//for modules held in atomistic meshes
	BError set_pointers_atomtransport(Atom_MeshCUDA* pMeshCUDA);

	__host__ void set_thermoelectric_mesh_flag(bool status) { set_gpu_value(is_thermoelectric_mesh, status); }
	__host__ void set_open_potential_flag(bool status) { set_gpu_value(is_open_potential, status); }

	//this evaluates the Poisson RHS when solving the Poisson equation on V
	__device__ cuBReal Poisson_RHS(int idx)
	{
		cuVEC_VC<cuBReal>& V = (pcuMesh ? *pcuMesh->pV : *pcuaMesh->pV);
		cuVEC_VC<cuBReal>& elC = (pcuMesh ? *pcuMesh->pelC : *pcuaMesh->pelC);

		if (!is_thermoelectric_mesh) {

			//no thermoelectric effect
			return -(V.grad_diri(idx) * elC.grad_sided(idx)) / elC[idx];
		}
		else {

			//include thermoelectric effect with Seebeck coefficient : delsq V = -S delsq T, obtained from div J = 0, where J = -sigma(grad V + S * grad T).
			//here we ignore gradients in sigma and S

			cuBReal Sc = (pcuMesh ? *pcuMesh->pSc : *pcuaMesh->pSc);
			cuBReal thermCond = (pcuMesh ? *pcuMesh->pthermCond : *pcuaMesh->pthermCond);
			if (pcuMesh) pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->pSc, Sc, *pcuMesh->pthermCond, thermCond);
			else pcuaMesh->update_parameters_ecoarse(idx, *pcuaMesh->pSc, Sc, *pcuaMesh->pthermCond, thermCond);

			cuVEC_VC<cuBReal>& Temp = (pcuMesh ? *pcuMesh->pTemp : *pcuaMesh->pTemp);

			//corresponding index in Temp
			int idx_temp = Temp.position_to_cellidx(V.cellidx_to_position(idx));

			return -Sc * Temp.delsq_robin(idx_temp, thermCond);
		}
	}

	//boundary differential of V for non-homogeneous Neumann boundary conditions (when thermoelectric effect is used)
	__device__ cuReal3 bdiff(int idx)
	{
		//Gradient of V normal to boundary is -S * grad T normal to boundary.
		//grad T normal to boundary obtained using Robin boundary condition as:

		//-K * grad T . n = heat flux normal to boundary = alpha(Tb - Ta), where alpha is Robin value, Ta ambient temperature and Tb temperature at boundary, K thermal conductivity
		//Thus grad V . n = S*alpha*(Tb-Ta) / K

		cuVEC_VC<cuBReal>& V = (pcuMesh ? *pcuMesh->pV : *pcuaMesh->pV);
		cuVEC_VC<cuBReal>& Temp = (pcuMesh ? *pcuMesh->pTemp : *pcuaMesh->pTemp);

		cuReal3 bdiff = cuReal3();
		cuReal3 shift = V.get_shift_to_emptycell(idx);

		if (!shift.IsNull()) {

			//corresponding index in Temp
			int idx_temp = Temp.position_to_cellidx(V.cellidx_to_position(idx));

			if (Temp.is_cmbnd(idx_temp)) {

				//at composite media boundaries (for Temp) cannot use Robin boundary for heat flux
				//instead we can approximate it using a sided differential in the cell just next to the boundary
				cuBReal Sc = (pcuMesh ? *pcuMesh->pSc : *pcuaMesh->pSc);
				if (pcuMesh) pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->pSc, Sc);
				else pcuaMesh->update_parameters_ecoarse(idx, *pcuaMesh->pSc, Sc);

				bdiff = Sc * Temp.grad_sided(idx_temp);
			}
			else {

				//boundary, not a cmbnd. Use grad V . n = S*alpha*(Tb-Ta) / K
				cuBReal Sc = (pcuMesh ? *pcuMesh->pSc : *pcuaMesh->pSc);
				cuBReal thermCond = (pcuMesh ? *pcuMesh->pthermCond : *pcuaMesh->pthermCond);
				if (pcuMesh) pcuMesh->update_parameters_ecoarse(idx, *pcuMesh->pSc, Sc, *pcuMesh->pthermCond, thermCond);
				else pcuaMesh->update_parameters_ecoarse(idx, *pcuaMesh->pSc, Sc, *pcuaMesh->pthermCond, thermCond);

				bdiff = Temp.get_robin_value(V.cellidx_to_position(idx), shift) * Sc / thermCond;
				//at negative boundary inverse sign since heat flux normal also has opposite sign
				if (shift.x < 0) bdiff.x *= -1;
				if (shift.y < 0) bdiff.y *= -1;
				if (shift.z < 0) bdiff.z *= -1;
			}
		}

		return bdiff;
	}

	//Charge transport only : V

	//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface
	//With thermoelectric effect this becomes:
	//1. No net current (e.g. no electrodes attached):
	//Jc = -sigma * grad V - sigma * Sc * grad T = a + b * grad V -> a = -sigma * Sc * grad T and b = -sigma taken at the interface
	//grad T normal to interface is given by Robin boundary condition : grad T . n = -alpha(Tb-Ta)/K
	//2. Net current generated (open potential condition):
	//Jc = sigma * Sc * grad T, so a = sigma * Sc * grad T, b = 0
	__device__ cuBReal a_func_pri(int cell1_idx, int cell2_idx, cuReal3 shift)
	{
		if (!is_thermoelectric_mesh) return 0.0;
		else {

			//include thermoelectric effect

			cuVEC_VC<cuBReal>& V = (pcuMesh ? *pcuMesh->pV : *pcuaMesh->pV);
			cuVEC_VC<cuBReal>& elC = (pcuMesh ? *pcuMesh->pelC : *pcuaMesh->pelC);
			cuVEC_VC<cuBReal>& Temp = (pcuMesh ? *pcuMesh->pTemp : *pcuaMesh->pTemp);

			cuBReal Sc = (pcuMesh ? *pcuMesh->pSc : *pcuaMesh->pSc);
			if (pcuMesh) pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pSc, Sc);
			else pcuaMesh->update_parameters_ecoarse(cell1_idx, *pcuaMesh->pSc, Sc);

			//normalized shift in normal direction to boundary: use * operator (dot product) with nshift to eliminate differentials orthogonal to the shift axis
			//do not use mod here as we need the temperature gradient to point in normal direction to boundary
			cuReal3 nshift = cu_normalize(shift);

			//corresponding index in Temp
			int idx_temp1 = Temp.position_to_cellidx(V.cellidx_to_position(cell1_idx));
			int idx_temp2 = Temp.position_to_cellidx(V.cellidx_to_position(cell2_idx));

			cuBReal T_grad1 = Temp.grad_sided(idx_temp1) * nshift;
			cuBReal T_grad2 = Temp.grad_sided(idx_temp2) * nshift;
			cuBReal T_grad = 1.5 * T_grad1 - 0.5 * T_grad2;

			//shift is from cell1 to cell2 so no need for minus sign adjustment
			return Sc * elC[cell1_idx] * T_grad;
		}
	}

	__device__ cuBReal a_func_sec(cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil)
	{
		if (!is_thermoelectric_mesh) return 0.0;
		else {

			//include thermoelectric effect

			cuVEC_VC<cuBReal>& V = (pcuMesh ? *pcuMesh->pV : *pcuaMesh->pV);
			cuVEC_VC<cuBReal>& elC = (pcuMesh ? *pcuMesh->pelC : *pcuaMesh->pelC);
			cuVEC_VC<cuBReal>& Temp = (pcuMesh ? *pcuMesh->pTemp : *pcuaMesh->pTemp);

			int cellm1_idx = V.position_to_cellidx(relpos_m1);

			//corresponding index in Temp
			int idx_temp1 = Temp.position_to_cellidx(relpos_m1);
			int idx_temp2 = Temp.position_to_cellidx(relpos_m1 + shift);

			cuBReal Sc = (pcuMesh ? *pcuMesh->pSc : *pcuaMesh->pSc);
			if (pcuMesh) pcuMesh->update_parameters_ecoarse(cellm1_idx, *pcuMesh->pSc, Sc);
			else pcuaMesh->update_parameters_ecoarse(cellm1_idx, *pcuaMesh->pSc, Sc);

			//normalized shift in normal direction to boundary: use * operator (dot product) with nshift to eliminate differentials orthogonal to the shift axis
			//do not use mod here as we need the temperature gradient to point in normal direction to boundary
			cuReal3 nshift = cu_normalize(shift);

			cuBReal T_grad1 = Temp.grad_sided(idx_temp1) * nshift;
			cuBReal T_grad2 = Temp.grad_sided(idx_temp2) * nshift;
			cuBReal T_grad = 1.5 * T_grad1 - 0.5 * T_grad2;

			//shift is from m1 to m2 cell, so use minus sign here if open potential mode
			if (is_open_potential) return -Sc * elC[cellm1_idx] * T_grad;
			else return Sc * elC[cellm1_idx] * T_grad;
		}
	}

	//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface	
	__device__ cuBReal b_func_pri(int cell1_idx, int cell2_idx)
	{
		if (is_thermoelectric_mesh && is_open_potential) return 0.0;
		else {

			cuVEC_VC<cuBReal>& elC = (pcuMesh ? *pcuMesh->pelC : *pcuaMesh->pelC);
			return -(1.5 * elC[cell1_idx] - 0.5 * elC[cell2_idx]);
		}
	}

	__device__ cuBReal b_func_sec(cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil)
	{
		if (is_thermoelectric_mesh && is_open_potential) return 0.0;
		else {

			cuVEC_VC<cuBReal>& elC = (pcuMesh ? *pcuMesh->pelC : *pcuaMesh->pelC);
			return -(1.5 * elC.weighted_average(relpos_m1, stencil) - 0.5 * elC.weighted_average(relpos_m1 + shift, stencil));
		}
	}

	//second order differential of V at cells either side of the boundary; delsq V = -grad V * grad elC / elC
	__device__ cuBReal diff2_pri(int cell1_idx, cuReal3 shift)
	{
		if (!is_thermoelectric_mesh) {

			//normalized, positive shift: use * operator (dot product) with nshift to eliminate differentials orthogonal to the shift axis
			cuReal3 nshift = cu_mod(cu_normalize(shift));

			cuVEC_VC<cuBReal>& V = (pcuMesh ? *pcuMesh->pV : *pcuaMesh->pV);
			cuVEC_VC<cuBReal>& elC = (pcuMesh ? *pcuMesh->pelC : *pcuaMesh->pelC);

			//no thermoelectric effect
			return -((V.grad_diri(cell1_idx) * nshift) * (elC.grad_sided(cell1_idx) * nshift)) / elC[cell1_idx];
		}
		else {

			//include thermoelectric effect with Seebeck coefficient

			cuVEC_VC<cuBReal>& V = (pcuMesh ? *pcuMesh->pV : *pcuaMesh->pV);
			cuVEC_VC<cuBReal>& Temp = (pcuMesh ? *pcuMesh->pTemp : *pcuaMesh->pTemp);

			cuBReal Sc = (pcuMesh ? *pcuMesh->pSc : *pcuaMesh->pSc);
			if (pcuMesh) pcuMesh->update_parameters_ecoarse(cell1_idx, *pcuMesh->pSc, Sc);
			else pcuaMesh->update_parameters_ecoarse(cell1_idx, *pcuaMesh->pSc, Sc);

			//corresponding index in Temp
			int idx_temp = Temp.position_to_cellidx(V.cellidx_to_position(cell1_idx));

			return -Sc * Temp.delsq_neu(idx_temp);
		}
	}

	__device__ cuBReal diff2_sec(cuReal3 relpos_m1, cuReal3 stencil, cuReal3 shift)
	{
		cuVEC_VC<cuBReal>& V = (pcuMesh ? *pcuMesh->pV : *pcuaMesh->pV);

		int cellm1_idx = V.position_to_cellidx(relpos_m1);
		return diff2_pri(cellm1_idx, shift);
	}
};

#endif

#endif
