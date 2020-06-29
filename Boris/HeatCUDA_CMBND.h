#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_HEAT

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "ManagedMeshCUDA.h"
#include "ManagedAtom_MeshCUDA.h"

#include "MeshParamsControlCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

class MeshCUDA;

//This is held as a cu_obj managed class in HeatCUDA modules
//It provides methods and access to mesh data for use in cuVEC_VC methods.
//The a_func, b_func and diff2_func methods are used to set CMBND conditions based on the continuity of a quantity and a flux.
//If V is the potential, then the flux is the function f(V) = a_func + b_func * V', where the V' differential direction is perpendicular to the interface.
class HeatCUDA_CMBND {

private:

	ManagedMeshCUDA* pcuMesh;
	ManagedAtom_MeshCUDA* pcuaMesh;

public:

	__host__ void construct_cu_obj(void) 
	{
		nullgpuptr(pcuMesh);
		nullgpuptr(pcuaMesh);
	}

	__host__ void destruct_cu_obj(void) {}

	//used in micromagnetic meshes
	BError set_pointers(MeshCUDA* pMeshCUDA);

	//used in atomistic meshes
	BError set_pointers(Atom_MeshCUDA* paMeshCUDA);

	//heat flux, f(T) = -K * grad T = a + b * grad T -> a = 0, b = -K
	__device__ cuBReal a_func_pri(int cell1_idx, int cell2_idx, cuReal3 shift)
	{
		return 0.0;
	}

	__device__ cuBReal a_func_sec(cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil)
	{
		return 0.0;
	}

	//heat flux, f(T) = -K * grad T = a + b * grad T -> a = 0, b = -K
	__device__ cuBReal b_func_pri(int cell1_idx, int cell2_idx)
	{
		//micromagnetic
		if (pcuMesh) {

			cuBReal thermCond = *pcuMesh->pthermCond;
			pcuMesh->update_parameters_tcoarse(cell1_idx, *pcuMesh->pthermCond, thermCond);

			return -1.0 * thermCond;
		}
		//atomistic
		else if (pcuaMesh) {

			cuBReal thermCond = *pcuaMesh->pthermCond;
			pcuaMesh->update_parameters_tcoarse(cell1_idx, *pcuaMesh->pthermCond, thermCond);

			return -1.0 * thermCond;
		}
		else return 0.0;
	}

	__device__ cuBReal b_func_sec(cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil)
	{
		//micromagnetic
		if (pcuMesh) {

			cuBReal thermCond = *pcuMesh->pthermCond;
			pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pthermCond, thermCond);

			return -1.0 * thermCond;
		}
		//atomistic
		else if (pcuaMesh) {

			cuBReal thermCond = *pcuaMesh->pthermCond;
			pcuaMesh->update_parameters_atposition(relpos_m1, *pcuaMesh->pthermCond, thermCond);

			return -1.0 * thermCond;
		}
		else return 0.0;
	}

	//second order differential of T at cells either side of the boundary; delsq T = -Jc^2 / K * elC - Q / K - many-temperature model coupling terms / K
	__device__ cuBReal diff2_pri(int cell1_idx, cuReal3 shift)
	{
		//micromagnetic
		if (pcuMesh) {

			cuVEC_VC<cuBReal>& Temp = *pcuMesh->pTemp;
			cuVEC_VC<cuBReal>& Temp_l = *pcuMesh->pTemp_l;

			cuVEC_VC<cuReal3>& E = *pcuMesh->pE;
			cuVEC_VC<cuBReal>& elC = *pcuMesh->pelC;

			cuBReal thermCond = *pcuMesh->pthermCond;

			if (E.linear_size() || cuIsNZ(pcuMesh->pQ->get0())) {

				pcuMesh->update_parameters_tcoarse(cell1_idx, *pcuMesh->pthermCond, thermCond);
			}
			else return 0.0;

			cuBReal value = 0.0;

			//Joule heating
			if (E.linear_size()) {

				int idx1_E = E.position_to_cellidx(Temp.cellidx_to_position(cell1_idx));
				value = -(elC[idx1_E] * E[idx1_E] * E[idx1_E]) / thermCond;
			}

			//heat source contribution if set
			if (cuIsNZ(pcuMesh->pQ->get0())) {

				cuBReal Q = *pcuMesh->pQ;
				pcuMesh->update_parameters_tcoarse(cell1_idx, *pcuMesh->pQ, Q);

				value -= Q / thermCond;
			}

			if (Temp_l.linear_size()) {

				cuBReal G_el = *pcuMesh->pG_e;
				pcuMesh->update_parameters_tcoarse(cell1_idx, *pcuMesh->pG_e, G_el);

				value += G_el * (Temp[cell1_idx] - Temp_l[cell1_idx]) / thermCond;
			}

			return value;
		}
		//atomistic
		else if (pcuaMesh) {

			cuVEC_VC<cuBReal>& Temp = *pcuaMesh->pTemp;
			cuVEC_VC<cuBReal>& Temp_l = *pcuaMesh->pTemp_l;

			cuVEC_VC<cuReal3>& E = *pcuaMesh->pE;
			cuVEC_VC<cuBReal>& elC = *pcuaMesh->pelC;

			cuBReal thermCond = *pcuaMesh->pthermCond;

			if (E.linear_size() || cuIsNZ(pcuaMesh->pQ->get0())) {

				pcuaMesh->update_parameters_tcoarse(cell1_idx, *pcuaMesh->pthermCond, thermCond);
			}
			else return 0.0;

			cuBReal value = 0.0;

			//Joule heating
			if (E.linear_size()) {

				int idx1_E = E.position_to_cellidx(Temp.cellidx_to_position(cell1_idx));
				value = -(elC[idx1_E] * E[idx1_E] * E[idx1_E]) / thermCond;
			}

			//heat source contribution if set
			if (cuIsNZ(pcuaMesh->pQ->get0())) {

				cuBReal Q = *pcuaMesh->pQ;
				pcuaMesh->update_parameters_tcoarse(cell1_idx, *pcuaMesh->pQ, Q);

				value -= Q / thermCond;
			}

			if (Temp_l.linear_size()) {

				cuBReal G_el = *pcuaMesh->pG_e;
				pcuaMesh->update_parameters_tcoarse(cell1_idx, *pcuaMesh->pG_e, G_el);

				value += G_el * (Temp[cell1_idx] - Temp_l[cell1_idx]) / thermCond;
			}

			return value;
		}
		else return 0.0;
	}

	__device__ cuBReal diff2_sec(cuReal3 relpos_m1, cuReal3 stencil, cuReal3 shift)
	{
		//micromagnetic
		if (pcuMesh) {

			cuVEC_VC<cuBReal>& Temp = *pcuMesh->pTemp;
			cuVEC_VC<cuBReal>& Temp_l = *pcuMesh->pTemp_l;

			cuVEC_VC<cuReal3>& E = *pcuMesh->pE;
			cuVEC_VC<cuBReal>& elC = *pcuMesh->pelC;

			cuBReal thermCond = *pcuMesh->pthermCond;

			if (E.linear_size() || cuIsNZ(pcuMesh->pQ->get0())) {

				pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pthermCond, thermCond);
			}
			else return 0.0;

			cuBReal value = 0.0;

			//Joule heating
			if (E.linear_size()) {

				int idx1_E = E.position_to_cellidx(relpos_m1);
				value = -(elC[idx1_E] * E[idx1_E] * E[idx1_E]) / thermCond;
			}

			//heat source contribution if set
			if (cuIsNZ(pcuMesh->pQ->get0())) {

				cuBReal Q = *pcuMesh->pQ;
				pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pQ, Q);

				value -= Q / thermCond;
			}

			if (Temp_l.linear_size()) {

				cuBReal G_el = *pcuMesh->pG_e;
				pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pG_e, G_el);

				value += G_el * (Temp.weighted_average(relpos_m1, stencil) - Temp_l.weighted_average(relpos_m1, stencil)) / thermCond;
			}

			return value;
		}
		//atomistic
		else if (pcuaMesh) {

			cuVEC_VC<cuBReal>& Temp = *pcuaMesh->pTemp;
			cuVEC_VC<cuBReal>& Temp_l = *pcuaMesh->pTemp_l;

			cuVEC_VC<cuReal3>& E = *pcuaMesh->pE;
			cuVEC_VC<cuBReal>& elC = *pcuaMesh->pelC;

			cuBReal thermCond = *pcuaMesh->pthermCond;

			if (E.linear_size() || cuIsNZ(pcuaMesh->pQ->get0())) {

				pcuaMesh->update_parameters_atposition(relpos_m1, *pcuaMesh->pthermCond, thermCond);
			}
			else return 0.0;

			cuBReal value = 0.0;

			//Joule heating
			if (E.linear_size()) {

				int idx1_E = E.position_to_cellidx(relpos_m1);
				value = -(elC[idx1_E] * E[idx1_E] * E[idx1_E]) / thermCond;
			}

			//heat source contribution if set
			if (cuIsNZ(pcuaMesh->pQ->get0())) {

				cuBReal Q = *pcuaMesh->pQ;
				pcuaMesh->update_parameters_atposition(relpos_m1, *pcuaMesh->pQ, Q);

				value -= Q / thermCond;
			}

			if (Temp_l.linear_size()) {

				cuBReal G_el = *pcuaMesh->pG_e;
				pcuaMesh->update_parameters_atposition(relpos_m1, *pcuaMesh->pG_e, G_el);

				value += G_el * (Temp.weighted_average(relpos_m1, stencil) - Temp_l.weighted_average(relpos_m1, stencil)) / thermCond;
			}

			return value;
		}
		else return 0.0;
	}
};

#endif

#endif

