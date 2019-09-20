#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_HEAT

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "ManagedMeshCUDA.h"

#include "MeshParamsControlCUDA.h"

class MeshCUDA;

//This is held as a cu_obj managed class in HeatCUDA modules
//It provides methods and access to mesh data for use in cuVEC_VC methods.
//The a_func, b_func and diff2_func methods are used to set CMBND conditions based on the continuity of a quantity and a flux.
//If V is the potential, then the flux is the function f(V) = a_func + b_func * V', where the V' differential direction is perpendicular to the interface.
class HeatCUDA_CMBND {

private:

	ManagedMeshCUDA* pcuMesh;

public:

	__host__ void construct_cu_obj(void) {}
	__host__ void destruct_cu_obj(void) {}

	BError set_pointers(MeshCUDA* pMeshCUDA);

	//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface	
	__device__ cuBReal a_func_pri(int cell1_idx, int cell2_idx, cuReal3 shift)
	{
		return 0.0;
	}

	__device__ cuBReal a_func_sec(cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil)
	{
		return 0.0;
	}

	//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface	
	__device__ cuBReal b_func_pri(int cell1_idx, int cell2_idx)
	{
		cuBReal thermCond = *pcuMesh->pthermCond;
		pcuMesh->update_parameters_tcoarse(cell1_idx, *pcuMesh->pthermCond, thermCond);

		return -1.0 * thermCond;
	}

	__device__ cuBReal b_func_sec(cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil)
	{
		cuBReal thermCond = *pcuMesh->pthermCond;
		pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pthermCond, thermCond);

		return -1.0 * thermCond;
	}

	//second order differential of V at cells either side of the boundary; delsq V = -grad V * grad elC / elC
	__device__ cuBReal diff2_pri(int cell1_idx)
	{
		cuVEC_VC<cuBReal>& Temp = *pcuMesh->pTemp;
		cuVEC_VC<cuReal3>& Jc = *pcuMesh->pJc;
		cuVEC_VC<cuBReal>& elC = *pcuMesh->pelC;

		cuBReal thermCond = *pcuMesh->pthermCond;

		if (Jc.linear_size() || cuIsNZ(pcuMesh->pQ->get0())) {

			pcuMesh->update_parameters_tcoarse(cell1_idx, *pcuMesh->pthermCond, thermCond);
		}
		else return 0.0;

		cuBReal value = 0.0;
		
		//Joule heating
		if (Jc.linear_size()) {

			int idx1_Jc = Jc.position_to_cellidx(Temp.cellidx_to_position(cell1_idx));
			value = -(Jc[idx1_Jc] * Jc[idx1_Jc]) / (elC[idx1_Jc] * thermCond);
		}

		//heat source contribution if set
		if (cuIsNZ(pcuMesh->pQ->get0())) {

			cuBReal Q = *pcuMesh->pQ;
			pcuMesh->update_parameters_tcoarse(cell1_idx, *pcuMesh->pQ, Q);

			value -= Q / thermCond;
		}

		return value;
	}

	__device__ cuBReal diff2_sec(cuReal3 relpos_m1, cuReal3 stencil)
	{
		cuVEC_VC<cuReal3>& Jc = *pcuMesh->pJc;
		cuVEC_VC<cuBReal>& elC = *pcuMesh->pelC;

		cuBReal thermCond = *pcuMesh->pthermCond;

		if (Jc.linear_size() || cuIsNZ(pcuMesh->pQ->get0())) {

			pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pthermCond, thermCond);
		}
		else return 0.0;

		cuBReal value = 0.0;

		//Joule heating
		if (Jc.linear_size()) {

			int idx1_Jc = Jc.position_to_cellidx(relpos_m1);
			value = -(Jc[idx1_Jc] * Jc[idx1_Jc]) / (elC[idx1_Jc] * thermCond);
		}
		
		//heat source contribution if set
		if (cuIsNZ(pcuMesh->pQ->get0())) {

			cuBReal Q = *pcuMesh->pQ;
			pcuMesh->update_parameters_atposition(relpos_m1, *pcuMesh->pQ, Q);

			value -= Q / thermCond;
		}

		return value;
	}
};

#endif

#endif

