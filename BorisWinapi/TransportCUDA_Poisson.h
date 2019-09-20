#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_TRANSPORT

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

class MeshCUDA;

//This is held as a cu_obj managed class in TransportCUDA modules
//It provides methods and access to mesh data for use in cuVEC_VC methods.
//The methods have fixed names, e.g. Poisson_RHS is used by Poisson solvers to evaluate the r.h.s. of the Poisson equation
//The a_func, b_func and diff2_func methods are used to set CMBND conditions based on the continuity of a quantity and a flux.
//If V is the potential, then the flux is the function f(V) = a_func + b_func * V', where the V' differential direction is perpendicular to the interface.
//This particular class is used for charge transport only.
class TransportCUDA_V_Funcs {

public:

	//these point the V and elC quantities held in *pMeshCUDA
	cuVEC_VC<cuBReal>* pV;
	cuVEC_VC <cuBReal>* pelC;

public:

	__host__ void construct_cu_obj(void) {}
	__host__ void destruct_cu_obj(void) {}

	BError set_pointers(MeshCUDA* pMeshCUDA);

	//this evaluates the Poisson RHS when solving the Poisson equation on V
	__device__ cuBReal Poisson_RHS(int idx)
	{
		cuVEC_VC<cuBReal>& V = *pV;
		cuVEC_VC<cuBReal>& elC = *pelC;

		return -(V.grad_diri(idx) * elC.grad_sided(idx)) / elC[idx];
	}

	//Charge transport only : V

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
		cuVEC_VC<cuBReal>& elC = *pelC;

		return -(1.5 * elC[cell1_idx] - 0.5 * elC[cell2_idx]);
	}

	__device__ cuBReal b_func_sec(cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil)
	{
		cuVEC_VC<cuBReal>& elC = *pelC;

		return -(1.5 * elC.weighted_average(relpos_m1, stencil) - 0.5 * elC.weighted_average(relpos_m1 + shift, stencil));
	}

	//second order differential of V at cells either side of the boundary; delsq V = -grad V * grad elC / elC
	__device__ cuBReal diff2_pri(int cell1_idx)
	{
		cuVEC_VC<cuBReal>& V = *pV;
		cuVEC_VC<cuBReal>& elC = *pelC;

		return -(V.grad_diri(cell1_idx) * elC.grad_sided(cell1_idx)) / elC[cell1_idx];
	}

	__device__ cuBReal diff2_sec(cuReal3 relpos_m1, cuReal3 stencil)
	{
		cuVEC_VC<cuBReal>& V = *pV;
		cuVEC_VC<cuBReal>& elC = *pelC;

		int cellm1_idx = V.position_to_cellidx(relpos_m1);

		return -(V.grad_diri(cellm1_idx) * elC.grad_sided(cellm1_idx)) / elC[cellm1_idx];
	}
};

#endif

#endif
