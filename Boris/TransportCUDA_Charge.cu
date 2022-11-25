#include "TransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "MeshCUDA.h"

#include "BorisCUDALib.cuh"

//-------------------Calculation Methods

void TransportCUDA::IterateChargeSolver_SOR(cu_obj<cuBReal>& damping, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value)
{
	//Note, TransportCUDA_V_Funcs covers both thermoelectric and no thermoelectric effect cases (but for thermoelectric effect must use the NNeu Poisson solver as nonhomogeneous Neumann boundary conditions are required for V)
	
	if (!is_thermoelectric_mesh) {

		//no thermoelectric effect
		pMeshCUDA->V()->IteratePoisson_SOR(pMeshCUDA->n_e.dim(), (TransportCUDA_V_Funcs&)poisson_V, damping, max_error, max_value);
	}
	else {

		//include thermoelectric effect
		pMeshCUDA->V()->IteratePoisson_NNeu_SOR(pMeshCUDA->n_e.dim(), (TransportCUDA_V_Funcs&)poisson_V, damping, max_error, max_value);
	}
}

#endif

#endif