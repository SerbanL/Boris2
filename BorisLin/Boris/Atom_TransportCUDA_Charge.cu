#include "Atom_TransportCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"

#include "BorisCUDALib.cuh"

//-------------------Calculation Methods

void Atom_TransportCUDA::IterateChargeSolver_SOR(cu_obj<cuBReal>& damping, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value)
{
	//Note, TransportCUDA_V_Funcs covers both thermoelectric and no thermoelectric effect cases (but for thermoelectric effect must use the NNeu Poisson solver as nonhomogeneous Neumann boundary conditions are required for V)

	if (!is_thermoelectric_mesh) {

		//no thermoelectric effect
		paMeshCUDA->V()->IteratePoisson_SOR(paMeshCUDA->n_e.dim(), (TransportCUDA_V_Funcs&)poisson_V, damping, max_error, max_value);
	}
	else {

		//include thermoelectric effect
		paMeshCUDA->V()->IteratePoisson_NNeu_SOR(paMeshCUDA->n_e.dim(), (TransportCUDA_V_Funcs&)poisson_V, damping, max_error, max_value);
	}
}

#endif

#endif