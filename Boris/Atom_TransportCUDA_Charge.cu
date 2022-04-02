#include "Atom_TransportCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"

#include "BorisCUDALib.cuh"

//-------------------Calculation Methods

void Atom_TransportCUDA::IterateChargeSolver_SOR(cu_obj<cuBReal>& damping, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value)
{
	paMeshCUDA->V()->IteratePoisson_SOR(paMeshCUDA->n_e.dim(), (TransportCUDA_V_Funcs&)poisson_V, damping, max_error, max_value);
}

#endif

#endif