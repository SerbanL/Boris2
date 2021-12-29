#include "stdafx.h"
#include "Atom_TransportCUDA_Poisson.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"

BError Atom_TransportCUDA_V_Funcs::set_pointers(Atom_MeshCUDA* paMeshCUDA)
{
	BError error(__FUNCTION__);

	//Mesh quantities
	
	if (set_gpu_value(pV, paMeshCUDA->V.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pelC, paMeshCUDA->elC.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	return error;
}

#endif

#endif