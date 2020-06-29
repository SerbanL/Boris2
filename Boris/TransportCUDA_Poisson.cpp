#include "stdafx.h"
#include "TransportCUDA_Poisson.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "MeshCUDA.h"

BError TransportCUDA_V_Funcs::set_pointers(MeshCUDA* pMeshCUDA)
{
	BError error(__FUNCTION__);

	//Mesh quantities
	
	if (set_gpu_value(pV, pMeshCUDA->V.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pelC, pMeshCUDA->elC.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	return error;
}

#endif

#endif