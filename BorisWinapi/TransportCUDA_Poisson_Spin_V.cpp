#include "stdafx.h"
#include "TransportCUDA_Poisson_Spin_V.h"

#if COMPILECUDA == 1

#ifdef MODULE_TRANSPORT

#include "MeshCUDA.h"

BError TransportCUDA_Spin_V_Funcs::set_pointers(MeshCUDA* pMeshCUDA, cu_obj<TransportCUDA_Spin_S_Funcs>& poisson_Spin_S)
{
	BError error(__FUNCTION__);

	//Mesh quantities

	if (set_gpu_value(pcuMesh, pMeshCUDA->cuMesh.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pPoisson_Spin_S, poisson_Spin_S.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	return error;
}

#endif

#endif