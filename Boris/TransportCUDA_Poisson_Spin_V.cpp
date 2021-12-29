#include "stdafx.h"
#include "TransportCUDA_Poisson_Spin_V.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "MeshCUDA.h"
#include "TransportCUDA.h"

BError TransportCUDA_Spin_V_Funcs::set_pointers(MeshCUDA* pMeshCUDA, TransportCUDA* pTransportCUDA)
{
	BError error(__FUNCTION__);

	//Mesh quantities

	if (set_gpu_value(pcuMesh, pMeshCUDA->cuMesh.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pdM_dt, pTransportCUDA->dM_dt.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pdelsq_V_fixed, pTransportCUDA->delsq_V_fixed.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(stsolve, pTransportCUDA->Get_STSolveType()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	return error;
}

#endif

#endif