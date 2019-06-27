#include "stdafx.h"
#include "HeatCUDA_CMBND.h"

#if COMPILECUDA == 1

#ifdef MODULE_HEAT

#include "MeshCUDA.h"

BError HeatCUDA_CMBND::set_pointers(MeshCUDA* pMeshCUDA)
{
	BError error(__FUNCTION__);

	if (set_gpu_value(pcuMesh, pMeshCUDA->cuMesh.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	return error;
}

#endif

#endif