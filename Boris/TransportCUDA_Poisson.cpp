#include "stdafx.h"
#include "TransportCUDA_Poisson.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "MeshCUDA.h"
#include "Atom_MeshCUDA.h"

//for modules held in micromagnetic meshes
BError TransportCUDA_V_Funcs::set_pointers_transport(MeshCUDA* pMeshCUDA)
{
	BError error(__FUNCTION__);

	//Mesh quantities
	
	if (set_gpu_value(pV, pMeshCUDA->V.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pelC, pMeshCUDA->elC.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	
	return error;
}

//for modules held in atomistic meshes
BError TransportCUDA_V_Funcs::set_pointers_atomtransport(Atom_MeshCUDA* paMeshCUDA)
{
	BError error(__FUNCTION__);

	//Mesh quantities

	if (set_gpu_value(pV, paMeshCUDA->V.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pelC, paMeshCUDA->elC.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	return error;
}

#endif

#endif