#include "stdafx.h"
#include "STransport.h"

#ifdef MODULE_COMPILATION_TRANSPORT

#include "SuperMesh.h"

//Calculate interface spin accumulation torque only in mesh with matching transport module
VEC<DBL3>& STransport::GetInterfacialSpinTorque(Transport* pMeshTrans)
{
	if (!pMeshTrans->PrepareDisplayVEC(pMeshTrans->pMesh->h))
		return pMeshTrans->displayVEC;

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		cu_obj<cuVEC<cuReal3>>& cudisplayVEC = GetInterfacialSpinTorqueCUDA(pMeshTrans);
		cudisplayVEC()->copy_to_cpuvec(pMeshTrans->displayVEC);

		return pMeshTrans->displayVEC;
	}
#endif

	//calculate interfacial spin torque in displayVEC from all contacts with matching mesh
	for (int idx1 = 0; idx1 < (int)CMBNDcontacts.size(); idx1++) {

		for (int idx2 = 0; idx2 < (int)CMBNDcontacts[idx1].size(); idx2++) {

			int idx_sec = CMBNDcontacts[idx1][idx2].mesh_idx.i;
			int idx_pri = CMBNDcontacts[idx1][idx2].mesh_idx.j;

			if (pTransport[idx_pri] == pMeshTrans) pTransport[idx_pri]->CalculateDisplaySAInterfaceTorque(pTransport[idx_sec], CMBNDcontacts[idx1][idx2]);
		}
	}

	return pMeshTrans->displayVEC;
}

#if COMPILECUDA == 1
//return interfacial spin torque in given mesh with matching transport module
cu_obj<cuVEC<cuReal3>>& STransport::GetInterfacialSpinTorqueCUDA(Transport* pMeshTrans)
{
	return dynamic_cast<STransportCUDA*>(pModuleCUDA)->GetInterfacialSpinTorque(dynamic_cast<TransportCUDA*>(pMeshTrans->pModuleCUDA));
}
#endif

#endif
