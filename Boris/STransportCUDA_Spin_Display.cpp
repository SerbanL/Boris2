#include "stdafx.h"
#include "STransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "STransport.h"
#include "SuperMesh.h"

//-------------------Getters

//return interfacial spin torque in given mesh with matching transport module
cu_obj<cuVEC<cuReal3>>& STransportCUDA::GetInterfacialSpinTorque(TransportBaseCUDA* pTransportBaseCUDA)
{
	if (!pTransportBaseCUDA->PrepareDisplayVEC(pTransportBaseCUDA->pMeshBaseCUDA->h))
		return pTransportBaseCUDA->displayVEC;

	//calculate and add Ts values in Ts_interf VECs
	for (int idx1 = 0; idx1 < (int)CMBNDcontacts.size(); idx1++) {

		for (int idx2 = 0; idx2 < (int)CMBNDcontacts[idx1].size(); idx2++) {

			int idx_sec = CMBNDcontacts[idx1][idx2].mesh_idx.i;
			int idx_pri = CMBNDcontacts[idx1][idx2].mesh_idx.j;

			if (pTransport[idx_pri] == pTransportBaseCUDA) {

				pTransport[idx_pri]->CalculateDisplaySAInterfaceTorque(pTransport[idx_sec], CMBNDcontactsCUDA[idx1][idx2], CMBNDcontacts[idx1][idx2].IsPrimaryTop());
			}
		}
	}

	return pTransportBaseCUDA->displayVEC;
}

#endif

#endif