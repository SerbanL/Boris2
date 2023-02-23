#include "STransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "BorisCUDALib.cuh"

void STransportCUDA::set_cmbnd_charge_transport(void)
{
	//calculate values at CMBND cells using boundary conditions
	for (int idx1 = 0; idx1 < (int)CMBNDcontacts.size(); idx1++) {

		for (int idx2 = 0; idx2 < (int)CMBNDcontacts[idx1].size(); idx2++) {

			int idx_sec = CMBNDcontacts[idx1][idx2].mesh_idx.i;
			int idx_pri = CMBNDcontacts[idx1][idx2].mesh_idx.j;
			size_t size = CMBNDcontacts[idx1][idx2].cells_box.size().dim();
			
			(*pV[idx_pri])()->set_cmbnd_continuous(
				size, *pV[idx_sec],
				(TransportCUDA_V_Funcs&)pTransport[idx_sec]->poisson_V,
				(TransportCUDA_V_Funcs&)pTransport[idx_pri]->poisson_V,
				(CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);
		}
	}
}

#endif

#endif