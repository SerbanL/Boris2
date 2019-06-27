#include "SHeatCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_HEAT

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"

//calculate and set values at composite media boundaries after all other cells have been computed and set
void SHeatCUDA::set_cmbnd_values(void)
{
	//calculate values at CMBND cells using boundary conditions
	for (int idx1 = 0; idx1 < (int)CMBNDcontacts.size(); idx1++) {

		for (int idx2 = 0; idx2 < (int)CMBNDcontacts[idx1].size(); idx2++) {

			int idx_sec = CMBNDcontacts[idx1][idx2].mesh_idx.i;
			int idx_pri = CMBNDcontacts[idx1][idx2].mesh_idx.j;
			size_t size = CMBNDcontacts[idx1][idx2].cells_box.size().dim();

			(*pTemp[idx_pri])()->set_cmbnd_continuous(size, *pTemp[idx_sec], (HeatCUDA_CMBND&)pHeat[idx_sec]->temp_cmbnd_funcs, (HeatCUDA_CMBND&)pHeat[idx_pri]->temp_cmbnd_funcs, (CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);
		}
	}
}

#endif

#endif