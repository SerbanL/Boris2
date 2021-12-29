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

			//1. Micromagnetic to micromagnetic meshes
			if (!pTransport[idx_pri]->pMeshBaseCUDA->is_atomistic() && !pTransport[idx_sec]->pMeshBaseCUDA->is_atomistic()) {

				(*pV[idx_pri])()->set_cmbnd_continuous(
					size, *pV[idx_sec], 
					(TransportCUDA_V_Funcs&)dynamic_cast<TransportCUDA*>(pTransport[idx_sec])->poisson_V, 
					(TransportCUDA_V_Funcs&)dynamic_cast<TransportCUDA*>(pTransport[idx_pri])->poisson_V, 
					(CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);
			}
			//2. Micromagnetic to atomistic meshes
			else if (!pTransport[idx_pri]->pMeshBaseCUDA->is_atomistic() && pTransport[idx_sec]->pMeshBaseCUDA->is_atomistic()) {

				(*pV[idx_pri])()->set_cmbnd_continuous(
					size, *pV[idx_sec],
					(Atom_TransportCUDA_V_Funcs&)dynamic_cast<Atom_TransportCUDA*>(pTransport[idx_sec])->poisson_V,
					(TransportCUDA_V_Funcs&)dynamic_cast<TransportCUDA*>(pTransport[idx_pri])->poisson_V,
					(CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);
			}
			//3. atomistic to micromagnetic meshes
			else if (pTransport[idx_pri]->pMeshBaseCUDA->is_atomistic() && !pTransport[idx_sec]->pMeshBaseCUDA->is_atomistic()) {

				(*pV[idx_pri])()->set_cmbnd_continuous(
					size, *pV[idx_sec],
					(TransportCUDA_V_Funcs&)dynamic_cast<TransportCUDA*>(pTransport[idx_sec])->poisson_V,
					(Atom_TransportCUDA_V_Funcs&)dynamic_cast<Atom_TransportCUDA*>(pTransport[idx_pri])->poisson_V,
					(CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);
			}
			//4. atomistic to atomistic meshes
			else if (pTransport[idx_pri]->pMeshBaseCUDA->is_atomistic() && pTransport[idx_sec]->pMeshBaseCUDA->is_atomistic()) {

				(*pV[idx_pri])()->set_cmbnd_continuous(
					size, *pV[idx_sec],
					(Atom_TransportCUDA_V_Funcs&)dynamic_cast<Atom_TransportCUDA*>(pTransport[idx_sec])->poisson_V, 
					(Atom_TransportCUDA_V_Funcs&)dynamic_cast<Atom_TransportCUDA*>(pTransport[idx_pri])->poisson_V,
					(CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);
			}
		}
	}
}

#endif

#endif