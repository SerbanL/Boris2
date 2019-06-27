#include "STransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_TRANSPORT

#include "BorisCUDALib.cuh"

#include "CUDAError.h"

#include "MeshCUDA.h"

void STransportCUDA::set_cmbnd_spin_transport_V(void)
{
	//calculate values at CMBND cells using boundary conditions
	for (int idx1 = 0; idx1 < (int)CMBNDcontacts.size(); idx1++) {

		for (int idx2 = 0; idx2 < (int)CMBNDcontacts[idx1].size(); idx2++) {

			int idx_sec = CMBNDcontacts[idx1][idx2].mesh_idx.i;
			int idx_pri = CMBNDcontacts[idx1][idx2].mesh_idx.j;
			size_t size = CMBNDcontacts[idx1][idx2].cells_box.size().dim();

			//use continuity of Jc and V across interface unless the interface is N-F type (normal metal - ferromagnetic) and the spin mixing conductance is not zero (i.e. continuous method disabled).

			//Is it an N-F contact?
			if ((pTransport[idx_pri]->pMeshCUDA->MComputation_Enabled() && !pTransport[idx_sec]->pMeshCUDA->Magnetisation_Enabled()) ||
				(pTransport[idx_sec]->pMeshCUDA->MComputation_Enabled() && !pTransport[idx_pri]->pMeshCUDA->Magnetisation_Enabled())) {

				//Yes we have an N-F contact. Is G interface enabled for this contact ? (remember top mesh sets G interface values)

				//if primary is top then we check GInterface in primary. If primary is bottom then we check GInterface in secondary as it must be top.
				if ((CMBNDcontacts[idx1][idx2].IsPrimaryTop() && pTransport[idx_pri]->pMeshCUDA->GInterface_Enabled()) ||
					(!CMBNDcontacts[idx1][idx2].IsPrimaryTop() && pTransport[idx_sec]->pMeshCUDA->GInterface_Enabled())) {

					//G interface method

					(*pV[idx_pri])()->set_cmbnd_continuousflux(
						size, *pV[idx_sec],
						(TransportCUDA_Spin_V_Funcs&)pTransport[idx_sec]->poisson_Spin_V, (TransportCUDA_Spin_V_Funcs&)pTransport[idx_pri]->poisson_Spin_V, (STransportCUDA_GInterf_V_Funcs&)gInterf_V,
						(CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);

					//next contact
					continue;
				}
			}

			//continuous interface method - the G interface method check above didn't pass so this is what we have left

			(*pV[idx_pri])()->set_cmbnd_continuous(
				size, *pV[idx_sec],
				(TransportCUDA_Spin_V_Funcs&)pTransport[idx_sec]->poisson_Spin_V, (TransportCUDA_Spin_V_Funcs&)pTransport[idx_pri]->poisson_Spin_V,
				(CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);
		}
	}
}

void STransportCUDA::set_cmbnd_spin_transport_S(void)
{
	//calculate values at CMBND cells using boundary conditions
	for (int idx1 = 0; idx1 < (int)CMBNDcontacts.size(); idx1++) {

		for (int idx2 = 0; idx2 < (int)CMBNDcontacts[idx1].size(); idx2++) {

			int idx_sec = CMBNDcontacts[idx1][idx2].mesh_idx.i;
			int idx_pri = CMBNDcontacts[idx1][idx2].mesh_idx.j;
			size_t size = CMBNDcontacts[idx1][idx2].cells_box.size().dim();

			//use continuity of Js and S across interface unless the interface is N-F type (normal metal - ferromagnetic) and the spin mixing conductance is not zero (i.e. continuous method disabled).

			//Is it an N-F contact?
			if ((pTransport[idx_pri]->pMeshCUDA->MComputation_Enabled() && !pTransport[idx_sec]->pMeshCUDA->Magnetisation_Enabled()) ||
				(pTransport[idx_sec]->pMeshCUDA->MComputation_Enabled() && !pTransport[idx_pri]->pMeshCUDA->Magnetisation_Enabled())) {

				//Yes we have an N-F contact. Is G interface enabled for this contact ? (remember top mesh sets G interface values)

				//if primary is top then we check GInterface in primary. If primary is bottom then we check GInterface in secondary as it must be top.
				if ((CMBNDcontacts[idx1][idx2].IsPrimaryTop() && pTransport[idx_pri]->pMeshCUDA->GInterface_Enabled()) ||
					(!CMBNDcontacts[idx1][idx2].IsPrimaryTop() && pTransport[idx_sec]->pMeshCUDA->GInterface_Enabled())) {

					//G interface method

					if (pTransport[idx_pri]->pMeshCUDA->MComputation_Enabled()) {

						//interface conductance method with F being the primary mesh

						(*pS[idx_pri])()->set_cmbnd_discontinuous(
							size, *pS[idx_sec],
							(TransportCUDA_Spin_S_Funcs&)pTransport[idx_sec]->poisson_Spin_S, (TransportCUDA_Spin_S_Funcs&)pTransport[idx_pri]->poisson_Spin_S, (STransportCUDA_GInterf_S_NF_Funcs&)gInterf_S_NF,
							(CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);
					}
					else {

						//interface conductance method with N being the primary mesh

						(*pS[idx_pri])()->set_cmbnd_discontinuous(
							size, *pS[idx_sec],
							(TransportCUDA_Spin_S_Funcs&)pTransport[idx_sec]->poisson_Spin_S, (TransportCUDA_Spin_S_Funcs&)pTransport[idx_pri]->poisson_Spin_S, (STransportCUDA_GInterf_S_FN_Funcs&)gInterf_S_FN,
							(CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);
					}

					//next contact
					continue;
				}
			}

			//continuous interface method - the G interface method check above didn't pass so this is what we have left

			(*pS[idx_pri])()->set_cmbnd_continuous(
				size, *pS[idx_sec],
				(TransportCUDA_Spin_S_Funcs&)pTransport[idx_sec]->poisson_Spin_S, (TransportCUDA_Spin_S_Funcs&)pTransport[idx_pri]->poisson_Spin_S,
				(CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);
		}
	}
}

#endif

#endif