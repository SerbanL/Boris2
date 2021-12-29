#include "STransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

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
			if ((pTransport[idx_pri]->Get_STSolveType() == STSOLVE_FERROMAGNETIC && pTransport[idx_sec]->Get_STSolveType() == STSOLVE_NORMALMETAL) ||
				(pTransport[idx_sec]->Get_STSolveType() == STSOLVE_FERROMAGNETIC && pTransport[idx_pri]->Get_STSolveType() == STSOLVE_NORMALMETAL)) {

				//Yes we have an N-F contact. Is G interface enabled for this contact ? (remember top mesh sets G interface values)

				//if primary is top then we check GInterface in primary. If primary is bottom then we check GInterface in secondary as it must be top.
				if ((CMBNDcontacts[idx1][idx2].IsPrimaryTop() && pTransport[idx_pri]->pMeshBaseCUDA->GInterface_Enabled()) ||
					(!CMBNDcontacts[idx1][idx2].IsPrimaryTop() && pTransport[idx_sec]->pMeshBaseCUDA->GInterface_Enabled())) {

					//G interface method
					
					//Micromagnetic to micromagnetic meshes
					//TO DO: atomistic to atomistic / micromagnetic meshes
					if (!pTransport[idx_pri]->pMeshBaseCUDA->is_atomistic() && !pTransport[idx_sec]->pMeshBaseCUDA->is_atomistic()) {

						(*pV[idx_pri])()->set_cmbnd_continuousflux(
							size, *pV[idx_sec],
							(TransportCUDA_Spin_V_Funcs&)dynamic_cast<TransportCUDA*>(pTransport[idx_sec])->poisson_Spin_V, (TransportCUDA_Spin_V_Funcs&)dynamic_cast<TransportCUDA*>(pTransport[idx_pri])->poisson_Spin_V, 
							(STransportCUDA_GInterf_V_Funcs&)gInterf_V,
							(CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);
					}
						
					//next contact
					continue;
				}
			}

			//continuous interface method - the G interface method check above didn't pass so this is what we have left
			
			//Micromagnetic to micromagnetic meshes
			//TO DO: atomistic to atomistic / micromagnetic meshes
			if (!pTransport[idx_pri]->pMeshBaseCUDA->is_atomistic() && !pTransport[idx_sec]->pMeshBaseCUDA->is_atomistic()) {

				(*pV[idx_pri])()->set_cmbnd_continuous(
					size, *pV[idx_sec],
					(TransportCUDA_Spin_V_Funcs&)dynamic_cast<TransportCUDA*>(pTransport[idx_sec])->poisson_Spin_V, (TransportCUDA_Spin_V_Funcs&)dynamic_cast<TransportCUDA*>(pTransport[idx_pri])->poisson_Spin_V,
					(CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);
			}
				
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

			if (pTransport[idx_pri]->Get_STSolveType() == STSOLVE_NONE || pTransport[idx_sec]->Get_STSolveType() == STSOLVE_NONE) continue;

			//use continuity of Js and S across interface unless the interface is N-F type (normal metal - ferromagnetic) and the spin mixing conductance is not zero (i.e. continuous method disabled).

			//Is it an N-F contact?
			if ((pTransport[idx_pri]->Get_STSolveType() == STSOLVE_FERROMAGNETIC && pTransport[idx_sec]->Get_STSolveType() == STSOLVE_NORMALMETAL) ||
				(pTransport[idx_sec]->Get_STSolveType() == STSOLVE_FERROMAGNETIC && pTransport[idx_pri]->Get_STSolveType() == STSOLVE_NORMALMETAL)) {

				//Yes we have an N-F contact. Is G interface enabled for this contact ? (remember top mesh sets G interface values)

				//if primary is top then we check GInterface in primary. If primary is bottom then we check GInterface in secondary as it must be top.
				if ((CMBNDcontacts[idx1][idx2].IsPrimaryTop() && pTransport[idx_pri]->pMeshBaseCUDA->GInterface_Enabled()) ||
					(!CMBNDcontacts[idx1][idx2].IsPrimaryTop() && pTransport[idx_sec]->pMeshBaseCUDA->GInterface_Enabled())) {

					//G interface method

					if (pTransport[idx_pri]->Get_STSolveType() == STSOLVE_FERROMAGNETIC) {

						//interface conductance method with F being the primary mesh
						
						//Micromagnetic to micromagnetic meshes
						//TO DO: atomistic to atomistic / micromagnetic meshes
						if (!pTransport[idx_pri]->pMeshBaseCUDA->is_atomistic() && !pTransport[idx_sec]->pMeshBaseCUDA->is_atomistic()) {

							(*pS[idx_pri])()->set_cmbnd_discontinuous(
								size, *pS[idx_sec],
								(TransportCUDA_Spin_S_Funcs&)dynamic_cast<TransportCUDA*>(pTransport[idx_sec])->poisson_Spin_S, (TransportCUDA_Spin_S_Funcs&)dynamic_cast<TransportCUDA*>(pTransport[idx_pri])->poisson_Spin_S, 
								(STransportCUDA_GInterf_S_NF_Funcs&)gInterf_S_NF,
								(CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);
						}
					}
					else {

						//interface conductance method with N being the primary mesh
						
						//Micromagnetic to micromagnetic meshes
						//TO DO: atomistic to atomistic / micromagnetic meshes
						if (!pTransport[idx_pri]->pMeshBaseCUDA->is_atomistic() && !pTransport[idx_sec]->pMeshBaseCUDA->is_atomistic()) {

							(*pS[idx_pri])()->set_cmbnd_discontinuous(
								size, *pS[idx_sec],
								(TransportCUDA_Spin_S_Funcs&)dynamic_cast<TransportCUDA*>(pTransport[idx_sec])->poisson_Spin_S, (TransportCUDA_Spin_S_Funcs&)dynamic_cast<TransportCUDA*>(pTransport[idx_pri])->poisson_Spin_S, 
								(STransportCUDA_GInterf_S_FN_Funcs&)gInterf_S_FN,
								(CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);
						}
					}

					//next contact
					continue;
				}
			}

			//continuous interface method - the G interface method check above didn't pass so this is what we have left

			//Micromagnetic to micromagnetic meshes
			//TO DO: atomistic to atomistic / micromagnetic meshes
			if (!pTransport[idx_pri]->pMeshBaseCUDA->is_atomistic() && !pTransport[idx_sec]->pMeshBaseCUDA->is_atomistic()) {

				(*pS[idx_pri])()->set_cmbnd_continuous(
					size, *pS[idx_sec],
					(TransportCUDA_Spin_S_Funcs&)dynamic_cast<TransportCUDA*>(pTransport[idx_sec])->poisson_Spin_S, (TransportCUDA_Spin_S_Funcs&)dynamic_cast<TransportCUDA*>(pTransport[idx_pri])->poisson_Spin_S,
					(CMBNDInfoCUDA&)CMBNDcontactsCUDA[idx1][idx2]);
			}
		}
	}
}

#endif

#endif