#include "stdafx.h"
#include "ExchangeBaseCUDA.h"
#include "ExchangeBase.h"
#include "Mesh_Ferromagnetic.h"
#include "Mesh_FerromagneticCUDA.h"

#if COMPILECUDA == 1

ExchangeBaseCUDA::ExchangeBaseCUDA(FMeshCUDA* pMeshCUDA_, ExchangeBase* pExchBase_)
{
	pMeshCUDA = pMeshCUDA_;
	pExchBase = pExchBase_;
}

BError ExchangeBaseCUDA::Initialize(void)
{
	BError error(CLASS_STR(ExchangeBaseCUDA));

	//initialize the cpu version of ExchangeBase
	//doing this will set cmbnd flags in M VECs which we can copy to gpu version here
	pExchBase->Initialize();

	//clear everything then rebuild
	CMBNDcontactsCUDA.clear();
	CMBNDcontacts.clear();
	pContactingManagedMeshes.clear();

	//copy managed cuda meshes to pContactingManagedMeshes using the pFMeshes initialized in ExchangeBase
	for (int idx = 0; idx < pExchBase->pFMeshes.size(); idx++) {

		pContactingManagedMeshes.push_back(&(pExchBase->pFMeshes[idx]->pMeshCUDA->cuMesh));
	}

	//set cmbnd flags
	if (!(pMeshCUDA->M)()->copyflags_from_cpuvec(pExchBase->pMesh->M)) error(BERROR_GPUERROR_CRIT);

	//copy CMBNDInfo for contacts
	for (int idx_contact = 0; idx_contact < pExchBase->CMBNDcontacts.size(); idx_contact++) {

		cu_obj<CMBNDInfoCUDA> contact;

		contact()->copy_from_CMBNDInfo<CMBNDInfo>(pExchBase->CMBNDcontacts[idx_contact]);

		CMBNDcontactsCUDA.push_back(contact);
		CMBNDcontacts.push_back(pExchBase->CMBNDcontacts[idx_contact]);
	}

	return error;
}

#endif
