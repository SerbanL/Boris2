#include "stdafx.h"
#include "AnisotropyCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_ANIUNI

#include "MeshCUDA.h"
#include "MeshDefs.h"
#include "DataDefs.h"

//--------------- UNIAXIAL

Anisotropy_UniaxialCUDA::Anisotropy_UniaxialCUDA(MeshCUDA* pMeshCUDA_)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
}

Anisotropy_UniaxialCUDA::~Anisotropy_UniaxialCUDA()
{}

BError Anisotropy_UniaxialCUDA::Initialize(void)
{
	BError error(CLASS_STR(Anisotropy_UniaxialCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect, 
		(MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_ANIUNI || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS),
		(MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_ANIUNI || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS),
		pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	return error;
}

BError Anisotropy_UniaxialCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Anisotropy_UniaxialCUDA));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 Anisotropy_UniaxialCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(pMeshCUDA->M, avRect);
}

#endif

#endif