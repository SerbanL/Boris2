#include "stdafx.h"
#include "AnisotropyCubiCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_ANICUBI

#include "MeshCUDA.h"
#include "MeshDefs.h"
#include "DataDefs.h"

//--------------- CUBIC

Anisotropy_CubicCUDA::Anisotropy_CubicCUDA(MeshCUDA* pMeshCUDA_)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
}

Anisotropy_CubicCUDA::~Anisotropy_CubicCUDA()
{}

BError Anisotropy_CubicCUDA::Initialize(void)
{
	BError error(CLASS_STR(Anisotropy_CubicCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect, 
		(MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_ANICUBI || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS),
		(MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_ANICUBI || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS),
		pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error) initialized = true;

	return error;
}

BError Anisotropy_CubicCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Anisotropy_CubicCUDA));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 Anisotropy_CubicCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(pMeshCUDA->M, avRect);
}

#endif

#endif