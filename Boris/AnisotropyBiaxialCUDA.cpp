#include "stdafx.h"
#include "AnisotropyBiaxialCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_ANIBI

#include "MeshCUDA.h"
#include "MeshDefs.h"
#include "DataDefs.h"

//--------------- UNIAXIAL

Anisotropy_BiaxialCUDA::Anisotropy_BiaxialCUDA(MeshCUDA* pMeshCUDA_)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
}

Anisotropy_BiaxialCUDA::~Anisotropy_BiaxialCUDA()
{}

BError Anisotropy_BiaxialCUDA::Initialize(void)
{
	BError error(CLASS_STR(Anisotropy_BiaxialCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect,
		(MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_ANIBI || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS),
		(MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_ANIBI || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS),
		pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	return error;
}

BError Anisotropy_BiaxialCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Anisotropy_BiaxialCUDA));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 Anisotropy_BiaxialCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(pMeshCUDA->M, avRect);
}

#endif

#endif