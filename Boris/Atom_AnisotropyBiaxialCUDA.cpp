#include "stdafx.h"
#include "Atom_AnisotropyBiaxialCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ANIBI) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"
#include "DataDefs.h"

//--------------- UNIAXIAL

Atom_Anisotropy_BiaxialCUDA::Atom_Anisotropy_BiaxialCUDA(Atom_MeshCUDA* paMeshCUDA_)
	: ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_Anisotropy_BiaxialCUDA::~Atom_Anisotropy_BiaxialCUDA()
{}

BError Atom_Anisotropy_BiaxialCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_Anisotropy_BiaxialCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect,
		(MOD_)paMeshCUDA->Get_Module_Heff_Display() == MOD_ANIBI || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS),
		(MOD_)paMeshCUDA->Get_Module_Energy_Display() == MOD_ANIBI || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS));
	if (!error)	initialized = true;

	return error;
}

BError Atom_Anisotropy_BiaxialCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Anisotropy_BiaxialCUDA));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 Atom_Anisotropy_BiaxialCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(paMeshCUDA->M1, avRect);
}

#endif

#endif