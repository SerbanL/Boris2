#include "stdafx.h"
#include "Atom_AnisotropyCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ANIUNI) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"
#include "DataDefs.h"

//--------------- UNIAXIAL

Atom_Anisotropy_UniaxialCUDA::Atom_Anisotropy_UniaxialCUDA(Atom_MeshCUDA* paMeshCUDA_)
	: ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_Anisotropy_UniaxialCUDA::~Atom_Anisotropy_UniaxialCUDA()
{}

BError Atom_Anisotropy_UniaxialCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_Anisotropy_UniaxialCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect, 
		(MOD_)paMeshCUDA->Get_Module_Heff_Display() == MOD_ANIUNI || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS),
		(MOD_)paMeshCUDA->Get_Module_Energy_Display() == MOD_ANIUNI || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS));
	if (!error)	initialized = true;

	return error;
}

BError Atom_Anisotropy_UniaxialCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Anisotropy_UniaxialCUDA));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 Atom_Anisotropy_UniaxialCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(paMeshCUDA->M1, avRect);
}

#endif

#endif