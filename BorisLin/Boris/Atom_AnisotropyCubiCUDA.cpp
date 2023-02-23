#include "stdafx.h"
#include "Atom_AnisotropyCubiCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ANICUBI) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"
#include "DataDefs.h"

//--------------- UNIAXIAL

Atom_Anisotropy_CubiCUDA::Atom_Anisotropy_CubiCUDA(Atom_MeshCUDA* paMeshCUDA_)
	: ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_Anisotropy_CubiCUDA::~Atom_Anisotropy_CubiCUDA()
{}

BError Atom_Anisotropy_CubiCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_Anisotropy_CubiCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect, 
		(MOD_)paMeshCUDA->Get_Module_Heff_Display() == MOD_ANICUBI || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS),
		(MOD_)paMeshCUDA->Get_Module_Energy_Display() == MOD_ANICUBI || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS));
	if (!error)	initialized = true;

	return error;
}

BError Atom_Anisotropy_CubiCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Anisotropy_CubiCUDA));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 Atom_Anisotropy_CubiCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(paMeshCUDA->M1, avRect);
}

#endif

#endif