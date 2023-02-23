#include "stdafx.h"
#include "Atom_SOTFieldCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_SOTFIELD) && ATOMISTIC == 1

#include "Atom_SOTField.h"
#include "Atom_MeshCUDA.h"
#include "MeshDefs.h"

Atom_SOTFieldCUDA::Atom_SOTFieldCUDA(Atom_MeshCUDA* paMeshCUDA_, Atom_SOTField* pHolderModule)
	: ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_SOTFieldCUDA::~Atom_SOTFieldCUDA()
{}

BError Atom_SOTFieldCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_SOTFieldCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect, 
		(MOD_)paMeshCUDA->Get_Module_Heff_Display() == MOD_SOTFIELD, 
		(MOD_)paMeshCUDA->Get_Module_Energy_Display() == MOD_SOTFIELD);
	if (!error)	initialized = true;

	//no energy density contribution here
	ZeroEnergy();

	return error;
}

BError Atom_SOTFieldCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_SOTFieldCUDA));

	Uninitialize();

	return error;
}

#endif

#endif
