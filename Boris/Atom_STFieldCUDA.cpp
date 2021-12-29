#include "stdafx.h"
#include "Atom_STFieldCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_STFIELD) && ATOMISTIC == 1

#include "Atom_STField.h"
#include "Atom_MeshCUDA.h"
#include "MeshDefs.h"

Atom_STFieldCUDA::Atom_STFieldCUDA(Atom_MeshCUDA* paMeshCUDA_, Atom_STField* pHolderModule)
	: ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_STFieldCUDA::~Atom_STFieldCUDA()
{}

BError Atom_STFieldCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_STFieldCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect, 
		(MOD_)paMeshCUDA->Get_Module_Heff_Display() == MOD_STFIELD, 
		(MOD_)paMeshCUDA->Get_Module_Energy_Display() == MOD_STFIELD);
	if (!error)	initialized = true;

	//no energy density contribution here
	ZeroEnergy();

	return error;
}

BError Atom_STFieldCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_STFieldCUDA));

	Uninitialize();

	return error;
}

#endif

#endif
