#include "stdafx.h"
#include "Atom_MOpticalCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_MOPTICAL) && ATOMISTIC == 1

#include "Atom_MOptical.h"
#include "Atom_MeshCUDA.h"
#include "DataDefs.h"

Atom_MOpticalCUDA::Atom_MOpticalCUDA(Atom_MeshCUDA* paMeshCUDA_) :
	ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_MOpticalCUDA::~Atom_MOpticalCUDA()
{}

BError Atom_MOpticalCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_MOpticalCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect, 
		(MOD_)paMeshCUDA->Get_Module_Heff_Display() == MOD_MOPTICAL || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_MOPTICAL),
		(MOD_)paMeshCUDA->Get_Module_Energy_Display() == MOD_MOPTICAL || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_MOPTICAL));
	if (!error)	initialized = true;

	return error;
}

BError Atom_MOpticalCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_MOpticalCUDA));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

		Uninitialize();
	}

	return error;
}

void Atom_MOpticalCUDA::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
}

#endif

#endif