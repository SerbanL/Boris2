#include "stdafx.h"
#include "Atom_Demag_NCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_DEMAG_N) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"
#include "DataDefs.h"

Atom_Demag_NCUDA::Atom_Demag_NCUDA(Atom_MeshCUDA* paMeshCUDA_) : 
	ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_Demag_NCUDA::~Atom_Demag_NCUDA()
{}

BError Atom_Demag_NCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_Demag_NCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect, 
		(MOD_)paMeshCUDA->Get_Module_Heff_Display() == MOD_DEMAG_N || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_DEMAG),
		(MOD_)paMeshCUDA->Get_Module_Energy_Display() == MOD_DEMAG_N || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_DEMAG));
	if (!error)	initialized = true;

	return error;
}

BError Atom_Demag_NCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Demag_NCUDA));

	Uninitialize();

	return error;
}

#endif

#endif
