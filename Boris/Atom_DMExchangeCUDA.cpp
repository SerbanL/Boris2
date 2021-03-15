#include "stdafx.h"
#include "Atom_DMExchangeCUDA.h"
#include "DataDefs.h"

#if defined(MODULE_COMPILATION_DMEXCHANGE) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"

#if COMPILECUDA == 1

#include "BorisLib.h"

Atom_DMExchangeCUDA::Atom_DMExchangeCUDA(Atom_MeshCUDA* paMeshCUDA_) :
	ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_DMExchangeCUDA::~Atom_DMExchangeCUDA()
{}

BError Atom_DMExchangeCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_DMExchangeCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect, 
		(MOD_)paMeshCUDA->Get_ActualModule_Heff_Display() == MOD_DMEXCHANGE || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_EXCH),
		(MOD_)paMeshCUDA->Get_ActualModule_Heff_Display() == MOD_DMEXCHANGE || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_EXCH));
	if (!error)	initialized = true;

	return error;
}

BError Atom_DMExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_DMExchangeCUDA));

	Uninitialize();

	return error;
}

#endif

#endif

