#include "stdafx.h"
#include "Atom_iDMExchangeCUDA.h"
#include "DataDefs.h"

#if defined(MODULE_COMPILATION_IDMEXCHANGE) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"

#if COMPILECUDA == 1

#include "BorisLib.h"

Atom_iDMExchangeCUDA::Atom_iDMExchangeCUDA(Atom_MeshCUDA* paMeshCUDA_) :
	ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_iDMExchangeCUDA::~Atom_iDMExchangeCUDA()
{}

BError Atom_iDMExchangeCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_iDMExchangeCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect, 
		(MOD_)paMeshCUDA->Get_ActualModule_Heff_Display() == MOD_IDMEXCHANGE || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_EXCH),
		(MOD_)paMeshCUDA->Get_ActualModule_Heff_Display() == MOD_IDMEXCHANGE || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_EXCH));
	if (!error)	initialized = true;

	return error;
}

BError Atom_iDMExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_DMExchangeCUDA));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 Atom_iDMExchangeCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(paMeshCUDA->M1, avRect);
}

#endif

#endif

