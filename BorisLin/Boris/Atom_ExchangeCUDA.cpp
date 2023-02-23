#include "stdafx.h"
#include "Atom_ExchangeCUDA.h"
#include "DataDefs.h"

#if defined(MODULE_COMPILATION_EXCHANGE) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"

#if COMPILECUDA == 1

#include "BorisLib.h"

Atom_ExchangeCUDA::Atom_ExchangeCUDA(Atom_MeshCUDA* paMeshCUDA_) :
	ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_ExchangeCUDA::~Atom_ExchangeCUDA()
{}

BError Atom_ExchangeCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_ExchangeCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect, 
		(MOD_)paMeshCUDA->Get_Module_Heff_Display() == MOD_EXCHANGE || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_EXCH),
		(MOD_)paMeshCUDA->Get_Module_Energy_Display() == MOD_EXCHANGE || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_EXCH));
	if (!error)	initialized = true;

	return error;
}

BError Atom_ExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_ExchangeCUDA));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 Atom_ExchangeCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(paMeshCUDA->M1, avRect);
}

#endif

#endif

