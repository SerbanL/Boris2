#include "stdafx.h"
#include "Atom_SurfExchangeCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_SURFEXCHANGE) && ATOMISTIC == 1

#include "Atom_SurfExchange.h"
#include "Atom_Mesh.h"
#include "Atom_MeshCUDA.h"
#include "DataDefs.h"

Atom_SurfExchangeCUDA::Atom_SurfExchangeCUDA(Atom_MeshCUDA* paMeshCUDA_, Atom_SurfExchange* paSurfExch_)
	: ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
	paSurfExch = paSurfExch_;
}

Atom_SurfExchangeCUDA::~Atom_SurfExchangeCUDA()
{}

BError Atom_SurfExchangeCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_SurfExchangeCUDA));

	//clear cu_arrs then rebuild them from information in SurfExchange module
	paMesh_Bot.clear();
	paMesh_Top.clear();

	//make sure information in SurfExchange module is up to date
	error = paSurfExch->Initialize();

	if (!error) {

		for (int idx = 0; idx < paSurfExch->paMesh_Bot.size(); idx++) {

			paMesh_Bot.push_back(paSurfExch->paMesh_Bot[idx]->paMeshCUDA->cuaMesh.get_managed_object());
		}

		for (int idx = 0; idx < paSurfExch->paMesh_Top.size(); idx++) {

			paMesh_Top.push_back(paSurfExch->paMesh_Top[idx]->paMeshCUDA->cuaMesh.get_managed_object());
		}

		//copy number of coupled cells to gpu memory : SurfExchange has been initialized so the count is correct
		coupled_cells.from_cpu(paSurfExch->coupled_cells);

		initialized = true;
	}

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect,
		(MOD_)paMeshCUDA->Get_Module_Heff_Display() == MOD_SURFEXCHANGE || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_SURFEXCH) || paMeshCUDA->IsOutputDataSet(DATA_T_SURFEXCH),
		(MOD_)paMeshCUDA->Get_Module_Energy_Display() == MOD_SURFEXCHANGE || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_SURFEXCH));
	if (error) initialized = false;

	if (initialized) set_Atom_SurfExchangeCUDA_pointers();

	return error;
}

BError Atom_SurfExchangeCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_SurfExchangeCUDA));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 Atom_SurfExchangeCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(paMeshCUDA->M1, avRect);
}

#endif

#endif

