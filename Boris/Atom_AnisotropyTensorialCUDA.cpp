#include "stdafx.h"
#include "Atom_AnisotropyTensorialCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ANITENS) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"
#include "DataDefs.h"

//--------------- UNIAXIAL

Atom_Anisotropy_TensorialCUDA::Atom_Anisotropy_TensorialCUDA(Atom_MeshCUDA* paMeshCUDA_)
	: ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_Anisotropy_TensorialCUDA::~Atom_Anisotropy_TensorialCUDA()
{}

BError Atom_Anisotropy_TensorialCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_Anisotropy_TensorialCUDA));

	if (!initialized) {

		if (!paMeshCUDA->Kt()->assign(cuSZ3(paMeshCUDA->get_tensorial_anisotropy().size(), 1, 1), cuReal4())) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		if (!paMeshCUDA->Kt()->copy_from_vector(paMeshCUDA->get_tensorial_anisotropy())) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		initialized = true;
	}

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect,
		(MOD_)paMeshCUDA->Get_Module_Heff_Display() == MOD_ANITENS || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS),
		(MOD_)paMeshCUDA->Get_Module_Energy_Display() == MOD_ANITENS || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_ANIS));
	if (error)	initialized = false;

	return error;
}

BError Atom_Anisotropy_TensorialCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Anisotropy_TensorialCUDA));

	Uninitialize();

	return error;
}

//-------------------Torque methods

cuReal3 Atom_Anisotropy_TensorialCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(paMeshCUDA->M1, avRect);
}

#endif

#endif