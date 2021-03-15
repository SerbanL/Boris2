#include "stdafx.h"
#include "ZeemanCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_ZEEMAN

#include "Zeeman.h"
#include "MeshCUDA.h"
#include "MeshDefs.h"
#include "DataDefs.h"

ZeemanCUDA::ZeemanCUDA(MeshCUDA* pMeshCUDA_, Zeeman* pZeeman_) :
	ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
	pZeeman = pZeeman_;

	//copy over any other data in holder module
	Ha.from_cpu(pZeeman->GetField());

	if (pZeeman->H_equation.is_set()) SetFieldEquation(pZeeman->H_equation.get_vector_fspec());
}

ZeemanCUDA::~ZeemanCUDA()
{}

BError ZeemanCUDA::Initialize(void)
{
	BError error(CLASS_STR(ZeemanCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect, 
		(MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_ZEEMAN || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_ZEE),
		(MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_ZEEMAN || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_ZEE),
		pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	return error;
}

BError ZeemanCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(ZeemanCUDA));

	//need this when we switch cuda mode
	if (!H_equation.is_set() && pZeeman->H_equation.is_set()) error = SetFieldEquation(pZeeman->H_equation.get_vector_fspec());

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

		Uninitialize();

		//update mesh dimensions in equation constants
		if (H_equation.is_set()) {

			error = SetFieldEquation(pZeeman->H_equation.get_vector_fspec());
		}

		Initialize();
	}

	return error;
}

void ZeemanCUDA::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		if (H_equation.is_set()) {

			SetFieldEquation(pZeeman->H_equation.get_vector_fspec());
		}
	}
	else if (cfgMessage == UPDATECONFIG_TEQUATION_CLEAR) {

		if (H_equation.is_set()) H_equation.clear();
	}
}

void ZeemanCUDA::SetField(cuReal3 Hxyz)
{
	if (H_equation.is_set()) H_equation.clear();

	Ha.from_cpu(Hxyz);

	//Note : this method is called from the CPU version, which makes any changes associated with non-zero Curie temperature and non-zero atomic moment, including in gpu memory, so not necessary to do anything else here.
	//This happens because the MatP parameters are updated, which automatically results in the corresponding MatPCUDA parameters being updated.
}

#endif

#endif