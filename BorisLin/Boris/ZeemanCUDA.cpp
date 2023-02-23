#include "stdafx.h"
#include "ZeemanCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_ZEEMAN

#include "Zeeman.h"
#include "MeshCUDA.h"
#include "MeshDefs.h"
#include "DataDefs.h"

#include "SuperMesh.h"

ZeemanCUDA::ZeemanCUDA(MeshCUDA* pMeshCUDA_, Zeeman* pZeeman_) :
	ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;
	pZeeman = pZeeman_;

	//copy over any other data in holder module
	Ha.from_cpu(pZeeman->GetField());

	pMeshCUDA->pHa = &Ha;

	if (pZeeman->H_equation.is_set()) SetFieldEquation(pZeeman->H_equation.get_vector_fspec());
}

ZeemanCUDA::~ZeemanCUDA()
{
	pMeshCUDA->pHa = nullptr;
}

//setup globalField transfer
BError ZeemanCUDA::InitializeGlobalField(void)
{
	BError error(__FUNCTION__);

	error = pZeeman->InitializeGlobalField();

	if (!error && pZeeman->globalField.linear_size()) {

		if (!globalField()->resize(pMeshCUDA->h, pMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		cu_arr<cuVEC<cuReal3>> pVal_from;
		cu_arr<cuVEC<cuReal3>> pVal_to;
		
		pVal_from.push_back((cuVEC<cuReal3>*&)pZeeman->pSMesh->GetGlobalFieldCUDA().get_managed_object());

		//Now copy mesh transfer object to cuda version
		if (!globalField()->copy_transfer_info(pVal_from, pVal_to, pZeeman->globalField.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		globalField()->transfer_in(pZeeman->globalField.linear_size(), pZeeman->globalField.size_transfer_in());
	}
	else globalField()->clear();

	return error;
}

BError ZeemanCUDA::Initialize(void)
{
	BError error(CLASS_STR(ZeemanCUDA));

	//If using Havec make sure size and resolution matches M
	if (Havec()->size_cpu().dim()) {
		if (!Havec()->resize((cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect)) {
			
			Havec()->clear();
			return error(BERROR_OUTOFGPUMEMORY_NCRIT);
			initialized = false;
		}
	}

	//if using global field, then initialize mesh transfer if needed
	error = InitializeGlobalField();

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)pMeshCUDA->h, (cuRect)pMeshCUDA->meshRect,
		(MOD_)pMeshCUDA->Get_Module_Heff_Display() == MOD_ZEEMAN || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_ZEE),
		(MOD_)pMeshCUDA->Get_Module_Energy_Display() == MOD_ZEEMAN || pMeshCUDA->IsOutputDataSet_withRect(DATA_E_ZEE),
		pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	if (initialized) set_ZeemanCUDA_pointers();

	return error;
}

BError ZeemanCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(ZeemanCUDA));

	//need this when we switch cuda mode
	if (!H_equation.is_set() && pZeeman->H_equation.is_set()) error = SetFieldEquation(pZeeman->H_equation.get_vector_fspec());

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_SMESH_GLOBALFIELD)) {

		Uninitialize();

		//update mesh dimensions in equation constants
		if (H_equation.is_set()) {

			error = SetFieldEquation(pZeeman->H_equation.get_vector_fspec());
		}

		//if global field not set, then also clear it here
		if (!pZeeman->globalField.linear_size()) globalField()->clear();
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
	if (Havec()->size_cpu().dim()) Havec()->clear();

	Ha.from_cpu(Hxyz);

	//Note : this method is called from the CPU version, which makes any changes associated with non-zero Curie temperature and non-zero atomic moment, including in gpu memory, so not necessary to do anything else here.
	//This happens because the MatP parameters are updated, which automatically results in the corresponding MatPCUDA parameters being updated.
}

BError ZeemanCUDA::SetFieldVEC(VEC<DBL3>& Havec_cpu)
{
	BError error(CLASS_STR(ZeemanCUDA));

	if (!Havec()->set_from_cpuvec(Havec_cpu)) error_on_create(BERROR_OUTOFGPUMEMORY_NCRIT);

	return error;
}

//-------------------Torque methods

cuReal3 ZeemanCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(pMeshCUDA->M, avRect);
}

#endif

#endif