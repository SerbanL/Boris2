#include "stdafx.h"
#include "Atom_ZeemanCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ZEEMAN) && ATOMISTIC == 1

#include "Atom_Zeeman.h"
#include "Atom_MeshCUDA.h"
#include "DataDefs.h"

#include "SuperMesh.h"

Atom_ZeemanCUDA::Atom_ZeemanCUDA(Atom_MeshCUDA* paMeshCUDA_, Atom_Zeeman* paZeeman_) :
	ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
	paZeeman = paZeeman_;

	//copy over any other data in holder module
	Ha.from_cpu(paZeeman->GetField());

	paMeshCUDA->pHa = &Ha;

	if (paZeeman->H_equation.is_set()) SetFieldEquation(paZeeman->H_equation.get_vector_fspec());
}

Atom_ZeemanCUDA::~Atom_ZeemanCUDA()
{
	paMeshCUDA->pHa = nullptr;
}

//setup globalField transfer
BError Atom_ZeemanCUDA::InitializeGlobalField(void)
{
	BError error(__FUNCTION__);

	error = paZeeman->InitializeGlobalField();

	if (!error && paZeeman->globalField.linear_size()) {

		if (!globalField()->resize(paMeshCUDA->h, paMeshCUDA->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		cu_arr<cuVEC<cuReal3>> pVal_from;
		cu_arr<cuVEC<cuReal3>> pVal_to;

		pVal_from.push_back((cuVEC<cuReal3>*&)paZeeman->pSMesh->GetGlobalFieldCUDA().get_managed_object());

		//Now copy mesh transfer object to cuda version
		if (!globalField()->copy_transfer_info(pVal_from, pVal_to, paZeeman->globalField.get_transfer())) return error(BERROR_OUTOFGPUMEMORY_CRIT);

		globalField()->transfer_in(paZeeman->globalField.linear_size(), paZeeman->globalField.size_transfer_in());
	}
	else globalField()->clear();

	return error;
}

BError Atom_ZeemanCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_ZeemanCUDA));

	//If using Havec make sure size and resolution matches M1
	if (Havec()->size_cpu().dim()) {
		if (!Havec()->resize((cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect)) {

			Havec()->clear();
			return error(BERROR_OUTOFGPUMEMORY_NCRIT);
			initialized = false;
		}
	}

	//if using global field, then initialize mesh transfer if needed
	error = InitializeGlobalField();

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect,
		(MOD_)paMeshCUDA->Get_Module_Heff_Display() == MOD_ZEEMAN || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_ZEE),
		(MOD_)paMeshCUDA->Get_Module_Energy_Display() == MOD_ZEEMAN || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_ZEE));
	if (!error)	initialized = true;

	if (initialized) set_Atom_ZeemanCUDA_pointers();

	return error;
}

BError Atom_ZeemanCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_ZeemanCUDA));

	//need this when we switch cuda mode
	if (!H_equation.is_set() && paZeeman->H_equation.is_set()) error = SetFieldEquation(paZeeman->H_equation.get_vector_fspec());

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_SMESH_GLOBALFIELD)) {

		Uninitialize();

		//update mesh dimensions in equation constants
		if (H_equation.is_set()) {

			error = SetFieldEquation(paZeeman->H_equation.get_vector_fspec());
		}

		//if global field not set, then also clear it here
		if (!paZeeman->globalField.linear_size()) globalField()->clear();
	}

	return error;
}

void Atom_ZeemanCUDA::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		if (H_equation.is_set()) {

			SetFieldEquation(paZeeman->H_equation.get_vector_fspec());
		}
	}
	else if (cfgMessage == UPDATECONFIG_TEQUATION_CLEAR) {

		if (H_equation.is_set()) H_equation.clear();
	}
}

void Atom_ZeemanCUDA::SetField(cuReal3 Hxyz)
{
	if (H_equation.is_set()) H_equation.clear();
	if (Havec()->size_cpu().dim()) Havec()->clear();

	Ha.from_cpu(Hxyz);
}

BError Atom_ZeemanCUDA::SetFieldVEC(VEC<DBL3>& Havec_cpu)
{
	BError error(CLASS_STR(Atom_ZeemanCUDA));

	if (!Havec()->set_from_cpuvec(Havec_cpu)) error_on_create(BERROR_OUTOFGPUMEMORY_NCRIT);

	return error;
}

//-------------------Torque methods

cuReal3 Atom_ZeemanCUDA::GetTorque(cuRect avRect)
{
	return CalculateTorque(paMeshCUDA->M1, avRect);
}

#endif

#endif