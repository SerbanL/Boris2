#include "stdafx.h"
#include "Atom_ZeemanCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ZEEMAN) && ATOMISTIC == 1

#include "Atom_Zeeman.h"
#include "Atom_MeshCUDA.h"
#include "DataDefs.h"

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

BError Atom_ZeemanCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_ZeemanCUDA));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		(cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect, 
		(MOD_)paMeshCUDA->Get_Module_Heff_Display() == MOD_ZEEMAN || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_ZEE),
		(MOD_)paMeshCUDA->Get_Module_Energy_Display() == MOD_ZEEMAN || paMeshCUDA->IsOutputDataSet_withRect(DATA_E_ZEE));
	if (!error)	initialized = true;

	//If using Havec make sure size and resolution matches M1
	if (Havec()->size_cpu().dim()) {
		if (!Havec()->resize((cuReal3)paMeshCUDA->h, (cuRect)paMeshCUDA->meshRect)) {

			Havec()->clear();
			return error(BERROR_OUTOFGPUMEMORY_NCRIT);
			initialized = false;
		}
	}

	if (initialized) set_Atom_ZeemanCUDA_pointers();

	return error;
}

BError Atom_ZeemanCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_ZeemanCUDA));

	//need this when we switch cuda mode
	if (!H_equation.is_set() && paZeeman->H_equation.is_set()) error = SetFieldEquation(paZeeman->H_equation.get_vector_fspec());

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

		Uninitialize();

		//update mesh dimensions in equation constants
		if (H_equation.is_set()) {

			error = SetFieldEquation(paZeeman->H_equation.get_vector_fspec());
		}
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