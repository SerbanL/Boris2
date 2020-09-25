#include "stdafx.h"
#include "Atom_ZeemanCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ZEEMAN) && ATOMISTIC == 1

#include "Atom_Zeeman.h"
#include "Atom_MeshCUDA.h"

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

	initialized = true;

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

		Initialize();
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

	Ha.from_cpu(Hxyz);
}

#endif

#endif