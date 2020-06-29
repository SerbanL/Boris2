#include "stdafx.h"
#include "Atom_AnisotropyCubiCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ANICUBI) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"

//--------------- UNIAXIAL

Atom_Anisotropy_CubiCUDA::Atom_Anisotropy_CubiCUDA(Atom_MeshCUDA* paMeshCUDA_)
	: ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_Anisotropy_CubiCUDA::~Atom_Anisotropy_CubiCUDA()
{}

BError Atom_Anisotropy_CubiCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_Anisotropy_CubiCUDA));

	initialized = true;

	return error;
}

BError Atom_Anisotropy_CubiCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Anisotropy_CubiCUDA));

	Uninitialize();

	Initialize();

	return error;
}

#endif

#endif