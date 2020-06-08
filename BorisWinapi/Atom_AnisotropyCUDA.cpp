#include "stdafx.h"
#include "Atom_AnisotropyCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_ANIUNI) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"

//--------------- UNIAXIAL

Atom_Anisotropy_UniaxialCUDA::Atom_Anisotropy_UniaxialCUDA(Atom_MeshCUDA* paMeshCUDA_)
	: ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_Anisotropy_UniaxialCUDA::~Atom_Anisotropy_UniaxialCUDA()
{}

BError Atom_Anisotropy_UniaxialCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_Anisotropy_UniaxialCUDA));

	initialized = true;

	return error;
}

BError Atom_Anisotropy_UniaxialCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Anisotropy_UniaxialCUDA));

	Uninitialize();

	Initialize();

	return error;
}

#endif

#endif