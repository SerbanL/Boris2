#include "stdafx.h"
#include "Atom_Demag_NCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_DEMAG_N) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"

Atom_Demag_NCUDA::Atom_Demag_NCUDA(Atom_MeshCUDA* paMeshCUDA_) : 
	ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_Demag_NCUDA::~Atom_Demag_NCUDA()
{}

BError Atom_Demag_NCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_Demag_NCUDA));

	initialized = true;

	return error;
}

BError Atom_Demag_NCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_Demag_NCUDA));

	Uninitialize();

	Initialize();

	return error;
}

#endif

#endif
