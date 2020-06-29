#include "stdafx.h"
#include "Atom_MOpticalCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_MOPTICAL) && ATOMISTIC == 1

#include "Atom_MOptical.h"
#include "Atom_MeshCUDA.h"

Atom_MOpticalCUDA::Atom_MOpticalCUDA(Atom_MeshCUDA* paMeshCUDA_) :
	ModulesCUDA()
{
	paMeshCUDA = paMeshCUDA_;
}

Atom_MOpticalCUDA::~Atom_MOpticalCUDA()
{}

BError Atom_MOpticalCUDA::Initialize(void)
{
	BError error(CLASS_STR(Atom_MOpticalCUDA));

	initialized = true;

	return error;
}

BError Atom_MOpticalCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_MOpticalCUDA));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

		Uninitialize();
		Initialize();
	}

	return error;
}

void Atom_MOpticalCUDA::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
}

#endif

#endif