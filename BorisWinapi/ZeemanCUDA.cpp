#include "stdafx.h"
#include "ZeemanCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_ZEEMAN

#include "Zeeman.h"
#include "Mesh_FerromagneticCUDA.h"

ZeemanCUDA::ZeemanCUDA(FMeshCUDA* pMeshCUDA_, Zeeman* pHolderModule)
	: ModulesCUDA()
{
	pMeshCUDA = pMeshCUDA_;

	//copy over any other data in holder module
	Ha.from_cpu(pHolderModule->GetField());
}

ZeemanCUDA::~ZeemanCUDA()
{}

BError ZeemanCUDA::Initialize(void)
{
	BError error(CLASS_STR(ZeemanCUDA));

	initialized = true;

	return error;
}

BError ZeemanCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(ZeemanCUDA));

	Uninitialize();

	Initialize();

	return error;
}

void ZeemanCUDA::SetField(cuReal3 Hxyz)
{
	Ha.from_cpu(Hxyz);

	//Note : this method is called from the CPU version, which makes any changes associated with non-zero Curie temperature and non-zero atomic moment, including in gpu memory, so not necessary to do anything else here.
	//This happens because the MatP parameters are updated, which automatically results in the corresponding MatPCUDA parameters being updated.
}

#endif

#endif