#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_OERSTED

#include "ModulesCUDA.h"

#include "ConvolutionCUDA.h"
#include "OerstedKernelCUDA.h"

class SuperMesh;
class Oersted;

class OerstedCUDA :
	public ModulesCUDA,
	public ConvolutionCUDA<OerstedCUDA, OerstedKernelCUDA>
{

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	SuperMesh * pSMesh;

	//the SDemag module (cpu version) holding this CUDA version
	Oersted* pOersted;

	//super-mesh magnetization values used for computing demag field on the super-mesh
	cu_obj<cuVEC<cuReal3>> sm_Vals;

public:

	OerstedCUDA(SuperMesh* pSMesh_, Oersted* pOersted_);
	~OerstedCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);

	//-------------------

	cu_obj<cuVEC<cuReal3>>& GetOerstedField(void) { return sm_Vals; }
};

#else

class OerstedCUDA
{
};

#endif

#endif


