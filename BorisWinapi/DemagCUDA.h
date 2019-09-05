#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_DEMAG

#include "ModulesCUDA.h"

#include "ConvolutionCUDA.h"
#include "DemagKernelCUDA.h"

class FMeshCUDA;
class Demag;

class DemagCUDA :
	public ModulesCUDA,	
	public ConvolutionCUDA<DemagKernelCUDA>
{

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	FMeshCUDA* pMeshCUDA;

	//point to cpu version of this module
	Demag *pDemag;

public:

	DemagCUDA(FMeshCUDA* pMeshCUDA_, Demag *pDemag_);
	~DemagCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	void UpdateField(void);

	//-------------------
};

#else

class DemagCUDA
{
};

#endif

#endif


