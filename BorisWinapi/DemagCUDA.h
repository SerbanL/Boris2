#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_DEMAG

#include "ModulesCUDA.h"

#include "ConvolutionCUDA.h"
#include "DemagKernelCUDA.h"

class MeshCUDA;
class Demag;

class DemagCUDA :
	public ModulesCUDA,	
	public ConvolutionCUDA<DemagKernelCUDA>
{

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

	//point to cpu version of this module
	Demag *pDemag;

	//The demag field computed separately : at certain steps in the ODE evaluation method we don't need to recalculate the demag field but can use a previous evaluation with an acceptable impact on the numerical error.
	//This mode needs to be enabled by the user, and can be much faster than the default mode. The default mode is to re-evaluate the demag field at every step.
	cu_obj<cuVEC<cuReal3>> Hdemag;

	//when using the evaluation speedup method we must ensure we have a previous Hdemag evaluation available
	bool Hdemag_calculated = false;

public:

	DemagCUDA(MeshCUDA* pMeshCUDA_, Demag *pDemag_);
	~DemagCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);
};

#else

class DemagCUDA
{
};

#endif

#endif


