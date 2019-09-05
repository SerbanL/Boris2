#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_DMEXCHANGE

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"
#include "ExchangeBaseCUDA.h"

class FMeshCUDA;
class DMExchange;

class DMExchangeCUDA :
	public ModulesCUDA,
	public ExchangeBaseCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	FMeshCUDA* pMeshCUDA;

public:

	DMExchangeCUDA(FMeshCUDA* pMeshCUDA_, DMExchange* pDMExchange);
	~DMExchangeCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	void UpdateField(void);

	//-------------------

};

#else

class DMExchangeCUDA
{
};

#endif

#endif

