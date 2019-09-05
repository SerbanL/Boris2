#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_IDMEXCHANGE

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"
#include "ExchangeBaseCUDA.h"

class FMeshCUDA;
class iDMExchange;

class iDMExchangeCUDA :
	public ModulesCUDA,
	public ExchangeBaseCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	FMeshCUDA* pMeshCUDA;

public:

	iDMExchangeCUDA(FMeshCUDA* pMeshCUDA_, iDMExchange* piDMExchange);
	~iDMExchangeCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	void UpdateField(void);

	//-------------------

};

#else

class iDMExchangeCUDA
{
};

#endif

#endif

