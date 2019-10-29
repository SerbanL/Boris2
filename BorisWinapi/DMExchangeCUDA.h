#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_DMEXCHANGE

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"
#include "ExchangeBaseCUDA.h"

class MeshCUDA;
class DMExchange;

class DMExchangeCUDA :
	public ModulesCUDA,
	public ExchangeBaseCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

public:

	DMExchangeCUDA(MeshCUDA* pMeshCUDA_, DMExchange* pDMExchange);
	~DMExchangeCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);

	void UpdateField(void);

	//-------------------

	//calculate exchange field at coupled cells in this mesh; accumulate energy density contribution in energy
	void CalculateExchangeCoupling(cu_obj<cuBReal>& energy);
};

#else

class DMExchangeCUDA
{
};

#endif

#endif

