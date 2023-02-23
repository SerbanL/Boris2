#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_VIDMEXCHANGE

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"
#include "ExchangeBaseCUDA.h"

class MeshCUDA;
class viDMExchange;

class viDMExchangeCUDA :
	public ModulesCUDA,
	public ExchangeBaseCUDA
{

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

public:

	viDMExchangeCUDA(MeshCUDA* pMeshCUDA_, viDMExchange* pviDMExchange);
	~viDMExchangeCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);

	//-------------------

	//calculate exchange field at coupled cells in this mesh; accumulate energy density contribution in energy : not used in this module
	void CalculateExchangeCoupling(cu_obj<cuBReal>& energy) {}

	//-------------------Torque methods

	cuReal3 GetTorque(cuRect avRect);

};

#else

class viDMExchangeCUDA
{
};

#endif

#endif