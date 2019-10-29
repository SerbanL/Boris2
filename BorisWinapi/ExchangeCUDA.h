#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_EXCHANGE

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"
#include "ExchangeBaseCUDA.h"

class MeshCUDA;
class Exch_6ngbr_Neu;

class Exch_6ngbr_NeuCUDA :
	public ModulesCUDA,
	public ExchangeBaseCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

public:

	Exch_6ngbr_NeuCUDA(MeshCUDA* pMeshCUDA_, Exch_6ngbr_Neu* pExch_6ngbr_Neu);
	~Exch_6ngbr_NeuCUDA();

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

class Exch_6ngbr_NeuCUDA
{
};

#endif

#endif

