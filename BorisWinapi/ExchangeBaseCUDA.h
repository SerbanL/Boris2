#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"
#include "ErrorHandler.h"

class FMeshCUDA;
class ManagedMeshCUDA;
class ExchangeBase;

class ExchangeBaseCUDA {

private:

	//CMBND contacts between this mesh and other ferromagnetic meshes (we do not require other ferromagnetic meshes to have an exchange module enabled, just this one).
	vector< cu_obj<CMBNDInfoCUDA> > CMBNDcontactsCUDA;
	//...and we also need a cpu-memory version
	vector<CMBNDInfoCUDA> CMBNDcontacts;

	//vector of pointers to all ferromagnetic meshes in managed mesh cuda form
	vector< cu_obj<ManagedMeshCUDA>* > pContactingManagedMeshes;

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module (the primary mesh)
	FMeshCUDA* pMeshCUDA;

	//pointer to cpu version of this (base for exchange-type module)
	ExchangeBase* pExchBase;

protected:

	//this is overloaded by inheriting Exchange-type modules. Need this to be virtual so if for any reason a base pointer is used, the overloaded method is called instead.
	BError Initialize(void);

	//calculate exchange field at coupled cells in this mesh; accumulate energy density contribution in energy
	void CalculateExchangeCoupling(cu_obj<cuReal>& energy);

	//protected constructor - this class should not be instantiated by itself, but only used as a base for an exchange-type module for purposes of code reuse
	ExchangeBaseCUDA(FMeshCUDA* pMeshCUDA_, ExchangeBase* pExchBase_);

	virtual ~ExchangeBaseCUDA() {}
};

#endif

