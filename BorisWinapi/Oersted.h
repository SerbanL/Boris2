#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Mesh;
class FMesh;
class SuperMesh;

#ifdef MODULE_OERSTED

#include "Convolution.h"
#include "OerstedKernel.h"

#if COMPILECUDA == 1
#include "OerstedCUDA.h"
#endif

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Full Demag using pre-calculated demag tensor

class Oersted : 
	public Modules,
	public Convolution<OerstedKernel>,
	public ProgramState<Oersted, tuple<>, tuple<>>
{

#if COMPILECUDA == 1
	friend OerstedCUDA;
#endif

private:

	SuperMesh * pSMesh;

	//super-mesh current density values used for computing Oersted field on the super-mesh
	VEC<DBL3> sm_Vals;

	//don't need to compute Oe field every iteration, only when a significant change in Jc occurs; but do need to compute it initially.
	bool oefield_computed = false;

public:

	Oersted(SuperMesh *pSMesh_);
	~Oersted() {}

	//-------------------Methods associated with saving/loading simulations

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Getters

	VEC<DBL3>& GetOerstedField(void) { return sm_Vals; }

#if COMPILECUDA == 1
	cu_obj<cuVEC<cuReal3>>& GetOerstedFieldCUDA(void) { return reinterpret_cast<OerstedCUDA*>(pModuleCUDA)->GetOerstedField(); }
#endif
};

#else

class Oersted :
	public Modules
{

private:

	//super-mesh current density values used for computing Oersted field on the super-mesh
	VEC<DBL3> sm_Vals;

#if COMPILECUDA == 1
	cu_obj<cuVEC<cuReal3>> sm_Vals_CUDA;
#endif

public:

	Oersted(SuperMesh *pSMesh_) {}
	~Oersted() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Getters

	VEC<DBL3>& GetOerstedField(void) { return sm_Vals; }

#if COMPILECUDA == 1
	cu_obj<cuVEC<cuReal3>>& GetOerstedFieldCUDA(void) { return sm_Vals_CUDA; }
#endif
};

#endif
