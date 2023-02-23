#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_VIDMEXCHANGE) && ATOMISTIC == 1

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class Atom_MeshCUDA;

//cannot include BorisLib.h (or even just VEC_VC.h and VEC.h) since it contains C++14 (ProgramState and Introspection in particular) code which nvcc doesn't compile.
//Since this header file inevitably gets included in .cu files (compiled with nvcc) then must forward declare VEC and VEC_VC here instead.
//To avoid incomplete type error then we must only declare pointers to these types in the class definition - otherwise how can the compiler know how much memory to allocate on object creation?
template <typename VType> class VEC;

class Atom_viDMExchangeCUDA :
	public ModulesCUDA
{

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	Atom_MeshCUDA* paMeshCUDA;

public:

	Atom_viDMExchangeCUDA(Atom_MeshCUDA* paMeshCUDA_);
	~Atom_viDMExchangeCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);

	//-------------------Torque methods

	cuReal3 GetTorque(cuRect avRect);
};

#else

class Atom_viDMExchangeCUDA
{
};

#endif

#endif
