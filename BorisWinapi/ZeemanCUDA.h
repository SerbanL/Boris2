#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_ZEEMAN

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class Zeeman;
class MeshCUDA;

class ZeemanCUDA :
	public ModulesCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

	//Applied field
	cu_obj<cuReal3> Ha;

public:

	ZeemanCUDA(MeshCUDA* pMeshCUDA_, Zeeman* pHolderModule);
	~ZeemanCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);

	void UpdateField(void);

	//-------------------

	void SetField(cuReal3 Hxyz);
	cuReal3 GetField(void) { return Ha.to_cpu(); }
};

#else

class ZeemanCUDA
{
};

#endif

#endif
