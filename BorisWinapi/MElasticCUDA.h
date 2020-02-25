#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_MELASTIC

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class MElastic;
class MeshCUDA;
class Mesh;

class MElasticCUDA :
	public ModulesCUDA
{
	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

	//pointer to cpu version of MeshCUDA
	Mesh *pMesh;

	//pointer to Zeeman holder module
	MElastic* pMElastic;

public:

	MElasticCUDA(Mesh* pMesh_, MElastic* pMElastic_);
	~MElasticCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	void UpdateField(void);

	//-------------------
	
	//copy all required mechanical VECs from their cpu versions
	BError copy_VECs_to_GPU(void);
};

#else

class MElasticCUDA
{
};

#endif

#endif
