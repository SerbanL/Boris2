#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_MELASTIC

#include <vector>

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

#include "MeshParamsCUDA.h"

class SuperMesh;
class SMElastic;
class MElasticCUDA;

class SMElasticCUDA :
	public ModulesCUDA
{
	friend SMElastic;

private:

	//pointer to cpu version of this super-mesh module
	SMElastic * pSMElastic;

	SuperMesh* pSMesh;

	//---------------------- CMBND data

	//CMBND contacts for all contacting elasticity meshes - these are ordered by first vector index; for each mesh there could be multiple contacting meshes and these are ordered by second vector index
	//we need a cu_obj managed copy of CMBNDcontacts from SMElastic so we can pass it to cuda kernels efficiently
	std::vector< std::vector< cu_obj<CMBNDInfoCUDA> > > CMBNDcontactsCUDA;
	//...and we also need a cpu-memory version, even though we can access it using pSMElastic - the problem is, we need the cpu data in .cu files where we cannot define SMElastic (as nvcc will then attempt to compile BorisLib)
	std::vector< std::vector<CMBNDInfoCUDA> > CMBNDcontacts;

	//list of all elastodynamics modules (same ordering as first vector in CMBNDcontacts)
	std::vector<MElasticCUDA*> pMElastic;

	//vector of pointers to all u_disp - need this to set cmbnd flags (same ordering as first vector in CMBNDcontacts)
	std::vector<cu_obj<cuVEC_VC<cuReal3>>*> pu_disp;

public:

	SMElasticCUDA(SuperMesh* pSMesh_, SMElastic* pSMElastic_);
	~SMElasticCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);
};

#else

class SMElasticCUDA
{
};

#endif

#endif
