#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_HEAT

#include <vector>
#include <tuple>

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

#include "HeatCUDA.h"
#include "MeshParamsCUDA.h"

using namespace std;

class SuperMesh;
class SHeat;

class SHeatCUDA :
	public ModulesCUDA
{

	friend SHeat;

private:

	//pointer to cpu version of this super-mesh module
	SHeat * pSHeat;

	SuperMesh* pSMesh;

	//---------------------- CMBND data

	//CMBND contacts for all contacting heat conduction meshes - these are ordered by first vector index; for each mesh there could be multiple contacting meshes and these are ordered by second vector index
	//CMBNDInfo describes the contact between 2 meshes, allowing calculation of values at cmbnd cells based on continuity of a potential and flux
	//we need a cu_obj managed copy of CMBNDcontacts from SHeat so we can pass it to cuda kernels efficiently
	vector< vector< cu_obj<CMBNDInfoCUDA> > > CMBNDcontactsCUDA;
	//...and we also need a cpu-memory version, even though we can access it using pSHeat - the problem is, we need the cpu data in .cu files where we cannot define SHeat (as nvcc will then attempt to compile BorisLib)
	vector< vector<CMBNDInfoCUDA> > CMBNDcontacts;

	//list of all transport modules in transport meshes (same ordering as first vector in CMBNDcontacts)
	vector<HeatCUDA*> pHeat;

	//vector of pointers to all V - need this to set cmbnd flags (same ordering as first vector in CMBNDcontacts)
	vector<cu_obj<cuVEC_VC<cuReal>>*> pTemp;

private:

	//calculate and set values at composite media boundaries after all other cells have been computed and set
	void set_cmbnd_values(void);

public:

	SHeatCUDA(SuperMesh* pSMesh_, SHeat* pSHeat_);
	~SHeatCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	void UpdateField(void);
};

#else

class SHeatCUDA
{
};

#endif

#endif
