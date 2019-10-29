#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_HEAT

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"
#include "HeatCUDA_CMBND.h"

class MeshCUDA;
class Mesh;
class SuperMesh;

class Heat;
class SHeatCUDA;

class HeatCUDA :
	public ModulesCUDA
{

	friend Heat;
	friend SHeatCUDA;

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA * pMeshCUDA;

	//pointer to cpuversion of MeshCUDA : *pMesh holds pMeshCUDA
	Mesh* pMesh;

	SuperMesh* pSMesh;

	Heat* pHeat;

	//evaluate heat equation and store result here. After this is done advance time for temperature based on values stored here.
	cu_arr<cuBReal> heatEq_RHS;

	//HeatCUDA_CMBND holds methods used in setting cmbnd conditions.
	//Pass the managed HeatCUDA_CMBND object to set_cmbnd_continuous in Temp
	cu_obj<HeatCUDA_CMBND> temp_cmbnd_funcs;

private:

	//-------------------Calculation Methods

	void IterateHeatEquation(cuBReal dT);

public:

	HeatCUDA(Mesh* pMesh_, SuperMesh* pSMesh_, Heat* pHeat_);
	~HeatCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);

	void UpdateField(void);

	//-------------------Setters

	//set Temp uniformly to base temperature
	void SetBaseTemperature(cuBReal Temperature);

	//set Temp non-uniformly as specified through the cT mesh parameter
	void SetBaseTemperature_Nonuniform(cuBReal Temperature);
};

#else

class HeatCUDA
{
};

#endif

#endif