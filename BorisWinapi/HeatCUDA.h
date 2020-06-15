#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_HEAT

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"
#include "HeatBaseCUDA.h"
#include "HeatCUDA_CMBND.h"

#include "Heat_Defs.h"

class MeshCUDA;
class Mesh;
class SuperMesh;

class Heat;
class SHeatCUDA;

class HeatCUDA :
	public ModulesCUDA,
	public HeatBaseCUDA
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

private:

	//-------------------Calculation Methods

	//1-temperature model
	void IterateHeatEquation_1TM(cuBReal dT);

	//2-temperature model
	void IterateHeatEquation_2TM(cuBReal dT);

public:

	HeatCUDA(Mesh* pMesh_, SuperMesh* pSMesh_, Heat* pHeat_);
	~HeatCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

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