#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#if defined(MODULE_HEAT) && ATOMISTIC == 1

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"
#include "HeatBaseCUDA.h"
#include "HeatCUDA_CMBND.h"

#include "Heat_Defs.h"

class Atom_MeshCUDA;
class Atom_Mesh;
class SuperMesh;

class Atom_Heat;
class SHeatCUDA;

class Atom_HeatCUDA :
	public ModulesCUDA,
	public HeatBaseCUDA
{

	friend Atom_Heat;
	friend SHeatCUDA;

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	Atom_MeshCUDA * paMeshCUDA;

	//pointer to cpuversion of MeshCUDA : *pMesh holds pMeshCUDA
	Atom_Mesh* paMesh;

	SuperMesh* pSMesh;

	Atom_Heat* paHeat;

private:

	//-------------------Calculation Methods

	//1-temperature model
	void IterateHeatEquation_1TM(cuBReal dT);

	//2-temperature model
	void IterateHeatEquation_2TM(cuBReal dT);

public:

	Atom_HeatCUDA(Atom_Mesh* paMesh_, SuperMesh* pSMesh_, Atom_Heat* paHeat_);
	~Atom_HeatCUDA();

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

class Atom_HeatCUDA
{
};

#endif

#endif
