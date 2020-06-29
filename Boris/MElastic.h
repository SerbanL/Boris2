#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Mesh;

#ifdef MODULE_COMPILATION_MELASTIC

class SuperMesh;

class MElastic :
	public Modules,
	public ProgramState<MElastic, tuple<DBL3>, tuple<>>
{

#if COMPILECUDA == 1
	friend class MElasticCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Mesh *pMesh;

	//Applied uniform stress vector (Cartesian coordinates)
	DBL3 Tsig;

	//pointer to supermesh
	SuperMesh* pSMesh;

private:

	//----------------------------------------------- Computational Helpers

	//compute strain form mechanical displacement
	void Calculate_Strain(void);

public:

	MElastic(Mesh *pMesh_);
	~MElastic();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//------------------- STRAIN GENERATION without SOLVER

	//Set/Get uniform stress
	void SetUniformStress(DBL3 Tsig_xyz);
	DBL3 GetUniformStress(void) { return Tsig; }

	//Set strain or displacement by loading from OVF2 files
	
	//Displacement : from this calculate strain tensor
	BError Load_Displacement_OVF2(string fileName);

	//Tensor : displacement is not calculated, as we only need the strain tensor to obtain the effective fields at runtime
	//For a strain tensor for cubic crystals, we only have 6 elements : diagonal and off-diagonal (symmetric).
	//These are specified using 2 separate OVF2 files containing vector data:
	//one for the xx, yy, zz elements (diagonal)
	//the other for the yz, xz, xy elements (off-diagonal, in this order)
	BError Load_Strain_OVF2(string fileName_Diag, string fileName_ODiag);
};

#else

class MElastic :
	public Modules
{

private:

private:

	//----------------------------------------------- Computational Helpers

	//compute strain form mechanical displacement
	void Calculate_Strain(void) {}

public:

	MElastic(Mesh *pMesh_) {}
	~MElastic() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------

	//Set/Get uniform stress
	void SetUniformStress(DBL3 Tsig_xyz) {}
	DBL3 GetUniformStress(void) { return DBL3(); }

	BError Load_Displacement_OVF2(string fileName) { return BError(); }
	BError Load_Strain_OVF2(string fileName_Diag, string fileName_ODiag) { return BError(); }
};

#endif