#pragma once

#include "BorisLib.h"
#include "Modules.h"

#if COMPILECUDA == 1
#include "SurfExchangeCUDA.h"
#endif

class Mesh;

#ifdef MODULE_COMPILATION_SURFEXCHANGE

//This SurfExchange module is only used in a ferromagnetic mesh; there is a SurfExchange_AFM module used in a two-sublattice model mesh

class SurfExchange :
	public Modules,
	public ProgramState<SurfExchange, std::tuple<>, std::tuple<>>
{

#if COMPILECUDA == 1
	friend SurfExchangeCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Mesh* pMesh;

	//magnetic meshes in surface exchange coupling with the mesh holding this module, top and bottom
	std::vector<Mesh*> pMesh_Bot, pMesh_Top;

	//number of coupled cells (either top or bottom) in this mesh
	int coupled_cells = 0;

public:

	SurfExchange(Mesh *pMesh_);
	~SurfExchange();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Energy methods

	//FM mesh
	double Get_EnergyChange(int spin_index, DBL3 Mnew);

	//-------------------Torque methods

	DBL3 GetTorque(Rect& avRect);
};

#else

class SurfExchange :
	public Modules
{

private:

public:

	SurfExchange(Mesh *pMesh_) {}
	~SurfExchange() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }
};

#endif