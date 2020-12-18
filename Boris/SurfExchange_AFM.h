#pragma once

#include "BorisLib.h"
#include "Modules.h"

#if COMPILECUDA == 1
#include "SurfExchangeCUDA_AFM.h"
#endif



class Mesh;

#ifdef MODULE_COMPILATION_SURFEXCHANGE

//This SurfExchange_AFM module is only used in a two-sublattice mesh; there is a SurfExchange module used in a ferromagnetic mesh

class SurfExchange_AFM :
	public Modules,
	public ProgramState<SurfExchange_AFM, std::tuple<>, std::tuple<>>
{

#if COMPILECUDA == 1
	friend SurfExchangeCUDA_AFM;
#endif

private:

	//pointer to mesh object holding this effective field module
	Mesh* pMesh;

	//magnetic meshes in surface exchange coupling with the mesh holding this module, top and bottom
	std::vector<Mesh*> pMesh_Bot, pMesh_Top;

	//number of coupled cells (either top or bottom) in this mesh
	int coupled_cells = 0;

public:

	SurfExchange_AFM(Mesh *pMesh_);
	~SurfExchange_AFM();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);
};

#else

class SurfExchange_AFM :
	public Modules
{

private:

public:

	SurfExchange_AFM(Mesh *pMesh_) {}
	~SurfExchange_AFM() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }
};

#endif