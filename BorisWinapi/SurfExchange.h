#pragma once

#include "BorisLib.h"
#include "Modules.h"

#if COMPILECUDA == 1
#include "SurfExchangeCUDA.h"
#endif

using namespace std;

class Mesh;
class FMesh;

#ifdef MODULE_SURFEXCHANGE

//Exchange modules can only be used in a ferromagnetic mesh

class SurfExchange :
	public Modules,
	public ProgramState<SurfExchange, tuple<>, tuple<>>
{

#if COMPILECUDA == 1
	friend SurfExchangeCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	FMesh* pMesh;

	//ferromagnetic meshes in surface exchange coupling with the mesh holding this module, top and bottom
	vector<FMesh*> pMesh_Bot, pMesh_Top;

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

	BError MakeCUDAModule(void);

	double UpdateField(void);
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

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC) { return BError(); }

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }
};

#endif