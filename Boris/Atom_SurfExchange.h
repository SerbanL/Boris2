#pragma once

#include "BorisLib.h"
#include "Modules.h"

#if COMPILECUDA == 1
#include "Atom_SurfExchangeCUDA.h"
#endif

class Atom_Mesh;
class Mesh;

#if defined(MODULE_COMPILATION_SURFEXCHANGE) && ATOMISTIC == 1

class Atom_SurfExchange :
	public Modules,
	public ProgramState<Atom_SurfExchange, std::tuple<>, std::tuple<>>
{

#if COMPILECUDA == 1
	friend Atom_SurfExchangeCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Atom_Mesh* paMesh;

	//magnetic meshes in surface exchange coupling with the mesh holding this module, top and bottom (atomistic meshes)
	std::vector<Atom_Mesh*> paMesh_Bot, paMesh_Top;

	//magnetic meshes in surface exchange coupling with the mesh holding this module, top and bottom (micromagnetic meshes)
	std::vector<Mesh*> pMesh_Bot, pMesh_Top;

public:

	Atom_SurfExchange(Atom_Mesh *paMesh_);
	~Atom_SurfExchange();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Torque methods

	DBL3 GetTorque(Rect& avRect);

	//-------------------Energy methods

	//For simple cubic mesh spin_index coincides with index in M1
	double Get_EnergyChange(int spin_index, DBL3 Mnew);
};

#else

class Atom_SurfExchange :
	public Modules
{

private:

public:

	Atom_SurfExchange(Atom_Mesh *paMesh_) {}
	~Atom_SurfExchange() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }
};

#endif