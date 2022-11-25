#pragma once

#include "BorisLib.h"
#include "Modules.h"

#if COMPILECUDA == 1
#include "STFieldCUDA.h"
#endif

class Mesh;
class Atom_Mesh;

#ifdef MODULE_COMPILATION_STFIELD

//STField module can only be used in a magnetic mesh

class STField :
	public Modules,
	public ProgramState<STField, std::tuple<>, std::tuple<>>
{

#if COMPILECUDA == 1
	friend STFieldCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Mesh * pMesh;

	//If STp is set to zero, these are the magnetic meshes which provide polarization vectors, top and bottom
	std::vector<Mesh*> pMesh_Bot, pMesh_Top;

	//as above for atomistic meshes top and bottom
	std::vector<Atom_Mesh*> paMesh_Bot, paMesh_Top;

private:

	//Mh is for "here", and Mc is what we're trying to couple to. xrel_h and yrel_h are relative to here, and zrel_c relative to coupling mesh
	bool check_cell_coupling(VEC_VC<DBL3>& Mh, VEC_VC<DBL3>& Mc, double xrel_h, double yrel_h, double zrel_c)
	{
		//relative coordinates to read value from coupled mesh (the one we're coupling to here)
		DBL3 cell_rel_pos = DBL3(xrel_h + Mh.rect.s.x - Mc.rect.s.x, yrel_h + Mh.rect.s.y - Mc.rect.s.y, zrel_c);
		//can't couple to an empty cell
		if (!Mc.rect.contains(cell_rel_pos + Mc.rect.s) || Mc.is_empty(cell_rel_pos)) return false;
		else return true;
	};

public:

	STField(Mesh *pMesh_);
	~STField();

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

class STField :
	public Modules
{

private:

public:

	STField(Mesh *pMesh_) {}
	~STField() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }
};

#endif