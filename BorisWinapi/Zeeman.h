#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Mesh;
class FMesh;

#ifdef MODULE_ZEEMAN

//Zeeman module can only be used in a ferromagnetic mesh

class Zeeman : 
	public Modules,
	public ProgramState<Zeeman, tuple<DBL3>, tuple<>>
{

private:

	//pointer to mesh object holding this effective field module
	FMesh *pMesh;

	//Applied field
	DBL3 Ha;

public:

	Zeeman(Mesh *pMesh_);
	~Zeeman();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	BError MakeCUDAModule(void);

	void UpdateField(void);

	//-------------------

	void SetField(DBL3 Hxyz);
	DBL3 GetField(void);
};

#else

class Zeeman :
	public Modules
{

private:

public:

	Zeeman(Mesh *pMesh_) {}
	~Zeeman() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC) { return BError(); }

	BError MakeCUDAModule(void) { return BError(); }

	void UpdateField(void) {}

	//-------------------

	void SetField(DBL3 Hxyz) {}
	DBL3 GetField(void) { return DBL3(); }
};

#endif