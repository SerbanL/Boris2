#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ZEEMAN) && ATOMISTIC == 1

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class Atom_Zeeman;
class Atom_MeshCUDA;

template <typename VType> class VEC;

class Atom_ZeemanCUDA :
	public ModulesCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	Atom_MeshCUDA* paMeshCUDA;

	//pointer to Zeeman holder module
	Atom_Zeeman* paZeeman;

	//Applied field
	cu_obj<cuReal3> Ha;

	//Applied field using user equation, thus allowing simultaneous spatial (x, y, z), stage time (t); base temperature (Tb) and stage step (Ss) are introduced as user constants.
	//A number of constants are always present : mesh dimensions in m (Lx, Ly, Lz)
	TEquationCUDA<cuBReal, cuBReal, cuBReal, cuBReal> H_equation;

	//Applied field but as a cuVEC (e.g. loaded from file), with same resolution as M.
	cu_obj<cuVEC<cuReal3>> Havec;

	//global field as obtained in this mesh
	//When a global field is set in supermesh, then globalField VEC is initialized here with mesh transfered values
	cu_obj<cuVEC<cuReal3>> globalField;

private:

	void set_Atom_ZeemanCUDA_pointers(void);

	//setup globalField transfer
	BError InitializeGlobalField(void);

public:

	Atom_ZeemanCUDA(Atom_MeshCUDA* paMeshCUDA_, Atom_Zeeman* paZeeman_);
	~Atom_ZeemanCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	void UpdateField(void);

	//-------------------

	void SetField(cuReal3 Hxyz);

	cu_obj<cuReal3>& GetField_cu_obj(void) { return Ha; }

	BError SetFieldEquation(const std::vector<std::vector< std::vector<EqComp::FSPEC> >>& fspec);

	BError SetFieldVEC(VEC<DBL3>& Havec_cpu);

	//-------------------Torque methods

	cuReal3 GetTorque(cuRect avRect);
};

#else

class Atom_ZeemanCUDA
{
};

#endif

#endif

