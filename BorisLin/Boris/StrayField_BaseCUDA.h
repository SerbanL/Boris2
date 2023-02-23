#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "BorisCUDALib.h"

class SuperMesh;
class DipoleMeshCUDA;
class ManagedMeshCUDA;

class StrayField_BaseCUDA
{

protected:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	SuperMesh* pSMesh;

	//strayField cuVEC taking values in this mesh
	cu_obj<cuVEC<cuReal3>> strayField;
	size_t strayField_size;

	//all currently set dipole meshes
	std::vector<DipoleMeshCUDA*> dipoles;
	//as above, but cu_obj managed versions
	cu_arr<ManagedMeshCUDA> managed_dipoles;

protected:

	StrayField_BaseCUDA(SuperMesh* pSMesh_);
	~StrayField_BaseCUDA();

	//call from Initialize before starting computations, so all dipole meshes can be collected
	void InitializeStrayField(void);

	//calculate stray field from all dipoles
	void CalculateStrayField(void);

	//check if stray field needs to be recalculated
	bool CheckRecalculateStrayField(void);
};

#else

class StrayField_BaseCUDA
{
};

#endif

#endif




