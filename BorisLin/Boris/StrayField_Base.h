#pragma once

#include "BorisLib.h"
#include "ErrorHandler.h"

#ifdef MODULE_COMPILATION_STRAYFIELD

class DipoleMesh;
class SuperMesh;

class StrayField_Base {

protected:

	SuperMesh *pSMesh;

	//strayField VEC taking values in this mesh
	VEC<DBL3> strayField;

	//all currently set dipole meshes
	std::vector<DipoleMesh*> dipoles;

protected:

	StrayField_Base(SuperMesh *pSMesh_);

	virtual ~StrayField_Base() {}

	//call from Initialize before starting computations, so all dipole meshes can be collected
	void InitializeStrayField(void);

	//calculate stray field from all dipoles
	void CalculateStrayField(void);

	//check if stray field needs to be recalculated
	bool CheckRecalculateStrayField(void);
};

#endif