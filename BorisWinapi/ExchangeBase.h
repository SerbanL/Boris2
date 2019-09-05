#pragma once

#include "BorisLib.h"
#include "ErrorHandler.h"

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1
#include "ExchangeBaseCUDA.h"
#endif

class SuperMesh;
class Mesh;
class FMesh;

using namespace std;

class ExchangeBase {

#if COMPILECUDA == 1
	friend ExchangeBaseCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	FMesh *pMesh;

	//CMBND contacts between this mesh and other ferromagnetic meshes (we do not require other ferromagnetic meshes to have an exchange module enabled, just this one).
	vector<CMBNDInfo> CMBNDcontacts;

	//vector of pointers to all M - CMBNDInfo has a data member INT2 mesh_idx; mesh_idx.secondary is an index in pM for a mesh in contact with this one
	vector<VEC_VC<DBL3>*> pM;

	//vector of pointers to all ferromagnetic meshes - same ordering as pM
	vector<FMesh*> pFMeshes;

protected:

	//this is overloaded by inheriting Exchange-type modules. Need this to be virtual so if for any reason a base pointer is used, the overloaded method is called instead.
	BError Initialize(void);

	//calculate exchange field at coupled cells in this mesh; accumulate energy density contribution in energy
	void CalculateExchangeCoupling(double& energy);

	//protected constructor - this class should not be instantiated by itself, but only used as a base for an exchange-type module for purposes of code reuse
	ExchangeBase(Mesh *pMesh_);

	virtual ~ExchangeBase() {}
};
