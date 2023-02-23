#pragma once

#include "BorisLib.h"
#include "ErrorHandler.h"

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1
#include "ExchangeBaseCUDA.h"
#endif

class SuperMesh;
class Mesh;

class ExchangeBase {

#if COMPILECUDA == 1
	friend ExchangeBaseCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Mesh *pMesh;

	//CMBND contacts between this mesh and other (anti)ferromagnetic meshes (we do not require other ferromagnetic meshes to have an exchange module enabled, just this one).
	std::vector<CMBNDInfo> CMBNDcontacts;

	//vector of pointers to all M - CMBNDInfo has a data member INT2 mesh_idx; mesh_idx.secondary is an index in pM for a mesh in contact with this one
	std::vector<VEC_VC<DBL3>*> pM;

	//vector of pointers to all ferromagnetic meshes - same ordering as pM
	std::vector<Mesh*> pMeshes;

protected:

protected:

	//this is overloaded by inheriting Exchange-type modules. Need this to be virtual so if for any reason a base pointer is used, the overloaded method is called instead.
	BError Initialize(void);

	//calculate exchange field at coupled cells in this mesh; accumulate energy density contribution in energy
	//the method only implements the coupling method, the actual computation must be supplied by the caller through the calculate_coupling function
	//this method takes the following inputs: cell1_idx, cell2_idx, relpos_m1, stencil, hshift_primary, Mesh_pri, Mesh_sec, and returns an energy density value to accumulate; hshift_primary is the normal direction cellsize (from primary to secondary)
	void CalculateExchangeCoupling(
		double& energy,
		std::function<double(int, int, DBL3, DBL3, DBL3, Mesh&, Mesh&)> calculate_coupling);

	//protected constructor - this class should not be instantiated by itself, but only used as a base for an exchange-type module for purposes of code reuse
	ExchangeBase(Mesh *pMesh_);

	virtual ~ExchangeBase() {}
};
