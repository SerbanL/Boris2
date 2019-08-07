#include "stdafx.h"
#include "Exchange.h"

#ifdef MODULE_EXCHANGE

#include "Mesh_Ferromagnetic.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "ExchangeCUDA.h"
#endif

/////////////////////////////////////////////////////////////////
//Exch_Const_6ngbr
//
//Hex_i = 2*A/(MU0*Ms^2) * delsq M_i       (delsq is the Laplace operator)
//
//delsq approximated using the 6 neighbor method (with Neumann boundary condition for the boundary cells): 
//
//delsq M_i = sum(j belongs to Ni) (M_j - M_i)/h^2, where Ni is the set of available neighbors for M_i : cells j that share a face with cell i
//
//Exchange energy density contribution from cell i is given by :
//
//Eex,i = -(MU0/2)*M_i.Hex_i = A/Ms^2 * sum(j belongs to Ni) M_i.(M_i-M_j)/h^2
//
// Using Neumann boundary conditions : at boundary, derivative normal to the boundary is zero. Thus to evaluate Laplacian we use:
//
// On inner cells (i, j, k) : M(i-1,j,k) + M(i+1,j,k) + M(i,j-1,k) + M(i,j+1,k) + M(i,j,k-1) + M(i,j,k+1) - 6 * M(i,j,k). Number of neighbors is 6 in this case.
// On boundary cells (including corners) : replace fictitious cell outside of mesh with value of opposite neighbor along same axis and use same formula as for inner cells.
// If both neighbors are missing along an axis (e.g. a 2D problem) then reduce number of neighbors used by 2 (i.e. decrease dimensionality of inner cell Laplacian formula)

Exch_6ngbr_Neu::Exch_6ngbr_Neu(Mesh *pMesh_) : 
	Modules(),
	ProgramStateNames(this, {}, {})
{
	pMesh = dynamic_cast<FMesh*>(pMesh_);

	error_on_create = UpdateConfiguration();

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

Exch_6ngbr_Neu::~Exch_6ngbr_Neu() 
{
}

BError Exch_6ngbr_Neu::Initialize(void) 
{
	BError error(CLASS_STR(Exch_6ngbr_Neu));

	initialized = true;

	return error;
}

BError Exch_6ngbr_Neu::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Exch_6ngbr_Neu));

	Uninitialize();

	Initialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration();
	}
#endif

	return error;
}

BError Exch_6ngbr_Neu::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Exch_6ngbr_Neu));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Exch_6ngbr_NeuCUDA(dynamic_cast<FMeshCUDA*>(pMesh->pMeshCUDA));
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Exch_6ngbr_Neu::UpdateField(void) 
{
	double energy = 0;

#pragma omp parallel for reduction(+:energy)
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			double Ms = pMesh->Ms;
			double A = pMesh->A;
			pMesh->update_parameters_mcoarse(idx, pMesh->A, A, pMesh->Ms, Ms);

			DBL3 Hexch = (2*A / (MU0*Ms*Ms)) * pMesh->M.delsq_neu(idx);

			pMesh->Heff[idx] += Hexch;

			energy += pMesh->M[idx] * Hexch;
		}
	}

	//average energy density, this is correct, see notes - e_ex ~= -(mu0/2) M.H as can be derived from the full expression for e_ex. It is not -mu0 M.H, which is the energy density for a dipole in a field !
	if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / (2 * (pMesh->M.get_nonempty_cells()));
	else energy = 0;

	this->energy = energy;

	return energy;
}

#endif
