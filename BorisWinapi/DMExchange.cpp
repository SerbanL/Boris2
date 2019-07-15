#include "stdafx.h"
#include "DMExchange.h"

#ifdef MODULE_DMEXCHANGE

#include "Mesh_Ferromagnetic.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "DMExchangeCUDA.h"
#endif

//For bulk Dzyaloshinskii-Moriya exchange, we have :
//
//Hdm,ex = -2D/(mu0*Ms) * curl m
//
//Hdm,ex adds to the usual isotropic direct exchange term, energy calculated using the total exchange field as usual.
//
//Here D is the DM exchange constant (J/m^2)

DMExchange::DMExchange(Mesh *pMesh_) :
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

DMExchange::~DMExchange()
{
}

BError DMExchange::Initialize(void)
{
	BError error(CLASS_STR(DMExchange));

	initialized = true;

	return error;
}

BError DMExchange::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(DMExchange));

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

BError DMExchange::MakeCUDAModule(void)
{
	BError error(CLASS_STR(DMExchange));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new DMExchangeCUDA(dynamic_cast<FMeshCUDA*>(pMesh->pMeshCUDA));
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

void DMExchange::UpdateField(void)
{
	double energy = 0;

	INT3 n = pMesh->n;

#pragma omp parallel for reduction(+:energy) 
	for (int idx = 0; idx < n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			double Ms = pMesh->Ms;
			double A = pMesh->A;
			double D = pMesh->D;
			pMesh->update_parameters_mcoarse(idx, pMesh->A, A, pMesh->D, D, pMesh->Ms, Ms);

			double Aconst = 2 * A / (MU0 * Ms * Ms);
			double Dconst = -2 * D / (MU0 * Ms * Ms);

			//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
			DBL3 bnd_dm_dx = (D / (2 * A * Ms)) * DBL3(0, -pMesh->M[idx].z, pMesh->M[idx].y);
			DBL3 bnd_dm_dy = (D / (2 * A * Ms)) * DBL3(pMesh->M[idx].z, 0, -pMesh->M[idx].x);
			DBL3 bnd_dm_dz = (D / (2 * A * Ms)) * DBL3(-pMesh->M[idx].y, pMesh->M[idx].x, 0);
			DBL33 bnd_nneu = DBL33(bnd_dm_dx, bnd_dm_dy, bnd_dm_dz);

			//direct exchange contribution
			DBL3 Hexch = Aconst * pMesh->M.delsq_nneu(idx, bnd_nneu);

			//Dzyaloshinskii-Moriya exchange contribution
				
			//Hdm, ex = -2D / (mu0*Ms) * curl m
			Hexch += Dconst * pMesh->M.curl_nneu(idx, bnd_nneu);
				
			pMesh->Heff[idx] += Hexch;

			energy += pMesh->M[idx] * Hexch;
		}
	}

	//average energy density, this is correct, see notes - e_ex ~= -(mu0/2) M.H as can be derived from the full expression for e_ex. It is not -mu0 M.H, which is the energy density for a dipole in a field !
	if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / (2 * (pMesh->M.get_nonempty_cells()));
	else energy = 0;

	this->energy = energy;
}

#endif