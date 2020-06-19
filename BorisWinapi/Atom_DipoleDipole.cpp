#include "stdafx.h"
#include "Atom_DipoleDipole.h"
#include "SuperMesh.h"

#if defined(MODULE_ATOM_DIPOLEDIPOLE) && ATOMISTIC == 1

#include "Atom_Mesh.h"

#if COMPILECUDA == 1
#include "Atom_DipoleDipoleCUDA.h"
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////

Atom_DipoleDipole::Atom_DipoleDipole(Atom_Mesh *paMesh_) :
	Modules(),
	Convolution<DipoleDipoleKernel>(paMesh_->n_dm, paMesh_->h_dm),
	ProgramStateNames(this, { VINFO(demag_pbc_images) }, {})
{
	paMesh = paMesh_;

	Uninitialize();

	error_on_create = Convolution_Error_on_Create();

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (paMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

Atom_DipoleDipole::~Atom_DipoleDipole()
{
	//when deleting the Demag module any pbc settings should no longer take effect in this mesh
	//thus must clear pbc flags in M1

	paMesh->M1.set_pbc(false, false, false);

	//same for the CUDA version if we are in cuda mode
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		paMesh->paMeshCUDA->M1()->copyflags_from_cpuvec(paMesh->M1);
	}
#endif
}

//Initialize mesh transfer from atomistic mesh to micromagnetic mesh for demag field computation
BError Atom_DipoleDipole::Initialize_Mesh_Transfer(void)
{
	BError error(CLASS_STR(Atom_DipoleDipole));

	using_macrocell = (paMesh->h_dm != paMesh->h);

	if (using_macrocell) {

		//M used only if using macrocell

		if (!M.resize(paMesh->h_dm, paMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!M.Initialize_MeshTransfer({ &paMesh->M1 }, {}, MESHTRANSFERTYPE_SUM)) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else M.clear();

	if (using_macrocell || paMesh->pSMesh->GetEvaluationSpeedup()) {

		//Hd used if using macrocell, or if using eval speedup

		if (!Hd.resize(paMesh->h_dm, paMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!Hd.Initialize_MeshTransfer({}, { &paMesh->Heff1 }, MESHTRANSFERTYPE_ENLARGED)) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else Hd.clear();

	return error;
}

BError Atom_DipoleDipole::Initialize(void)
{
	BError error(CLASS_STR(Atom_DipoleDipole));

	using_macrocell = (paMesh->h_dm != paMesh->h);

	if (!initialized) {

		error = Initialize_Mesh_Transfer();

		//only include self demag for macrocells. if not using macrocell (i.e. h_dm same as h) then there's no self demag term
		//Remember the kernel is pre-multiplied by muB, so the correct field results since the stored moments are in units of muB.
		error = Calculate_DipoleDipole_Kernels(using_macrocell);

		if (!error) initialized = true;
	}

	if (using_macrocell) {

		//need to calculate non-empty cells here so we don't waste time during computations (M is a VEC, not a VEC_VC, which means non-empty cells need to be calculated on every call)
		M.transfer_in();
		non_empty_volume = M.get_nonempty_cells() * M.h.dim();
	}
	else {

		non_empty_volume = paMesh->M1.get_nonempty_cells() * paMesh->M1.h.dim();
	}

	Hd_calculated = false;

	return error;
}

BError Atom_DipoleDipole::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_DipoleDipole));

	Uninitialize();

	//must enforce conditions:
	//1) h_dm has to have an equal integer number of h cells included in all 3 dimensions
	//2) h_dm still has to divide meshRect into n_dm

	//n_dm, h_dm combination is already correct as we know, thus we only check condition 1) here
	//if condition 1) not verified then set h_dm equal to h (and n_dm equal to n) : it is up to the user to make sure h_dm is set correctly.
	//this is an advanced user setting so it's a reasonable expectation.

	INT3 p_integer = paMesh->h_dm / paMesh->h;
	DBL3 p_double = paMesh->h_dm / paMesh->h;

	if (p_double != DBL3(p_integer) || !IsE(p_double.x, p_double.y) || !IsE(p_double.x, p_double.z) || !IsE(p_double.y, p_double.z)) {

		paMesh->h_dm = paMesh->h;
		paMesh->n_dm = paMesh->n;
	}

	//only need to uninitialize if n or h have changed, or pbc settings have changed
	if (!CheckDimensions(paMesh->n_dm, paMesh->h_dm, demag_pbc_images) || cfgMessage == UPDATECONFIG_MESHCHANGE) {

		Uninitialize();

		//Set convolution dimensions for embedded multiplication and required PBC conditions
		error = SetDimensions(paMesh->n_dm, paMesh->h_dm, true, demag_pbc_images);

		Hd.clear();
		M.clear();
	}

	Hd_calculated = false;

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

//Set PBC
BError Atom_DipoleDipole::Set_PBC(INT3 demag_pbc_images_)
{
	BError error(__FUNCTION__);

	demag_pbc_images = demag_pbc_images_;

	paMesh->Set_Magnetic_PBC(demag_pbc_images);

	//update will be needed if pbc settings have changed
	error = UpdateConfiguration(UPDATECONFIG_DEMAG_CONVCHANGE);

	return error;
}

BError Atom_DipoleDipole::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Atom_Demag));

#if COMPILECUDA == 1

	if (paMesh->paMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Atom_DipoleDipoleCUDA(paMesh->paMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Atom_DipoleDipole::UpdateField(void)
{
	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////// EVAL SPEEDUP /////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (paMesh->pSMesh->GetEvaluationSpeedup()) {

		//use evaluation speedup method (calculate Hdemag only when required, otherwise just update effective field with previously calculated Hdemag)

		int update_type = paMesh->pSMesh->Check_Step_Update();

		if (update_type != EVALSPEEDUPSTEP_SKIP || !Hd_calculated) {

			//calculate field required

			if (using_macrocell) {

				//transfer magnetic moments to macrocell mesh
				M.transfer_in();

				//convolute and get "energy" value
				Convolute(M, Hd, true);
				//energy not calculated in macrocell mode : would need to correct for use of self demag term in macrocell
				energy = 0.0;
			}
			else {

				//not using macrocell so get moments directly from mesh

				//convolute and get "energy" value
				energy = Convolute(paMesh->M1, Hd, true);

				//finish off energy value
				if (non_empty_volume) energy *= -MUB_MU0 / (2 * non_empty_volume);
				else energy = 0;
			}

			Hd_calculated = true;
		}

		//transfer dipole-dipole field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
		Hd.transfer_out();
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////// NO SPEEDUP //////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else {

		//don't use evaluation speedup

		if (using_macrocell) {

			//transfer magnetic moments to macrocell mesh
			M.transfer_in();

			//convolute and get "energy" value
			Convolute(M, Hd, true);
			//energy not calculated in macrocell mode : would need to correct for use of self demag term in macrocell
			energy = 0.0;

			//transfer dipole-dipole field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
			Hd.transfer_out();
		}
		else {

			//not using macrocell so get moments directly from mesh

			//convolute and get "energy" value
			energy = Convolute(paMesh->M1, paMesh->Heff1, false);

			//finish off energy value
			if (non_empty_volume) energy *= -MUB_MU0 / (2 * non_empty_volume);
			else energy = 0;
		}
	}

	return energy;
}

#endif


