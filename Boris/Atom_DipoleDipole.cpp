#include "stdafx.h"
#include "Atom_DipoleDipole.h"
#include "SuperMesh.h"

#if defined(MODULE_COMPILATION_ATOM_DIPOLEDIPOLE) && ATOMISTIC == 1

#include "Atom_Mesh.h"

#if COMPILECUDA == 1
#include "Atom_DipoleDipoleCUDA.h"
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////

Atom_DipoleDipole::Atom_DipoleDipole(Atom_Mesh *paMesh_) :
	Modules(),
	Convolution<Atom_DipoleDipole, DipoleDipoleKernel>(paMesh_->n_dm, paMesh_->h_dm),
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

	paMesh->M1.set_pbc(0, 0, 0);

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

		if (!Hd.resize(paMesh->h_dm, paMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!Hd.Initialize_MeshTransfer({}, { &paMesh->Heff1 }, MESHTRANSFERTYPE_ENLARGED)) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		M.clear();
		Hd.clear();
	}

	if (paMesh->pSMesh->GetEvaluationSpeedup()) {

		//make sure to allocate memory for Hdemag if we need it
		if (paMesh->pSMesh->GetEvaluationSpeedup() >= 3) { if (!Hdemag3.resize(paMesh->h_dm, paMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT); }
		else Hdemag3.clear();

		if (paMesh->pSMesh->GetEvaluationSpeedup() >= 2) { if (!Hdemag2.resize(paMesh->h_dm, paMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT); }
		else Hdemag2.clear();

		if (paMesh->pSMesh->GetEvaluationSpeedup() >= 1) { if (!Hdemag.resize(paMesh->h_dm, paMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT); }
		else Hdemag.clear();

		if (Hdemag.linear_size()) if (!Hdemag.Initialize_MeshTransfer({}, { &paMesh->Heff1 }, MESHTRANSFERTYPE_ENLARGED)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (Hdemag2.linear_size()) if (!Hdemag2.Initialize_MeshTransfer({}, { &paMesh->Heff1 }, MESHTRANSFERTYPE_ENLARGED)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (Hdemag3.linear_size()) if (!Hdemag3.Initialize_MeshTransfer({}, { &paMesh->Heff1 }, MESHTRANSFERTYPE_ENLARGED)) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		Hdemag.clear();
		Hdemag2.clear();
		Hdemag3.clear();
	}

	num_Hdemag_saved = 0;
	   
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

		if (using_macrocell) selfDemagCoeff = DipoleDipoleTFunc().SelfDemag_PBC(paMesh->h_dm, paMesh->n_dm, demag_pbc_images);

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

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		paMesh->h, paMesh->meshRect, 
		(MOD_)paMesh->Get_Module_Heff_Display() == MOD_ATOM_DIPOLEDIPOLE || paMesh->IsOutputDataSet_withRect(DATA_E_DEMAG),
		(MOD_)paMesh->Get_Module_Energy_Display() == MOD_ATOM_DIPOLEDIPOLE || paMesh->IsOutputDataSet_withRect(DATA_E_DEMAG));
	if (error)	initialized = false;

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

		Hdemag.clear();
		Hdemag2.clear();
		Hdemag3.clear();
	}

	num_Hdemag_saved = 0;

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
	BError error(CLASS_STR(Atom_DipoleDipole));

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
	//transfer magnetic moments to macrocell mesh
	if (using_macrocell) M.transfer_in();

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////// NO SPEEDUP //////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (!paMesh->pSMesh->GetEvaluationSpeedup() || (num_Hdemag_saved < paMesh->pSMesh->GetEvaluationSpeedup() && !paMesh->pSMesh->Check_Step_Update())) {

		//don't use evaluation speedup

		if (using_macrocell) {

			//convolute and get "energy" value
			if (Module_Heff.linear_size()) Convolute(M, Hd, true, &Module_Heff, &Module_energy);
			else Convolute(M, Hd, true);
			//energy not calculated in macrocell mode : would need to correct for use of self demag term in macrocell
			energy = 0.0;

			//transfer dipole-dipole field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
			Hd.transfer_out();
		}
		else {

			//not using macrocell so get moments directly from mesh

			//convolute and get "energy" value
			if (Module_Heff.linear_size()) energy = Convolute(paMesh->M1, paMesh->Heff1, false, &Module_Heff, &Module_energy);
			else energy = Convolute(paMesh->M1, paMesh->Heff1, false);

			//finish off energy value
			if (non_empty_volume) energy *= -MUB_MU0 / (2 * non_empty_volume);
			else energy = 0;
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////// EVAL SPEEDUP /////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else {

		//use evaluation speedup method (Hdemag will have memory allocated - this was done in the Initialize method)

		//update if required by ODE solver or if we don't have enough previous evaluations saved to extrapolate
		if (paMesh->pSMesh->Check_Step_Update() || num_Hdemag_saved < paMesh->pSMesh->GetEvaluationSpeedup()) {

			VEC<DBL3>* pHdemag;

			if (num_Hdemag_saved < paMesh->pSMesh->GetEvaluationSpeedup()) {

				//don't have enough evaluations, so save next one
				switch (num_Hdemag_saved)
				{
				case 0:
					pHdemag = &Hdemag;
					time_demag1 = paMesh->pSMesh->Get_EvalStep_Time();
					break;
				case 1:
					pHdemag = &Hdemag2;
					time_demag2 = paMesh->pSMesh->Get_EvalStep_Time();
					break;
				case 2:
					pHdemag = &Hdemag3;
					time_demag3 = paMesh->pSMesh->Get_EvalStep_Time();
					break;
				}

				num_Hdemag_saved++;
			}
			else {

				//have enough evaluations saved, so just cycle between them now

				//QUADRATIC
				if (paMesh->pSMesh->GetEvaluationSpeedup() == 3) {

					//1, 2, 3 -> next is 1
					if (time_demag3 > time_demag2 && time_demag2 > time_demag1) {

						pHdemag = &Hdemag;
						time_demag1 = paMesh->pSMesh->Get_EvalStep_Time();
					}
					//2, 3, 1 -> next is 2
					else if (time_demag3 > time_demag2 && time_demag1 > time_demag2) {

						pHdemag = &Hdemag2;
						time_demag2 = paMesh->pSMesh->Get_EvalStep_Time();
					}
					//3, 1, 2 -> next is 3, leading to 1, 2, 3 again
					else {

						pHdemag = &Hdemag3;
						time_demag3 = paMesh->pSMesh->Get_EvalStep_Time();
					}
				}
				//LINEAR
				else if (paMesh->pSMesh->GetEvaluationSpeedup() == 2) {

					//1, 2 -> next is 1
					if (time_demag2 > time_demag1) {

						pHdemag = &Hdemag;
						time_demag1 = paMesh->pSMesh->Get_EvalStep_Time();
					}
					//2, 1 -> next is 2, leading to 1, 2 again
					else {

						pHdemag = &Hdemag2;
						time_demag2 = paMesh->pSMesh->Get_EvalStep_Time();
					}
				}
				//STEP
				else {

					pHdemag = &Hdemag;
				}
			}

			//convolute and get "energy" value
			if (using_macrocell) {

				//convolute and get "energy" value
				if (Module_Heff.linear_size()) Convolute(M, *pHdemag, true, &Module_Heff, &Module_energy);
				else Convolute(M, *pHdemag, true);
				//energy not calculated in macrocell mode : would need to correct for use of self demag term in macrocell
				energy = 0.0;
			}
			else {

				//not using macrocell so get moments directly from mesh

				//convolute and get "energy" value
				if (Module_Heff.linear_size()) energy = Convolute(paMesh->M1, *pHdemag, true, &Module_Heff, &Module_energy);
				else energy = Convolute(paMesh->M1, *pHdemag, true);

				//finish off energy value
				if (non_empty_volume) energy *= -MUB_MU0 / (2 * non_empty_volume);
				else energy = 0;
			}

			//transfer demagnetising field to atomistic mesh effective field : all atomistic cells within the larger micromagnetic cell receive the same field
			pHdemag->transfer_out();

			//subtract self contribution if using macrocell
			if (using_macrocell) {
				
				#pragma omp parallel for
				for (int idx = 0; idx < pHdemag->linear_size(); idx++) {

					//subtract self contribution: we'll add in again for the new moment, so it least the self contribution is exact
					(*pHdemag)[idx] -= (selfDemagCoeff & M[idx]);
				}
			}
		}
		else {

			//not required to update, and we have enough previous evaluations: use previous Hdemag saves to extrapolate for current evaluation

			double a1 = 1.0, a2 = 0.0, a3 = 0.0;
			double time = paMesh->pSMesh->Get_EvalStep_Time();

			//QUADRATIC
			if (paMesh->pSMesh->GetEvaluationSpeedup() == 3) {

				if (time_demag2 != time_demag1 && time_demag2 != time_demag3 && time_demag1 != time_demag3) {

					a1 = (time - time_demag2) * (time - time_demag3) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3));
					a2 = (time - time_demag1) * (time - time_demag3) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3));
					a3 = (time - time_demag1) * (time - time_demag2) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2));
				}

				//construct effective field approximation
				if (using_macrocell) {

					#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

					 	Hd[idx] = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + (selfDemagCoeff & M[idx]);
					}

					//add to Heff in the atomistic mesh
					Hd.transfer_out();
				}
				else {

					#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						paMesh->Heff1[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3;
					}
				}
			}
			//LINEAR
			else if (paMesh->pSMesh->GetEvaluationSpeedup() == 2) {

				if (time_demag2 != time_demag1) {

					a1 = (time - time_demag2) / (time_demag1 - time_demag2);
					a2 = (time - time_demag1) / (time_demag2 - time_demag1);
				}

				//construct effective field approximation
				if (using_macrocell) {

					#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						Hd[idx] = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + (selfDemagCoeff & M[idx]);
					}

					//add to Heff in the atomistic mesh
					Hd.transfer_out();
				}
				else {

					#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						paMesh->Heff1[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2;
					}
				}
			}
			//STEP
			else {

				//construct effective field approximation
				if (using_macrocell) {

					#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						Hd[idx] = Hdemag[idx] + (selfDemagCoeff & M[idx]);
					}

					//add to Heff in the atomistic mesh
					Hd.transfer_out();
				}
				else {

					#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						paMesh->Heff1[idx] += Hdemag[idx];
					}
				}
			}
		}
	}

	return energy;
}

#endif


