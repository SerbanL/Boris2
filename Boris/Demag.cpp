#include "stdafx.h"
#include "Demag.h"
#include "SuperMesh.h"

#ifdef MODULE_COMPILATION_DEMAG

#include "SimScheduleDefs.h"

#include "Mesh.h"

#if COMPILECUDA == 1
#include "DemagCUDA.h"
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////

Demag::Demag(Mesh *pMesh_) : 
	Modules(),
	Convolution<Demag, DemagKernel>(pMesh_->GetMeshSize(), pMesh_->GetMeshCellsize()),
	ProgramStateNames(this, {VINFO(demag_pbc_images)}, {})
{
	pMesh = pMesh_;

	Uninitialize();

	error_on_create = Convolution_Error_on_Create();

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

Demag::~Demag()
{
	//when deleting the Demag module any pbc settings should no longer take effect in this mesh
	//thus must clear pbc flags in M

	pMesh->M.set_pbc(0, 0, 0);
	pMesh->M2.set_pbc(0, 0, 0);

	//same for the CUDA version if we are in cuda mode
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pMesh->pMeshCUDA->M()->copyflags_from_cpuvec(pMesh->M);
	}
#endif
}

BError Demag::Initialize(void) 
{	
	BError error(CLASS_STR(Demag));

	if (!initialized) {
		
		error = Calculate_Demag_Kernels();

		selfDemagCoeff = DemagTFunc().SelfDemag_PBC(pMesh->h, pMesh->n, demag_pbc_images);

		if (!error) initialized = true;
	}

	//make sure to allocate memory for Hdemag if we need it
	if (pMesh->pSMesh->GetEvaluationSpeedup() >= 6) { if (!Hdemag6.resize(pMesh->h, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT); }
	else Hdemag6.clear();

	if (pMesh->pSMesh->GetEvaluationSpeedup() >= 5) { if (!Hdemag5.resize(pMesh->h, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT); }
	else Hdemag5.clear();

	if (pMesh->pSMesh->GetEvaluationSpeedup() >= 4) { if (!Hdemag4.resize(pMesh->h, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT); }
	else Hdemag4.clear();

	if (pMesh->pSMesh->GetEvaluationSpeedup() >= 3) { if (!Hdemag3.resize(pMesh->h, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT); }
	else Hdemag3.clear();

	if (pMesh->pSMesh->GetEvaluationSpeedup() >= 2) { if (!Hdemag2.resize(pMesh->h, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT); }
	else Hdemag2.clear();

	if (pMesh->pSMesh->GetEvaluationSpeedup() >= 1) { if (!Hdemag.resize(pMesh->h, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT); }
	else Hdemag.clear();

	num_Hdemag_saved = 0;

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect, 
		(MOD_)pMesh->Get_Module_Heff_Display() == MOD_DEMAG || pMesh->IsOutputDataSet_withRect(DATA_E_DEMAG) || pMesh->IsStageSet(SS_MONTECARLO),
		(MOD_)pMesh->Get_Module_Energy_Display() == MOD_DEMAG || pMesh->IsOutputDataSet_withRect(DATA_E_DEMAG) || pMesh->IsStageSet(SS_MONTECARLO));
	if (error) initialized = false;

	//if a Monte Carlo stage is set then we need to compute fields
	if (pMesh->IsStageSet(SS_MONTECARLO)) pMesh->Set_Force_MonteCarlo_ComputeFields(true);

	return error;
}

BError Demag::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Demag));

	//only need to uninitialize if n or h have changed, or pbc settings have changed
	if (!CheckDimensions(pMesh->n, pMesh->h, demag_pbc_images) || cfgMessage == UPDATECONFIG_DEMAG_CONVCHANGE) {
		
		Uninitialize();

		//Set convolution dimensions for embedded multiplication and required PBC conditions
		error = SetDimensions(pMesh->n, pMesh->h, true, demag_pbc_images);

		//if memory needs to be allocated for Hdemag, it will be done through Initialize 
		Hdemag.clear();
		Hdemag2.clear();
		Hdemag3.clear();
		Hdemag4.clear();
		Hdemag5.clear();
		Hdemag6.clear();
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
BError Demag::Set_PBC(INT3 demag_pbc_images_)
{
	BError error(__FUNCTION__);

	demag_pbc_images = demag_pbc_images_;

	pMesh->Set_Magnetic_PBC(demag_pbc_images);

	//update will be needed if pbc settings have changed
	error = UpdateConfiguration(UPDATECONFIG_DEMAG_CONVCHANGE);

	return error;
}

BError Demag::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Demag));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new DemagCUDA(pMesh->pMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Demag::UpdateField(void) 
{
	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////// NO SPEEDUP //////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (!pMesh->pSMesh->GetEvaluationSpeedup() || (num_Hdemag_saved < pMesh->pSMesh->GetEvaluationSpeedup() && !pMesh->pSMesh->Check_Step_Update())) {

		//don't use evaluation speedup, so no need to use Hdemag (this won't have memory allocated anyway) - or else we are using speedup but don't yet have enough previous evaluations at steps where we should be extrapolating

		//convolute and get "energy" value
		if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

			if (Module_Heff.linear_size()) energy = Convolute_AveragedInputs_DuplicatedOutputs(pMesh->M, pMesh->M2, pMesh->Heff, pMesh->Heff2, false, &Module_Heff, &Module_energy);
			else energy = Convolute_AveragedInputs_DuplicatedOutputs(pMesh->M, pMesh->M2, pMesh->Heff, pMesh->Heff2, false);
		}
		else {

			if (Module_Heff.linear_size()) energy = Convolute(pMesh->M, pMesh->Heff, false, &Module_Heff, &Module_energy);
			else energy = Convolute(pMesh->M, pMesh->Heff, false);
		}

		//finish off energy value
		if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / (2 * pMesh->M.get_nonempty_cells());
		else energy = 0;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////// EVAL SPEEDUP /////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	
	else {

		//use evaluation speedup method (Hdemag will have memory allocated - this was done in the Initialize method)

		//update if required by ODE solver or if we don't have enough previous evaluations saved to extrapolate
		if (pMesh->pSMesh->Check_Step_Update() || num_Hdemag_saved < pMesh->pSMesh->GetEvaluationSpeedup()) {

			VEC<DBL3>* pHdemag;

			if (num_Hdemag_saved < pMesh->pSMesh->GetEvaluationSpeedup()) {

				//don't have enough evaluations, so save next one
				switch (num_Hdemag_saved)
				{
				case 0:
					pHdemag = &Hdemag;
					time_demag1 = pMesh->pSMesh->Get_EvalStep_Time();
					break;
				case 1:
					pHdemag = &Hdemag2;
					time_demag2 = pMesh->pSMesh->Get_EvalStep_Time();
					break;
				case 2:
					pHdemag = &Hdemag3;
					time_demag3 = pMesh->pSMesh->Get_EvalStep_Time();
					break;
				case 3:
					pHdemag = &Hdemag4;
					time_demag4 = pMesh->pSMesh->Get_EvalStep_Time();
					break;
				case 4:
					pHdemag = &Hdemag5;
					time_demag5 = pMesh->pSMesh->Get_EvalStep_Time();
					break;
				case 5:
					pHdemag = &Hdemag6;
					time_demag6 = pMesh->pSMesh->Get_EvalStep_Time();
					break;
				}

				num_Hdemag_saved++;
			}
			else {

				//have enough evaluations saved, so just cycle between them now
				
				//QUINTIC
				if (pMesh->pSMesh->GetEvaluationSpeedup() == 6) {

					//1, 2, 3, 4, 5, 6 -> next is 1
					if (time_demag6 > time_demag5 && time_demag5 > time_demag4 && time_demag4 > time_demag3 && time_demag3 > time_demag2 && time_demag2 > time_demag1) {

						pHdemag = &Hdemag;
						time_demag1 = pMesh->pSMesh->Get_EvalStep_Time();
					}
					//2, 3, 4, 5, 6, 1 -> next is 2
					else if (time_demag1 > time_demag2) {

						pHdemag = &Hdemag2;
						time_demag2 = pMesh->pSMesh->Get_EvalStep_Time();
					}
					//3, 4, 5, 6, 1, 2 -> next is 3
					else if (time_demag2 > time_demag3) {

						pHdemag = &Hdemag3;
						time_demag3 = pMesh->pSMesh->Get_EvalStep_Time();
					}
					//4, 5, 6, 1, 2, 3 -> next is 4
					else if (time_demag3 > time_demag4) {

						pHdemag = &Hdemag4;
						time_demag4 = pMesh->pSMesh->Get_EvalStep_Time();
					}
					//5, 6, 1, 2, 3, 4 -> next is 5
					else if (time_demag4 > time_demag5) {

						pHdemag = &Hdemag5;
						time_demag5 = pMesh->pSMesh->Get_EvalStep_Time();
					}
					else {

						pHdemag = &Hdemag6;
						time_demag6 = pMesh->pSMesh->Get_EvalStep_Time();
					}
				}
				//QUARTIC
				else if (pMesh->pSMesh->GetEvaluationSpeedup() == 5) {

					//1, 2, 3, 4, 5 -> next is 1
					if (time_demag5 > time_demag4 && time_demag4 > time_demag3 && time_demag3 > time_demag2 && time_demag2 > time_demag1) {

						pHdemag = &Hdemag;
						time_demag1 = pMesh->pSMesh->Get_EvalStep_Time();
					}
					//2, 3, 4, 5, 1 -> next is 2
					else if (time_demag1 > time_demag2) {

						pHdemag = &Hdemag2;
						time_demag2 = pMesh->pSMesh->Get_EvalStep_Time();
					}
					//3, 4, 5, 1, 2 -> next is 3
					else if (time_demag2 > time_demag3) {

						pHdemag = &Hdemag3;
						time_demag3 = pMesh->pSMesh->Get_EvalStep_Time();
					}
					//4, 5, 1, 2, 3 -> next is 4
					else if (time_demag3 > time_demag4) {

						pHdemag = &Hdemag4;
						time_demag4 = pMesh->pSMesh->Get_EvalStep_Time();
					}
					else {

						pHdemag = &Hdemag5;
						time_demag5 = pMesh->pSMesh->Get_EvalStep_Time();
					}
				}
				//CUBIC
				else if (pMesh->pSMesh->GetEvaluationSpeedup() == 4) {

					//1, 2, 3, 4 -> next is 1
					if (time_demag4 > time_demag3 && time_demag3 > time_demag2 && time_demag2 > time_demag1) {

						pHdemag = &Hdemag;
						time_demag1 = pMesh->pSMesh->Get_EvalStep_Time();
					}
					//2, 3, 4, 1 -> next is 2
					else if (time_demag1 > time_demag2) {

						pHdemag = &Hdemag2;
						time_demag2 = pMesh->pSMesh->Get_EvalStep_Time();
					}
					//3, 4, 1, 2 -> next is 3
					else if (time_demag2 > time_demag3) {

						pHdemag = &Hdemag3;
						time_demag3 = pMesh->pSMesh->Get_EvalStep_Time();
					}
					else {

						pHdemag = &Hdemag4;
						time_demag4 = pMesh->pSMesh->Get_EvalStep_Time();
					}
				}
				//QUADRATIC
				else if (pMesh->pSMesh->GetEvaluationSpeedup() == 3) {

					//1, 2, 3 -> next is 1
					if (time_demag3 > time_demag2 && time_demag2 > time_demag1) {

						pHdemag = &Hdemag;
						time_demag1 = pMesh->pSMesh->Get_EvalStep_Time();
					}
					//2, 3, 1 -> next is 2
					else if (time_demag1 > time_demag2) {

						pHdemag = &Hdemag2;
						time_demag2 = pMesh->pSMesh->Get_EvalStep_Time();
					}
					//3, 1, 2 -> next is 3, leading to 1, 2, 3 again
					else {

						pHdemag = &Hdemag3;
						time_demag3 = pMesh->pSMesh->Get_EvalStep_Time();
					}
				}
				//LINEAR
				else if (pMesh->pSMesh->GetEvaluationSpeedup() == 2) {

					//1, 2 -> next is 1
					if (time_demag2 > time_demag1) {

						pHdemag = &Hdemag;
						time_demag1 = pMesh->pSMesh->Get_EvalStep_Time();
					}
					//2, 1 -> next is 2, leading to 1, 2 again
					else {

						pHdemag = &Hdemag2;
						time_demag2 = pMesh->pSMesh->Get_EvalStep_Time();
					}
				}
				//STEP
				else {

					pHdemag = &Hdemag;
				}
			}
			
			//do evaluation
			if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

				if (Module_Heff.linear_size()) energy = Convolute_AveragedInputs(pMesh->M, pMesh->M2, *pHdemag, true, &Module_Heff, &Module_energy);
				else energy = Convolute_AveragedInputs(pMesh->M, pMesh->M2, *pHdemag, true);
			}
			else {

				if (Module_Heff.linear_size()) energy = Convolute(pMesh->M, *pHdemag, true, &Module_Heff, &Module_energy);
				else energy = Convolute(pMesh->M, *pHdemag, true);
			}

			//finish off energy value
			if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / (2 * pMesh->M.get_nonempty_cells());
			else energy = 0;

			if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

				//add contribution to Heff and Heff2
#pragma omp parallel for
				for (int idx = 0; idx < pHdemag->linear_size(); idx++) {

					pMesh->Heff[idx] += (*pHdemag)[idx];
					pMesh->Heff2[idx] += (*pHdemag)[idx];
					//subtract self demag contribution: we'll add in again for the new magnetization, so it least the self demag is exact
					(*pHdemag)[idx] -= (selfDemagCoeff & (pMesh->M[idx] + pMesh->M2[idx]) / 2);
				}
			}
			else {

				//add contribution to Heff
#pragma omp parallel for
				for (int idx = 0; idx < pHdemag->linear_size(); idx++) {

					pMesh->Heff[idx] += (*pHdemag)[idx];
					//subtract self demag contribution: we'll add in again for the new magnetization, so it least the self demag is exact
					(*pHdemag)[idx] -= (selfDemagCoeff & pMesh->M[idx]);
				}
			}
		}
		else {

			//not required to update, and we have enough previous evaluations: use previous Hdemag saves to extrapolate for current evaluation

			double a1 = 1.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, a5 = 0.0, a6 = 0.0;
			double time = pMesh->pSMesh->Get_EvalStep_Time();

			//QUINTIC
			if (pMesh->pSMesh->GetEvaluationSpeedup() == 6) {

				a1 = (time - time_demag2) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) * (time - time_demag6) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3) * (time_demag1 - time_demag4) * (time_demag1 - time_demag5) * (time_demag1 - time_demag6));
				a2 = (time - time_demag1) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) * (time - time_demag6) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3) * (time_demag2 - time_demag4) * (time_demag2 - time_demag5) * (time_demag2 - time_demag6));
				a3 = (time - time_demag1) * (time - time_demag2) * (time - time_demag4) * (time - time_demag5) * (time - time_demag6) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2) * (time_demag3 - time_demag4) * (time_demag3 - time_demag5) * (time_demag3 - time_demag6));
				a4 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag5) * (time - time_demag6) / ((time_demag4 - time_demag1) * (time_demag4 - time_demag2) * (time_demag4 - time_demag3) * (time_demag4 - time_demag5) * (time_demag4 - time_demag6));
				a5 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag4) * (time - time_demag6) / ((time_demag5 - time_demag1) * (time_demag5 - time_demag2) * (time_demag5 - time_demag3) * (time_demag5 - time_demag4) * (time_demag5 - time_demag6));
				a6 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) / ((time_demag6 - time_demag1) * (time_demag6 - time_demag2) * (time_demag6 - time_demag3) * (time_demag6 - time_demag4) * (time_demag6 - time_demag5));

				if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					//add contribution to Heff and Heff2
#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						DBL3 Hdemag_value = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + Hdemag6[idx] * a6 + (selfDemagCoeff & (pMesh->M[idx] + pMesh->M2[idx]) / 2);
						pMesh->Heff[idx] += Hdemag_value;
						pMesh->Heff2[idx] += Hdemag_value;
					}
				}
				else {

					//add contribution to Heff
#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						pMesh->Heff[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + Hdemag6[idx] * a6 + (selfDemagCoeff & pMesh->M[idx]);
					}
				}
			}
			//QUARTIC
			else if (pMesh->pSMesh->GetEvaluationSpeedup() == 5) {

				a1 = (time - time_demag2) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3) * (time_demag1 - time_demag4) * (time_demag1 - time_demag5));
				a2 = (time - time_demag1) * (time - time_demag3) * (time - time_demag4) * (time - time_demag5) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3) * (time_demag2 - time_demag4) * (time_demag2 - time_demag5));
				a3 = (time - time_demag1) * (time - time_demag2) * (time - time_demag4) * (time - time_demag5) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2) * (time_demag3 - time_demag4) * (time_demag3 - time_demag5));
				a4 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag5) / ((time_demag4 - time_demag1) * (time_demag4 - time_demag2) * (time_demag4 - time_demag3) * (time_demag4 - time_demag5));
				a5 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) * (time - time_demag4) / ((time_demag5 - time_demag1) * (time_demag5 - time_demag2) * (time_demag5 - time_demag3) * (time_demag5 - time_demag4));

				if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					//add contribution to Heff and Heff2
#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						DBL3 Hdemag_value = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + (selfDemagCoeff & (pMesh->M[idx] + pMesh->M2[idx]) / 2);
						pMesh->Heff[idx] += Hdemag_value;
						pMesh->Heff2[idx] += Hdemag_value;
					}
				}
				else {

					//add contribution to Heff
#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						pMesh->Heff[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + (selfDemagCoeff & pMesh->M[idx]);
					}
				}
			}
			//CUBIC
			else if (pMesh->pSMesh->GetEvaluationSpeedup() == 4) {

				a1 = (time - time_demag2) * (time - time_demag3) * (time - time_demag4) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3) * (time_demag1 - time_demag4));
				a2 = (time - time_demag1) * (time - time_demag3) * (time - time_demag4) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3) * (time_demag2 - time_demag4));
				a3 = (time - time_demag1) * (time - time_demag2) * (time - time_demag4) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2) * (time_demag3 - time_demag4));
				a4 = (time - time_demag1) * (time - time_demag2) * (time - time_demag3) / ((time_demag4 - time_demag1) * (time_demag4 - time_demag2) * (time_demag4 - time_demag3));

				if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					//add contribution to Heff and Heff2
#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						DBL3 Hdemag_value = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + (selfDemagCoeff & (pMesh->M[idx] + pMesh->M2[idx]) / 2);
						pMesh->Heff[idx] += Hdemag_value;
						pMesh->Heff2[idx] += Hdemag_value;
					}
				}
				else {

					//add contribution to Heff
#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						pMesh->Heff[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + (selfDemagCoeff & pMesh->M[idx]);
					}
				}
			}
			//QUADRATIC
			else if (pMesh->pSMesh->GetEvaluationSpeedup() == 3) {

				a1 = (time - time_demag2) * (time - time_demag3) / ((time_demag1 - time_demag2) * (time_demag1 - time_demag3));
				a2 = (time - time_demag1) * (time - time_demag3) / ((time_demag2 - time_demag1) * (time_demag2 - time_demag3));
				a3 = (time - time_demag1) * (time - time_demag2) / ((time_demag3 - time_demag1) * (time_demag3 - time_demag2));

				if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					//add contribution to Heff and Heff2
#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						DBL3 Hdemag_value = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + (selfDemagCoeff & (pMesh->M[idx] + pMesh->M2[idx]) / 2);
						pMesh->Heff[idx] += Hdemag_value;
						pMesh->Heff2[idx] += Hdemag_value;
					}
				}
				else {

					//add contribution to Heff
#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						pMesh->Heff[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + (selfDemagCoeff & pMesh->M[idx]);
					}
				}
			}
			//LINEAR
			else if (pMesh->pSMesh->GetEvaluationSpeedup() == 2) {

				a1 = (time - time_demag2) / (time_demag1 - time_demag2);
				a2 = (time - time_demag1) / (time_demag2 - time_demag1);

				if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					//add contribution to Heff and Heff2
#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						DBL3 Hdemag_value = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + (selfDemagCoeff & (pMesh->M[idx] + pMesh->M2[idx]) / 2);
						pMesh->Heff[idx] += Hdemag_value;
						pMesh->Heff2[idx] += Hdemag_value;
					}
				}
				else {

					//add contribution to Heff
#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						pMesh->Heff[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2 + (selfDemagCoeff & pMesh->M[idx]);
					}
				}
			}
			//STEP
			else {

				if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					//add contribution to Heff and Heff2
#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						DBL3 Hdemag_value = Hdemag[idx] + (selfDemagCoeff & (pMesh->M[idx] + pMesh->M2[idx]) / 2);
						pMesh->Heff[idx] += Hdemag_value;
						pMesh->Heff2[idx] += Hdemag_value;
					}
				}
				else {

					//add contribution to Heff
#pragma omp parallel for
					for (int idx = 0; idx < Hdemag.linear_size(); idx++) {

						pMesh->Heff[idx] += Hdemag[idx] + (selfDemagCoeff & pMesh->M[idx]);
					}
				}
			}
		}
	}

	return energy;
}

//-------------------Energy methods

//FM mesh
double Demag::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//Module_Heff needs to be calculated (done during a Monte Carlo simulation, where this method would be used)
	if (Module_Heff.linear_size()) {

		//do not divide by 2 as we are not double-counting here
		if (Mnew != DBL3()) return -pMesh->h.dim() * MU0 * Module_Heff[pMesh->M.cellidx_to_position(spin_index)] * (Mnew - pMesh->M[spin_index]);
		else return -pMesh->h.dim() * MU0 * Module_Heff[pMesh->M.cellidx_to_position(spin_index)] * pMesh->M[spin_index];
	}
	else return 0.0;
}

//AFM mesh
DBL2 Demag::Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B)
{
	//Module_Heff needs to be calculated (done during a Monte Carlo simulation, where this method would be used)
	if (Module_Heff.linear_size() && Module_Heff2.linear_size()) {

		DBL3 M = (pMesh->M[spin_index] + pMesh->M2[spin_index]) / 2;
		DBL3 Mnew = (Mnew_A + Mnew_B) / 2;

		double energy_ = 0.0;

		//do not divide by 2 as we are not double-counting here
		if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

			energy_ = -pMesh->h.dim() * MU0 * Module_Heff[pMesh->M.cellidx_to_position(spin_index)] * (Mnew - M);
		}
		else {

			energy_ = -pMesh->h.dim() * MU0 * Module_Heff[pMesh->M.cellidx_to_position(spin_index)] * M;
		}

		return DBL2(energy_, energy_);
	}
	else return DBL2();
}

#endif



