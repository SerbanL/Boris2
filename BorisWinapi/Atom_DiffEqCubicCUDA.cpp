#include "stdafx.h"
#include "Atom_DiffEqCubicCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_DiffEqCubic.h"

#include "Atom_Mesh_Cubic.h"
#include "Atom_Mesh_CubicCUDA.h"

Atom_DifferentialEquationCubicCUDA::Atom_DifferentialEquationCubicCUDA(Atom_DifferentialEquation *pameshODE) :
	Atom_DifferentialEquationCUDA(pameshODE)
{
	error_on_create = AllocateMemory(true);

	cuaDiffEq()->set_pointers(this);
	SetODEMethodPointers();
}

Atom_DifferentialEquationCubicCUDA::~Atom_DifferentialEquationCubicCUDA()
{
	//If called_from_destructor is true then do not attempt to transfer data to cpu where this is held in the derived class of DifferentialEquation
	//This derived class will have already destructed so attempting to copy over data to it is an invalid action and can crash the program.
	//This can happen when a mesh is deleted with CUDA switched on
	//It's also possible this destructor was called simply due to switching CUDA off, in which case called_from_destructor should be false
	CleanupMemory(!pameshODE->called_from_destructor);
}

//---------------------------------------- SET-UP METHODS  : DiffEqCUDA.cpp and DiffEqCUDA.cu

//allocate memory depending on set evaluation method - also cleans up previously allocated memory by calling CleanupMemory();
BError Atom_DifferentialEquationCubicCUDA::AllocateMemory(bool copy_from_cpu)
{
	BError error(CLASS_STR(Atom_DifferentialEquationCubicCUDA));

	//first make sure everything not needed is cleared
	CleanupMemory();

	if (!sM1()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	//copy values from cpu : it's possible the user switches to CUDA during a simulation and expects it to continue as normal
	//some evaluation methods need saved data (sM1, sEval... etc.) to carry over.
	else if (copy_from_cpu) sM1()->copy_from_cpuvec(pameshODE->sM1);

	switch (pameshODE->evalMethod) {

	case EVAL_EULER:
		break;

	case EVAL_AHEUN:
	case EVAL_TEULER:
		break;

	case EVAL_RK4:
		if (!sEval0()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval0()->copy_from_cpuvec(pameshODE->sEval0);

		if (!sEval1()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval1()->copy_from_cpuvec(pameshODE->sEval1);

		if (!sEval2()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval2()->copy_from_cpuvec(pameshODE->sEval2);
		break;

	case EVAL_ABM:
		if (!sEval0()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval0()->copy_from_cpuvec(pameshODE->sEval0);

		if (!sEval1()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval1()->copy_from_cpuvec(pameshODE->sEval1);
		break;

	case EVAL_RK23:
		if (!sEval0()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval0()->copy_from_cpuvec(pameshODE->sEval0);

		if (!sEval1()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval1()->copy_from_cpuvec(pameshODE->sEval1);

		if (!sEval2()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval2()->copy_from_cpuvec(pameshODE->sEval2);
		break;

	case EVAL_RKF:
		if (!sEval0()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if(copy_from_cpu) sEval0()->copy_from_cpuvec(pameshODE->sEval0);

		if (!sEval1()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval1()->copy_from_cpuvec(pameshODE->sEval1);

		if (!sEval2()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval2()->copy_from_cpuvec(pameshODE->sEval2);

		if (!sEval3()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval3()->copy_from_cpuvec(pameshODE->sEval3);

		if (!sEval4()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval4()->copy_from_cpuvec(pameshODE->sEval4);
		break;

	case EVAL_RKCK:
		if (!sEval0()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval0()->copy_from_cpuvec(pameshODE->sEval0);

		if (!sEval1()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval1()->copy_from_cpuvec(pameshODE->sEval1);

		if (!sEval2()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval2()->copy_from_cpuvec(pameshODE->sEval2);

		if (!sEval3()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval3()->copy_from_cpuvec(pameshODE->sEval3);

		if (!sEval4()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval4()->copy_from_cpuvec(pameshODE->sEval4);
		break;

	case EVAL_RKDP:
		if (!sEval0()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval0()->copy_from_cpuvec(pameshODE->sEval0);

		if (!sEval1()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval1()->copy_from_cpuvec(pameshODE->sEval1);

		if (!sEval2()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval2()->copy_from_cpuvec(pameshODE->sEval2);

		if (!sEval3()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval3()->copy_from_cpuvec(pameshODE->sEval3);

		if (!sEval4()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval4()->copy_from_cpuvec(pameshODE->sEval4);

		if (!sEval5()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval5()->copy_from_cpuvec(pameshODE->sEval5);
		break;

	case EVAL_SD:
		if (!sEval0()->resize((cuSZ3)paMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) sEval0()->copy_from_cpuvec(pameshODE->sEval0);
		break;
	}

	//For stochastic equations must also allocate memory for thermal VECs

	bool prng_used = false;

	switch (pameshODE->setODE) {

	case ODE_SLLG:
	case ODE_SLLGSTT:
	case ODE_SLLGSA:
		if (!H_Thermal()->resize((cuReal3)paMesh->h, (cuRect)paMesh->meshRect)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
		else if (copy_from_cpu) H_Thermal()->copy_from_cpuvec(pameshODE->H_Thermal);
		prng_used = true;
		break;
	}

	if (prng_used) {

		//initialize the pseudo-random number generator with a seed and memory size - recommended use kernel size divided by 128
		if (prng()->initialize(GetTickCount(), paMesh->n.dim() / 128) != cudaSuccess) error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}

	return error;
}

//deallocate memory before re-allocating it (depending on evaluation method previously allocated memory might not be used again, so need clean-up before)
void Atom_DifferentialEquationCubicCUDA::CleanupMemory(bool copy_to_cpu)
{
	//copy values to cpu before erasing : it's possible the user switches CUDA off during a simulation and expects it to continue as normal
	//some evaluation methods need saved data (sM1, sEval... etc.) to carry over.
	//CleanupMemory will be called by destructor in this case, so before destroying gpu data copy it over to cpu if possible
	//CleanupMemory may also be called in other circumstances, in particular from the cpu version of CleanupMemory, after having cleaned cpu vecs, thus in this case the copy methods will not run.

	//Only clear vectors not used for current evaluation method
	if (copy_to_cpu && sM1()->size_cpu() == pameshODE->sM1.size()) sM1()->copy_to_cpuvec(pameshODE->sM1);
	sM1()->clear();

	if (
		pameshODE->evalMethod != EVAL_RK4 &&
		pameshODE->evalMethod != EVAL_ABM &&
		pameshODE->evalMethod != EVAL_RK23 &&
		pameshODE->evalMethod != EVAL_RKF &&
		pameshODE->evalMethod != EVAL_RKCK &&
		pameshODE->evalMethod != EVAL_RKDP &&
		pameshODE->evalMethod != EVAL_SD) {

		if (copy_to_cpu && sEval0()->size_cpu() == pameshODE->sEval0.size()) sEval0()->copy_to_cpuvec(pameshODE->sEval0);
		sEval0()->clear();
	}

	if (pameshODE->evalMethod != EVAL_RK4 &&
		pameshODE->evalMethod != EVAL_ABM &&
		pameshODE->evalMethod != EVAL_RK23 &&
		pameshODE->evalMethod != EVAL_RKF &&
		pameshODE->evalMethod != EVAL_RKCK &&
		pameshODE->evalMethod != EVAL_RKDP) {

		if (copy_to_cpu && sEval1()->size_cpu() == pameshODE->sEval1.size()) sEval1()->copy_to_cpuvec(pameshODE->sEval1);
		sEval1()->clear();
	}

	if (pameshODE->evalMethod != EVAL_RK4 &&
		pameshODE->evalMethod != EVAL_RK23 &&
		pameshODE->evalMethod != EVAL_RKF &&
		pameshODE->evalMethod != EVAL_RKCK &&
		pameshODE->evalMethod != EVAL_RKDP) {

		if (copy_to_cpu && sEval2()->size_cpu() == pameshODE->sEval2.size()) sEval2()->copy_to_cpuvec(pameshODE->sEval2);
		sEval2()->clear();
	}

	if (pameshODE->evalMethod != EVAL_RK4 &&
		pameshODE->evalMethod != EVAL_RKF &&
		pameshODE->evalMethod != EVAL_RKCK &&
		pameshODE->evalMethod != EVAL_RKDP) {

		if (copy_to_cpu && sEval3()->size_cpu() == pameshODE->sEval3.size()) sEval3()->copy_to_cpuvec(pameshODE->sEval3);
		sEval3()->clear();

		if (copy_to_cpu && sEval4()->size_cpu() == pameshODE->sEval4.size()) sEval4()->copy_to_cpuvec(pameshODE->sEval4);
		sEval4()->clear();
	}

	if (pameshODE->evalMethod != EVAL_RKDP) {

		if (copy_to_cpu && sEval5()->size_cpu() == pameshODE->sEval5.size()) sEval5()->copy_to_cpuvec(pameshODE->sEval5);
		sEval5()->clear();
	}

	//For thermal vecs only clear if not used for current set ODE
	if (pameshODE->setODE != ODE_SLLG &&
		pameshODE->setODE != ODE_SLLGSTT &&
		pameshODE->setODE != ODE_SLLGSA) {

		if (copy_to_cpu && H_Thermal()->size_cpu() == pameshODE->H_Thermal.size()) H_Thermal()->copy_to_cpuvec(pameshODE->H_Thermal);
		H_Thermal()->clear();
	}

	prng()->clear();
}


BError Atom_DifferentialEquationCubicCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_DifferentialEquationCubicCUDA));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHDELETED)) {

		//if a mesh is deleted then a Atom_DifferentialEquationCUDA object can be deleted
		//this results in deletion of static data in Atom_ODECommonCUDA
		//whilst the static data is remade by UpdateConfiguration in Atom_ODECommonCUDA following this, our ManagedAtom_DiffEqCubicCUDA object now has pointers which are not linked correctly, so need to update them
		cuaDiffEq()->set_pointers(this);
	}

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_ODE_SOLVER)) {

		error = AllocateMemory();
	}

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_ODE_MOVEMESH)) {

		if (!error) {

			//set skip cells flags for moving mesh if enabled
			if (pameshODE->moving_mesh) {

				Rect mesh_rect = paMesh->GetMeshRect();

				DBL3 end_size = mesh_rect.size() & DBL3(MOVEMESH_ENDRATIO, 1, 1);

				Rect end_rect_left = Rect(mesh_rect.s, mesh_rect.s + end_size);
				Rect end_rect_right = Rect(mesh_rect.e - end_size, mesh_rect.e);

				paMeshCUDA->M1()->set_skipcells((cuRect)end_rect_left);
				paMeshCUDA->M1()->set_skipcells((cuRect)end_rect_right);
			}
			else {

				paMeshCUDA->M1()->clear_skipcells();
			}
		}
	}

	return error;
}

#endif
#endif