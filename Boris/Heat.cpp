#include "stdafx.h"
#include "Heat.h"

#ifdef MODULE_COMPILATION_HEAT

#include "Mesh.h"
#include "SuperMesh.h"

Heat::Heat(Mesh *pMesh_) :
	Modules(),
	pMesh(pMesh_),
	pSMesh(pMesh_->pSMesh),
	HeatBase(pMesh_),
	ProgramStateNames(this, 
		{ 
			VINFO(tmtype),
			VINFO(T_ambient), VINFO(alpha_boundary), 
			VINFO(insulate_px), VINFO(insulate_nx), VINFO(insulate_py), VINFO(insulate_ny), VINFO(insulate_pz), VINFO(insulate_nz),
			VINFO(Q_equation)
		}, {})
{
	//this needs to go here, not in HeatBase, as it calls UpdateConfiguration before all pointers have been set
	error_on_create = Set_TMType();

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

Heat::~Heat()
{
	//free memory used for heat solver
	pMesh->Temp.clear();
	pMesh->Temp_l.clear();
}

//-------------------Abstract base class method implementations

BError Heat::Initialize(void)
{
	BError error(CLASS_STR(Heat));

	initialized = true;

	SetRobinBoundaryConditions();

	return error;
}

BError Heat::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Heat));

	Uninitialize();

	bool success = true;

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MODULEADDED, UPDATECONFIG_MODULEDELETED, UPDATECONFIG_HEAT_MODELTYPE)) {

		//make sure the cellsize divides the mesh rectangle
		pMesh->n_t = round(pMesh->meshRect / pMesh->h_t);
		if (pMesh->n_t.x < 2) pMesh->n_t.x = 2;
		if (pMesh->n_t.y < 2) pMesh->n_t.y = 2;
		if (pMesh->n_t.z < 2) pMesh->n_t.z = 2;
		pMesh->h_t = pMesh->meshRect / pMesh->n_t;

		//update mesh dimensions in equation constants
		if (Q_equation.is_set()) {

			DBL3 meshDim = pMesh->GetMeshDimensions();

			Q_equation.set_constant("Lx", meshDim.x, false);
			Q_equation.set_constant("Ly", meshDim.y, false);
			Q_equation.set_constant("Lz", meshDim.z, false);
			Q_equation.remake_equation();
		}

		//make sure correct memory is assigned for heat solver quantities

		//first clear any VECs which are not needed
		switch (tmtype) {

		case TMTYPE_1TM:
			//single temperature only
			pMesh->Temp_l.clear();
			break;

		case TMTYPE_2TM:
			//itinerant electrons temperature <-> lattice temperature
			break;
		}

		//if we've just changed the model type, make sure any new temperature VECs start from a reasonable setting - same as the primary Temp VEC.
		bool initialize_Temp_l = !(pMesh->Temp_l.linear_size());

		//now allocate/resize required memory
		if (pMesh->M.linear_size()) {

			//in a ferromagnet the magnetization sets the shape only on initialization. If already initialized then shape already set.
			if (pMesh->Temp.linear_size()) {

				success = pMesh->Temp.resize(pMesh->h_t, pMesh->meshRect);

				//lattice temperature for many-temperature models
				if (tmtype == TMTYPE_2TM) success &= pMesh->Temp_l.resize(pMesh->h_t, pMesh->meshRect);
			}
			else {

				success = pMesh->Temp.assign(pMesh->h_t, pMesh->meshRect, pMesh->base_temperature, pMesh->M);

				//lattice temperature for many-temperature models
				if (tmtype == TMTYPE_2TM) success &= pMesh->Temp_l.assign(pMesh->h_t, pMesh->meshRect, pMesh->base_temperature, pMesh->Temp);
			}

		}
		else if (pMesh->elC.linear_size()) {

			//in a normal metal the electrical conductivity sets the shape

			//before doing this must make sure elC was itself set in the Transport module (it could be this module is being updated before the Transport module)
			//Transport module is guaranteed to be set otherwise elC would have zero size - it does mean Transport has UpdateConfiguration called twice but it doesn't matter.
			if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

				error = pMesh->CallModuleMethod(&Transport::UpdateConfiguration, cfgMessage);
				if (error) return error;
			}
			else if (pMesh->IsModuleSet(MOD_TMR)) {

				error = pMesh->CallModuleMethod(&TMR::UpdateConfiguration, cfgMessage);
				if (error) return error;
			}

			if (success) {

				if (pMesh->Temp.linear_size()) {

					success = pMesh->Temp.resize(pMesh->h_t, pMesh->meshRect, pMesh->elC);

					//lattice temperature for many-temperature models
					if (tmtype == TMTYPE_2TM) success &= pMesh->Temp_l.resize(pMesh->h_t, pMesh->meshRect, pMesh->Temp);
				}
				else {

					success = pMesh->Temp.assign(pMesh->h_t, pMesh->meshRect, pMesh->base_temperature, pMesh->elC);

					//lattice temperature for many-temperature models
					if (tmtype == TMTYPE_2TM) success &= pMesh->Temp_l.assign(pMesh->h_t, pMesh->meshRect, pMesh->base_temperature, pMesh->Temp);
				}
			}
		}
		else {

			//in an insulator (or conductor with Transport module not set) the shape is set directly in Temp
			if (pMesh->Temp.linear_size()) {

				success = pMesh->Temp.resize(pMesh->h_t, pMesh->meshRect);

				//lattice temperature for many-temperature models
				if (tmtype == TMTYPE_2TM) success &= pMesh->Temp_l.resize(pMesh->h_t, pMesh->meshRect);
			}
			else {

				success = pMesh->Temp.assign(pMesh->h_t, pMesh->meshRect, pMesh->base_temperature);
			}
		}

		//if we've just changed the model type, make sure any new temperature VECs start from a reasonable setting - same as the primary Temp VEC.
		if (initialize_Temp_l && pMesh->Temp_l.linear_size()) pMesh->Temp_l.copy_values(pMesh->Temp);

		//allocate memory for the heatEq_RHS auxiliary vector
		if (success) success = malloc_vector(heatEq_RHS, pMesh->n_t.dim(), 0.0);
	}

	if (!success) return error(BERROR_OUTOFMEMORY_CRIT);

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

void Heat::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		UpdateTEquationUserConstants(false);
	}
	else if (cfgMessage == UPDATECONFIG_TEQUATION_CLEAR) {

		if (Q_equation.is_set()) Q_equation.clear();
	}

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pModuleCUDA->UpdateConfiguration_Values(cfgMessage);
	}
#endif
}

BError Heat::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Heat));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new HeatCUDA(pMesh, pSMesh, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

#endif