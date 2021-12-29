#include "stdafx.h"
#include "Mesh_Ferromagnetic.h"

#ifdef MESH_COMPILATION_FERROMAGNETIC

#include "SuperMesh.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

FMesh::FMesh(SuperMesh *pSMesh_) :
	Mesh(MESH_FERROMAGNETIC, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId), 
			VINFO(displayedPhysicalQuantity), VINFO(displayedBackgroundPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), VINFO(Module_Heff_Display), VINFO(Module_Energy_Display),
			VINFO(meshRect), VINFO(n), VINFO(h), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(n_m), VINFO(h_m), VINFO(n_s), VINFO(h_s), VINFO(link_stochastic),
			VINFO(M), 
			VINFO(V), VINFO(E), VINFO(S), VINFO(elC), 
			VINFO(Temp), VINFO(Temp_l), 
			VINFO(u_disp), VINFO(strain_diag), VINFO(strain_odiag),
			VINFO(pMod), 
			VINFO(exclude_from_multiconvdemag),
			//Members in this derived class
			VINFO(move_mesh_trigger), VINFO(skyShift), VINFO(exchange_couple_to_meshes),
			VINFO(mc_cone_angledeg), VINFO(mc_acceptance_rate), VINFO(mc_parallel), VINFO(mc_disabled), VINFO(mc_constrain), VINFO(cmc_n),
			//Material Parameters
			VINFO(grel), VINFO(alpha), VINFO(Ms), VINFO(Nxy), 
			VINFO(A), VINFO(D), VINFO(D_dir), VINFO(J1), VINFO(J2),
			VINFO(K1), VINFO(K2), VINFO(K3), VINFO(mcanis_ea1), VINFO(mcanis_ea2), VINFO(mcanis_ea3),
			VINFO(Kt),
			VINFO(susrel), VINFO(susprel), VINFO(cHA), VINFO(cHmo),
			VINFO(elecCond), VINFO(amrPercentage), VINFO(P), VINFO(beta), VINFO(De), VINFO(n_density), 
			VINFO(SHA), VINFO(flSOT), VINFO(STq), VINFO(STa), VINFO(STp), 
			VINFO(betaD), VINFO(l_sf), VINFO(l_ex), VINFO(l_ph), VINFO(Gi), VINFO(Gmix),
			VINFO(ts_eff), VINFO(tsi_eff), VINFO(pump_eff), VINFO(cpump_eff), VINFO(the_eff),
			VINFO(base_temperature), VINFO(T_equation), VINFO(T_Curie), VINFO(T_Curie_material), VINFO(atomic_moment), 
			VINFO(density), VINFO(MEc), VINFO(Ym), VINFO(Pr),
			VINFO(thermCond), VINFO(shc), VINFO(shc_e), VINFO(G_e), VINFO(cT), VINFO(Q),
			
			//OBSOLETE - must keep to allow older simulation files to load if they have these defined
			VINFO(lambda), VINFO(cTsig), VINFO(dTsig)
		},
		{
			//Modules Implementations
			IINFO(Demag_N), IINFO(Demag), IINFO(SDemag_Demag),
			IINFO(Exch_6ngbr_Neu), IINFO(DMExchange), IINFO(iDMExchange), IINFO(viDMExchange), IINFO(SurfExchange),
			IINFO(Zeeman), IINFO(MOptical), IINFO(MElastic), IINFO(Roughness),
			IINFO(Anisotropy_Uniaxial), IINFO(Anisotropy_Cubic), IINFO(Anisotropy_Biaxial), IINFO(Anisotropy_Tensorial),
			IINFO(Transport), IINFO(Heat),
			IINFO(SOTField), IINFO(STField)
		}),
	meshODE(this)
{}

FMesh::FMesh(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_) :
	Mesh(MESH_FERROMAGNETIC, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId),
			VINFO(displayedPhysicalQuantity), VINFO(displayedBackgroundPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), VINFO(Module_Heff_Display), VINFO(Module_Energy_Display),
			VINFO(meshRect), VINFO(n), VINFO(h), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(n_m), VINFO(h_m), VINFO(n_s), VINFO(h_s), VINFO(link_stochastic),
			VINFO(M), 
			VINFO(V), VINFO(E), VINFO(S), VINFO(elC), 
			VINFO(Temp), VINFO(Temp_l), 
			VINFO(u_disp), VINFO(strain_diag), VINFO(strain_odiag),
			VINFO(pMod),
			VINFO(exclude_from_multiconvdemag),
			//Members in this derived class
			VINFO(move_mesh_trigger), VINFO(skyShift), VINFO(exchange_couple_to_meshes),
			VINFO(mc_cone_angledeg), VINFO(mc_acceptance_rate), VINFO(mc_parallel), VINFO(mc_disabled), VINFO(mc_constrain), VINFO(cmc_n),
			//Material Parameters
			VINFO(grel), VINFO(alpha), VINFO(Ms), VINFO(Nxy),
			VINFO(A), VINFO(D), VINFO(D_dir), VINFO(J1), VINFO(J2),
			VINFO(K1), VINFO(K2), VINFO(K3), VINFO(mcanis_ea1), VINFO(mcanis_ea2), VINFO(mcanis_ea3),
			VINFO(Kt),
			VINFO(susrel), VINFO(susprel), VINFO(cHA), VINFO(cHmo),
			VINFO(elecCond), VINFO(amrPercentage), VINFO(P), VINFO(beta), VINFO(De), VINFO(n_density), 
			VINFO(SHA), VINFO(flSOT), VINFO(STq), VINFO(STa), VINFO(STp), 
			VINFO(betaD), VINFO(l_sf), VINFO(l_ex), VINFO(l_ph), VINFO(Gi), VINFO(Gmix),
			VINFO(ts_eff), VINFO(tsi_eff), VINFO(pump_eff), VINFO(cpump_eff), VINFO(the_eff),
			VINFO(base_temperature), VINFO(T_equation), VINFO(T_Curie), VINFO(T_Curie_material), VINFO(atomic_moment),
			VINFO(density), VINFO(MEc), VINFO(Ym), VINFO(Pr),
			VINFO(thermCond), VINFO(shc), VINFO(shc_e), VINFO(G_e), VINFO(cT), VINFO(Q),
			
			//OBSOLETE - must keep to allow older simulation files to load if they have these defined
			VINFO(lambda), VINFO(cTsig), VINFO(dTsig)
		},
		{
			//Modules Implementations
			IINFO(Demag_N), IINFO(Demag), IINFO(SDemag_Demag),
			IINFO(Exch_6ngbr_Neu), IINFO(DMExchange), IINFO(iDMExchange), IINFO(viDMExchange), IINFO(SurfExchange),
			IINFO(Zeeman), IINFO(MOptical), IINFO(MElastic), IINFO(Roughness),
			IINFO(Anisotropy_Uniaxial), IINFO(Anisotropy_Cubic), IINFO(Anisotropy_Biaxial), IINFO(Anisotropy_Tensorial),
			IINFO(Transport), IINFO(Heat),
			IINFO(SOTField), IINFO(STField)
		}),
	meshODE(this)
{
	//default settings
	displayedPhysicalQuantity = MESHDISPLAY_MAGNETIZATION;

	meshRect = meshRect_;
	
	h = h_;
	h_e = h_;
	h_t = h_;
	h_m = h_;
	h_s = h_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//default modules configuration
	if (!error_on_create) error_on_create = AddModule(MOD_DEMAG);
	if (!error_on_create) error_on_create = AddModule(MOD_EXCHANGE);
	if (!error_on_create) error_on_create = AddModule(MOD_ZEEMAN);

	//--------------------------
	
	//If cuda is enabled we also need to make the cuda mesh version
	if (cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

FMesh::~FMesh()
{
}

void FMesh::RepairObjectState(void) 
{
	//at this point Heff is empty and must not be since MComputation_Enabled will report wrong
	Heff.assign(h, meshRect, DBL3(0, 0, 0));

	//also need to re-calculate special functions used in text formulas as they might not correspond to the loaded parameters now
	
	//calculate scaling function (Curie Weiss law)
	if (T_Curie > 0) {

		DBL3 Ha = CallModuleMethod(&ZeemanBase::GetField);
		pCurieWeiss->Initialize_CurieWeiss(MU0 * (MUB / BOLTZMANN) * atomic_moment * Ha.norm(), T_Curie);
		pLongRelSus->Initialize_LongitudinalRelSusceptibility(pCurieWeiss->get_data(), atomic_moment, T_Curie);
		pAlpha1->Initialize_Alpha1(DBL2(1.0, 0.0), 0.0);
	}
	else {

		pCurieWeiss->Initialize_CurieWeiss(0.0, 1.0);
		pLongRelSus->Initialize_LongitudinalRelSusceptibility(pCurieWeiss->get_data(), atomic_moment, 1.0);
		pAlpha1->Initialize_Alpha1(DBL2(1.0, 0.0), 0.0);
	}

	//Make sure Zeeman module is always the first one in the list : Zeeman module sets Heff (if Zeeman module disabled then PrepareIteration clears Heff)
	if (IsModuleSet(MOD_ZEEMAN)) {

		int idxZeeman = pMod.get_index_from_ID(MOD_ZEEMAN);
		if (idxZeeman != 0) pMod.move(idxZeeman);
	}
}

//----------------------------------- IMPORTANT CONTROL METHODS

//call when a configuration change has occurred - some objects might need to be updated accordingly
BError FMesh::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(std::string(CLASS_STR(FMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

	///////////////////////////////////////////////////////
	//Mesh specific configuration
	///////////////////////////////////////////////////////

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

		n = round(meshRect / h);
		if (n.x == 0) n.x = 1;
		if (n.y == 0) n.y = 1;
		if (n.z == 0) n.z = 1;
		h = meshRect / n;

		//resize arrays held in Mesh - if changing mesh to new size then use resize : this maps values from old mesh to new mesh size (keeping magnitudes). Otherwise assign magnetization along x in whole mesh.
		if (M.linear_size()) {

			if (!M.resize(h, meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
		}
		else if (!M.assign(h, meshRect, DBL3(-Ms, 0, 0))) return error(BERROR_OUTOFMEMORY_CRIT);

		if (!Heff.assign(h, meshRect, DBL3(0, 0, 0))) return error(BERROR_OUTOFMEMORY_CRIT);

		//update material parameters spatial dependence as cellsize and rectangle could have changed
		if (!error && !update_meshparam_var()) error(BERROR_OUTOFMEMORY_NCRIT);

		//update any text equations used in mesh parameters (dependence on mesh dimensions possible)
		if (!error) update_all_meshparam_equations();
	}

	//erase any unused skyrmion trackers in this mesh
	skyShift.UpdateConfiguration(saveDataList);

	if (cfgMessage == UPDATECONFIG_PARAMVALUECHANGED_MLENGTH) {

		//renormalization of magnetization length done in meshODE. Here we need to check Ms has not been set to zero for non-empty cells if a spatial variation is set.
		//This can happen e.g. when shape_setparam is used for the first time.
		Ms.fix_s_scaling_zeros(M, Ms.get0());
	}

	//------------------------ CUDA UpdateConfiguration if set
	
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (!error) error = pMeshCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	//------------------------
	
	///////////////////////////////////////////////////////
	//Update configuration for mesh modules
	///////////////////////////////////////////////////////

	for (int idx = 0; idx < (int)pMod.size(); idx++) {
		
		if (pMod[idx] && !error) {

			error = pMod[idx]->UpdateConfiguration(cfgMessage);
		}
	}

	///////////////////////////////////////////////////////
	//Update configuration for mesh ode solver
	///////////////////////////////////////////////////////

	if (!error) error = meshODE.UpdateConfiguration(cfgMessage);

	return error;
}

void FMesh::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	///////////////////////////////////////////////////////
	//Update configuration in this mesh
	///////////////////////////////////////////////////////

	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		//Update text equations other than those used in mesh parameters
		UpdateTEquationUserConstants();

		//update any text equations used in mesh parameters
		update_all_meshparam_equations();
	}
	else if (cfgMessage == UPDATECONFIG_TEQUATION_CLEAR) {

		if (T_equation.is_set()) T_equation.clear();
	}

	///////////////////////////////////////////////////////
	//Update configuration for mesh modules
	///////////////////////////////////////////////////////

	for (int idx = 0; idx < (int)pMod.size(); idx++) {

		if (pMod[idx]) {

			pMod[idx]->UpdateConfiguration_Values(cfgMessage);
		}
	}

	///////////////////////////////////////////////////////
	//Update configuration for mesh ode solver
	///////////////////////////////////////////////////////

	meshODE.UpdateConfiguration_Values(cfgMessage);

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pMeshCUDA) {

		pMeshCUDA->UpdateConfiguration_Values(cfgMessage);
	}
#endif
}

BError FMesh::SwitchCUDAState(bool cudaState)
{
	BError error(std::string(CLASS_STR(FMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

#if COMPILECUDA == 1

	//--------------------------------------------

	//are we switching to cuda?
	if (cudaState) {

		if (!pMeshBaseCUDA) {

			//then make MeshCUDA object, copying over currently held cpu data
			pMeshBaseCUDA = new FMeshCUDA(this);
			pMeshCUDA = dynamic_cast<MeshCUDA*>(pMeshBaseCUDA);

			error = pMeshBaseCUDA->Error_On_Create();
			if (!error) error = pMeshCUDA->cuMesh()->set_pointers(pMeshCUDA);
		}
	}
	else {

		//delete MeshCUDA object and null
		if (pMeshBaseCUDA) delete pMeshBaseCUDA;
		pMeshBaseCUDA = nullptr;
		pMeshCUDA = nullptr;
	}

	//--------------------------------------------

	//SwitchCUDA state for differential equation
	if (!error) error = meshODE.SwitchCUDAState(cudaState);

	//--------------------------------------------

	//SwitchCUDA state for all active modules in this mesh
	for (int idx = 0; idx < (int)pMod.size(); idx++) {

		if (pMod[idx] && !error) error = pMod[idx]->SwitchCUDAState(cudaState);
	}

#endif

	return error;
}

//couple this mesh to touching dipoles by setting skip cells as required : used for domain wall moving mesh algorithms
void FMesh::CoupleToDipoles(bool status)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) pMeshCUDA->M()->copy_to_cpuvec(M);
#endif

	//check all meshes to see if there are any dipole meshes touching this mesh
	for (int idx = 0; idx < pSMesh->size(); idx++) {

		if ((*pSMesh)[idx]->GetMeshType() == MESH_DIPOLE) {

			Rect mesh_intersection = (*pSMesh)[idx]->meshRect.get_intersection(meshRect);

			//for this dipole to "touch" this mesh, the rects intersection must be exactly a plane
			if (mesh_intersection.IsPlane()) {

				//y-z plane : x is the perpendicular direction
				if (IsZ(mesh_intersection.s.x - mesh_intersection.e.x)) {

					//is the selected mesh on the positive or negative side of the boundary? Set flag to use. Also adjust mesh_intersection so it contains all the cells to set flags for.
					if (IsZ(meshRect.s.x - mesh_intersection.s.x)) {

						mesh_intersection.e.x += h.x;
					}
					else {

						mesh_intersection.s.x -= h.x;
					}
				}
				//x-z plane : y is the perpendicular direction
				else if (IsZ(mesh_intersection.s.y - mesh_intersection.e.y)) {

					//is the selected mesh on the positive or negative side of the boundary? Set flag to use. Also adjust mesh_intersection so it contains all the cells to set flags for.
					if (IsZ(meshRect.s.y - mesh_intersection.s.y)) {

						mesh_intersection.e.y += h.y;
					}
					else {

						mesh_intersection.s.y -= h.y;
					}
				}
				//x-y plane : z is the perpendicular direction
				else if (IsZ(mesh_intersection.s.z - mesh_intersection.e.z)) {

					//is the selected mesh on the positive or negative side of the boundary? Set flag to use. Also adjust mesh_intersection so it contains all the cells to set flags for.
					if (IsZ(meshRect.s.z - mesh_intersection.s.z)) {

						mesh_intersection.e.z += h.z;
					}
					else {

						mesh_intersection.s.z -= h.z;
					}
				}

				//found a touching dipole - mark all cells in the intersection as skip cells
				M.set_skipcells(mesh_intersection, status);

				if (status) {

					//set interface cells to have magnetization direction along the touching dipole direction
					DBL3 Mdipole_direction = dynamic_cast<Mesh*>((*pSMesh)[idx])->GetAverageMagnetization().normalized();

					Box box = M.box_from_rect_max(mesh_intersection);

					for (int i = box.s.i; i < box.e.i; i++) {
						for (int j = box.s.j; j < box.e.j; j++) {
							for (int k = box.s.k; k < box.e.k; k++) {

								double Mag = M[INT3(i, j, k)].norm();
								M[INT3(i, j, k)] = Mag * Mdipole_direction;
							}
						}
					}
				}
			}
		}
	}

#if COMPILECUDA == 1
	if (pMeshCUDA) pMeshCUDA->M()->copy_from_cpuvec(M);
#endif
}

double FMesh::CheckMoveMesh(void)
{

#if COMPILECUDA == 1
	if (pMeshCUDA) return pMeshCUDA->CheckMoveMesh(meshODE.MoveMeshAntisymmetric(), meshODE.MoveMeshThreshold());
#endif

	//move mesh algorithm applied to systems containing domain walls in order to simulate domain wall movement.
	//In this case the magnetization at the ends of the mesh must be fixed (magnetization evolver not applied to the ends).
	//Find the average magnetization projected along the absolute direction of the magnetization at the ends.
	//For a domain wall at the centre this value should be close to zero. If it exceeds a threshold then move mesh by a cellsize along +x or -x.

	if (!move_mesh_trigger) return 0;

	//the ends should not be completely empty, and must have a constant magnetization direction
	int cells_fixed = (int)ceil_epsilon((double)n.x * MOVEMESH_ENDRATIO);

	DBL3 M_end = M.average_omp(Box(0, 0, 0, cells_fixed, n.y, n.z));
	DBL3 M_av_left = M.average_omp(Box(cells_fixed, 0, 0, n.x / 2, n.y, n.z));
	DBL3 M_av_right = M.average_omp(Box(n.x / 2, 0, 0, n.x - cells_fixed, n.y, n.z));

	double direction = (2 * double(meshODE.MoveMeshAntisymmetric()) - 1);

	DBL3 M_av = M_av_right + direction * M_av_left;

	if (GetMagnitude(M_end) * GetMagnitude(M_av)) {

		double Mratio = (M_end * M_av) / (GetMagnitude(M_end) * GetMagnitude(M_av));

		if (Mratio > meshODE.MoveMeshThreshold()) return -h.x * direction;
		if (Mratio < -meshODE.MoveMeshThreshold()) return h.x * direction;
	}

	return 0;
}

//----------------------------------- OVERLOAD MESH VIRTUAL METHODS

void FMesh::SetCurieTemperature(double Tc, bool set_default_dependences)
{
	//Curie temperature is a constant in mesh parameter equations, so update them
	if (Tc != T_Curie) update_all_meshparam_equations();

	if (Tc > 0) {

		T_Curie = Tc;

		//calculate scaling function (Curie Weiss law)
		DBL3 Ha = CallModuleMethod(&ZeemanBase::GetField);
		pCurieWeiss->Initialize_CurieWeiss(MU0 * (MUB / BOLTZMANN) * atomic_moment * Ha.norm(), T_Curie);
		pLongRelSus->Initialize_LongitudinalRelSusceptibility(pCurieWeiss->get_data(), atomic_moment, T_Curie);

		pAlpha1->Initialize_Alpha1(DBL2(1.0, 0.0), 0.0);

		if (set_default_dependences) {

			//set default temperature dependences
			Ms.set_t_scaling_equation(std::string("me(T/Tc)"), userConstants, T_Curie, base_temperature);
			P.set_t_scaling_equation(std::string("me(T/Tc)"), userConstants, T_Curie, base_temperature);
			A.set_t_scaling_equation(std::string("me(T/Tc)^2"), userConstants, T_Curie, base_temperature);
			D.set_t_scaling_equation(std::string("me(T/Tc)^2"), userConstants, T_Curie, base_temperature);
			K1.set_t_scaling_equation(std::string("me(T/Tc)^3"), userConstants, T_Curie, base_temperature);
			K2.set_t_scaling_equation(std::string("me(T/Tc)^3"), userConstants, T_Curie, base_temperature);

			susrel.set_t_scaling_equation(std::string("chi(T/Tc)"), userConstants, T_Curie, base_temperature);

			alpha.set_t_scaling_equation(std::string("alpha1(T/Tc)"), userConstants, T_Curie, base_temperature);
		}

		//make sure to also update them - this method can be called during a simulation, e.g. if field changes.
		Ms.update(base_temperature);
		P.update(base_temperature);
		A.update(base_temperature);
		D.update(base_temperature);
		K1.update(base_temperature);
		K2.update(base_temperature);
		susrel.update(base_temperature);
		alpha.update(base_temperature);
	}
	else {

		//turn it off
		T_Curie = 0.0;

		if (set_default_dependences) {

			//reset temperature dependencies for affected parameters
			Ms.clear_t_scaling();
			P.clear_t_scaling();
			A.clear_t_scaling();
			D.clear_t_scaling();
			K1.clear_t_scaling();
			K2.clear_t_scaling();
			susrel.clear_t_scaling();
			alpha.clear_t_scaling();
		}
	}

#if COMPILECUDA == 1
	if (pMeshCUDA) {

		//T_Curie changed : sync with cuda version
		pMeshCUDA->T_Curie.from_cpu((cuBReal)T_Curie);

		pMeshCUDA->set_special_functions_data();
	}
#endif
}

void FMesh::SetAtomicMoment(DBL2 atomic_moment_ub)
{
	atomic_moment = atomic_moment_ub.i;

	//setting atomic_moment will affect the temperature dependence of me (normalised equilibrium magnetization), so some parameter temperature dependencies must change - calling SetCurieTemperature with the current Curie temperature will do this
	SetCurieTemperature(T_Curie, false);
}

#endif