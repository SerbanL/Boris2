#include "stdafx.h"
#include "Atom_Mesh_Cubic.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "SuperMesh.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

Atom_Mesh_Cubic::Atom_Mesh_Cubic(SuperMesh *pSMesh_) :
	Atom_Mesh(MESH_ATOM_CUBIC, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId), 
			VINFO(displayedPhysicalQuantity), VINFO(displayedBackgroundPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), VINFO(Module_Heff_Display), VINFO(Module_Energy_Display),
			VINFO(meshRect), VINFO(n), VINFO(h), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(n_m), VINFO(h_m), VINFO(n_dm), VINFO(h_dm),
			VINFO(M1), 
			VINFO(V), VINFO(E), VINFO(S), VINFO(elC),
			VINFO(Temp), VINFO(Temp_l),
			VINFO(pMod), 
			VINFO(exclude_from_multiconvdemag),
			//Members in this derived class
			VINFO(move_mesh_trigger), VINFO(skyShift), VINFO(exchange_couple_to_meshes),
			VINFO(mc_cone_angledeg), VINFO(mc_acceptance_rate), VINFO(mc_parallel), VINFO(mc_disabled), VINFO(mc_constrain), VINFO(cmc_n),
			//Material Parameters
			VINFO(grel), VINFO(alpha), VINFO(mu_s), VINFO(Nxy),
			VINFO(J), VINFO(D), VINFO(D_dir), VINFO(Js), VINFO(Js2),
			VINFO(K1), VINFO(K2),  VINFO(K3), VINFO(mcanis_ea1), VINFO(mcanis_ea2), VINFO(mcanis_ea3),
			VINFO(Kt),
			VINFO(cHA), VINFO(cHmo),
			VINFO(s_eff),
			VINFO(elecCond), VINFO(amrPercentage), VINFO(P), VINFO(beta), VINFO(De), VINFO(n_density),
			VINFO(SHA), VINFO(flSOT), VINFO(STq), VINFO(STa), VINFO(STp),
			VINFO(betaD), VINFO(l_sf), VINFO(l_ex), VINFO(l_ph), VINFO(Gi), VINFO(Gmix),
			VINFO(ts_eff), VINFO(tsi_eff), VINFO(pump_eff), VINFO(cpump_eff), VINFO(the_eff),
			VINFO(base_temperature), VINFO(T_equation), 
			VINFO(density),
			VINFO(thermCond), VINFO(shc), VINFO(shc_e), VINFO(G_e), VINFO(cT), VINFO(Q)
		},
		{
			//Modules Implementations
			IINFO(Atom_Demag_N), IINFO(Atom_Demag), IINFO(Atom_DipoleDipole), 
			IINFO(Atom_Zeeman), IINFO(Atom_MOptical),
			IINFO(Atom_Exchange), IINFO(Atom_DMExchange), IINFO(Atom_iDMExchange), IINFO(Atom_viDMExchange),IINFO(Atom_SurfExchange),
			IINFO(Atom_Anisotropy_Uniaxial), IINFO(Atom_Anisotropy_Cubic), IINFO(Atom_Anisotropy_Biaxial), IINFO(Atom_Anisotropy_Tensorial),
			IINFO(Atom_Transport),
			IINFO(Atom_Heat),
			IINFO(Atom_SOTField), IINFO(Atom_STField),
			IINFO(StrayField_AtomMesh)
		}),
	meshODE(this)
{
}

Atom_Mesh_Cubic::Atom_Mesh_Cubic(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_) :
	Atom_Mesh(MESH_ATOM_CUBIC, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId),
			VINFO(displayedPhysicalQuantity), VINFO(displayedBackgroundPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), VINFO(Module_Heff_Display), VINFO(Module_Energy_Display),
			VINFO(meshRect), VINFO(n), VINFO(h), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(n_m), VINFO(h_m), VINFO(n_dm), VINFO(h_dm),
			VINFO(M1), 
			VINFO(V), VINFO(E), VINFO(S), VINFO(elC),
			VINFO(Temp), VINFO(Temp_l),
			VINFO(pMod),
			VINFO(exclude_from_multiconvdemag),
			//Members in this derived class
			VINFO(move_mesh_trigger), VINFO(skyShift), VINFO(exchange_couple_to_meshes),
			VINFO(mc_cone_angledeg), VINFO(mc_acceptance_rate), VINFO(mc_parallel), VINFO(mc_disabled), VINFO(mc_constrain), VINFO(cmc_n),
			//Material Parameters
			VINFO(grel), VINFO(alpha), VINFO(mu_s), VINFO(Nxy),
			VINFO(J), VINFO(D), VINFO(D_dir), VINFO(Js), VINFO(Js2),
			VINFO(K1), VINFO(K2),  VINFO(K3), VINFO(mcanis_ea1), VINFO(mcanis_ea2), VINFO(mcanis_ea3),
			VINFO(Kt),
			VINFO(cHA), VINFO(cHmo),
			VINFO(s_eff),
			VINFO(elecCond), VINFO(amrPercentage), VINFO(P), VINFO(beta), VINFO(De), VINFO(n_density),
			VINFO(SHA), VINFO(flSOT), VINFO(STq), VINFO(STa), VINFO(STp),
			VINFO(betaD), VINFO(l_sf), VINFO(l_ex), VINFO(l_ph), VINFO(Gi), VINFO(Gmix),
			VINFO(ts_eff), VINFO(tsi_eff), VINFO(pump_eff), VINFO(cpump_eff), VINFO(the_eff),
			VINFO(base_temperature), VINFO(T_equation),
			VINFO(density),
			VINFO(thermCond), VINFO(shc), VINFO(shc_e), VINFO(G_e), VINFO(cT), VINFO(Q)
		},
		{
			//Modules Implementations
			IINFO(Atom_Demag_N), IINFO(Atom_Demag), IINFO(Atom_DipoleDipole),
			IINFO(Atom_Zeeman), IINFO(Atom_MOptical),
			IINFO(Atom_Exchange), IINFO(Atom_DMExchange), IINFO(Atom_iDMExchange), IINFO(Atom_viDMExchange), IINFO(Atom_SurfExchange),
			IINFO(Atom_Anisotropy_Uniaxial), IINFO(Atom_Anisotropy_Cubic), IINFO(Atom_Anisotropy_Biaxial), IINFO(Atom_Anisotropy_Tensorial),
			IINFO(Atom_Transport),
			IINFO(Atom_Heat),
			IINFO(Atom_SOTField), IINFO(Atom_STField),
			IINFO(StrayField_AtomMesh)
		}),
	meshODE(this)
{
	//default settings
	displayedPhysicalQuantity = MESHDISPLAY_MOMENT;

	meshRect = meshRect_;
	
	h = h_;
	h_e = h_;
	h_t = h_;
	h_m = h_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//default modules configuration
	if (!error_on_create) error_on_create = AddModule(MOD_ATOM_DIPOLEDIPOLE);
	if (!error_on_create) error_on_create = AddModule(MOD_EXCHANGE);
	if (!error_on_create) error_on_create = AddModule(MOD_ZEEMAN);

	//--------------------------
	
	//If cuda is enabled we also need to make the cuda mesh version
	if (cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

Atom_Mesh_Cubic::~Atom_Mesh_Cubic()
{
}

void Atom_Mesh_Cubic::RepairObjectState(void)
{
	//at this point Heff1 is empty and must not be since MComputation_Enabled will report wrong
	Heff1.assign(h, meshRect, DBL3(0, 0, 0));

	//Make sure Zeeman module is always the first one in the list : Zeeman module sets Heff (if Zeeman module disabled then PrepareIteration clears Heff)
	if (IsModuleSet(MOD_ZEEMAN)) {

		int idxZeeman = pMod.get_index_from_ID(MOD_ZEEMAN);
		if (idxZeeman != 0) pMod.move(idxZeeman);
	}
}

//----------------------------------- IMPORTANT CONTROL METHODS

//call when a configuration change has occurred - some objects might need to be updated accordingly
BError Atom_Mesh_Cubic::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(std::string(CLASS_STR(Atom_Mesh_Cubic)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

	///////////////////////////////////////////////////////
	//Mesh specific configuration
	///////////////////////////////////////////////////////

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE)) {

		n = round(meshRect / h);
		if (n.x == 0) n.x = 1;
		if (n.y == 0) n.y = 1;
		if (n.z == 0) n.z = 1;
		h = meshRect / n;

		n_dm = round(meshRect / h_dm);
		if (n_dm.x == 0) n_dm.x = 1;
		if (n_dm.y == 0) n_dm.y = 1;
		if (n_dm.z == 0) n_dm.z = 1;
		h_dm = meshRect / n_dm;

		//resize arrays held in Mesh - if changing mesh to new size then use resize : this maps values from old mesh to new mesh size (keeping magnitudes). Otherwise assign magnetization along x in whole mesh.
		if (M1.linear_size()) {

			if (!M1.resize(h, meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
		}
		else if (!M1.assign(h, meshRect, DBL3(-mu_s, 0, 0))) return error(BERROR_OUTOFMEMORY_CRIT);

		if (!Heff1.assign(h, meshRect, DBL3(0, 0, 0))) return error(BERROR_OUTOFMEMORY_CRIT);

		//update material parameters spatial dependence as cellsize and rectangle could have changed
		if (!error && !update_meshparam_var()) error(BERROR_OUTOFMEMORY_NCRIT);

		//update any text equations used in mesh parameters (dependence on mesh dimensions possible)
		if (!error) update_all_meshparam_equations();
	}

	if (cfgMessage == UPDATECONFIG_PARAMVALUECHANGED_MLENGTH) {

		//renormalization of magnetization length done in meshODE. Here we need to check mu_s has not been set to zero for non-empty cells if a spatial variation is set.
		//This can happen e.g. when shape_setparam is used for the first time.
		mu_s.fix_s_scaling_zeros(M1, mu_s.get0());
	}

	//------------------------ CUDA UpdateConfiguration if set
	
#if COMPILECUDA == 1
	if (paMeshCUDA) {

		if (!error) error = paMeshCUDA->UpdateConfiguration(cfgMessage);
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

void Atom_Mesh_Cubic::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
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
	if (paMeshCUDA) {

		paMeshCUDA->UpdateConfiguration_Values(cfgMessage);
	}
#endif
}

BError Atom_Mesh_Cubic::SwitchCUDAState(bool cudaState)
{
	BError error(std::string(CLASS_STR(Atom_Mesh_Cubic)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

#if COMPILECUDA == 1

	//--------------------------------------------

	//are we switching to cuda?
	if (cudaState) {

		if (!pMeshBaseCUDA) {

			//then make MeshCUDA object, copying over currently held cpu data
			pMeshBaseCUDA = new Atom_Mesh_CubicCUDA(this);
			paMeshCUDA = dynamic_cast<Atom_MeshCUDA*>(pMeshBaseCUDA);

			error = pMeshBaseCUDA->Error_On_Create();
			if (!error) error = paMeshCUDA->cuaMesh()->set_pointers(paMeshCUDA);
		}
	}
	else {

		//delete MeshCUDA object and null
		if (pMeshBaseCUDA) delete pMeshBaseCUDA;
		pMeshBaseCUDA = nullptr;
		paMeshCUDA = nullptr;
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
void Atom_Mesh_Cubic::CoupleToDipoles(bool status)
{
#if COMPILECUDA == 1
	if (paMeshCUDA) paMeshCUDA->M1()->copy_to_cpuvec(M1);
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
				M1.set_skipcells(mesh_intersection, status);

				if (status) {

					//set interface cells to have magnetization direction along the touching dipole direction
					DBL3 Mdipole_direction = dynamic_cast<Mesh*>((*pSMesh)[idx])->GetAverageMagnetization().normalized();

					Box box = M1.box_from_rect_max(mesh_intersection);

					for (int i = box.s.i; i < box.e.i; i++) {
						for (int j = box.s.j; j < box.e.j; j++) {
							for (int k = box.s.k; k < box.e.k; k++) {

								double Mag = M1[INT3(i, j, k)].norm();
								M1[INT3(i, j, k)] = Mag * Mdipole_direction;
							}
						}
					}
				}
			}
		}
	}

#if COMPILECUDA == 1
	if (paMeshCUDA) paMeshCUDA->M1()->copy_from_cpuvec(M1);
#endif
}

double Atom_Mesh_Cubic::CheckMoveMesh(void)
{
#if COMPILECUDA == 1
	if (paMeshCUDA) return paMeshCUDA->CheckMoveMesh(meshODE.MoveMeshAntisymmetric(), meshODE.MoveMeshThreshold());
#endif

	//move mesh algorithm applied to systems containing domain walls in order to simulate domain wall movement.
	//In this case the moments at the ends of the mesh must be fixed (ODE solver not applied to the ends).
	//Find the average moment projected along the absolute direction of the moment at the ends.
	//For a domain wall at the centre this value should be close to zero. If it exceeds a threshold then move mesh by a unit cellsize along +x or -x.

	if (!move_mesh_trigger) return 0;

	//the ends should not be completely empty, and must have a constant magnetization direction
	int cells_fixed = (int)ceil_epsilon((double)n.x * MOVEMESH_ENDRATIO);

	DBL3 M_end = M1.average_omp(Box(0, 0, 0, cells_fixed, n.y, n.z));
	DBL3 M_av_left = M1.average_omp(Box(cells_fixed, 0, 0, n.x / 2, n.y, n.z));
	DBL3 M_av_right = M1.average_omp(Box(n.x / 2, 0, 0, n.x - cells_fixed, n.y, n.z));

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

//----------------------------------- OTHER CONTROL METHODS : implement pure virtual Atom_Mesh methods

#endif