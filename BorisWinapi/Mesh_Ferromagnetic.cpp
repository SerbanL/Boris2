#include "stdafx.h"
#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

FMesh::FMesh(SuperMesh *pSMesh_) :
	Mesh(MESH_FERROMAGNETIC, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId), VINFO(displayedPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), VINFO(meshRect), VINFO(n), VINFO(h), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(M), VINFO(V), VINFO(S), VINFO(elC), VINFO(Temp), VINFO(pMod),
			//Members in this derived class
			VINFO(move_mesh_trigger), VINFO(skyShift), VINFO(exchange_couple_to_meshes),
			//Material Parameters
			VINFO(grel), VINFO(alpha), VINFO(Ms), VINFO(Nxy), VINFO(A), VINFO(D), VINFO(J1), VINFO(J2), VINFO(K1), VINFO(K2), VINFO(mcanis_ea1), VINFO(mcanis_ea2), VINFO(susrel), VINFO(susprel), VINFO(cHA),
			VINFO(elecCond), VINFO(amrPercentage), VINFO(P), VINFO(beta), VINFO(De), VINFO(SHA), VINFO(flSOT), VINFO(betaD), VINFO(l_sf), VINFO(l_ex), VINFO(l_ph), VINFO(Gi), VINFO(Gmix), VINFO(ts_eff), VINFO(tsi_eff), VINFO(pump_eff),
			VINFO(base_temperature), VINFO(T_Curie), VINFO(T_Curie_material), VINFO(atomic_moment), VINFO(thermCond), VINFO(density), VINFO(shc), VINFO(cT), VINFO(Q)
		},
		{
			//Modules Implementations
			IINFO(Demag_N), IINFO(Demag), IINFO(SDemag_Demag),
			IINFO(Exch_6ngbr_Neu), IINFO(DMExchange), IINFO(iDMExchange), IINFO(SurfExchange),
			IINFO(Zeeman),
			IINFO(Anisotropy_Uniaxial), IINFO(Anisotropy_Cubic),
			IINFO(Transport), IINFO(Heat),
			IINFO(SOTField), 
			IINFO(Roughness)
		}),
	meshODE(this)
{}

FMesh::FMesh(Rect meshRect_, DBL3 h_, SuperMesh *pSMesh_) :
	Mesh(MESH_FERROMAGNETIC, pSMesh_),
	ProgramStateNames(this,
		{
			//Mesh members
			VINFO(meshType), VINFO(meshIdCounter), VINFO(meshId), VINFO(displayedPhysicalQuantity), VINFO(vec3rep), VINFO(displayedParamVar), VINFO(meshRect), VINFO(n), VINFO(h), VINFO(n_e), VINFO(h_e), VINFO(n_t), VINFO(h_t), VINFO(M), VINFO(V), VINFO(S), VINFO(elC), VINFO(Temp), VINFO(pMod),
			//Members in this derived class
			VINFO(move_mesh_trigger), VINFO(skyShift), VINFO(exchange_couple_to_meshes),
			//Material Parameters
			VINFO(grel), VINFO(alpha), VINFO(Ms), VINFO(Nxy), VINFO(A), VINFO(D), VINFO(J1), VINFO(J2), VINFO(K1), VINFO(K2), VINFO(mcanis_ea1), VINFO(mcanis_ea2), VINFO(susrel), VINFO(susprel), VINFO(cHA),
			VINFO(elecCond), VINFO(amrPercentage), VINFO(P), VINFO(beta), VINFO(De), VINFO(SHA), VINFO(flSOT), VINFO(betaD), VINFO(l_sf), VINFO(l_ex), VINFO(l_ph), VINFO(Gi), VINFO(Gmix), VINFO(ts_eff), VINFO(tsi_eff), VINFO(pump_eff),
			VINFO(base_temperature), VINFO(T_Curie), VINFO(T_Curie_material), VINFO(atomic_moment), VINFO(thermCond), VINFO(density), VINFO(shc), VINFO(cT), VINFO(Q)
		},
		{
			//Modules Implementations
			IINFO(Demag_N), IINFO(Demag), IINFO(SDemag_Demag),
			IINFO(Exch_6ngbr_Neu), IINFO(DMExchange), IINFO(iDMExchange), IINFO(SurfExchange),
			IINFO(Zeeman),
			IINFO(Anisotropy_Uniaxial), IINFO(Anisotropy_Cubic),
			IINFO(Transport), IINFO(Heat),
			IINFO(SOTField),
			IINFO(Roughness)
		}),
	meshODE(this)
{
	//default settings
	displayedPhysicalQuantity = MESHDISPLAY_MAGNETIZATION;

	meshRect = meshRect_;
	
	h = h_;
	h_e = h_;
	h_t = h_;

	error_on_create = UpdateConfiguration();

	//when creating a new ferromagnetic mesh set the default Curie temperature for permalloy, including temperature dependence for parameters which depend on it
	if (!error_on_create) SetCurieTemperature(T_Curie);

	//default modules configuration
	if (!error_on_create) error_on_create = AddModule(MOD_DEMAG);
	if (!error_on_create) error_on_create = AddModule(MOD_EXCHANGE6NGBR);
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
}

//----------------------------------- IMPORTANT CONTROL METHODS

//call when the mesh dimensions have changed - sets every quantity to the right dimensions
BError FMesh::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(string(CLASS_STR(FMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

	//get number of cells in each dimension : first divide the mesh rectangle with rounding to get an integer number of cells, then adjust cellsize so h * n gives the rectangle dimensions
	
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

	//erase any unused skyrmion trackers in this mesh
	skyShift.UpdateConfiguration(saveDataList);

	//------------------------ CUDA UpdateConfiguration if set
	
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (!error) error = pMeshCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	//------------------------
	
	//change mesh dimensions in all currently set effective field modules
	for (int idx = 0; idx < (int)pMod.size(); idx++) {
		
		if (pMod[idx] && !error) {

			error = pMod[idx]->UpdateConfiguration(cfgMessage);
		}
	}

	if (!error) error = meshODE.UpdateConfiguration(cfgMessage);
	
	//update material parameters spatial dependence as cellsize and rectangle could have changed
	if (!error && !update_meshparam_var()) error(BERROR_OUTOFMEMORY_NCRIT);

	return error;
}

void FMesh::CoupleToDipoles(bool status)
{
#if COMPILECUDA == 1
	if(pMeshCUDA) pMeshCUDA->M()->copy_to_cpuvec(M);
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

					//set interface cells to have magnetisation direction along the touching dipole direction
					DBL3 Mdipole_direction = (*pSMesh)[idx]->GetAverageMagnetisation().normalized();

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

BError FMesh::SwitchCUDAState(bool cudaState)
{
	BError error(string(CLASS_STR(FMesh)) + "(" + (*pSMesh).key_from_meshId(meshId) + ")");

#if COMPILECUDA == 1

	//--------------------------------------------

	//are we switching to cuda?
	if (cudaState) {

		//delete MeshCUDA object and null (just in case)
		if (pMeshCUDA) delete pMeshCUDA;
		pMeshCUDA = nullptr;

		//then make MeshCUDA object, copying over currently held cpu data
		pMeshCUDA = new FMeshCUDA(this);
		error = pMeshCUDA->Error_On_Create();
		if(!error) error = pMeshCUDA->cuMesh()->set_pointers(pMeshCUDA);
	}
	else {

		//delete MeshCUDA object and null
		if (pMeshCUDA) delete pMeshCUDA;
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

double FMesh::CheckMoveMesh(void)
{

#if COMPILECUDA == 1
	if (pMeshCUDA) return reinterpret_cast<FMeshCUDA*>(pMeshCUDA)->CheckMoveMesh(meshODE.MoveMeshAntisymmetric(), meshODE.MoveMeshThreshold());
#endif

	//move mesh algorithm applied to systems containing domain walls in order to simulate domain wall movement.
	//In this case the magnetisation at the ends of the mesh must be fixed (magnetisation evolver not applied to the ends).
	//Find the average magnetisation projected along the absolute direction of the magnetisation at the ends.
	//For a domain wall at the centre this value should be close to zero. If it exceeds a threshold then move mesh by a cellsize along +x or -x.

	if (!move_mesh_trigger) return 0;

	//the ends should not be completely empty, and must have a constant magnetisation direction
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

void FMesh::SetCurieTemperature(double Tc)
{
	if (Tc >= 1) {

		T_Curie = Tc;

		//also recalculate parameters which have a default temperature dependence based on Tc. Note, if these have been assigned a custom temperature dependence it will be overwritten.

		//1. damping
		//2. Ms -> me
		//3. A -> me^2, D -> me^2
		//4. susrel
		//5. K1 and K2 anisotropy constants -> me^3
		//6. P -> me^2

		//1. alpha = alpha0 * (1 - T/3Tc) up to Tc then 2 * alpha0 * T / 3 Tc above Tc -> for this reason use an array, not a formula

		//-------------------------------

		//2. me = B(me * 3Tc/T + mu*mu0*Ha/kBT), where B(x) = coth(x) - 1/x, Ha applied field, mu is the atomic moment, Tc Curie temperature, kB Boltzmann constant, T temperature
		//Solve me by finding root of equation F(x) = B(cx + d) - x, where x = me, c = 3Tc/T, d = mu*mu0*Ha/kBT - Newton Raphson works very well
		//Note me is the scaling for Ms0 (0 temperature Ms value)

		//3. Scaling for A is me^2

		//4. suspar (parallel susceptiblity) is given by dMe/dHa, where Me = me * Ms0, Ms0 being the zero temperature saturation magnetization
		//Thus suspar = mu0 Ms0 d0 B'(c me + d) / (1 - cB'(c me + d)), where d0 = mu/kBT, so d = d0 * mu0 * Ha
		//
		//We don't want a Ms0 dependence (what happens when Ms0 value is changed??), so we define susrel = suspar / mu0*Ms0 - this is the MatP we'll calculate here. In the LLB equation convert back to suspar to use for the longitudinal relaxation term.

		//5. K1 and K2 scale as me^3

		//6. P scales as me

		//This is B'(x) = 1 - coth^2(x) + 1/x^2
		auto Bdiff = [](double x) -> double {

			double coth = (1 + exp(-2 * x)) / (1 - exp(-2 * x));
			return ((1 - coth * coth) + 1 / (x*x));
		};

		DBL3 Ha = CallModuleMethod(&Zeeman::GetField);

		std::vector<double> t_scaling_me, t_scaling_me2, t_scaling_me3, t_scaling_sus, t_scaling_dam;

		t_scaling_me.assign(MAX_TEMPERATURE + 1, 1.0);
		t_scaling_me2.assign(MAX_TEMPERATURE + 1, 1.0);
		t_scaling_me3.assign(MAX_TEMPERATURE + 1, 1.0);
		t_scaling_sus.assign(MAX_TEMPERATURE + 1, 0.0);
		t_scaling_dam.assign(MAX_TEMPERATURE + 1, 1.0);

		int chunk = (int)floor_epsilon(MAX_TEMPERATURE / OmpThreads);

		#pragma omp parallel for
		for (int thread = 0; thread < OmpThreads; thread++) {

			double me = 1.0;

			for (int t_value = 1 + thread * chunk; t_value <= (thread != OmpThreads - 1 ? (thread + 1) * chunk : (int)MAX_TEMPERATURE); t_value++) {

				double c = 3 * T_Curie / t_value;
				double d0 = (MUB / BOLTZMANN) * atomic_moment / t_value;
				double d = MU0 * d0 * Ha.norm();

				//This is F(x) = B(cx + d) - x
				std::function<double(double)> F = [&](double x) -> double {

					double arg = c * x + d;
					return ((1 + exp(-2 * arg)) / (1 - exp(-2 * arg))) - 1 / arg - x;
				};

				//This is F'(x) = B'(cx + d) * c - 1
				std::function<double(double)> Fdiff = [&](double x) -> double {

					double arg = c * x + d;
					double coth = (1 + exp(-2 * arg)) / (1 - exp(-2 * arg));
					return (((1 - coth * coth) + 1 / (arg*arg)) * c - 1);
				};

				//solve for me at temperature = t_value : Ms0 scaling
				//next starting value for the root finder is me - this makes the NewtonRaphson algorithm extremely efficient (typically a single iteration is required to keep within 1e-6 accuracy!)
				me = Root_NewtonRaphson(F, Fdiff, me, 1e-6);

				//Ms and P scaling
				t_scaling_me[t_value] = me;
				
				//A scaling
				t_scaling_me2[t_value] = me*me;
				
				//K1 and K2 scaling
				t_scaling_me3[t_value] = me * me * me;
				
				//susrel scaling
				t_scaling_sus[t_value] = d0 * Bdiff(c * me + d) / (1.0 - c * Bdiff(c * me + d));
				
				//alpha scaling
				if (t_value < T_Curie) {

					t_scaling_dam[t_value] = 1.0 - (double)t_value / (3.0 * T_Curie);
				}
				else {

					t_scaling_dam[t_value] = 2.0 * (double)t_value / (3.0 * T_Curie);
				}
			}
		}

		//set scaling arrays now
		Ms.set_precalculated_scaling_array(t_scaling_me);
		P.set_precalculated_scaling_array(t_scaling_me);
		A.set_precalculated_scaling_array(t_scaling_me2);
		D.set_precalculated_scaling_array(t_scaling_me2);
		K1.set_precalculated_scaling_array(t_scaling_me3);
		K2.set_precalculated_scaling_array(t_scaling_me3);
		susrel.set_precalculated_scaling_array(t_scaling_sus);
		alpha.set_precalculated_scaling_array(t_scaling_dam);

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

#if COMPILECUDA == 1
	if (pMeshCUDA) {

		//T_Curie changed : sync with cuda version
		pMeshCUDA->T_Curie.from_cpu((cuBReal)T_Curie);
	}
#endif
}

void FMesh::SetAtomicMoment(double atomic_moment_ub)
{
	atomic_moment = atomic_moment_ub;

	//setting atomic_moment will affect the temperature dependence of me (normalised equilibrium magnetisation), so some parameter temperature dependencies must change - calling SetCurieTemperature with the current Curie temperature will do this
	SetCurieTemperature(T_Curie);
}