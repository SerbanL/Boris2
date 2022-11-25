#include "stdafx.h"
#include "MElastic.h"

#ifdef MODULE_COMPILATION_MELASTIC

#include "Mesh_Ferromagnetic.h"
#include "MeshParamsControl.h"

#include "SuperMesh.h"

#include "MElastic_Boundaries.h"

#if COMPILECUDA == 1
#include "MElasticCUDA.h"
#endif

MElastic::MElastic(Mesh *pMesh_) :
	Sd_equation({ "x", "y", "z", "t" }), Sod_equation({ "x", "y", "z", "t" }),
	Modules(),
	ProgramStateNames(this, { 
	VINFO(Tsig),
	VINFO(vx), VINFO(vy), VINFO(vz), VINFO(sdd), VINFO(sxy), VINFO(sxz), VINFO(syz), 
	VINFO(Sd_equation), VINFO(Sod_equation) }, {})
{
	pMesh = pMesh_;
	pSMesh = pMesh->pSMesh;

	Tsig = DBL3(0, 0, 0);

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

MElastic::~MElastic()
{
	//free memory used for heat solver
	pMesh->u_disp.clear();
	pMesh->strain_diag.clear();
	pMesh->strain_odiag.clear();
}

BError MElastic::Initialize(void)
{
	BError error(CLASS_STR(MElastic));
	
	if (!initialized) {
		
		//Must have at least one fixed surface defined.
		if (!pSMEl->fixed_u_surfaces.size()) return error(BERROR_INCORRECTCONFIG);

		//build fixed surfaces for this mesh
		fixed_u_surfaces.clear();

		for (int idx = 0; idx < pSMEl->fixed_u_surfaces.size(); idx++) {

			Rect intersection = pSMEl->fixed_u_surfaces[idx].get_intersection(pMesh->meshRect);
			if (!intersection.IsNull() && intersection.IsPlane()) fixed_u_surfaces.push_back(intersection);
		}

		//build stress surfaces for this mesh
		external_stress_surfaces.clear();

		std::vector<std::pair<std::string, double>> constants(pMesh->userConstants.size());
		for (int idx = 0; idx < pMesh->userConstants.size(); idx++) {

			constants[idx] = { pMesh->userConstants.get_key_from_index(idx), pMesh->userConstants[idx] };
		}

		DBL3 meshDim = pMesh->GetMeshDimensions();
		constants.push_back({ "Lx", meshDim.x });
		constants.push_back({ "Ly", meshDim.y });
		constants.push_back({ "Lz", meshDim.z });

		for (int idx = 0; idx < pSMEl->stress_surfaces_rect.size(); idx++) {

			Rect intersection = pSMEl->stress_surfaces_rect[idx].get_intersection(pMesh->meshRect);
			if (!intersection.IsNull() && intersection.IsPlane()) {
				
				external_stress_surfaces.push_back(MElastic_Boundary());
				external_stress_surfaces.back().setup_surface(pMesh->u_disp, intersection);
				if (!external_stress_surfaces.back().setup_equation_stimulus(pSMEl->stress_surfaces_equations[idx], constants)) return error(BERROR_INCORRECTSTRING);
			}
		}

		//setup MElastic_Boundary entries for CMBND interfaces - this is done in SMelastic (initialized after all MElastic modules have initialized)

		//fixed_u_surfaces rectangles must not intersect with any external_stress_surfaces rectangles
		for (auto stress_surf : external_stress_surfaces) {
			for (auto fixed_rect : fixed_u_surfaces) {

				Rect intersection = fixed_rect.get_intersection(stress_surf.get_surface());
				//if there's any intersection it cannot be a plane (line and point are still fine)
				if (intersection.IsPlane()) return error(BERROR_INCORRECTCONFIG);
			}
		}
		
		//set Dirichlet conditions for u_disp (zero, i.e. fixed, or zero displacement, points)
		pMesh->u_disp.clear_dirichlet_flags();
		for (auto fixed_rect : fixed_u_surfaces) pMesh->u_disp.set_dirichlet_conditions(fixed_rect, DBL3());

		//set Dirichlet conditions for strain_diag (external force)
		pMesh->strain_diag.clear_dirichlet_flags();
		for (auto external_stress : external_stress_surfaces) pMesh->strain_diag.set_dirichlet_conditions(external_stress.get_surface(), DBL3());
	}

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect,
		(MOD_)pMesh->Get_Module_Heff_Display() == MOD_MELASTIC || pMesh->IsOutputDataSet_withRect(DATA_E_MELASTIC),
		(MOD_)pMesh->Get_Module_Energy_Display() == MOD_MELASTIC || pMesh->IsOutputDataSet_withRect(DATA_E_MELASTIC),
		pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	return error;
}

BError MElastic::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(MElastic));

	Uninitialize();

	bool success = true;

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE)) {

		//make sure the cellsize divides the mesh rectangle
		pMesh->n_m = round(pMesh->meshRect / pMesh->h_m);
		if (pMesh->n_m.x < 2) pMesh->n_m.x = 2;
		if (pMesh->n_m.y < 2) pMesh->n_m.y = 2;
		if (pMesh->n_m.z < 2) pMesh->n_m.z = 2;
		pMesh->h_m = pMesh->meshRect / pMesh->n_m;

		//make sure correct memory is assigned for heat solver quantities

		if (pMesh->M.linear_size()) {

			//in a ferromagnet the magnetization sets the shape only on initialization. If already initialized then shape already set.
			if (pMesh->u_disp.linear_size()) {

				success = pMesh->u_disp.resize(pMesh->h_m, pMesh->meshRect);
			}
			else {

				success = pMesh->u_disp.assign(pMesh->h_m, pMesh->meshRect, DBL3(), pMesh->M);
			}
		}
		else {

			//the shape is set directly
			if (pMesh->u_disp.linear_size()) {

				success = pMesh->u_disp.resize(pMesh->h_m, pMesh->meshRect);
			}
			else {

				success = pMesh->u_disp.assign(pMesh->h_m, pMesh->meshRect, DBL3());
			}
		}

		//strain tensor - set empty cells using information in u_disp
		if (pMesh->strain_diag.linear_size()) {

			success &= pMesh->strain_diag.resize(pMesh->h_m, pMesh->meshRect, pMesh->u_disp);
			success &= pMesh->strain_odiag.resize(pMesh->h_m, pMesh->meshRect, pMesh->u_disp);
		}
		else {

			success &= pMesh->strain_diag.assign(pMesh->h_m, pMesh->meshRect, DBL3(), pMesh->u_disp);
			success &= pMesh->strain_odiag.assign(pMesh->h_m, pMesh->meshRect, DBL3(), pMesh->u_disp);
		}

		//correct size for FDTD data
		success &= vx.resize(SZ3(pMesh->n_m.x, pMesh->n_m.y + 1, pMesh->n_m.z + 1));
		success &= vy.resize(SZ3(pMesh->n_m.x + 1, pMesh->n_m.y, pMesh->n_m.z + 1));
		success &= vz.resize(SZ3(pMesh->n_m.x + 1, pMesh->n_m.y + 1, pMesh->n_m.z));

		success &= sdd.resize(SZ3(pMesh->n_m.x + 1, pMesh->n_m.y + 1, pMesh->n_m.z + 1));
		success &= sxy.resize(SZ3(pMesh->n_m.x, pMesh->n_m.y, pMesh->n_m.z + 1));
		success &= sxz.resize(SZ3(pMesh->n_m.x, pMesh->n_m.y + 1, pMesh->n_m.z));
		success &= syz.resize(SZ3(pMesh->n_m.x + 1, pMesh->n_m.y, pMesh->n_m.z));

		//update mesh dimensions in equation constants
		if (Sd_equation.is_set() || Sod_equation.is_set()) {

			DBL3 meshDim = pMesh->GetMeshDimensions();

			if (Sd_equation.is_set()) {

				Sd_equation.set_constant("Lx", meshDim.x, false);
				Sd_equation.set_constant("Ly", meshDim.y, false);
				Sd_equation.set_constant("Lz", meshDim.z, false);
				Sd_equation.remake_equation();
			}

			if (Sod_equation.is_set()) {

				Sod_equation.set_constant("Lx", meshDim.x, false);
				Sod_equation.set_constant("Ly", meshDim.y, false);
				Sod_equation.set_constant("Lz", meshDim.z, false);
				Sod_equation.remake_equation();
			}
		}
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

void MElastic::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		//if this affects external stress surfaces, equations will need to be remade
		Uninitialize();

		UpdateTEquationUserConstants(false);
	}

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pModuleCUDA->UpdateConfiguration_Values(cfgMessage);
	}
#endif
}

BError MElastic::MakeCUDAModule(void)
{
	BError error(CLASS_STR(MElastic));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new MElasticCUDA(pMesh, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double MElastic::UpdateField(void)
{
	if (Sd_equation.is_set_vector() || Sod_equation.is_set_vector()) {

		//strain specified using a formula
		Set_Strain_From_Formula();
	}

	return 0.0;
}

//-------------------Energy methods

//FM mesh
double MElastic::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (pMesh->M.is_not_empty(spin_index)) {

		double Ms = pMesh->Ms;
		DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
		DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
		DBL3 mcanis_ea3 = pMesh->mcanis_ea3;
		DBL2 MEc = pMesh->MEc;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->Ms, Ms, pMesh->MEc, MEc, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2, pMesh->mcanis_ea3, mcanis_ea3);

		DBL3 position = pMesh->M.cellidx_to_position(spin_index);
		//xx, yy, zz
		DBL3 Sd = pMesh->strain_diag[position];
		//yz, xz, xy
		DBL3 Sod = pMesh->strain_odiag[position];

		//normalised magnetization
		//Magneto-elastic term here applicable for a cubic crystal. We use the mcanis_ea1 and mcanis_ea2 axes to fix the cubic lattice orientation, thus rotate the m, Sd and Sod vectors.

		Sd = DBL3(Sd * mcanis_ea1, Sd * mcanis_ea2, Sd * mcanis_ea3);
		Sod = DBL3(Sod * mcanis_ea1, Sod * mcanis_ea2, Sod * mcanis_ea3);

		auto Get_Energy = [&](DBL3 M) -> double
		{
			DBL3 m = DBL3(M * mcanis_ea1, M * mcanis_ea2, M * mcanis_ea3) / Ms;

			DBL3 Hmel_1 = (-2.0 * MEc.i / (MU0 * Ms)) * DBL3(
				m.x*Sd.x*mcanis_ea1.x + m.y*Sd.y*mcanis_ea2.x + m.z*Sd.z*mcanis_ea3.x,
				m.x*Sd.x*mcanis_ea1.y + m.y*Sd.y*mcanis_ea2.y + m.z*Sd.z*mcanis_ea3.y,
				m.x*Sd.x*mcanis_ea1.z + m.y*Sd.y*mcanis_ea2.z + m.z*Sd.z*mcanis_ea3.z);

			DBL3 Hmel_2 = (-2.0 * MEc.j / (MU0 * Ms)) * DBL3(
				Sod.z * (mcanis_ea1.x*m.y + mcanis_ea2.x*m.x) + Sod.y * (mcanis_ea1.x*m.z + mcanis_ea3.x*m.x) + Sod.x * (mcanis_ea2.x*m.z + mcanis_ea3.x*m.y),
				Sod.z * (mcanis_ea1.y*m.y + mcanis_ea2.y*m.x) + Sod.y * (mcanis_ea1.y*m.z + mcanis_ea3.y*m.x) + Sod.x * (mcanis_ea2.y*m.z + mcanis_ea3.y*m.y),
				Sod.z * (mcanis_ea1.z*m.y + mcanis_ea2.z*m.x) + Sod.y * (mcanis_ea1.z*m.z + mcanis_ea3.z*m.x) + Sod.x * (mcanis_ea2.z*m.z + mcanis_ea3.z*m.y));

			return -MU0 * M * (Hmel_1 + Hmel_2) / 2;
		};
		
		if (Mnew != DBL3()) return pMesh->h.dim() * (Get_Energy(Mnew) - Get_Energy(pMesh->M[spin_index]));
		else return pMesh->h.dim() *  Get_Energy(pMesh->M[spin_index]);
	}
	else return 0.0;
}

//AFM mesh
DBL2 MElastic::Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B)
{
	if (pMesh->M.is_not_empty(spin_index)) {

		DBL2 Ms_AFM = pMesh->Ms_AFM;
		DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
		DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
		DBL3 mcanis_ea3 = pMesh->mcanis_ea3;
		DBL2 MEc = pMesh->MEc;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->Ms_AFM, Ms_AFM, pMesh->MEc, MEc, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2, pMesh->mcanis_ea3, mcanis_ea3);

		DBL3 position = pMesh->M.cellidx_to_position(spin_index);
		//xx, yy, zz
		DBL3 Sd = pMesh->strain_diag[position];
		//yz, xz, xy
		DBL3 Sod = pMesh->strain_odiag[position];

		Sd = DBL3(Sd * mcanis_ea1, Sd * mcanis_ea2, Sd * mcanis_ea3);
		Sod = DBL3(Sod * mcanis_ea1, Sod * mcanis_ea2, Sod * mcanis_ea3);

		auto Get_Energy = [&](DBL3 M, DBL3 M2) -> DBL2 {

			//normalised magnetization
			//Magneto-elastic term here applicable for a cubic crystal. We use the mcanis_ea1 and mcanis_ea2 axes to fix the cubic lattice orientation, thus rotate the m, Sd and Sod vectors.

			DBL3 mA = DBL3(M * mcanis_ea1, M * mcanis_ea2, M * mcanis_ea3) / Ms_AFM.i;
			DBL3 mB = DBL3(M2 * mcanis_ea1, M2 * mcanis_ea2, M2 * mcanis_ea3) / Ms_AFM.j;

			DBL3 Hmel_1_A = (-2.0 * MEc.i / (MU0 * Ms_AFM.i)) * DBL3(
				mA.x*Sd.x*mcanis_ea1.x + mA.y*Sd.y*mcanis_ea2.x + mA.z*Sd.z*mcanis_ea3.x,
				mA.x*Sd.x*mcanis_ea1.y + mA.y*Sd.y*mcanis_ea2.y + mA.z*Sd.z*mcanis_ea3.y,
				mA.x*Sd.x*mcanis_ea1.z + mA.y*Sd.y*mcanis_ea2.z + mA.z*Sd.z*mcanis_ea3.z);

			DBL3 Hmel_2_A = (-2.0 * MEc.j / (MU0 * Ms_AFM.i)) * DBL3(
				Sod.z * (mcanis_ea1.x*mA.y + mcanis_ea2.x*mA.x) + Sod.y * (mcanis_ea1.x*mA.z + mcanis_ea3.x*mA.x) + Sod.x * (mcanis_ea2.x*mA.z + mcanis_ea3.x*mA.y),
				Sod.z * (mcanis_ea1.y*mA.y + mcanis_ea2.y*mA.x) + Sod.y * (mcanis_ea1.y*mA.z + mcanis_ea3.y*mA.x) + Sod.x * (mcanis_ea2.y*mA.z + mcanis_ea3.y*mA.y),
				Sod.z * (mcanis_ea1.z*mA.y + mcanis_ea2.z*mA.x) + Sod.y * (mcanis_ea1.z*mA.z + mcanis_ea3.z*mA.x) + Sod.x * (mcanis_ea2.z*mA.z + mcanis_ea3.z*mA.y));

			DBL3 Hmel_1_B = (-2.0 * MEc.i / (MU0 * Ms_AFM.j)) * DBL3(
				mB.x*Sd.x*mcanis_ea1.x + mB.y*Sd.y*mcanis_ea2.x + mB.z*Sd.z*mcanis_ea3.x,
				mB.x*Sd.x*mcanis_ea1.y + mB.y*Sd.y*mcanis_ea2.y + mB.z*Sd.z*mcanis_ea3.y,
				mB.x*Sd.x*mcanis_ea1.z + mB.y*Sd.y*mcanis_ea2.z + mB.z*Sd.z*mcanis_ea3.z);

			DBL3 Hmel_2_B = (-2.0 * MEc.j / (MU0 * Ms_AFM.j)) * DBL3(
				Sod.z * (mcanis_ea1.x*mB.y + mcanis_ea2.x*mB.x) + Sod.y * (mcanis_ea1.x*mB.z + mcanis_ea3.x*mB.x) + Sod.x * (mcanis_ea2.x*mB.z + mcanis_ea3.x*mB.y),
				Sod.z * (mcanis_ea1.y*mB.y + mcanis_ea2.y*mB.x) + Sod.y * (mcanis_ea1.y*mB.z + mcanis_ea3.y*mB.x) + Sod.x * (mcanis_ea2.y*mB.z + mcanis_ea3.y*mB.y),
				Sod.z * (mcanis_ea1.z*mB.y + mcanis_ea2.z*mB.x) + Sod.y * (mcanis_ea1.z*mB.z + mcanis_ea3.z*mB.x) + Sod.x * (mcanis_ea2.z*mB.z + mcanis_ea3.z*mB.y));

			return DBL2(-MU0 * M * (Hmel_1_A + Hmel_2_A) / 2, -MU0 * M2 * (Hmel_1_B + Hmel_2_B) / 2);
		};

		if (Mnew_A != DBL3() && Mnew_B != DBL3()) return pMesh->h.dim() * (Get_Energy(Mnew_A, Mnew_B) - Get_Energy(pMesh->M[spin_index], pMesh->M2[spin_index]));
		else return pMesh->h.dim() *  Get_Energy(pMesh->M[spin_index], pMesh->M2[spin_index]);
	}
	else return DBL2();
}

//------------------- Configuration

//reset stress-strain solver to initial values (zero velocity, displacement and stress)
void MElastic::Reset_ElSolver(void)
{
	vx.set(0.0); vy.set(0.0); vz.set(0.0);
	sdd.set(DBL3());
	sxy.set(0.0); sxz.set(0.0); syz.set(0.0);
	
	pMesh->u_disp.set(DBL3());
	pMesh->strain_diag.set(DBL3());
	pMesh->strain_odiag.set(DBL3());

#if COMPILECUDA == 1
	if (pModuleCUDA) dynamic_cast<MElasticCUDA*>(pModuleCUDA)->Reset_ElSolver();
#endif
}

//set diagonal and shear strain text equations
BError MElastic::Set_Sd_Equation(std::string text_equation)
{
	BError error(CLASS_STR(MElastic));

	Sd_equation.clear();

	DBL3 meshDim = pMesh->GetMeshDimensions();

	//important to set user constants first : if these are required then the make_from_string call will fail. This may mean the user constants are not set when we expect them to.
	UpdateTEquationUserConstants(false);

	if (!Sd_equation.make_from_string(text_equation, { {"Lx", meshDim.x}, {"Ly", meshDim.y}, {"Lz", meshDim.z} })) return error(BERROR_INCORRECTSTRING);

#if COMPILECUDA == 1
	if (pModuleCUDA) dynamic_cast<MElasticCUDA*>(pModuleCUDA)->Set_Sd_Equation(Sd_equation.get_vector_fspec());
#endif

	return error;
}

BError MElastic::Set_Sod_Equation(std::string text_equation)
{
	BError error(CLASS_STR(MElastic));

	Sod_equation.clear();

	DBL3 meshDim = pMesh->GetMeshDimensions();

	//important to set user constants first : if these are required then the make_from_string call will fail. This may mean the user constants are not set when we expect them to.
	UpdateTEquationUserConstants(false);

	if (!Sod_equation.make_from_string(text_equation, { {"Lx", meshDim.x}, {"Ly", meshDim.y}, {"Lz", meshDim.z} })) return error(BERROR_INCORRECTSTRING);

#if COMPILECUDA == 1
	if (pModuleCUDA) dynamic_cast<MElasticCUDA*>(pModuleCUDA)->Set_Sod_Equation(Sod_equation.get_vector_fspec());
#endif

	return error;
}

//clear text equations
void MElastic::Clear_Sd_Sod_Equations(void)
{
	Sd_equation.clear();
	Sod_equation.clear();

#if COMPILECUDA == 1
	if (pModuleCUDA) dynamic_cast<MElasticCUDA*>(pModuleCUDA)->Clear_Sd_Sod_Equations();
#endif
}

//Update TEquation object with user constants values
void MElastic::UpdateTEquationUserConstants(bool makeCuda)
{
	if (pMesh->userConstants.size()) {

		std::vector<std::pair<std::string, double>> constants(pMesh->userConstants.size());
		for (int idx = 0; idx < pMesh->userConstants.size(); idx++) {

			constants[idx] = { pMesh->userConstants.get_key_from_index(idx), pMesh->userConstants[idx] };
		}

		Sd_equation.set_constants(constants);
		Sod_equation.set_constants(constants);

		//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
		if (pModuleCUDA && makeCuda) {

			dynamic_cast<MElasticCUDA*>(pModuleCUDA)->Set_Sd_Equation(Sd_equation.get_vector_fspec());
			dynamic_cast<MElasticCUDA*>(pModuleCUDA)->Set_Sod_Equation(Sod_equation.get_vector_fspec());
		}
#endif
	}
}

#endif
