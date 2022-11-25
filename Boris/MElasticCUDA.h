#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_MELASTIC

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

#include "MElastic_BoundariesCUDA.h"

class SMElasticCUDA;
class MElastic;
class MeshCUDA;
class Mesh;

class MElasticCUDA :
	public ModulesCUDA
{
	friend SMElasticCUDA;
	friend MElastic;

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

	//pointer to cpu version of MeshCUDA
	Mesh *pMesh;

	//pointer to MElastic holder module
	MElastic* pMElastic;

private:

	//---------------------- FDTD scheme velocity-stress representation. Need separate components as they are staggered.

	//velocity x : mid x edges, size nx * (ny + 1) * (nz + 1)
	cu_obj<cuVEC<cuBReal>> vx;

	//velocity y : mid y edges, size (nx + 1) * ny * (nz + 1)
	cu_obj<cuVEC<cuBReal>> vy;

	//velocity z : mid z edges, size (nx + 1) * (ny + 1) * nz
	cu_obj<cuVEC<cuBReal>> vz;

	//diagonal stress components : vertices, size (nx + 1)*(ny + 1)*(nz + 1)
	cu_obj<cuVEC<cuReal3>> sdd;

	//off-diagonal stress sigma_xy : mid xy faces, size nx*ny*(nz + 1)
	cu_obj<cuVEC<cuBReal>> sxy;

	//off-diagonal stress sigma_xz : mid xz faces, size nx*(ny + 1)*nz
	cu_obj<cuVEC<cuBReal>> sxz;

	//off-diagonal stress sigma_yz : mid yz faces, size (nx + 1)*ny*nz
	cu_obj<cuVEC<cuBReal>> syz;

	//----------------------

	//corresponds to MElastic::external_stress_surfaces
	std::vector<cu_obj<MElastic_BoundaryCUDA>> external_stress_surfaces;

	//same information as in external_stress_surfaces, but stored in a cu_arr so we can pass it whole into a CUDA kernel
	cu_arr<MElastic_BoundaryCUDA> external_stress_surfaces_arr;

	//with cuda switched on this will hold the text equation object (TEquationCUDA manages GPU memory but is held in CPU memory itself, so cannot place it directly in MElastic_BoundaryCUDA)
	//vector size same as MElastic(CUDA)::external_stress_surfaces
	std::vector<TEquationCUDA<cuBReal, cuBReal, cuBReal>> Fext_equationCUDA;

	//----------------------

	//Strain using user equation, thus allowing simultaneous spatial (x, y, z), and stage time (t) dependence.
	//A number of constants are always present : mesh dimensions in m (Lx, Ly, Lz)
	//When text equation set, then elastodynamics solver is disabled.
	//diagonal
	TEquationCUDA<cuBReal, cuBReal, cuBReal, cuBReal> Sd_equation;
	//off-diagonal
	TEquationCUDA<cuBReal, cuBReal, cuBReal, cuBReal> Sod_equation;

private:

	//----------------------------------------------- Auxiliary

	//Run-time auxiliary to set strain directly from user supplied text formulas
	void Set_Strain_From_Formula(void);

	//update velocity for dT time increment (also updating displacement)
	void Iterate_Elastic_Solver_Velocity(double dT);
	//update stress for dT time increment
	void Iterate_Elastic_Solver_Stress(double dT);
	//update strain from stress
	void Calculate_Strain_From_Stress(void);

	//compute magnetoelastic effective field to use in magnetization equation. Accumulate energy.
	void Calculate_MElastic_Field(void);

	//---------------------------------------------- CMBND

	void make_velocity_continuous(
		size_t size, int axis,
		cu_obj<CMBNDInfoCUDA>& contact,
		cu_obj<cuVEC<cuBReal>>& vx_sec, cu_obj<cuVEC<cuBReal>>& vy_sec, cu_obj<cuVEC<cuBReal>>& vz_sec, cu_obj<cuVEC_VC<cuReal3>>& u_disp_sec);

	void make_stress_continuous(
		size_t size, int axis,
		cu_obj<CMBNDInfoCUDA>& contact,
		cu_obj<cuVEC<cuReal3>>& sdd_sec, cu_obj<cuVEC<cuBReal>>& sxy_sec, cu_obj<cuVEC<cuBReal>>& sxz_sec, cu_obj<cuVEC<cuBReal>>& syz_sec,
		cu_obj<cuVEC_VC<cuReal3>>& u_disp_sec);

public:

	MElasticCUDA(Mesh* pMesh_, MElastic* pMElastic_);
	~MElasticCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	void UpdateField(void);

	//------------------- Configuration

	//reset stress-strain solver to initial values (zero velocity, displacement and stress)
	void Reset_ElSolver(void);

	//set diagonal and shear strain text equations
	BError Set_Sd_Equation(std::vector<std::vector< std::vector<EqComp::FSPEC> >> fspec);
	BError Set_Sod_Equation(std::vector<std::vector< std::vector<EqComp::FSPEC> >> fspec);
	//clear text equations
	void Clear_Sd_Sod_Equations(void);

	//-------------------
	
	//copy all required mechanical VECs from their cpu versions
	BError copy_VECs_to_GPU(void);
};

#else

class MElasticCUDA
{
};

#endif

#endif