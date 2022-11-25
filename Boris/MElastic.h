#pragma once

#include "BorisLib.h"
#include "Modules.h"

class Mesh;

#ifdef MODULE_COMPILATION_MELASTIC

class SuperMesh;
class SMElastic;
class MElastic_Boundary;

class MElastic :
	public Modules,
	public ProgramState<MElastic, std::tuple<
	DBL3, 
	VEC<double>, VEC<double>, VEC<double>, VEC<DBL3>, VEC<double>, VEC<double>, VEC<double>,
	TEquation<double, double, double, double>, TEquation<double, double, double, double>>, std::tuple<>>
{
	friend SMElastic;

#if COMPILECUDA == 1
	friend class MElasticCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Mesh *pMesh;

	//pointer to supermesh
	SuperMesh* pSMesh;

	//pointer to SMElastic supermesh module; this is set from SMElastic in its UpdateConfiguration method
	//need it to access some SMElastic properties
	SMElastic* pSMEl = nullptr;

	//----------------------

	//Applied uniform stress vector (Cartesian coordinates) - used only for no solve mode (el_dT zero)
	DBL3 Tsig;

	//---------------------- FDTD scheme velocity-stress representation. Need separate components as they are staggered.

	//velocity x : mid x edges, size nx * (ny + 1) * (nz + 1)
	VEC<double> vx;

	//velocity y : mid y edges, size (nx + 1) * ny * (nz + 1)
	VEC<double> vy;

	//velocity z : mid z edges, size (nx + 1) * (ny + 1) * nz
	VEC<double> vz;

	//diagonal stress components : vertices, size (nx + 1)*(ny + 1)*(nz + 1)
	VEC<DBL3> sdd;

	//off-diagonal stress sigma_xy : mid xy faces, size nx*ny*(nz + 1)
	VEC<double> sxy;

	//off-diagonal stress sigma_xz : mid xz faces, size nx*(ny + 1)*nz
	VEC<double> sxz;

	//off-diagonal stress sigma_yz : mid yz faces, size (nx + 1)*ny*nz
	VEC<double> syz;

	//----------------------

	//Fixed surfaces (plane rectangles in absolute coordinates on mesh surfaces where we set zero Dirichlet value for displacement and velocity)
	std::vector<Rect> fixed_u_surfaces;

	//Surfaces (plane rectangles in absolute coordinates on mesh surfaces) on which we define external forces. Must not overlap with fixed_u_surfaces.
	//The external forces are defined using a vector equation (or a constant force) with 3 components for Fx, Fy, Fz.
	//The equation depends on x, y, t variables : x, y are positions relative to corresponding plane rectangle. t is the time (simulation time)
	//These are used to set non-homogeneous Neumann boundary conditions
	std::vector<MElastic_Boundary> external_stress_surfaces;

	//----------------------

	//Strain using user equation, thus allowing simultaneous spatial (x, y, z), and stage time (t) dependence.
	//A number of constants are always present : mesh dimensions in m (Lx, Ly, Lz)
	//When text equation set, then elastodynamics solver is disabled.
	//diagonal
	TEquation<double, double, double, double> Sd_equation;
	//off-diagonal
	TEquation<double, double, double, double> Sod_equation;

private:

	//----------------------------------------------- Auxiliary

	//Run-time auxiliary to set strain directly from user supplied text formulas
	void Set_Strain_From_Formula(void);

	//Update TEquation object with user constants values
	void UpdateTEquationUserConstants(bool makeCuda = true);

	//----------------------------------------------- Computational Helpers

	//update velocity for dT time increment (also updating displacement)
	void Iterate_Elastic_Solver_Velocity(double dT);
	//update stress for dT time increment
	void Iterate_Elastic_Solver_Stress(double dT);
	//update strain from stress
	void Calculate_Strain_From_Stress(void);

	//compute magnetoelastic effective field to use in magnetization equation; return energy
	double Calculate_MElastic_Field(void);

	//compute strain form mechanical displacement (only used in the no solve mode, i.e. with externally supplied mechanical displacement, or in constant and uniform stress mode).
	void Calculate_Strain(void);

	//---------------------------------------------- CMBND

	void make_velocity_continuous(
		CMBNDInfo& contact,
		VEC<double>& vx_sec, VEC<double>& vy_sec, VEC<double>& vz_sec, VEC_VC<DBL3>& u_disp_sec);

	void make_stress_continuous(
		CMBNDInfo& contact,
		VEC<DBL3>& sdd_sec, VEC<double>& sxy_sec, VEC<double>& sxz_sec, VEC<double>& syz_sec,
		VEC_VC<DBL3>& u_disp_sec);

public:

	MElastic(Mesh *pMesh_);
	~MElastic();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Energy methods

	//FM Mesh
	double Get_EnergyChange(int spin_index, DBL3 Mnew);

	//AFM Mesh
	DBL2 Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B);

	//------------------- STRAIN GENERATION without SOLVER

	//Set/Get uniform stress
	void SetUniformStress(DBL3 Tsig_xyz);
	DBL3 GetUniformStress(void) { return Tsig; }

	//Set strain or displacement by loading from OVF2 files
	
	//Displacement : from this calculate strain tensor
	BError Load_Displacement_OVF2(std::string fileName);

	//Tensor : displacement is not calculated, as we only need the strain tensor to obtain the effective fields at runtime
	//For a strain tensor for cubic crystals, we only have 6 elements : diagonal and off-diagonal (symmetric).
	//These are specified using 2 separate OVF2 files containing vector data:
	//one for the xx, yy, zz elements (diagonal)
	//the other for the yz, xz, xy elements (off-diagonal, in this order)
	BError Load_Strain_OVF2(std::string fileName_Diag, std::string fileName_ODiag);

	//------------------- Configuration

	//reset stress-strain solver to initial values (zero velocity, displacement and stress)
	void Reset_ElSolver(void);

	//set diagonal and shear strain text equations
	BError Set_Sd_Equation(std::string text_equation);
	BError Set_Sod_Equation(std::string text_equation);
	//clear text equations
	void Clear_Sd_Sod_Equations(void);

	std::string Get_Sd_Equation(void) { return Sd_equation.show_equation(); }
	std::string Get_Sod_Equation(void) { return Sod_equation.show_equation(); }
};

#else

class MElastic :
	public Modules
{

private:

private:

	//----------------------------------------------- Computational Helpers

	//compute strain form mechanical displacement
	void Calculate_Strain(void) {}

public:

	MElastic(Mesh *pMesh_) {}
	~MElastic() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------

	//Set/Get uniform stress
	void SetUniformStress(DBL3 Tsig_xyz) {}
	DBL3 GetUniformStress(void) { return DBL3(); }

	BError Load_Displacement_OVF2(std::string fileName) { return BError(); }
	BError Load_Strain_OVF2(std::string fileName_Diag, std::string fileName_ODiag) { return BError(); }

	//------------------- Configuration

	//reset stress-strain solver to initial values (zero velocity, displacement and stress)
	void Reset_ElSolver(void) {}

	//set diagonal and shear strain text equations
	BError Set_Sd_Equation(std::string text_equation) { return BError(); }
	BError Set_Sod_Equation(std::string text_equation) { return BError(); }
	//clear text equations
	void Clear_Sd_Sod_Equations(void) {}

	std::string Get_Sd_Equation(void) { return ""; }
	std::string Get_Sod_Equation(void) { return ""; }
};

#endif