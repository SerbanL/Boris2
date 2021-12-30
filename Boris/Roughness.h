#pragma once

#include "BorisLib.h"
#include "Modules.h"

class Mesh;

#ifdef MODULE_COMPILATION_ROUGHNESS

#include "Convolution.h"
#include "RoughnessKernel.h"

#include "Boris_Enums_Defs.h"

#if COMPILECUDA == 1
#include "RoughnessCUDA.h"
#endif

////////////////////////////////////////////////////////////////////////////////////////////////

class Roughness :
	public Modules,
	public Convolution<Roughness, RoughnessKernel<Roughness>>,
	public ProgramState<Roughness, std::tuple<INT3, VEC_VC<double>>, std::tuple<>>
{
	friend RoughnessKernel<Roughness>;

#if COMPILECUDA == 1
	friend RoughnessCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Mesh* pMesh;

	//cellsize refine factors : if pMesh->n is the number of cells in the F mesh, then refine & pMesh->n is the number of cells in the refined mesh for Roughness computations
	//Similarly the cellsize of the refined mesh is given by pMesh->h / refine
	INT3 refine = INT3(1);

	//multiplicative functions used to obtain roughness field in the coarse mesh
	VEC<DBL3> Fmul_rough;
	VEC<DBL3> Fomul_rough;

	//magnetization shape on the fine mesh, i.e. a refined version of pMesh->M, but we only need to keep track of empty / non-empty cells (i.e. the shape)
	//When loading an image file, if Roughness module enabled then load it into Mshape_fine first, then set pMesh->M from Mfine such that
	//the coarse magnetic body includes the fine one (i.e. if a coarse cell has any non-empty fine cells then set them as non-empty)
	VEC_VC<double> Mshape_fine;

	//flag used by RoughnessKernel to choose which kernel multiplication to use during a convolution
	bool off_diagonal_convolution = false;

private:

	//compute Fmul_rough (diagonal) or Fomul_rough (off-diagonal) after the Roughness Kernel has been computed
	BError Initialize_Fmul(VEC<DBL3>& Fout, bool off_diagonal);

	//set shape of M from Mshape_fine : any coarse cells which have at least one fine cell non-empty must be set non-empty
	void Set_Coarse_M_Inclusion(void);
	
	//set values in Mshape_fine so it display nicely with color coding along z axis, indicating roughness
	void Calculate_Mshape_fine(void);

public:

	Roughness(Mesh *pMesh_);
	~Roughness() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void);

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Energy methods

	//FM mesh
	double Get_EnergyChange(int spin_index, DBL3 Mnew);

	//AFM mesh
	DBL2 Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B);

	//-------------------Setters

	//change cellsize refinement
	void set_refine(INT3 refine_) { refine = refine_; UpdateConfiguration(UPDATECONFIG_ROUGHNESS_CHANGE); }

	//called from Mesh::change_mesh_shape : apply shape to Mshape_fine, then set pMesh->M from it
	template  <typename Lambda, typename ... PType>
	BError change_mesh_shape(Lambda& run_this, PType& ... params)
	{
		BError error(__FUNCTION__);

		Uninitialize();

		//apply shape to Mshape_fine
		error = run_this(Mshape_fine, 1.0, params...);
		
		if (!error) {

			//set shape of M from Mshape_fine
			Set_Coarse_M_Inclusion();
			
			Calculate_Mshape_fine();
		}

		return error;
	}

	//copy roughness from another Roughness module (e.g. from a different mesh - the values, including shape, and discretisation are copied)
	BError copy_roughness(Roughness* pRough_copy_this);

	//clear roughness: set fine shape to coarse shape
	void clear_roughness(void);

	//-------------------Getters

	INT3 get_refine(void) { return refine; }

	//Get Mshape_fine for display
	VEC_VC<double>& GetRoughness(void) { return Mshape_fine; }
};

#else

class Roughness :
	public Modules
{

private:

	VEC_VC<double> Mshape_fine;

private:

public:

	Roughness(Mesh *pMesh_) {}
	~Roughness() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Setters

	//change cellsize refinement
	void set_refine(INT3 refine_) {}

	//called from Mesh::change_mesh_shape : apply shape to Mshape_fine, then set pMesh->M from it
	template  <typename Lambda, typename ... PType>
	BError change_mesh_shape(Lambda& run_this, PType& ... params) { return BError(); }

	//copy roughness from another Roughness module (e.g. from a different mesh - the values, including shape, and discretisation are copied)
	BError copy_roughness(Roughness* pRough_copy_this) { return BError(); }

	//clear roughness: set fine shape to coarse shape
	void clear_roughness(void) {}

	//-------------------Getters

	INT3 get_refine(void) { return INT3(); }

	//Get Mshape_fine for display
	VEC_VC<double>& GetRoughness(void) { return Mshape_fine; }
};

#endif