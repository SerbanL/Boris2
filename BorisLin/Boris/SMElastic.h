#pragma once

#include "BorisLib.h"
#include "Modules.h"

class SuperMesh;
class MElastic;

#ifdef MODULE_COMPILATION_MELASTIC

#if COMPILECUDA == 1
#include "SMElasticCUDA.h"
#endif

class SMElastic :
	public Modules,
	public ProgramState<SMElastic, std::tuple<double, bool, vector_lut<Rect>, vector_lut<Rect>, std::vector<std::string>>, std::tuple<>>
{

#if COMPILECUDA == 1
	friend SMElasticCUDA;
#endif

	friend MElastic;

private:

	//pointer to supermesh
	SuperMesh* pSMesh;

	//---------------------- CMBND data

	//CMBND contacts for all contacting elastic meshes - these are ordered by first vector index; for each mesh there could be multiple contacting meshes and these are ordered by second vector index
	std::vector< std::vector<CMBNDInfo> > CMBNDcontacts;

	//list of all MElastic modules in meshes (same ordering as first vector in CMBNDcontacts)
	std::vector<MElastic*> pMElastic;

	//vector of pointers to all u_disp - need this to set cmbnd flags (same ordering as first vector in CMBNDcontacts)
	std::vector<VEC_VC<DBL3>*> pu_disp;

	//----------------------

	//Fixed surfaces (plane rectangles in absolute coordinates on mesh surfaces where we set zero Dirichlet value for displacement and velocity)
	vector_lut<Rect> fixed_u_surfaces;

	//external stress surfaces (plane rectangles in absolute coordinates on mesh surfaces where we have an external stimulus)
	vector_lut<Rect> stress_surfaces_rect;
	//external stimulus specified using a text vector equation (same size as stress_surfaces_rect)
	std::vector<std::string> stress_surfaces_equations;

	//----------------------

	//time step for the elastic solver equations - if in a magnetic mesh must always be smaller or equal to dT (the magnetization equation time-step)
	double el_dT = 0.0;

	//if true, then set el_dT same as the dT used in the magnetic equation
	bool linked_el_dT = true;

	//save the last magnetic dT used: when advancing the elastic solver this is the time we need to advance by. 
	//Update magnetic_dT after each elastic solve advance (in case an adaptive time-step method is used for the magnetic part).
	double magnetic_dT;

public:

	SMElastic(SuperMesh *pSMesh_);
	~SMElastic() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Getters

	double get_el_dT(void) { return el_dT; }
	bool get_linked_el_dT(void) { return linked_el_dT; }

	//-------------------Setters

	void set_el_dT(double dT) { el_dT = dT; linked_el_dT = false; }
	void set_linked_el_dT(bool flag);

	//------------------- Fixed and Stress Surfaces

	//--------- Fixed

	BError Add_Fixed_Surface(Rect surface_rect);

	int Get_Num_Fixed_Surfaces(void) { return fixed_u_surfaces.size(); }
	Rect Get_Fixed_Surface(int idx) { if (GoodIdx(fixed_u_surfaces.size(), idx)) return fixed_u_surfaces[idx]; else return Rect(); }

	void Del_Fixed_Surface(int idx) { if (GoodIdx(fixed_u_surfaces.size(), idx)) fixed_u_surfaces.erase(idx); }

	void Edit_Fixed_Surface(int idx, Rect rect) { if (GoodIdx(fixed_u_surfaces.size(), idx)) fixed_u_surfaces[idx] = rect; }

	//get fixed surface index for surfgace with given minor Id (major Id is 0)
	int Get_Fixed_Surface_Index(int minorId) { return fixed_u_surfaces.get_index_from_id(INT2(0, minorId)); }

	//get the minor id of fixed surface with given index
	int Get_Fixed_Surface_id(int index) { if (GoodIdx(fixed_u_surfaces.size(), index)) return fixed_u_surfaces.get_id_from_index(index).minor; else return -1; }

	//--------- Stress

	BError Add_Stress_Surface(Rect surface_rect, std::string equation);

	int Get_Num_Stress_Surfaces(void) { return stress_surfaces_rect.size(); }
	Rect Get_Stress_Surface(int idx) { if (GoodIdx(stress_surfaces_rect.size(), idx)) return stress_surfaces_rect[idx]; else return Rect(); }
	std::string Get_Stress_Surface_Equation(int idx) { if (GoodIdx(stress_surfaces_equations.size(), idx)) return stress_surfaces_equations[idx]; else return ""; }

	void Del_Stress_Surface(int idx)
	{
		if (GoodIdx(stress_surfaces_rect.size(), idx)) {

			stress_surfaces_rect.erase(idx);
			stress_surfaces_equations.erase(stress_surfaces_equations.begin() + idx);
		}
	}

	void Edit_Stress_Surface_Rectangle(int idx, Rect rect) { if (GoodIdx(stress_surfaces_rect.size(), idx)) stress_surfaces_rect[idx] = rect; }
	void Edit_Stress_Surface_Equation(int idx, std::string equation) { if (GoodIdx(stress_surfaces_equations.size(), idx)) stress_surfaces_equations[idx] = equation; }

	//get stress surface index for surfgace with given minor Id (major Id is 0)
	int Get_Stress_Surface_Index(int minorId) { return stress_surfaces_rect.get_index_from_id(INT2(0, minorId)); }

	//get the minor id of stress surface with given index
	int Get_Stress_Surface_id(int index) { if (GoodIdx(stress_surfaces_rect.size(), index)) return stress_surfaces_rect.get_id_from_index(index).minor; else return -1; }
};

#else

class SMElastic :
	public Modules
{

private:

private:

public:

	SMElastic(SuperMesh *pSMesh_) {}
	~SMElastic() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Getters

	double get_el_dT(void) { return 0.0; }

	//-------------------Setters

	void set_el_dT(double dT) {}

	//------------------- Fixed and Stress Surfaces

	BError Add_Fixed_Surface(Rect surface_rect) { return BError(); }

	int Get_Num_Fixed_Surfaces(void) { return 0; }
	Rect Get_Fixed_Surface(int idx) { return Rect(); }
	void Del_Fixed_Surface(int idx) {}

	BError Add_Stress_Surface(Rect surface_rect, std::string equation) { return BError(); }

	int Get_Num_Stress_Surfaces(void) { return 0; }
	Rect Get_Stress_Surface(int idx) { return Rect; }
	std::string Get_Stress_Surface_Equation(int idx) { return ""; }
	void Del_Stress_Surface(int idx) {}
};

#endif