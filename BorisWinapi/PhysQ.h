#pragma once

#include "SimSharedData.h"

#include "BorisLib.h"

#include "PhysQDefs.h"

using namespace std;

//------------------------------- The physical quantities collected from meshes - to be computed into a physical quantity representation

class PhysQ :
	public SimulationSharedData
{

	//represent all possible quantites to display : VEC or VEC_VC with their associated template type (the stored type)
	enum VECTYPE_ { VECTYPE_VOID = 0, VECTYPE_VEC, VECTYPE_VEC_VC };
	enum VECSTYPE_ { VECSTYPE_NONE = 0, VECSTYPE_FLOAT, VECSTYPE_DOUBLE, VECSTYPE_FLT3, VECSTYPE_DBL3 };

private:

	//VEC or VEC_VC quantity : vtype specifies the exact type we need to cast this to
	void *pVEC = nullptr;

	//second vec to display if in dual display mode : must be same type as pVEC and same dimensions
	void *pVEC2 = nullptr;

	//representation of type of *pVEC
	VECTYPE_ vtype = VECTYPE_VOID;
	//and its stored type
	VECSTYPE_ vstype = VECSTYPE_NONE;

	//what is the displayed quantity? (this will be a value from MESHDISPLAY_ enum)
	unsigned displayedType;

	//Rectangle of physical quantity (VEC) : dimensions * cellsize. Still keep this in case both pvecQ and pscaQ are null but we want to display a mesh frame.
	Rect rect;

	//the cellsize
	DBL3 h;

	//number of cells
	SZ3 n;

	//the type of representation to use for vectorial quantities
	VEC3REP_ vec3rep = VEC3REP_FULL;

	//there should be one mesh which is in focus - this sets conversion factors for display and centers view on the focused mesh
	bool meshInFocus = false;

	//the mesh name from which this PhysQ was collected
	string meshName;

public:

	//------------------------------------------- Constructors

	//disabled mesh display
	PhysQ(Rect rect_, DBL3 h_, unsigned displayedType_);

	//constructors for various displayed types
	PhysQ(VEC<DBL3> *pVEC_, unsigned displayedType_, VEC3REP_ vec3rep_);
	PhysQ(VEC<FLT3> *pVEC_, unsigned displayedType_, VEC3REP_ vec3rep_);

	PhysQ(VEC<DBL3> *pVEC_, VEC<DBL3> *pVEC2_, unsigned displayedType_, VEC3REP_ vec3rep_);
	PhysQ(VEC<FLT3> *pVEC_, VEC<FLT3> *pVEC2_, unsigned displayedType_, VEC3REP_ vec3rep_);
	
	PhysQ(VEC<double> *pVEC_, unsigned displayedType_);
	PhysQ(VEC<float> *pVEC_, unsigned displayedType_);
	
	PhysQ(VEC_VC<DBL3> *pVEC_, unsigned displayedType_, VEC3REP_ vec3rep_);
	PhysQ(VEC_VC<FLT3> *pVEC_, unsigned displayedType_, VEC3REP_ vec3rep_);
	
	PhysQ(VEC_VC<DBL3> *pVEC_, VEC_VC<DBL3> *pVEC2_, unsigned displayedType_, VEC3REP_ vec3rep_);
	PhysQ(VEC_VC<FLT3> *pVEC_, VEC_VC<FLT3> *pVEC2_, unsigned displayedType_, VEC3REP_ vec3rep_);
	
	PhysQ(VEC_VC<double> *pVEC_, unsigned displayedType_);
	PhysQ(VEC_VC<float> *pVEC_, unsigned displayedType_);

	PhysQ(const PhysQ& copyThis) { *this = copyThis; }
	PhysQ& operator=(const PhysQ& copyThis);

	~PhysQ() {}

	//------------------------------------------- Access

	//get stored type casted to correct quantity for external use
	VEC<DBL3>* get_vec_dbl3(void) { return static_cast<VEC<DBL3>*>(pVEC); }
	VEC_VC<DBL3>* get_vec_vc_dbl3(void) { return static_cast<VEC_VC<DBL3>*>(pVEC); }

	VEC<FLT3>* get_vec_flt3(void) { return static_cast<VEC<FLT3>*>(pVEC); }
	VEC_VC<FLT3>* get_vec_vc_flt3(void) { return static_cast<VEC_VC<FLT3>*>(pVEC); }

	VEC<DBL3>* get2_vec_dbl3(void) { return static_cast<VEC<DBL3>*>(pVEC2); }
	VEC_VC<DBL3>* get2_vec_vc_dbl3(void) { return static_cast<VEC_VC<DBL3>*>(pVEC2); }

	VEC<FLT3>* get2_vec_flt3(void) { return static_cast<VEC<FLT3>*>(pVEC2); }
	VEC_VC<FLT3>* get2_vec_vc_flt3(void) { return static_cast<VEC_VC<FLT3>*>(pVEC2); }

	VEC<double>* get_vec_double(void) { return static_cast<VEC<double>*>(pVEC); }
	VEC_VC<double>* get_vec_vc_double(void) { return static_cast<VEC_VC<double>*>(pVEC); }

	VEC<float>* get_vec_float(void) { return static_cast<VEC<float>*>(pVEC); }
	VEC_VC<float>* get_vec_vc_float(void) { return static_cast<VEC_VC<float>*>(pVEC); }

	//get value at coordinate, where we know if the stored type is vectorial or scalar (e.g. DBL3 vs double)
	DBL3 get_vec_point(DBL3 coordinate);
	double get_sca_point(DBL3 coordinate);

	//------------------------------------------- Properties

	//set meshInFocus flag as well as mesh name
	PhysQ& set_focus(bool status, string meshName_) { meshInFocus = status; meshName = meshName_; return *this; }
	bool is_focused(void) { return meshInFocus; }

	string get_name(void) { return meshName; }
	string get_unit(void) { return meshDisplay_unit(displayedType); }
	string get_type_name(void) { return displayHandles(displayedType); }
	int get_type(void) { return displayedType; }

	VEC3REP_ get_vec3rep(void) { return vec3rep; }

	SZ3 size(void);

	Rect rectangle(void) { return rect; }
	DBL3 cellsize(void) { return h; }

	//stored types enquiries
	bool is_vectorial(void) { return (vstype == VECSTYPE_FLT3 || vstype == VECSTYPE_DBL3); }
	
	bool is_vec_vc(void) { return (vtype == VECTYPE_VEC_VC); }

	bool is_double_precision(void) { return (vstype == VECSTYPE_DBL3 || vstype == VECSTYPE_DOUBLE); }

	bool is_empty(void) { return (pVEC == nullptr || !n.dim()); }

	bool is_dual(void) { return pVEC2 != nullptr; }
};