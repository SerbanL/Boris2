#pragma once

#include "BorisLib.h"

#include "MaterialParameter.h"
#include "ErrorHandler.h"
#include "Boris_Enums_Defs.h"
#include "MeshDefs.h"
#include "OVF2_Handlers.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	
//		Abstract base class for mesh parameters, applicable for both micromagnetic and atomistic meshes
//		This is inherited by MeshBase so that public MeshParams methods can be called using MeshBase pointer
//		It is also inherited by MeshParams (and Atom_MeshParams), as they implement this interface
//
//		Note the implementation is done through the virtual inheritance mechanism, i.e.:
//		MeshParamsBase -(virtual)-> MeshBase -> (Atom_)Mesh -> FinalMesh
//		MeshParamsBase -(virtual)-> (Atom_)MeshParams -> (Atom_)Mesh -> FinalMesh
//
//		This technique allows the pure virtual methods here to be visible at MeshBase level, and yet
//		the whole structure can still be compiled as FinalMesh is now a proper implementation,
//		since the pure virtual methods will appear to it as implemented in MeshParams.
//		Without the virtual inheritance the pure virtual methods inherited in MeshBase would
//		not be implemented, so the structure wouldn't compile.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MeshParamsBase
{

private:

protected:

	//Why is a void* needed? It's because we want to access methods in implementations of this abstract base class, and cannot use the much better CRTP mechanism.
	//CRTP would require MeshBase to be templated also, otherwise it's not a useful abstract base class for both micromagnetic and atomistic meshes, and then we cannot store pointers to it in a vector.
	//One workaround is the void*, requiring casting to the correct implementation using a switch based on meshType : see run_on_param_switch method
	//This feels like a very ugly solution, and goes against many design principles, but can't help feel it's the right way to go lacking other solutions - didn't take this design decision lightly!
	//The reason is really to do with maintainability, not so much code reuse per se : if you need to add new params manipulation methods, or tweak an old one, you can simply do it just once here, not twice in MeshParams, and in Atom_MeshParams. 
	//Also what if in the future you might decide to add a different base type of mesh apart from micromagnetic or atomistic ones? The code wouldn't grow gracefuly and would become a nightmare to maintain.
	//In some ways it's even safer : what if you forget to make the same change everywhere? 
	MESHTYPE_ meshBaseType;
	void* pImplementation;

	//list of all mesh parameters storing the units, handles and parameter type : used for console info display
	vector_key_lut<MeshParamDescriptor> meshParams;

public:

	//-------------------------------- LIST ALL MESH PARAMETERS HERE (which really need to be here, and not in the implementation)

	//the mesh base temperature (K)
	double base_temperature = 0.0;

	//text equation for base temperature dependence, allowing time dependence; stage step (Ss) is introduced as user constant.
	TEquation<double> T_equation;

	//Heat source stimulus in heat equation. Ideally used with a spatial variation. (W//m3)
	MatP<double, double> Q = 0.0;

	//set temperature spatial variation coefficient (unitless) - used with temperature settings in a simulation schedule only, not with console command directly
	MatP<double, double> cT = 1.0;

private:

	//-------------------------Parameter control

	//casts pImplementation to the correct implementation based on meshType, then calls the run_this code through the run_on_param templated function defined in the implementation
	//This is exactly the problem: run_on_param is templated, so cannot be made pure virtual as the compiler needs to know what it is statically.
	template <typename RType, typename Lambda, typename ... PType>
	RType run_on_param_switch(PARAM_ paramID, Lambda& run_this, PType& ... run_this_args);

protected:

	//-------------------------Setters

	//called in implementation constructor : (Atom_)MeshParams
	//remember a class inherited through virtual inheritance is constructed first from the class with highest inheritance
	//this means it will call the void constructor, unless you specifically ask it there to call another constructor
	//calling the base class (virtually inherited) constructor anywhere else is wasted code
	//Calling the base constructor with a parameter at the class with highest inheritance is really ugly, and unsafe (what if you forget to do it when you define a new implementation? - there are loads of them : FMesh, AFMesh, FCC, BCC, HCP, etc. etc.)
	//Much better to do it in (Atom_)MeshParams once and forget about it.
	template <typename Implementation>
	void set_meshparamsbase_implementation(Implementation* pImplementation_, MESHTYPE_ meshBaseType_)
	{
		pImplementation = pImplementation_;
		meshBaseType = meshBaseType_;
	}

public:

	//------------------------CTOR/DTOR

	MeshParamsBase(void) :
		T_equation({ "t" })
	{}
	
	virtual ~MeshParamsBase() {}

	//-------------------------Getters

	//return number of mesh parameters
	int get_num_meshparams(void) { return meshParams.size(); }

	//get id of indexed mesh parameter (this is the value from PARAM_ enum)
	int get_meshparam_id(int index) { return meshParams.get_ID_from_index(index); }
	int get_meshparam_id(std::string paramHandle) { return meshParams.get_ID_from_key(paramHandle); }

	//get handle of indexed mesh parameter
	std::string get_meshparam_handle(int index) { return meshParams.get_key_from_index(index); }
	std::string get_meshparam_handle(PARAM_ paramID) { return meshParams.get_key_from_ID(paramID); }

	//get unit of indexed mesh parameter
	std::string get_meshparam_unit(int index) { return meshParams[index].unit; }
	std::string get_meshparam_unit(PARAM_ paramID) { return meshParams(paramID).unit; }

	bool contains_param(PARAM_ paramID) { return meshParams.is_ID_set(paramID); }
	bool contains_param(std::string paramHandle) { return meshParams.has_key(paramHandle); }

	//check if parameter changes magnetization length
	bool param_changes_mlength(PARAM_ paramID) { return (paramID == PARAM_MS || paramID == PARAM_MS_AFM || paramID == PARAM_ATOM_SC_MUS); }

	//get value of indexed mesh parameter as a std::string (with unit)
	std::string get_meshparam_value(int index);
	std::string get_meshparam_value(PARAM_ paramID);

	//get value of indexed mesh parameter as a std::string (without unit)
	std::string get_meshparam_value_sci(int index);
	std::string get_meshparam_value_sci(PARAM_ paramID);

	PARAMTYPE_ get_meshparam_type(PARAM_ paramID) { return meshParams(paramID).get_type(); }

	//returns a std::string describing the set temperature dependence ("none", "array" or set equation : "name parameters...") 
	std::string get_paraminfo_string(PARAM_ paramID);

	//returns a std::string describing the set spatial dependence with any parameters
	std::string get_paramvarinfo_string(PARAM_ paramID);

	//check if the given parameter has a temperature dependence set
	bool is_paramtemp_set(PARAM_ paramID);
	//check if the given parameters has a temperature dependence specified using a text equation
	bool is_paramtempequation_set(PARAM_ paramID);
	//check if the given parameter has a  spatial variation set
	bool is_paramvar_set(PARAM_ paramID);
	//check if the spatial dependence is set using a text equation
	bool is_paramvarequation_set(PARAM_ paramID);
	//check if scaling array is scalar (or else vectorial)
	bool is_paramvar_scalar(PARAM_ paramID);
	//check if the given parameter has a  temperature dependence or a spatial variation set
	bool is_param_nonconst(PARAM_ paramID);

	//get mesh parameter temperature scaling up to max_temperature : return a vector from 0K up to and including max_temperature with scaling coefficients
	bool get_meshparam_tempscaling(PARAM_ paramID, double max_temperature, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z);

	//is this param hidden or can we display it?
	bool is_param_hidden(PARAM_ paramID) { return meshParams(paramID).hidden; }

	//-------------------------Spatial scaling VEC get / calculate methods

	//get reference to mesh parameter spatial scaling VEC
	//use void* since the scaling is templated. Caller must recast it correctly - see is_paramvar_scalar method
	void* get_meshparam_s_scaling(PARAM_ paramID);

	//get value of mesh parameter spatial scaling coefficient at given position (and time)
	Any get_meshparam_s_scaling_value(PARAM_ paramID, DBL3 rel_pos, double stime);

	//calculate spatial variation into the provided VECs - intended to be used when the spatial variation is set using a text equation
	void calculate_meshparam_s_scaling(PARAM_ paramID, VEC<double>& displayVEC_SCA, double stime);
	void calculate_meshparam_s_scaling(PARAM_ paramID, VEC<DBL3>& displayVEC_VEC, double stime);

	//-------------------------Setters : value and temperature dependence

	//set value from std::string for named parameter (units allowed in std::string)
	void set_meshparam_value(PARAM_ paramID, std::string value_text);

	//clear mesh parameter temperature dependence
	void clear_meshparam_temp(PARAM_ paramID);

	//set mesh parameter array scaling
	bool set_meshparam_tscaling_array(PARAM_ paramID, std::vector<double>& temp, std::vector<double>& scaling_x, std::vector<double>& scaling_y, std::vector<double>& scaling_z);

	//set temperature dependence info std::string for console display purposes
	void set_meshparam_tscaling_info(PARAM_ paramID, std::string info_text);

	//-------------------------Setters : spatial variation

	//set the mesh parameter spatial variation equation with given user constants
	void set_meshparam_s_equation(PARAM_ paramID, std::string& equationText, vector_key<double>& userConstants, DBL3 meshDimensions);

	//clear mesh parameter spatial variation (all if paramID == PARAM_ALL)
	void clear_meshparam_variation(PARAM_ paramID);

	//update mesh parameter spatial variation (e.g. cellsize or rectangle could have changed)
	bool update_meshparam_var(PARAM_ paramID, DBL3 h, Rect rect);

	//set parameter spatial variation using a given generator and arguments (arguments passed as a std::string to be interpreted and converted using ToNum)
	BError set_meshparam_var(PARAM_ paramID, MATPVAR_ generatorID, DBL3 h, Rect rect, std::string generatorArgs, std::function<std::vector<unsigned char>(std::string, INT2)>& bitmap_loader);

	//set parameter spatial variation using a shape : set value in given shape only
	BError set_meshparam_shape(PARAM_ paramID, DBL3 h, Rect rect, std::vector<MeshShape> shapes, std::string value_text);

	//-------------------------Setters/Updaters : text equations

	//set the mesh parameter temperature equation with given user constants
	virtual void set_meshparam_t_equation(PARAM_ paramID, std::string& equationText, vector_key<double>& userConstants) = 0;

	//update text equations for mesh parameters with user constants, mesh dimensions, Curie temperature, base temperature
	virtual bool update_meshparam_equations(PARAM_ paramID, vector_key<double>& userConstants, DBL3 meshDimensions) = 0;

	//-------------------------General Updaters

	//call this to update given parameter output value to current base_temperature
	void update_parameters(PARAM_ paramID = PARAM_ALL);

	//-------------------------Special setters/getters

	virtual std::string get_tensorial_anisotropy_string(void) { return ""; }
};