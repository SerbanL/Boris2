#pragma once

#include "BorisLib.h"

using namespace std;

enum MOD_;
enum PARAM_;
enum MESHDISPLAY_;

class DatumConfig;

//Container for Simulation class data (Simulation inherits from this) to which other lower ranked objects need access - inherit from this (lower ranked means an object which is a member of Simulation, or a member of a member of Simulation, etc.)
//All members are static, built in the SimulationSharedData constructor
//!!NOTE!!! data in SimulationSharedData could be built in the Simulation constructor, but this may cause problems due to object initialization order. 
//e.g. if Simulation has an object which depends on the data in SimulationSharedData, it will be initialized before the data held here, potentially causing bugs. This can be avoided by declaring that object as a pointer in Simulation and
//making it only after the data here has been initialized, or initializing the data in SimulationSharedData rather than Simulation - this seems preferable as it's foolproof.
class SimulationSharedData {

protected:

	//available pre-set formulas to use for material parameters temperature dependence
	//key: name of formula to appear in console, lut: a MATPFORM_ entry, stored int : number of coefficients required for the formula
	static vector_key_lut<int> formula_descriptor;

	//available pre-set spatial variation generators to use for material parameters spatial dependence
	//key: name of generator to appear in console, lut: a MATPVAR_ entry, stored string : default parameters for generator stored as text, usable with ToNum
	static vector_key_lut<string> vargenerator_descriptor;

	//unit of displayed mesh quantity
	static vector_lut<string> meshDisplay_unit;
	
	//quantities to display in meshes
	static vector_lut<string> displayHandles;

	//available modules for each mesh type
	static vector_lut< vector<MOD_> > modules_for_meshtype;

	//allowed quantities to display for each mesh type
	static vector_lut< vector<MESHDISPLAY_> > meshAllowedDisplay;

	//available material parameters for each mesh type
	static vector_lut< vector<PARAM_> > params_for_meshtype;

	//modules which cannot be available simultaneously : exclusive modules.
	static exclusions<MOD_> exclusiveModules;

	//modules that run on the super-mesh will have counterparts in individual meshes, which should not run
	static exclusions<MOD_> superMeshExclusiveModules;

	//some modules can only be active if there is a companion super-mesh module active also
	static exclusions<MOD_> superMeshCompanionModules;

	//list of data to output during a simulation (e.g. to output file and/or to processing arrays)
	static vector_lut<DatumConfig> saveDataList;

	//modifier for mesh shaping functions
	//all general shaping functions go through Mesh::change_mesh_shape to apply the shape
	//Note, there are specific exceptions e.g. grain generators which only work on M)
	//by default this flag is false, meaning the shape is carried over to all relevant primary quantities in the mesh, e.g. M, elC, Temp (if enabled).
	//if this flag is set to true, then Mesh::change_mesh_shape only applies the shape to the quantity currently focused in the display, e.g. only M, only elC, or only Temp
	//this allows setting different shapes for these quantities in the same mesh
	static bool shape_change_individual;

	//if set to true, the transport solver is iterated only at the end of a step or stage
	//this is useful for static problems, e.g. MR loops, where it's more efficient to have the magnetization state solved separately.
	static bool static_transport_solver;

	//enable all cuda computations?
	static bool cudaEnabled;

	//free and total memory for gpu and cpu - update these as memory is allocated / freed so we can keep track of it in console interactive objects without having to interrogate at every refresh
	static size_t gpuMemFree_MB;
	static size_t gpuMemTotal_MB;
	static size_t cpuMemFree_MB;
	static size_t cpuMemTotal_MB;

public:

	//initialize all data here, not in the derived Simulation class - see NOTE above
	SimulationSharedData(bool called_from_Simulation = false);

	virtual ~SimulationSharedData() {}
};