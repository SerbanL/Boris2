#pragma once

#include "BorisLib.h"

#include "ModulesDefs.h"
#include "ParametersDefs.h"
#include "MeshDefs.h"

#include "SimulationData.h"

class StageConfig;

//Container for Simulation class data (Simulation inherits from this) to which other lower ranked objects need access - inherit from this (lower ranked means an object which is a member of Simulation, or a member of a member of Simulation, etc.)
//All members are static, built in the SimulationSharedData constructor
//!!NOTE!!! data in SimulationSharedData could be built in the Simulation constructor, but this may cause problems due to object initialization order. 
//e.g. if Simulation has an object which depends on the data in SimulationSharedData, it will be initialized before the data held here, potentially causing bugs. This can be avoided by declaring that object as a pointer in Simulation and
//making it only after the data here has been initialized, or initializing the data in SimulationSharedData rather than Simulation - this seems preferable as it's foolproof.
class SimulationSharedData {

protected:

	//material temperature dependence types : stored std::string : text to appear in console, lut : entry from MATPTDEP_ enum
	static vector_lut<std::string> temperature_dependence_type;

	//available pre-set spatial variation generators to use for material parameters spatial dependence
	//key: name of generator to appear in console, lut: a MATPVAR_ entry, stored std::string : default parameters for generator stored as text, usable with ToNum
	static vector_key_lut<std::string> vargenerator_descriptor;

	//unit of displayed mesh quantity
	static vector_lut<std::string> meshDisplay_unit;
	
	//quantities to display in meshes
	static vector_lut<std::string> displayHandles;

	//available modules for each mesh type
	static vector_lut< std::vector<MOD_> > modules_for_meshtype;
	//subset of modules_for_meshtype with modules only suitable for effective field display
	static vector_lut< std::vector<MOD_> > displaymodules_for_meshtype;

	//allowed quantities to display for each mesh type
	static vector_lut< std::vector<MESHDISPLAY_> > meshAllowedDisplay;

	//available material parameters for each mesh type
	static vector_lut< std::vector<PARAM_> > params_for_meshtype;

	//entries from PARAM_, specifying if temperature dependence (first) and spatial variation (second) are enabled.
	static vector_lut<std::pair<bool, bool>> params_enabled_props;

	//modules which cannot be available simultaneously : exclusive modules.
	static exclusions<MOD_> exclusiveModules;

	//modules that run on the super-mesh will have counterparts in individual meshes, which should not run
	static exclusions<MOD_> superMeshExclusiveModules;

	//some modules can only be active if there is a companion super-mesh module active also
	static exclusions<MOD_> superMeshCompanionModules;

	//list of data to output during a simulation (e.g. to output file and/or to processing arrays)
	static vector_lut<DatumConfig> saveDataList;

	//simulation stages describing the simulation schedule
	static vector_lut<StageConfig> simStages;

	//modifier for mesh shaping functions
	//all general shaping functions go through Mesh::change_mesh_shape to apply the shape
	//Note, there are specific exceptions e.g. grain generators which only work on M)
	//by default this flag is false, meaning the shape is carried over to all relevant primary quantities in the mesh, e.g. M, elC, Temp (if enabled).
	//if this flag is set to true, then Mesh::change_mesh_shape only applies the shape to the quantity currently focused in the display, e.g. only M, only elC, or only Temp
	//this allows setting different shapes for these quantities in the same mesh
	static bool shape_change_individual;

	//display transparency for foreground (primary display - major) and background (minor); useful if you want to display both background and foreground and see them both
	//only takes effect for dual display : foreground has some transparency by default. Values from 0 (fully transparent) to 1 (opaque).
	static DBL2 displayTransparency;

	//minimum and maximum values for mesh display; if both set to 0 0 then ignore limits.
	static DBL2 displayThresholds;

	//set which component to trigger display thresholds on when using vector quantities (x, y, z, or magnitude)
	static int displayThresholdTrigger;

	//working directory
	static std::string directory;

	//free and total memory for gpu and cpu - update these as memory is allocated / freed so we can keep track of it in console interactive objects without having to interrogate at every refresh
	static size_t gpuMemFree_MB;
	static size_t gpuMemTotal_MB;
	static size_t cpuMemFree_MB;
	static size_t cpuMemTotal_MB;

	//collect available cuda device major versions here, indexed from 0
	static std::vector<std::pair<int, std::string>> cudaDeviceVersions;

public:

	//if set to true, the transport solver is iterated only at the end of a step or stage
	//this is useful for static problems, e.g. MR loops, where it's more efficient to have the magnetization state solved separately.
	static bool static_transport_solver;

	//if set to true, the transport solver iteration is entirely disabled, so all transport quantities remain constant.
	static bool disabled_transport_solver;

	//enable all cuda computations?
	static bool cudaEnabled;

	//there may be more than 1 cuda devices, numbered from 0 up. This does not switch cuda on, but indicates which devices is used when cuda is switched on
	//always check the current program cuda architecture matches the device cuda compute version when switching on
	static int cudaDeviceSelect;

	//simulation schedule indexes : stage (main index in simStages) and step (sub-index in stage)
	static INT2 stage_step;
	//run a single simulation stage? (default is false, but CMD_RUNSTAGE can set this true, and when simulation stops this flag is set back to false)
	static bool single_stage_run;

	//constants defined by user at runtime to be used with TEquation objects; the key is the user constant name
	static vector_key<double> userConstants;

public:

	//initialize all data here, not in the derived Simulation class - see NOTE above
	SimulationSharedData(bool called_from_Simulation = false);

	virtual ~SimulationSharedData() {}
};