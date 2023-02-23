#pragma once

#include "PhysQ.h"

#include "CompileFlags.h"
#if GRAPHICS == 1

//------------------------------------------- Stored component in Physical Quantity Representation

struct PhysQRepComponent :
	public ProgramState<PhysQRepComponent,
	std::tuple<bool, double, int, bool, std::string, std::string, std::string>, std::tuple<>>
{
	//scale the dimension of drawn objects to the largest value (modulus) in the physical quantity?
	bool scale_to_magnitude = true;

	//exponent of scaling factor when scale_to_magnitude is used (1 results in linear scaling. Looks better with lower values when physical quantities have large extremes).
	double exponent = 0.2;

	//object selector to draw the transformBatch with
	int obSelector = 0;

	//the calculated batch of transforms used in rendering - its rect is the same as the PhysQ rect, its cellsize depends on the level of detail used to generate this PhysQRepComponent
	VEC<CBObjectTransform> transformBatch;
	std::vector<bool> emptyCell;
	//if there are many cells to draw then use transformBatch_render instead: the problem is transformBatch can also have many cells marked not to be rendered, and we will be looping over them in the drawing routine, wasting iterations.
	//instead save only the cells marked to be rendered to this vector
	std::vector<CBObjectTransform> transformBatch_render;

	//drawe mesh outline using this calculated transform applied to CDO_CUBEFRAME
	CBObjectTransform meshFrame = CBObjectTransform();
	bool drawFrame = true;

	//unit of PhysQ values form which this PhysQRepComponent was generated
	std::string unit;

	//the displayed type name (e.g. magnetization)
	std::string typeName;

	//the mesh name from which the PhysQ which generated this PhysQRepComponent was collected
	std::string meshName;

	//the type of representation to use for vectorial quantities (specified when a PhysQ is returned, so value set here before calculating a representation for it)
	VEC3REP_ vec3rep;

	PhysQRepComponent(void) :
		ProgramStateNames(this,
			{
				VINFO(scale_to_magnitude), VINFO(exponent), VINFO(obSelector), VINFO(drawFrame), VINFO(unit), VINFO(typeName), VINFO(meshName)
			}, {})
	{}

	PhysQRepComponent(const PhysQRepComponent& copyThis) :
		ProgramStateNames(this,
			{
				VINFO(scale_to_magnitude), VINFO(exponent), VINFO(obSelector), VINFO(drawFrame), VINFO(unit), VINFO(typeName), VINFO(meshName)
			}, {})
	{
		scale_to_magnitude = copyThis.scale_to_magnitude;
		exponent = copyThis.exponent;

		obSelector = copyThis.obSelector;

		transformBatch = copyThis.transformBatch;
		transformBatch_render = copyThis.transformBatch_render;
		emptyCell = copyThis.emptyCell;

		meshFrame = copyThis.meshFrame;
		drawFrame = copyThis.drawFrame;

		unit = copyThis.unit;
		typeName = copyThis.typeName;
		meshName = copyThis.meshName;
	}

	void RepairObjectState(void) 
	{ 
		//this must be recalculated (will be done on first call to update the physical quantity representation)
		transformBatch.clear();
		transformBatch_render.clear();
		emptyCell.clear();
	}			
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//collection of graphical display settings
struct PhysQRepSettings {

	//1. rectangle of mesh in focus - the 3D coordinate system origin is placed at the center of this rect (when converted to logical units)
	Rect focusRect;

	//2. conversion factor from meters to logical units : multiply mesh dimensions by this to obtain logical dimensions
	double m_to_l;

	//3. level of displayed detail (sets the number of display elements per logical unit)
	double detail_level;

	//4. Eye position (camX, camY, camZ) - use D3D::SetCameraPosition
	float camX, camY, camZ;

	//5. Up position (camUx, camUy, camUz) - use D3D::SetCameraUp
	float camUx, camUy, camUz;

	//6. View shift (view_shiftX, view_shiftY, view_shiftZ) - use D3D::Set3DOriginPixelPosition
	float view_shiftX, view_shiftY, view_shiftZ;

	//7. FOV
	float fov;

	//---------------------------------------------------------------------------------------------

	//calculate settings for default view (these are the settings set when the mesh focus changes)
	PhysQRepSettings(std::vector<PhysQ>& physQ, UINT wndWidth, D2D1_RECT_F spaceRect);

	PhysQRepSettings(
		Rect focusRect_, double m_to_l_, double detail_level_,
		FLT3 camPos, FLT3 camUp, FLT3 viewShift, float fov_)
	{
		focusRect = focusRect_;
		m_to_l = m_to_l_;
		detail_level = detail_level_;
		camX = camPos.x; camY = camPos.y; camZ = camPos.z;
		camUx = camUp.x; camUy = camUp.y; camUz = camUp.z;
		view_shiftX = viewShift.x;  view_shiftY = viewShift.y; view_shiftZ = viewShift.z;
		fov = fov_;
	}

	PhysQRepSettings(const PhysQRepSettings& copyThis) { *this = copyThis; }

	PhysQRepSettings& operator=(const PhysQRepSettings& copyThis)
	{
		focusRect = copyThis.focusRect;
		m_to_l = copyThis.m_to_l;
		detail_level = copyThis.detail_level;
		camX = copyThis.camX; camY = copyThis.camY; camZ = copyThis.camZ;
		camUx = copyThis.camUx; camUy = copyThis.camUy; camUz = copyThis.camUz;
		view_shiftX = copyThis.view_shiftX;  view_shiftY = copyThis.view_shiftY; view_shiftZ = copyThis.view_shiftZ;
		fov = copyThis.fov;

		return *this;
	}

	PhysQRepSettings operator*(double mul)
	{
		return PhysQRepSettings(
			this->focusRect * mul, this->m_to_l * mul, this->detail_level * mul,
			FLT3(this->camX, this->camY, this->camZ) * mul,
			FLT3(this->camUx, this->camUy, this->camUz) * mul,
			FLT3(this->view_shiftX, this->view_shiftY, this->view_shiftZ) * mul,
			this->fov * mul);
	}

	PhysQRepSettings operator+(const PhysQRepSettings& rhs)
	{
		return PhysQRepSettings(
			this->focusRect + rhs.focusRect, this->m_to_l + rhs.m_to_l, this->detail_level + rhs.detail_level,
			FLT3(this->camX + rhs.camX, this->camY + rhs.camY, this->camZ + rhs.camZ),
			FLT3(this->camUx + rhs.camUx, this->camUy + rhs.camUy, this->camUz + rhs.camUz),
			FLT3(this->view_shiftX + rhs.view_shiftX, this->view_shiftY + rhs.view_shiftY, this->view_shiftZ + rhs.view_shiftZ),
			this->fov + rhs.fov);
	}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//------------------------------- The Physical Quantity Representation

class PhysQRep :
	public GraphicalObject,
	public ProgramState<PhysQRep,
	std::tuple<double, Rect, double, int, int, int, std::string, std::string, DBL2, std::string, std::vector<PhysQRepComponent>>, std::tuple<>>
{
private:

	OmpReduction<double> omp_reduction, omp_reduction2;

	//conversion factor from meters to logical units : multiply mesh dimensions by this to obtain logical dimensions
	double m_to_l = 1e9;

	//rectangle of mesh in focus - the 3D coordinate system origin is placed at the center of this rect (when converted to logical units)
	Rect focusRect;

	//level of displayed detail. This is the side of a cubic mesh cell. The number of displayed elements in each mesh is obtained by dividing its rectangle by this cubic cell (to the nearest integer in each dimension).
	double detail_level = 5e-9;

	//for vector quantity display switch to using simpler elements if number of displayed cells exceeds this amount; value of zero disables this mode
	int renderspeedup1_cells = NUMDRAWNCELLS_RENDERSPEEDUP1;
	//for vector and scalar quantity display switch to not drawing completely surrounded cells if number of displayed cells exceeds this amount; value of zero disables this mode
	int renderspeedup2_cells = NUMDRAWNCELLS_RENDERSPEEDUP2;
	//for vector quantity display only draw cells on a checkerboard pattern
	int renderspeedup3_cells = NUMDRAWNCELLS_RENDERSPEEDUP3;

	//the name of mesh in focus
	std::string meshName_focused;

	//the name of displayed type for mesh in focus
	std::string typeName_focused;

	//minimum and maximum values in focused mesh (if vectorial quantities these are magnitudes)
	DBL2 minmax_focusedmeshValues = DBL2();

	//unit of minmax values
	std::string unit;

	//computed physical representations for each mesh with displayed quantities
	std::vector<PhysQRepComponent> physQRep;

	//vector used for obtaining info when mouse hovers over displayed mesh - declared here so we don't have to allocate them every time; contains distance, mesh index, intersection point
	std::vector<std::tuple<double, int, DBL3>> mouseInfo;

	//highlight cell mouse is hovering over
	bool highlightCell = false;
	CBObjectTransform highlightedCell = CBObjectTransform();

private:

	//calculate representation for a component and return minimum and maximum values (either of components, or magnitude), as well as maximum magnitude value - for VEC<DBL3> / VEC_VC<DBL3>
	template <typename VECType>
	DBL2 CalculateRepresentation_VEC(VECType* pQ, PhysQRepComponent& physQRepComponent, double maxTransparency, DBL2 displayThresholds, int displayThresholdTrigger);

	//dual VEC display mode
	template <typename VECType>
	DBL2 CalculateRepresentation_VEC(VECType* pQ, VECType* pQ2, PhysQRepComponent& physQRepComponent, double maxTransparency, DBL2 displayThresholds, int displayThresholdTrigger);

	//after global minimum and maximum found use this to make adjustments - scaling changes
	void AdjustMagnitude_VEC(PhysQRepComponent& physQRepComponent, DBL2 minmax);

	//calculate representation for a component and return minimum and maximum values - for VEC<double> / VEC_VC<double>
	template <typename VECType>
	DBL2 CalculateRepresentation_SCA(VECType* pQ, PhysQRepComponent& physQRepComponent, double maxTransparency, DBL2 displayThresholds);

	//after global minimum and maximum found use this to make adjustments - coloring changes
	void AdjustMagnitude_SCA(PhysQRepComponent& physQRepComponent, DBL2 minmax);

	//get alpha transparency value from given object position in world space
	float get_alpha_value(DBL3 position);

public:

	//------------------------------------------- Constructor

	PhysQRep(void) :
		GraphicalObject(),
		ProgramStateNames(this, { 
		VINFO(m_to_l), VINFO(focusRect), 
		VINFO(detail_level), VINFO(renderspeedup1_cells), VINFO(renderspeedup2_cells), VINFO(renderspeedup3_cells), 
		VINFO(meshName_focused), VINFO(typeName_focused), VINFO(minmax_focusedmeshValues), VINFO(unit), 
		VINFO(physQRep) }, {})
	{
	}

	void RepairObjectState(void) { highlightCell = false; highlightedCell = CBObjectTransform(); }

	//------------------------------------------- Properties

	size_t size(void) { return physQRep.size(); }
	
	std::string get_focused_meshName(void) { return meshName_focused; }
	std::string get_focused_typeName(void) { return typeName_focused; }

	std::string get_min_value_string(void) { return ToString(minmax_focusedmeshValues.i, unit); }
	std::string get_max_value_string(void) { return ToString(minmax_focusedmeshValues.j, unit); }

	double get_detail_level(void) { return detail_level; }

	INT3 get_display_renderthresholds(void) { return INT3(renderspeedup1_cells, renderspeedup2_cells, renderspeedup3_cells); }

	//------------------------------------------- indexing

	PhysQRepComponent& operator[](const int& index) { return physQRep[index]; }

	//------------------------------------------- Settings

	void adjust_detail_level(double multiplier) 
	{ 
		if(detail_level*multiplier > 0) detail_level *= multiplier;

		if (detail_level < DETAILVALUEMIN) detail_level = DETAILVALUEMIN;
		if (detail_level > DETAILVALUEMAX) detail_level = DETAILVALUEMAX;
	}

	void set_detail_level(double detail_level_)
	{
		detail_level = detail_level_;

		if (detail_level < DETAILVALUEMIN) detail_level = DETAILVALUEMIN;
		if (detail_level > DETAILVALUEMAX) detail_level = DETAILVALUEMAX;
	}

	void set_display_renderthresholds(INT3 renderthresholds) 
	{ 
		renderspeedup1_cells = renderthresholds.i;
		renderspeedup2_cells = renderthresholds.j;
		renderspeedup3_cells = renderthresholds.k;
	}

	//get current display settings
	PhysQRepSettings get_current_settings(void);

	//------------------------------------------- Calculations

	//Calculate default settings and from physQ calculate physQRep
	void CalculateRepresentation_AutoSettings(std::vector<PhysQ> physQ, D2D1_RECT_F spaceRect);

	//Set new settings, then calculate representation using these new settings
	void CalculateRepresentation_NewSettings(std::vector<PhysQ> physQ, PhysQRepSettings newSettings);

	//Calculate physical quantity representation at current settings, without changing them
	void CalculateRepresentation(std::vector<PhysQ> physQ);

	//------------------------------------------- Drawing

	//calculate transparency (objects closer to camera become more transparent)
	void SetAlphaBlend(void);

	//draw the precomputed physical quantities representations
	void DrawPhysQRep(void);

	//------------------------------------------- Value reading

	//from mouse coordinates obtain mesh coordinates as seen on screen, together with value at the position. If mouse is not on a mesh point return false.
	bool GetMouseInfo(INT2 mouse, std::string* pmeshName, std::string* ptypeName, DBL3* pMeshPosition = nullptr, double* pMeshValue = nullptr, std::string* pValueUnit = nullptr);
};

#endif