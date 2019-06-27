#pragma once

#include "BorisLib.h"

#include "MaterialParameterFormulas.h"

#include "ErrorHandler.h"

#include "Boris_Enums_Defs.h"

#if COMPILECUDA == 1
#include "MaterialParameterCUDA.h"
#endif

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////
//Class holding a single material parameter with any associated temperature dependence and spatial variation.
//This is designed to be used just as the type stored (PType) would be used on its own in mathematical formulas.
//
//PType is the parameter type
//SType is the type used for spatial variation.
//e.g. most parameters just use a scalar for the spatial variation, but there are special cases where we need a vector, e.g. DBL3 for anisotropy symmetry axes, specifying rotations for spatial variation.
//
//This is achieved by updating the parameter output value (current_value) if the temperature changes.
//1. For uniform temperature just call update(Temperature) whenever the uniform base temperature is changed
//2. For non-uniform temperature the value must be updated in every computational cell before being used, using current_value = get(Temperature). 
//This should be done by an external variadic template function if there are multiple parameters to update.

template <typename PType, typename SType>
class MatP :
	public MatPFormula,
	public ProgramState<MatP<PType, SType>, tuple<PType, PType, vector<double>, VEC<SType>, int, string, int, vector<double>>, tuple<>>
{

private:

#if COMPILECUDA == 1
	//when CUDA enabled this will point to the cu_obj managed MatPCUDA object which mirrors this cpu memory parameter in gpu memory.
	//Pointer assigned and nulled in MeshParamsCUDA ctor and dtor. Use it to update the MatPCUDA parameter accordingly when this parameter changes (if not nullptr).
	//To use it, cast it to MatPCUDA<conversion(PType)>. The conversion is done to a cuReal variant, e.g. double/float to cuReal, DBL3/FLT3 to cuReal3 etc.
	void *p_cu_obj_mpcuda = nullptr;
#endif

	//the parameter value at 0K : when constructing it, this is the value to set. This is also the value displayed/changed in the console. Temperaure dependence set separately.
	PType value_at_0K;

	//this is the value that is actually read out and is obtained from value_at_0K depending on the set temperature
	PType current_value;

	//array describing temperature scaling of value_at_0K. The index in this vector is the temperature value in K, i.e. values are set at 1K increments starting at 0K
	vector<double> t_scaling;

	//VEC with spatial dependence scaling of value_at_0K. The VEC rect must match the rect of the mesh to which this MatP applies.
	VEC<SType> s_scaling;

	//value from MATPVAR_ enum
	int s_scaling_type = MATPVAR_NONE;

	//spatial scaling setting info : text with ";" separators. First field is the scaling set type (e.g. none, custom, random, jagged, etc.); The other fields are parameters for spatial generators.
	string s_scaling_info = "none";

	//function object pointing to selected formula from MatPFormula. No coefficients passed, only the temperature, these are set in MatPFormula already. Always return a scaling factor as a double.
	function<double(const MatPFormula&, double)> t_scaling_formula;

private:

	//---------Output value update : CUDA helpers

#if COMPILECUDA == 1
	//update just the value at 0K and current value in the corresponding MatPCUDA
	void update_cuda_value(void);

	//full update of corresponding MatPCUDA
	void update_cuda_object(void);
#endif

	//---------Set spatial dependence : generator methods

	//custom, set from an image file : coefficients vary from 0 : black to 1 : white
	bool set_s_custom(DBL3 h, Rect rect, double offset, double scale, string fileName, function<vector<BYTE>(string, INT2)>& bitmap_loader);

	//random points in given range
	bool set_s_random(DBL3 h, Rect rect, DBL2 range, int seed);

	//random points in given range at given square spacing in xy plane : in between fill using bilinear interpolation
	bool set_s_jagged(DBL3 h, Rect rect, DBL2 range, double spacing, int seed);

	//generate circular defects with a tanh radial profile with values in the given range, diameter range and average spacing (prng instantiated with given seed). The defect positioning is random.
	bool set_s_defects(DBL3 h, Rect rect, DBL2 range, DBL2 diameter_range, double spacing, int seed);

	//generate line faults in the given range length, orientation length(degrees azimuthal) and average spacing (prng instantiated with given seed).
	bool set_s_faults(DBL3 h, Rect rect, DBL2 range, DBL2 length_range, DBL2 orientation_range, double spacing, int seed);

	//generate voronoi 2d tessellation with values in the given range, and average spacing (prng instantiated with given seed). The values range is used to generate coefficient multipliers which are constant in each voronoi cell.
	bool set_s_voronoi2d(DBL3 h, Rect rect, DBL2 range, double spacing, int seed);

	//generate voronoi 3d tessellation with values in the given range, and average spacing (prng instantiated with given seed). The values range is used to generate coefficient multipliers which are constant in each voronoi cell.
	bool set_s_voronoi3d(DBL3 h, Rect rect, DBL2 range, double spacing, int seed);

	//generate voronoi 2d tessellation with average spacing. Set coefficient values randomly in the given range only at the Voronoi cell boundaries (prng instantiated with given seed).
	bool set_s_voronoiboundary2d(DBL3 h, Rect rect, DBL2 range, double spacing, int seed);
	
	//generate voronoi 3d tessellation with average spacing. Set coefficient values randomly in the given range only at the Voronoi cell boundaries (prng instantiated with given seed).
	bool set_s_voronoiboundary3d(DBL3 h, Rect rect, DBL2 range, double spacing, int seed);

	//generate voronoi 2d tessellation with average spacing. This method is applicable only to DBL3 PType, where a rotation operation is applied, fixed in each Voronoi cell. 
	//The rotation uses the values for polar (theta) and azimuthal (phi) angles specified in given ranges. prng instantiated with given seed.
	bool set_s_voronoirotation2d(DBL3 h, Rect rect, DBL2 theta, DBL2 phi, double spacing, int seed);

	//generate voronoi 3d tessellation with average spacing. This method is applicable only to DBL3 PType, where a rotation operation is applied, fixed in each Voronoi cell. 
	//The rotation uses the values for polar (theta) and azimuthal (phi) angles specified in given ranges. prng instantiated with given seed.
	bool set_s_voronoirotation3d(DBL3 h, Rect rect, DBL2 theta, DBL2 phi, double spacing, int seed);

public:

	//---------Constructors

	MatP(void) :
		ProgramStateNames(this, { VINFO(value_at_0K), VINFO(current_value), VINFO(t_scaling), VINFO(s_scaling), VINFO(s_scaling_type), VINFO(s_scaling_info), VINFO(formula_selector), VINFO(coeff) }, {})
	{
		value_at_0K = PType();
		current_value = value_at_0K;
		t_scaling_formula = get_set_function();
	}

	MatP(PType value) :
		value_at_0K(value),
		ProgramStateNames(this, { VINFO(value_at_0K), VINFO(current_value), VINFO(t_scaling), VINFO(s_scaling), VINFO(s_scaling_type), VINFO(s_scaling_info), VINFO(formula_selector), VINFO(coeff) }, {})

	{
		current_value = value_at_0K;	//default to 0K value
		t_scaling_formula = get_set_function();
	}

	//An assignment from a MatP to another must always have matched template parameters. 
	//Conversion when parameters are mismatched should not be used, but we need to define the conversion operator to use the copy_parameters method since the compiler sees such a conversion in the compilation tree as theoretically possible (so throws error)
	template <typename PType_, typename SType_>
	operator MatP<PType_, SType_>&() { return *(new MatP<PType_, SType_>()); }

	//---------Assignment operator - copy a MatP to another of same type

	MatP<PType, SType>& operator=(const MatP<PType, SType>& copy_this);

	//---------Implement ProgramState method

	//must set t_scaling_formula now depending on matp_formula_selector
	void RepairObjectState(void) { t_scaling_formula = get_set_function((MATPFORM_)formula_selector); }

	//---------Output value update

	//update output value for given temperature (use if base temperature changes), including updating of CUDA value if set
	void update(double Temperature = 0);

	//set current value for given temperature, but do not update CUDA value
	void set_current(double Temperature);

	//---------Output value get

	//get value at given temperature and no spatial scaling, but do not update output (meant for use with non-uniform temperature; for uniform temperature it's faster to just read the value through the conversion operator)
	//Use with temperature dependence : YES, spatial variation : NO
	PType get(double Temperature = 0);

	//get value (at base temperature) with spatial scaling (must be set so check before!).
	//Use with temperature dependence : NO, spatial variation : YES
	PType get(const DBL3& position);

	//get value with spatial scaling (must be set so check before!) and temperature dependence
	//Use with temperature dependence : YES, spatial variation : YES
	PType get(const DBL3& position, double Temperature);

	//get 0K value
	PType get0(void) const { return value_at_0K; }

	//get current output value (at base temperature with no spatial scaling)
	PType get_current(void) const { return current_value; }

	//---------Set scaling formula

	//set scaling formula using enum entry
	void set_scaling_formula(MATPFORM_ formula_selector_, vector<double> coefficients);

	//---------Set scaling array (temperature)

	//clear any temperature dependence
	void clear_t_scaling(void);

	//calculate t_scaling from given temperature dependence of scaling coefficients
	bool set_scaling_array(vector<double>& temp_arr, vector<double>& scaling_arr);

	//set scaling array directly, where the value indexes must correspond to temperature values (e.g. t_scaling[0] is for temperature 0 K, etc.).
	void set_precalculated_scaling_array(vector<double>& scaling_arr);

	//---------Set spatial dependence

	//clear spatial scaling
	void clear_s_scaling(void);

	//update spatial scaling VEC for cellsize and rectangle, depending on currently set spatial scaling
	bool update_s_scaling(DBL3 h, Rect rect);

	//set parameter spatial variation using a given generator and arguments (arguments passed as a string to be interpreted and converted using ToNum)
	BError set_s_scaling(DBL3 h, Rect rect, MATPVAR_ generatorID, string generatorArgs, function<vector<BYTE>(string, INT2)>& bitmap_loader);

	//---------Get temperature dependence

	//get temperature scaling coefficients as an array from 0K up to and including max_temperature
	vector<double> get_temperature_scaling(double max_temperature) const;

	//---------Read value

	//use the material parameter just like a PType type to read from it
	operator PType() const { return current_value; }

	//---------Set value

	//set value at 0K, and make current value equal to it
	MatP<PType, SType>& operator=(PType set_value);

	//---------Get Info (intended for console display)

	//returns a string describing the set temperature dependence ("none", "array" or set formula : "name parameters...") 
	string get_info_string(void) const;

	//returns a string describing the set spatial dependence ("none", "array" or set formula : "name parameters...") 
	string get_varinfo_string(void) const { return s_scaling_info; }

	//does it have a temperature dependence?
	bool is_tdep(void) const { return (formula_selector != MATPFORM_NONE || t_scaling.size()); }

	//does it have a spacial dependence?
	bool is_sdep(void) const { return (s_scaling_type != MATPVAR_NONE); }

	//---------Comparison operators

	//Not currently needed

	//---------Arithmetic operators

	//the return type might not simply be PType! (e.g. PType is a DBL3, so '*' for DBL3 is a scalar product which results in a double). 
	//Need to tell the compiler exactly how to determine the return type. Similar problems may occur for the other arithmetic operators

	//multiply with a MatP - '^' is used to signify vector product for VAL3
	template <typename PType_, typename SType_>
	auto operator^(const MatP<PType_, SType_>& rhs) const -> decltype(declval<PType_>() ^ declval<PType_>())
	{
		return current_value ^ rhs.current_value;
	}

	//multiply with a MatP - '&' is used to signify component-by-component product for VAL3
	template <typename PType_, typename SType_>
	auto operator&(const MatP<PType_, SType_>& rhs) const -> decltype(declval<PType_>() & declval<PType_>())
	{
		return current_value & rhs.current_value;
	}

	//multiply with a MatP
	auto operator*(const MatP& rhs) const -> decltype(declval<PType>() * declval<PType>())
	{
		return current_value * rhs.current_value;
	}

	//multiply with a value on the RHS
	template <typename PType_, std::enable_if_t<!std::is_same<PType_, MatP<PType, SType>>::value>* = nullptr>
	auto operator*(const PType_& rhs) const -> decltype(declval<PType>() * declval<PType_>())
	{
		return current_value * rhs;
	}

	//multiply with a value on the LHS
	template <typename _PType, std::enable_if_t<!std::is_same<_PType, MatP<PType, SType>>::value>* = nullptr>
	friend auto operator*(const _PType& lhs, const MatP<PType, SType>& rhs) -> decltype(declval<_PType>() * declval<PType>())
	{
		return lhs * rhs.current_value;
	}

	//divide by a MatP
	auto operator/(const MatP& rhs) const -> decltype(declval<PType>() / declval<PType>())
	{
		return current_value / rhs.current_value;
	}

	//division with a value on the RHS
	template <typename PType_, std::enable_if_t<!std::is_same<PType_, MatP<PType, SType>>::value>* = nullptr>
	auto operator/(const PType_& rhs) const -> decltype(declval<PType>() / declval<PType_>())
	{
		return current_value / rhs;
	}

	//division with a value on the LHS
	template <typename _PType, std::enable_if_t<!std::is_same<_PType, MatP<PType, SType>>::value>* = nullptr>
	friend auto operator/(const _PType& lhs, const MatP<PType, SType>& rhs) -> decltype(declval<_PType>() / declval<PType>())
	{
		return lhs / rhs.current_value;
	}

	//addition with a MatP
	auto operator+(const MatP& rhs) const -> decltype(declval<PType>() + declval<PType>())
	{
		return current_value + rhs.current_value;
	}

	//addition with a value on the RHS
	template <typename PType_, std::enable_if_t<!std::is_same<PType_, MatP<PType, SType>>::value>* = nullptr>
	auto operator+(const PType_& rhs) const -> decltype(declval<PType>() + declval<PType_>())
	{
		return current_value + rhs;
	}

	//addition with a value on the LHS
	template <typename _PType, std::enable_if_t<!std::is_same<_PType, MatP<PType, SType>>::value>* = nullptr>
	friend auto operator+(const _PType& lhs, const MatP<PType, SType>& rhs) -> decltype(declval<_PType>() + declval<PType>())
	{
		return lhs + rhs.current_value;
	}

	//difference with a MatP
	auto operator-(const MatP& rhs) const -> decltype(declval<PType>() - declval<PType>())
	{
		return current_value - rhs.current_value;
	}

	//difference with a value on the RHS
	template <typename PType_, std::enable_if_t<!std::is_same<PType_, MatP<PType, SType>>::value>* = nullptr>
	auto operator-(const PType_& rhs) const -> decltype(declval<PType>() - declval<PType_>())
	{
		return current_value - rhs;
	}

	//difference with a value on the LHS
	template <typename _PType, std::enable_if_t<!std::is_same<_PType, MatP<PType, SType>>::value>* = nullptr>
	friend auto operator-(const _PType& lhs, const MatP<PType, SType>& rhs) -> decltype(declval<_PType>() - declval<PType>())
	{
		return lhs - rhs.current_value;
	}

	//---------Special access methods

	vector<double>& t_scaling_ref(void) { return t_scaling; }

	VEC<SType>& s_scaling_ref(void) { return s_scaling; }

#if COMPILECUDA == 1
	template <typename cu_obj_MatPCUDA_RType>
	void set_p_cu_obj_mpcuda(cu_obj_MatPCUDA_RType* ptr) { p_cu_obj_mpcuda = ptr; }
	void null_p_cu_obj_mpcuda(void) { p_cu_obj_mpcuda = nullptr; }
#endif
};

//---------Assignment operator - copy a MatP to another of same type

template <typename PType, typename SType>
MatP<PType, SType>& MatP<PType, SType>::operator=(const MatP<PType, SType>& copy_this)
{
	//copy value at 0K:
	value_at_0K = copy_this.value_at_0K;

	//copy current value:
	current_value = copy_this.current_value;

	//copy temperature scaling array:
	clear_vector(t_scaling);
	t_scaling = copy_this.t_scaling;

	//copy MatPFormula base
	MatPFormula::operator=(copy_this);

	//now set scaling formula to match
	t_scaling_formula = get_function();

	//copy any spatial scaling
	
	//s_scaling may be empty but copy_this.s_scaling not; in this case the following command doesn't have an effect, but this copy operator should be followed by an update function which will generate the required s_scaling (update_s_scaling called)
	//Note, this copy is still needed here in case s_scaling and copy_this.s_scaling match in dimensions but have different values - update_s_scaling will skip if it sees matching dimensions.
	s_scaling.copy_values(copy_this.s_scaling);

	//s_scaling info so the update will generate the required s_scaling
	s_scaling_type = copy_this.s_scaling_type;
	s_scaling_info = copy_this.s_scaling_info;

	//update cuda object also if set
#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif

	return *this;
}

//---------- TEMPERATURE ONLY

//get value at given temperature and no spatial scaling, but do not update output (meant for use with non-uniform temperature; for uniform temperature it's faster to just read the value through the conversion operator)
//Use with temperature dependence : YES, spatial variation : NO
template <typename PType, typename SType>
PType MatP<PType, SType>::get(double Temperature)
{
	if (formula_selector != MATPFORM_NONE) {

		//use pre-set formula
		return value_at_0K * t_scaling_formula(*this, Temperature);
	}
	else if (t_scaling.size()) {

		//use custom temperature scaling
		int index = (int)floor_epsilon(Temperature);
		if (index + 1 < (int)t_scaling.size() && index >= 0) {

			//use linear interpolation for temperature in range index to index + 1
			return (value_at_0K * (t_scaling[index] * (double(index + 1) - Temperature) + t_scaling[index + 1] * (Temperature - double(index))));
		}
		//for temperatures higher than the loaded array just use the last scaling value
		else return (value_at_0K * t_scaling.back());
	}
	//no temperature dependence set
	else return value_at_0K;
}

//---------- SPATIAL ONLY

//get value (at base temperature) with spatial scaling (must be set so check before!)
//Use with temperature dependence : NO, spatial variation : YES
template <>
inline double MatP<double, double>::get(const DBL3& position)
{
	return current_value * s_scaling[position];
}

template <>
inline DBL2 MatP<DBL2, double>::get(const DBL3& position)
{
	return current_value * s_scaling[position];
}

template <>
inline DBL3 MatP<DBL3, double>::get(const DBL3& position)
{
	return current_value * s_scaling[position];
}

template <>
inline DBL3 MatP<DBL3, DBL3>::get(const DBL3& position)
{
	//rotation not scalar multiplication
	return rotate_polar(current_value, s_scaling[position]);
}

//---------- TEMPERATURE and SPATIAL

//get value with spatial scaling (must be set so check before!) and temperature dependence
//Use with temperature dependence : YES, spatial variation : YES
template <>
inline double MatP<double, double>::get(const DBL3& position, double Temperature)
{
	if (formula_selector != MATPFORM_NONE) {

		//use pre-set formula
		return (value_at_0K * t_scaling_formula(*this, Temperature)) * s_scaling[position];
	}
	else if (t_scaling.size()) {

		//use custom temperature scaling
		int index = (int)floor_epsilon(Temperature);
		if (index + 1 < (int)t_scaling.size() && index >= 0) {

			//use linear interpolation for temperature in range index to index + 1
			return (value_at_0K * (t_scaling[index] * (double(index + 1) - Temperature) + t_scaling[index + 1] * (Temperature - double(index)))) * s_scaling[position];
		}
		//for temperatures higher than the loaded array just use the last scaling value
		else return (value_at_0K * t_scaling.back()) * s_scaling[position];
	}
	//no temperature dependence set
	else return value_at_0K * s_scaling[position];
}

template <>
inline DBL2 MatP<DBL2, double>::get(const DBL3& position, double Temperature)
{
	if (formula_selector != MATPFORM_NONE) {

		//use pre-set formula
		return (value_at_0K * t_scaling_formula(*this, Temperature)) * s_scaling[position];
	}
	else if (t_scaling.size()) {

		//use custom temperature scaling
		int index = (int)floor_epsilon(Temperature);
		if (index + 1 < (int)t_scaling.size() && index >= 0) {

			//use linear interpolation for temperature in range index to index + 1
			return (value_at_0K * (t_scaling[index] * (double(index + 1) - Temperature) + t_scaling[index + 1] * (Temperature - double(index)))) * s_scaling[position];
		}
		//for temperatures higher than the loaded array just use the last scaling value
		else return (value_at_0K * t_scaling.back()) * s_scaling[position];
	}
	//no temperature dependence set
	else return value_at_0K * s_scaling[position];
}

template <>
inline DBL3 MatP<DBL3, double>::get(const DBL3& position, double Temperature)
{
	if (formula_selector != MATPFORM_NONE) {

		//use pre-set formula
		return (value_at_0K * t_scaling_formula(*this, Temperature)) * s_scaling[position];
	}
	else if (t_scaling.size()) {

		//use custom temperature scaling
		int index = (int)floor_epsilon(Temperature);
		if (index + 1 < (int)t_scaling.size() && index >= 0) {

			//use linear interpolation for temperature in range index to index + 1
			return (value_at_0K * (t_scaling[index] * (double(index + 1) - Temperature) + t_scaling[index + 1] * (Temperature - double(index)))) * s_scaling[position];
		}
		//for temperatures higher than the loaded array just use the last scaling value
		else return (value_at_0K * t_scaling.back()) * s_scaling[position];
	}
	//no temperature dependence set
	else return value_at_0K * s_scaling[position];
}

template <>
inline DBL3 MatP<DBL3, DBL3>::get(const DBL3& position, double Temperature)
{
	if (formula_selector != MATPFORM_NONE) {

		//use pre-set formula
		return rotate_polar(value_at_0K * t_scaling_formula(*this, Temperature), s_scaling[position]);
	}
	else if (t_scaling.size()) {

		//use custom temperature scaling
		int index = (int)floor_epsilon(Temperature);
		if (index + 1 < (int)t_scaling.size() && index >= 0) {

			//use linear interpolation for temperature in range index to index + 1
			return rotate_polar(value_at_0K * (t_scaling[index] * (double(index + 1) - Temperature) + t_scaling[index + 1] * (Temperature - double(index))), s_scaling[position]);
		}
		//for temperatures higher than the loaded array just use the last scaling value
		else return rotate_polar(value_at_0K * t_scaling.back(), s_scaling[position]);
	}
	//no temperature dependence set
	else return rotate_polar(value_at_0K, s_scaling[position]);
}

//---------Output value update

template <typename PType, typename SType>
void MatP<PType, SType>::update(double Temperature)
{
	current_value = get(Temperature);

#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_value();
#endif
}

//set current value for given temperature, but do not update CUDA value
template <typename PType, typename SType>
void MatP<PType, SType>::set_current(double Temperature)
{
	if (formula_selector != MATPFORM_NONE) {

		//use pre-set formula
		current_value = value_at_0K * t_scaling_formula(*this, Temperature);
	}
	else if (t_scaling.size()) {

		//use custom temperature scaling
		int index = (int)floor_epsilon(Temperature);
		if (index + 1 < (int)t_scaling.size() && index >= 0) {

			//use linear interpolation for temperature in range index to index + 1
			current_value = (value_at_0K * (t_scaling[index] * (double(index + 1) - Temperature) + t_scaling[index + 1] * (Temperature - double(index))));
		}
		//for temperatures higher than the loaded array just use the last scaling value
		else current_value = (value_at_0K * t_scaling.back());
	}
	//no temperature dependence set
	else current_value = value_at_0K;
}

//---------Set value

template <typename PType, typename SType>
MatP<PType, SType>& MatP<PType, SType>::operator=(PType set_value)
{
	value_at_0K = set_value;
	current_value = value_at_0K;	//default to 0K value when using this

#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_value();
#endif

	return *this;
}

//---------Get Info

template <typename PType, typename SType>
string MatP<PType, SType>::get_info_string(void) const
{
	if (formula_selector == MATPFORM_NONE) {

		if (t_scaling.size()) return "array";
		else return "none";
	}
	else return get_formula_info();
}

//---------Set scaling formula

template <typename PType, typename SType>
void MatP<PType, SType>::set_scaling_formula(MATPFORM_ formula_selector_, vector<double> coefficients)
{
	t_scaling_formula = get_set_function(formula_selector_);

	set_formula_coeffs(coefficients);

	//if setting the none formula (no analytical temperature dependence) then also clear any array temperature dependence : i.e. set no temperature dependence
	if (formula_selector == MATPFORM_NONE) t_scaling.clear();

	update();

#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif
}

//---------Set scaling array

//clear any temperature dependence
template <typename PType, typename SType>
void MatP<PType, SType>::clear_t_scaling(void)
{
	//reset formula
	t_scaling_formula = get_set_function(MATPFORM_NONE);

	//clear scaling array
	t_scaling.clear();

	//update output value to value at 0K
	update();

#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif
}

//calculate t_scaling from given temperature dependence of scaling coefficients
template <typename PType, typename SType>
bool MatP<PType, SType>::set_scaling_array(vector<double>& temp_arr, vector<double>& scaling_arr)
{
	//vectors must have same size and must have at least 2 points
	if (temp_arr.size() < 2 || temp_arr.size() != scaling_arr.size()) return false;

	//make sure the temperature values are in increasing order
	quicksort(temp_arr, scaling_arr);

	//temperature values cannot be negative
	if (temp_arr[0] < 0) return false;

	//make sure temperature values are not repeated - values must be strictly increasing
	delete_repeats(temp_arr, scaling_arr);

	//set tscaling size : covers temperature from 0K up to maximum possible temperature in 1K increments.
	double max_temp = minimum(temp_arr.back(), MAX_TEMPERATURE);
	t_scaling.resize((size_t)(floor_epsilon(max_temp) + 1));

	int temp_arr_idx = 0;

	//first temperature point : 0K
	double temp_a = 0;
	//next temperature point - will use interpolation to fill in any integer temperature values in between; this may coincide with temp_a
	double temp_b = temp_arr[0];

	//first scaling coefficient (at 0K) : extrapolate from first 2 points in the arrays
	double scaling_a = interpolate(DBL2(temp_arr[0], scaling_arr[0]), DBL2(temp_arr[1], scaling_arr[1]), 0);
	//next scaling coefficient
	double scaling_b = scaling_arr[0];

	//we know what the first value must be:
	t_scaling[temp_a] = scaling_a;

	for (int t_value = 1; t_value < (int)t_scaling.size(); t_value++) {

		//have we finished this [Ta, Tb) interval? Use a while since the temperature points in temp_arr may be finely spaced (less than 1K intervals)
		while ((double)t_value >= temp_b) {
		
			temp_arr_idx++;

			//have we reached the end? i.e. was temp_b the last temperature value in temp_arr? This can only happen if temp_b is an integer value, in which case t_value = temp_b now.
			if (temp_arr_idx >= (int)temp_arr.size()) {

				t_scaling[t_value] = scaling_arr.back();
				break;
			}

			//next [Ta, Tb) interval
			temp_a = temp_b;
			temp_b = temp_arr[temp_arr_idx];

			scaling_a = scaling_b;
			scaling_b = scaling_arr[temp_arr_idx];
		}

		//t_value is now in the interval [temp_a, temp_b). Set t_scaling value at t_value index using interpolation
		t_scaling[t_value] = interpolate(DBL2(temp_a, scaling_a), DBL2(temp_b, scaling_b), t_value);
	}

	//no scaling formula - using array instead
	t_scaling_formula = get_set_function(MATPFORM_NONE);
	update();

#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif
	
	return true;
}

template <typename PType, typename SType>
void MatP<PType, SType>::set_precalculated_scaling_array(vector<double>& scaling_arr)
{
	t_scaling = scaling_arr;

	//no scaling formula - using array instead
	t_scaling_formula = get_set_function(MATPFORM_NONE);
	update();

#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif
}

//---------Set spatial dependence

//clear spatial scaling
template <typename PType, typename SType>
void MatP<PType, SType>::clear_s_scaling(void)
{
	s_scaling.clear();
	s_scaling_type = MATPVAR_NONE;
	s_scaling_info = "none";

	//update cuda object also if set
#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif
}

//update spatial scaling VEC for cellsize and rectangle, depending on currently set spatial scaling
template <typename PType, typename SType>
bool MatP<PType, SType>::update_s_scaling(DBL3 h, Rect rect)
{
	//nothing to do
	if (s_scaling_type == MATPVAR_NONE) return true;

	//there is a spatial scaling set
	if (s_scaling_type == MATPVAR_CUSTOM) {

		//for custom spatial scaling just stretch to new sizes, but keep it 2D
		INT3 cells = round(rect / h);
		if (s_scaling.linear_size()) s_scaling.resize(rect / SZ3(cells.x, cells.y, 1), rect);
		return true;
	}
	else {

		//if cellsize and rect match then nothing to do
		if ((s_scaling.h == h) && (s_scaling.rect == rect)) return true;

		//fields in s_scaling_info first contains the generator name, then its args - all separated using "; " or ", "
		vector<string> fields = split(s_scaling_info, { ", ", "; " });

		//get rid of first field, then recombine into a string using " " as separators
		BError error = set_s_scaling(h, rect, (MATPVAR_)s_scaling_type, combine(subvec(fields, 1), " "), function<vector<BYTE>(string, INT2)>());

		if (!error) return true;
	}

	return false;
}

//set parameter spatial variation using a given generator and arguments (arguments passed as a string to be interpreted and converted using ToNum)
template <typename PType, typename SType>
BError MatP<PType, SType>::set_s_scaling(DBL3 h, Rect rect, MATPVAR_ generatorID, string generatorArgs, function<vector<BYTE>(string, INT2)>& bitmap_loader)
{
	BError error(__FUNCTION__);

	switch (generatorID) {

	case MATPVAR_CUSTOM:
	{
		vector<string> fields = split(generatorArgs, " ");

		if (fields.size() < 3) return error(BERROR_INCORRECTVALUE);

		double offset = ToNum(fields[0]);
		double scale = ToNum(fields[1]);
		string fileName = combine(subvec(fields, 2), " ");

		if (!set_s_custom(h, rect, offset, scale, fileName, bitmap_loader)) error(BERROR_OUTOFMEMORY_NCRIT);
	}
	break;

	case MATPVAR_RANDOM:
	{
		vector<string> fields = split(generatorArgs, " ");

		if (fields.size() != 3) return error(BERROR_INCORRECTVALUE);

		DBL2 range = ToNum(fields[0] + " " + fields[1]);
		int seed = ToNum(fields[2]);

		if (!set_s_random(h, rect, range, seed)) error(BERROR_OUTOFMEMORY_NCRIT);
	}
	break;

	case MATPVAR_JAGGED:
	{
		vector<string> fields = split(generatorArgs, " ");

		if (fields.size() != 4) return error(BERROR_INCORRECTVALUE);

		DBL2 range = ToNum(fields[0] + " " + fields[1]);
		double spacing = ToNum(fields[2], "m");
		int seed = ToNum(fields[3]);

		if (!set_s_jagged(h, rect, range, spacing, seed)) error(BERROR_OUTOFMEMORY_NCRIT);
	}
	break;

	case MATPVAR_DEFECTS:
	{
		vector<string> fields = split(generatorArgs, " ");

		if (fields.size() != 6) return error(BERROR_INCORRECTVALUE);

		DBL2 range = ToNum(fields[0] + " " + fields[1]);
		DBL2 diameter_range = ToNum(fields[2] + " " + fields[3], "m");
		double spacing = ToNum(fields[4], "m");
		int seed = ToNum(fields[5]);

		if (!set_s_defects(h, rect, range, diameter_range, spacing, seed)) error(BERROR_OUTOFMEMORY_NCRIT);
	}
	break;

	case MATPVAR_FAULTS:
	{
		vector<string> fields = split(generatorArgs, " ");

		if (fields.size() != 8) return error(BERROR_INCORRECTVALUE);

		DBL2 range = ToNum(fields[0] + " " + fields[1]);
		DBL2 length_range = ToNum(fields[2] + " " + fields[3], "m");
		DBL2 orientation_range = ToNum(fields[4] + " " + fields[5]);
		double spacing = ToNum(fields[6], "m");
		int seed = ToNum(fields[7]);

		if (!set_s_faults(h, rect, range, length_range, orientation_range, spacing, seed)) error(BERROR_OUTOFMEMORY_NCRIT);
	}
	break;

	case MATPVAR_VORONOI2D:
	{
		vector<string> fields = split(generatorArgs, " ");

		if (fields.size() != 4) return error(BERROR_INCORRECTVALUE);

		DBL2 range = ToNum(fields[0] + " " + fields[1]);
		double spacing = ToNum(fields[2], "m");
		int seed = ToNum(fields[3]);

		if (!set_s_voronoi2d(h, rect, range, spacing, seed)) error(BERROR_OUTOFMEMORY_NCRIT);
	}
	break;

	case MATPVAR_VORONOI3D:
	{
		vector<string> fields = split(generatorArgs, " ");

		if (fields.size() != 4) return error(BERROR_INCORRECTVALUE);

		DBL2 range = ToNum(fields[0] + " " + fields[1]);
		double spacing = ToNum(fields[2], "m");
		int seed = ToNum(fields[3]);

		if (!set_s_voronoi3d(h, rect, range, spacing, seed)) error(BERROR_OUTOFMEMORY_NCRIT);
	}
	break;

	case MATPVAR_VORONOIBND2D:
	{
		vector<string> fields = split(generatorArgs, " ");

		if (fields.size() != 4) return error(BERROR_INCORRECTVALUE);

		DBL2 range = ToNum(fields[0] + " " + fields[1]);
		double spacing = ToNum(fields[2], "m");
		int seed = ToNum(fields[3]);

		if (!set_s_voronoiboundary2d(h, rect, range, spacing, seed)) error(BERROR_OUTOFMEMORY_NCRIT);
	}
	break;

	case MATPVAR_VORONOIBND3D:
	{
		vector<string> fields = split(generatorArgs, " ");

		if (fields.size() != 4) return error(BERROR_INCORRECTVALUE);

		DBL2 range = ToNum(fields[0] + " " + fields[1]);
		double spacing = ToNum(fields[2], "m");
		int seed = ToNum(fields[3]);

		if (!set_s_voronoiboundary3d(h, rect, range, spacing, seed)) error(BERROR_OUTOFMEMORY_NCRIT);
	}
	break;

	case MATPVAR_VORONOIROT2D:
	{
		vector<string> fields = split(generatorArgs, " ");

		if (fields.size() != 6) return error(BERROR_INCORRECTVALUE);

		DBL2 range_theta = ToNum(fields[0] + " " + fields[1]);
		DBL2 range_phi = ToNum(fields[2] + " " + fields[3]);
		double spacing = ToNum(fields[4], "m");
		int seed = ToNum(fields[5]);

		if (!set_s_voronoirotation2d(h, rect, range_theta, range_phi, spacing, seed)) error(BERROR_OUTOFMEMORY_NCRIT);
	}
	break;

	case MATPVAR_VORONOIROT3D:
	{
		vector<string> fields = split(generatorArgs, " ");

		if (fields.size() != 6) return error(BERROR_INCORRECTVALUE);

		DBL2 range_theta = ToNum(fields[0] + " " + fields[1]);
		DBL2 range_phi = ToNum(fields[2] + " " + fields[3]);
		double spacing = ToNum(fields[4], "m");
		int seed = ToNum(fields[5]);

		if (!set_s_voronoirotation3d(h, rect, range_theta, range_phi, spacing, seed)) error(BERROR_OUTOFMEMORY_NCRIT);
	}
	break;
	}

	return error;
}

//---------Set spatial dependence : generator methods

//custom, set from an image file : coefficients vary from 0 : black to 1 : white
template <typename PType, typename SType>
bool MatP<PType, SType>::set_s_custom(DBL3 h, Rect rect, double offset, double scale, string fileName, function<vector<BYTE>(string, INT2)>& bitmap_loader)
{
	if (GetFileTermination(fileName) != ".png")
		fileName += ".png";
	
	INT3 cells = round(rect / h);
	bool success = s_scaling.generate_custom_2D(SZ3(cells.x, cells.y, 1), rect, offset, scale, bitmap_loader(fileName, INT2(cells.x, cells.y)));

	//don't display path or termination in paramvar desription
	ExtractFilenameDirectory(fileName);
	ExtractFilenameTermination(fileName);

	s_scaling_type = MATPVAR_CUSTOM;
	s_scaling_info = "custom; " + ToString(offset) + "; " + ToString(scale) + "; " + fileName;

	//update cuda object also if set
#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif

	//allow setting a wrong filename without causing error messages
	return true;
}

//random points in given range
template <typename PType, typename SType>
bool MatP<PType, SType>::set_s_random(DBL3 h, Rect rect, DBL2 range, int seed)
{
	bool success = s_scaling.generate_random(h, rect, range, seed);

	if (success) {

		s_scaling_type = MATPVAR_RANDOM;
		s_scaling_info = "random; " + ToString(range) + "; " + ToString(seed);
	}
	else clear_s_scaling();

	//update cuda object also if set
#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif

	return success;
}

//random points in given range at given square spacing in xy plane : in between fill using bilinear interpolation
template <typename PType, typename SType>
bool MatP<PType, SType>::set_s_jagged(DBL3 h, Rect rect, DBL2 range, double spacing, int seed)
{
	bool success = s_scaling.generate_jagged(h, rect, range, spacing, seed);

	if (success) {

		s_scaling_type = MATPVAR_JAGGED;
		s_scaling_info = "jagged; " + ToString(range) + "; " + ToString(spacing, "m") + "; " + ToString(seed);
	}
	else clear_s_scaling();

	//update cuda object also if set
#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif

	return success;
}

//generate circular defects with a tanh radial profile with values in the given range, diameter range and average spacing(prng instantiated with given seed).The defect positioning is random.
template <typename PType, typename SType>
bool MatP<PType, SType>::set_s_defects(DBL3 h, Rect rect, DBL2 range, DBL2 diameter_range, double spacing, int seed)
{
	bool success = s_scaling.generate_defects(h, rect, range, 1.0, diameter_range, spacing, seed);

	if (success) {

		s_scaling_type = MATPVAR_DEFECTS;
		s_scaling_info = "defects; " + ToString(range) + "; " + ToString(diameter_range, "m") + "; " + ToString(spacing, "m") + "; " + ToString(seed);
	}
	else clear_s_scaling();

	//update cuda object also if set
#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif

	return success;
}

//generate circular defects with a tanh radial profile with values in the given range, diameter range and average spacing(prng instantiated with given seed).The defect positioning is random.
template <typename PType, typename SType>
bool MatP<PType, SType>::set_s_faults(DBL3 h, Rect rect, DBL2 range, DBL2 length_range, DBL2 orientation_range, double spacing, int seed)
{
	bool success = s_scaling.generate_faults(h, rect, range, 1.0, length_range, orientation_range, spacing, seed);

	if (success) {

		s_scaling_type = MATPVAR_FAULTS;
		s_scaling_info = "faults; " + ToString(range) + "; " + ToString(length_range, "m") + "; " + ToString(orientation_range) + "; " + ToString(spacing, "m") + "; " + ToString(seed);
	}
	else clear_s_scaling();

	//update cuda object also if set
#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif

	return success;
}

//generate voronoi 2d tessellation with values in the given range, and average spacing (prng instantiated with given seed). The values range is used to generate coefficient multipliers which are constant in each voronoi cell.
template <typename PType, typename SType>
bool MatP<PType, SType>::set_s_voronoi2d(DBL3 h, Rect rect, DBL2 range, double spacing, int seed)
{
	bool success = s_scaling.generate_voronoi2d(h, rect, range, spacing, seed);

	if (success) {

		s_scaling_type = MATPVAR_VORONOI2D;
		s_scaling_info = "vor2D; " + ToString(range) + "; " + ToString(spacing, "m") + "; " + ToString(seed);
	}
	else clear_s_scaling();

	//update cuda object also if set
#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif

	return success;
}

//generate voronoi 3d tessellation with values in the given range, and average spacing (prng instantiated with given seed). The values range is used to generate coefficient multipliers which are constant in each voronoi cell.
template <typename PType, typename SType>
bool MatP<PType, SType>::set_s_voronoi3d(DBL3 h, Rect rect, DBL2 range, double spacing, int seed)
{
	bool success = s_scaling.generate_voronoi3d(h, rect, range, spacing, seed);

	if (success) {

		s_scaling_type = MATPVAR_VORONOI3D;
		s_scaling_info = "vor3D; " + ToString(range) + "; " + ToString(spacing, "m") + "; " + ToString(seed);
	}
	else clear_s_scaling();

	//update cuda object also if set
#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif

	return success;
}

//generate voronoi 2d tessellation with average spacing. Set coefficient values randomly in the given range only at the Voronoi cell boundaries (prng instantiated with given seed).
template <typename PType, typename SType>
bool MatP<PType, SType>::set_s_voronoiboundary2d(DBL3 h, Rect rect, DBL2 range, double spacing, int seed)
{
	bool success = s_scaling.generate_voronoiboundary2d(h, rect, range, SType(1), spacing, seed);

	if (success) {

		s_scaling_type = MATPVAR_VORONOIBND2D;
		s_scaling_info = "vorbnd2D; " + ToString(range) + "; " + ToString(spacing, "m") + "; " + ToString(seed);
	}
	else clear_s_scaling();

	//update cuda object also if set
#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif

	return success;
}

//generate voronoi 3d tessellation with average spacing. Set coefficient values randomly in the given range only at the Voronoi cell boundaries (prng instantiated with given seed).
template <typename PType, typename SType>
bool MatP<PType, SType>::set_s_voronoiboundary3d(DBL3 h, Rect rect, DBL2 range, double spacing, int seed)
{
	bool success = s_scaling.generate_voronoiboundary3d(h, rect, range, SType(1), spacing, seed);

	if (success) {

		s_scaling_type = MATPVAR_VORONOIBND3D;
		s_scaling_info = "vorbnd3D; " + ToString(range) + "; " + ToString(spacing, "m") + "; " + ToString(seed);
	}
	else clear_s_scaling();

	//update cuda object also if set
#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif

	return success;
}

//generate voronoi 2d tessellation with average spacing. This method is applicable only to DBL3 PType, where a rotation operation is applied, fixed in each Voronoi cell. 
//The rotation uses the values for polar (theta) and azimuthal (phi) angles specified in given ranges. prng instantiated with given seed.
template <typename PType, typename SType>
bool MatP<PType, SType>::set_s_voronoirotation2d(DBL3 h, Rect rect, DBL2 theta, DBL2 phi, double spacing, int seed)
{
	bool success = s_scaling.generate_voronoirotation2d(h, rect, theta, phi, spacing, seed);

	if (success) {

		s_scaling_type = MATPVAR_VORONOIROT2D;
		s_scaling_info = "vorrot2D; " + ToString(theta) + "; " + ToString(phi) + "; " + ToString(spacing, "m") + "; " + ToString(seed);
	}
	else clear_s_scaling();

	//update cuda object also if set
#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif

	return success;
}

//generate voronoi 3d tessellation with average spacing. This method is applicable only to DBL3 PType, where a rotation operation is applied, fixed in each Voronoi cell. 
//The rotation uses the values for polar (theta) and azimuthal (phi) angles specified in given ranges. prng instantiated with given seed.
template <typename PType, typename SType>
bool MatP<PType, SType>::set_s_voronoirotation3d(DBL3 h, Rect rect, DBL2 theta, DBL2 phi, double spacing, int seed)
{
	bool success = s_scaling.generate_voronoirotation3d(h, rect, theta, phi, spacing, seed);

	if (success) {

		s_scaling_type = MATPVAR_VORONOIROT3D;
		s_scaling_info = "vorrot3D; " + ToString(theta) + "; " + ToString(phi) + "; " + ToString(spacing, "m") + "; " + ToString(seed);
	}
	else clear_s_scaling();

	//update cuda object also if set
#if COMPILECUDA == 1
	if (p_cu_obj_mpcuda) update_cuda_object();
#endif

	return success;
}

//---------Get temperature dependence

//get temperature scaling coefficients as an array from 0K up to and including max_temperature
template <typename PType, typename SType>
vector<double> MatP<PType, SType>::get_temperature_scaling(double max_temperature) const
{
	vector<double> get_scaling;

	if (max_temperature < 0 || max_temperature > MAX_TEMPERATURE) return get_scaling;

	get_scaling.resize((int)floor_epsilon(max_temperature) + 1);

	for (int t_value = 0; t_value <= max_temperature; t_value++) {

		if (formula_selector != MATPFORM_NONE) {

			//use pre-set formula
			get_scaling[t_value] = t_scaling_formula(*this, (double)t_value);
		}
		else if (t_value < t_scaling.size()) {

			get_scaling[t_value] = t_scaling[t_value];
		}
		else if (t_scaling.size()) {

			get_scaling[t_value] = t_scaling.back();
		}
		else get_scaling[t_value] = 1.0;
	}

	return get_scaling;
}