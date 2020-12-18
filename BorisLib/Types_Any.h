#pragma once

// Store any simple type in an Any, out of the possible cases treated in MatchType_CallFunc method. This is not templated, so an Any can be stored in a std::vector and change type stored at run-time.
//
// Usage:
//
// Any a(5);			//store an integer
// Any b = 5;			//store an integer
// b = 6.6;				//double is converted to integer (will hold value 6)
// int c = b;			//c will hold value 6
// b = std::string("test");	//conversion not possible, b still stores a 6
//
// Any can be used with stream << operators and std::stringstream >> operator. e.g.:
// std::cout << b << std::endl;	//output 6
//
// Any d(5.5);			//d stores a double
// d = b;				//d now stores an integer
//
// Any e(&c);			//e stores a reference to the integer c (only use in scopes smaller or same as the variable pointed to). 
//						//Typical use case: function with parameter pack variables passed through reference. Make:
						//std::vector<Any> vec_any = { Any(&variables)... }; //Use vec_any as needed: changes to stored values will propagate to calling variables. 
// e = 7;				//both c and e store value 7
//
// Any f(0.0);
// f.convert_string("8.1987km", "m");	//stores value 8198.7. Similar to >> std::stringstream operator for conversion, but allows use of units.
//
// f.convert_to_string("m");	//similar to << std::stringstream operator for conversion, but allows use of units. Will output 8.1987km as a std::string.
//
// Any g;
// g.convert_string_set_type("8000", btypeinfo<double>().name());	//stores 8000 as a double - use this when type of Any has not been set yet, or to overwrite set type.
// g.clear();			//clear stored type and value
// g.IsNull();			//anything stored?
// g.get_type();		//get btypeinfo name

#include <string>
#include <vector>

#include "Types_Conversion.h"
#include "Types_VAL.h"
#include "Types_ReIm.h"
#include "Types_Rect.h"
#include "Types_Sequences.h"
#include "Types_Info.h"

//------------------------------- Parameters for function calls (also act as tags for tag dispatching)

//NEWVALUE
struct NewValue_Params {

	void *pValue;

	NewValue_Params(void *pValue_) { pValue = pValue_; }
};

//SETVALUE
template <typename Type>
struct SetValue_Params {

	Type value;

	SetValue_Params(Type value_) { value = value_; }
};

template <>
struct SetValue_Params<void> {

};

//FROMSTRING
struct ConvertFromString_Params {

	std::string text;
	std::string unit;

	ConvertFromString_Params(const std::string& text_, const std::string& unit_)
		: text(text_), unit(unit_)
	{}
};

//FROMSTRING_COMPONENTS
struct ConvertFromStringComponents_Params {

	std::vector<std::string> components;
	std::string unit;

	ConvertFromStringComponents_Params(const std::vector<std::string>& components_, const std::string& unit_)
		: components(components_), unit(unit_)
	{}
};

//TOSTRING
struct ConvertToString_Params {

	std::string unit;

	ConvertToString_Params(std::string unit_) { unit = unit_; }
};

//CONVERTTYPE
struct ConvertType_Params {

};

//DELETETYPE
struct DeleteType_Params {

};

class Any {

private:

	//the actual value stored in pValue
	void* pValue = nullptr;

	//has this object allocated pValue, or is it just a pointer to somewhere else? 
	//Whenever new operator is used here set allocated = true, and only delete here if allocated == true.
	//This allows for Any to become a reference to a pointer (constructed through the pointer parameter constructor) - only use this in special circumstances in a restricted scope where it can be guaranteed to original pointer will not delete before this Any goes out of scope.
	bool allocated = false;

	//the type of the value stored in pValue, encoded as a std::string obtained from typeid(...).name;
	std::string type_name;

private:

	//------------------------------- Overloaded functions selected through tag dispatching (where the tag also acts as useful input parameter). RType is the return type (must be selected correctly).

	//NEWVALUE
	template <typename RType, typename Type> RType RunThisMethod(NewValue_Params param)
	{
		pValue = new Type(*reinterpret_cast<Type*>(param.pValue));
		allocated = true;
	}

	//SETVALUE
	template <typename RType, typename Type> RType Convert_RType_to_Type(SetValue_Params<RType> param, std::true_type)
	{
		//can convert RType to Type : do it
		RType value = *reinterpret_cast<RType*>(&param.value);

		*reinterpret_cast<Type*>(pValue) = static_cast<Type>(value);
		return value;
	}

	template <typename RType, typename Type> RType Convert_RType_to_Type(...) { return RType(); }	//cannot convert : do nothing

	template <typename RType, typename Type> RType RunThisMethod(SetValue_Params<RType> param)
	{
		return Convert_RType_to_Type<RType, Type>(param, std::is_convertible<RType, Type>{});
	}

	//FROMSTRING
	template <typename RType, typename Type> RType RunThisMethod(ConvertFromString_Params param)
	{
		//set value from std::string to Type as obtained from type_name
		Type value = ToNum(param.text, param.unit);

		if (!pValue) {
			pValue = new Type(value);
			allocated = true;
		}
		else {

			*reinterpret_cast<Type*>(pValue) = value;
		}
	}

	//FROMSTRING_COMPONENTS
	template <typename RType, typename Type> RType RunThisMethod(ConvertFromStringComponents_Params param)
	{
		//find list of possible number of parameters accepted in Type constructors (all convertible from std::string using ToNum)
		std::vector<int> ctor_parameters;
		constructor_parameters<Type>::get_constructors(ctor_parameters);

		//for strings restrict to 1 parameter only (constructor_parameters<std::string>::get_constructors() also returns one with 3 parameters)
		if (std::is_same<Type, std::string>::value) ctor_parameters.assign(1, 1);

		//number of std::string components available to use
		int available_parameters = (int)param.components.size();

		//start from largest number of parameters and convert as soon as a constructor is found which takes the number of available parameters - return number of parameters used
		for (auto it = ctor_parameters.rbegin(); it != ctor_parameters.rend(); ++it) {

			if (available_parameters >= *it) {

				//found one - convert
				std::string convert_text;
				for (int idx = 0; idx < *it; ) {

					convert_text += param.components[idx];
					if (++idx != *it) convert_text += " ";
				}

				//convert and return number of std::string components used (ToNum uses std::stringstream >> operator, which mirrors the available constructors in each case) - see Types.h
				Type value = ToNum(convert_text, param.unit);
				*reinterpret_cast<Type*>(pValue) = value;

				return *it;
			}
		}

		//no suitable constructor found
		return 0;
	}

	//TOSTRING
	template <typename RType, typename Type> RType RunThisMethod(ConvertToString_Params param)
	{
		return ToString(*reinterpret_cast<Type*>(pValue), param.unit);
	}

	//CONVERTTYPE
	template <typename RType, typename Type> RType Convert_Type_to_RType(std::true_type)
	{
		Type value = *reinterpret_cast<Type*>(pValue);
		return (RType)value;
	}

	template <typename RType, typename Type> RType Convert_Type_to_RType(std::false_type) { return RType(); }

	template <typename RType, typename Type> RType RunThisMethod(ConvertType_Params param)
	{
		return Convert_Type_to_RType<RType, Type>(std::is_convertible<Type, RType>());
	}

	//DELETETYPE
	template <typename RType, typename Type> RType RunThisMethod(DeleteType_Params param)
	{
		if (pValue && allocated) {

			delete reinterpret_cast<Type*>(pValue);
			pValue = nullptr;
			allocated = false;
		}

		return RType();
	}

	//Type selector and dispatcher based on the PType tag (which also serves as useful input parameter in most cases so no need to optimize by static_cast<void> it away).
	//The type is decided by type_name. Add all variable types which can be handled by Any objects here.
	template <typename RType, typename PType> RType MatchType_CallFunc(PType param) {

		if (type_name == btype_info<bool>().name()) {

			return RunThisMethod<RType, bool>(param);
		}

		else if (type_name == btype_info<char>().name()) {

			return RunThisMethod<RType, char>(param);
		}

		else if (type_name == btype_info<int>().name()) {

			return RunThisMethod<RType, int>(param);
		}

		else if (type_name == btype_info<float>().name()) {

			return RunThisMethod<RType, float>(param);
		}

		else if (type_name == btype_info<double>().name()) {

			return RunThisMethod<RType, double>(param);
		}

		else if (type_name == btype_info<std::string>().name()) {

			return RunThisMethod<RType, std::string>(param);
		}

		else if (type_name == btype_info<INT2>().name()) {

			return RunThisMethod<RType, INT2>(param);
		}

		else if (type_name == btype_info<FLT2>().name()) {

			return RunThisMethod<RType, FLT2>(param);
		}

		else if (type_name == btype_info<DBL2>().name()) {

			return RunThisMethod<RType, DBL2>(param);
		}

		else if (type_name == btype_info<INT3>().name()) {

			return RunThisMethod<RType, INT3>(param);
		}

		else if (type_name == btype_info<FLT3>().name()) {

			return RunThisMethod<RType, FLT3>(param);
		}

		else if (type_name == btype_info<DBL3>().name()) {

			return RunThisMethod<RType, DBL3>(param);
		}

		else if (type_name == btype_info<INT4>().name()) {

			return RunThisMethod<RType, INT4>(param);
		}

		else if (type_name == btype_info<FLT4>().name()) {

			return RunThisMethod<RType, FLT4>(param);
		}

		else if (type_name == btype_info<DBL4>().name()) {

			return RunThisMethod<RType, DBL4>(param);
		}

		else if (type_name == btype_info<Box>().name()) {

			return RunThisMethod<RType, Box>(param);
		}

		else if (type_name == btype_info<Rect>().name()) {

			return RunThisMethod<RType, Rect>(param);
		}

		else if (type_name == btype_info<SEQ>().name()) {

			return RunThisMethod<RType, SEQ>(param);
		}

		else if (type_name == btype_info<SEQ3>().name()) {

			return RunThisMethod<RType, SEQ3>(param);
		}

		else if (type_name == btype_info<SEQP>().name()) {

			return RunThisMethod<RType, SEQP>(param);
		}

		else if (type_name == btype_info<COSSEQ>().name()) {

			return RunThisMethod<RType, COSSEQ>(param);
		}

		else if (type_name == btype_info<COSSEQ3>().name()) {

			return RunThisMethod<RType, COSSEQ3>(param);
		}

		else if (type_name == btype_info<SINOSC>().name()) {

			return RunThisMethod<RType, SINOSC>(param);
		}

		else if (type_name == btype_info<COSOSC>().name()) {

			return RunThisMethod<RType, COSOSC>(param);
		}

		else if (type_name == btype_info<SINOSC3>().name()) {

			return RunThisMethod<RType, SINOSC3>(param);
		}

		else if (type_name == btype_info<COSOSC3>().name()) {

			return RunThisMethod<RType, COSOSC3>(param);
		}

		else if (type_name == btype_info<StringSequence>().name()) {

			return RunThisMethod<RType, StringSequence>(param);
		}

		else if (type_name == btype_info<FILESEQ>().name()) {

			return RunThisMethod<RType, FILESEQ>(param);
		}

		else if (type_name == btype_info<FILESEQ3>().name()) {

			return RunThisMethod<RType, FILESEQ3>(param);
		}

		return RType();
	}

	//------------------------------- AUX FUNCTIONS

	//copy the value and type in copyThis by allocating new memory - results in independent object
	void CopyAny(const Any& copyThis)
	{
		if (pValue && allocated)
			MatchType_CallFunc<void, DeleteType_Params>(DeleteType_Params());

		//copy type name
		type_name = copyThis.type_name;

		//now match the type name with the actual type and allocate memory and set value
		MatchType_CallFunc<void, NewValue_Params>(NewValue_Params(copyThis.pValue));
	}

public:

	//------------------------------- CONSTRUCTORS

	//empty constructor : null object
	Any(void)
	{
		pValue = nullptr;
		allocated = false;
		type_name = "";
	}

	//value constructor - makes an independent Any object
	template <typename Type>
	Any(Type value)
	{
		pValue = new Type(value);
		allocated = true;

		type_name = btype_info<Type>().name();
	}

	//pointer constructor - makes a pointer Any object, i.e. pValue points to passed value
	template <typename Type> Any(Type *value)
	{
		pValue = value;
		allocated = false;

		type_name = btype_info<Type>().name();
	}

	//copy constructor - make an exact copy. This means: independent if object to copy from is independent (i.e. has allocated == true), else keep it as a dependent object - this means the value set in pValue is a pointer to some value allocated elsewhere.
	Any(const Any& copyThis)
	{
		if (copyThis.allocated) {

			//make new Any object, independent of copyThis
			CopyAny(copyThis);
		}
		else {

			//copyThis points to a value allocated elsewhere. Keep this object the same.
			pValue = copyThis.pValue;
			type_name = copyThis.type_name;
		}
	}

	~Any()
	{
		if (pValue && allocated)
			MatchType_CallFunc<void, DeleteType_Params>(DeleteType_Params());
	}

	//------------------------------- READ VALUE FROM THIS

	//conversion to external type, e.g. double value = a; //here a is of type Any. The value stored in a is converted to the type of value (double in this case). If conversion not possible then value is set to its default, here RType()
	template <typename RType>
	operator RType()
	{
		if (pValue)
			return MatchType_CallFunc<RType, ConvertType_Params>(ConvertType_Params());
		else return RType();
	}

	//------------------------------- WRITE VALUE TO THIS

	//set value, e.g. Any a; a = 5;
	//If type of a is already set, the external value is converted to this type. e.g. if Any a(5); a = 6.5; results in a value of 6 stored, not 6.5, since a has been set to store a type integer.
	template <typename RType>
	Any& operator=(const RType& value)
	{
		if (pValue) {

			MatchType_CallFunc<RType, SetValue_Params<RType>>(SetValue_Params<RType>(value));
		}
		else {

			pValue = new RType(value);
			allocated = true;

			type_name = btype_info<RType>().name();
		}

		return *this;
	}

	//assignment operator - make an exact copy. This means: independent if object to copy from is independent (i.e. has allocated == true), else keep it as a dependent object - this means the value set in pValue is a pointer to some value allocated elsewhere.
	Any& operator=(const Any& copyThis)
	{
		if (copyThis.allocated) {

			//make new Any object, independent of copyThis
			CopyAny(copyThis);
		}
		else {

			if (pValue && allocated)
				MatchType_CallFunc<void, DeleteType_Params>(DeleteType_Params());

			//copyThis points to a value allocated elsewhere. Keep this object the same.
			allocated = false;
			pValue = copyThis.pValue;
			type_name = copyThis.type_name;
		}

		return *this;
	}

	//------------------------------- CONVERSION FROM/TO STRING OR STREAMS

	std::string convert_to_string(const std::string& unit = "")
	{
		if (pValue)
			return MatchType_CallFunc<std::string, ConvertToString_Params>(ConvertToString_Params(unit));
		else return "";
	}

	//this allows conversion to std::string using << via a std::stringstream (and also needed by ToString method). 
	//Reason for using const_cast:
	//const Any& rhs is required (ToString receives a const Type& rhs). convert_to_string() cannot be made const since it calls MatchType_CallFunc which cannot be made const (it must change values when handling other type of calls)
	//Thus cannot call convert_to_string using a const Any object (rhs) => must const_cast to get rid of the const. So, why not const_cast<Any>(rhs) ? Due to the copy constructor (not fully clear on this last one).
	friend std::ostream& operator<<(std::ostream &os, const Any& rhs) { os << const_cast<Any*>(&rhs)->convert_to_string(); return os; }

	//convert from std::string using a std::stringstream (streamable in)
	friend Any& operator>>(const std::stringstream &ss, Any &rhs)
	{
		rhs.convert_string(ss.str());
		return rhs;
	}

	//convert std::string to currently set value type. You can also use ToNum, but need to specify the type, e.g. : Any a = (DBL3)ToNum("1.0"); If the type is already set but not sure what type then a.convert_string(...); can be used
	void convert_string(const std::string& text, const std::string& unit = "")
	{
		//if this is an empty Any then cannot use this : need to know what type to convert to - use convert_string_set_type instead
		if (!text.length() || !pValue) return;

		MatchType_CallFunc<void, ConvertFromString_Params>(ConvertFromString_Params(text, unit));
	}

	//convert std::string to currently set value type. You can also use ToNum, but need to specify the type, e.g. : Any a = (DBL3)ToNum("1.0"); If the type is already set but not sure what type then a.convert_string(...); can be used
	void convert_string_set_type(const std::string& text, const std::string& new_type)
	{
		if (!text.length() || !new_type.length()) return;

		//set type name
		type_name = new_type;

		//delete currently held value so we can set a new type from the std::string and specified type name (if type name not correct then this Any will in effect be cleared)
		if (pValue && allocated)
			MatchType_CallFunc<void, DeleteType_Params>(DeleteType_Params());

		pValue = nullptr;
		allocated = false;

		MatchType_CallFunc<void, ConvertFromString_Params>(ConvertFromString_Params(text, ""));
	}

	//convert to currently set value type by using as many components as possible from the start of the std::vector. e.g. for a simple type (int, double, etc.) a single component is used, but DBL3 uses 3, Rect and Box use 6.
	//Note the constructors of multi-component types allow use of fewer components (e.g. DBL3 can be initialized using 1 component) - also convert in these cases if number of components in the std::vector is less than the maximum possible.
	//return number of components used
	int convert_string(const std::vector<std::string>& components, const std::string& unit = "")
	{
		if (!components.size() || !pValue) return 0;

		return MatchType_CallFunc<int, ConvertFromStringComponents_Params>(ConvertFromStringComponents_Params(components, unit));
	}

	//------------------------------- COMPARISON OPERATORS

	// >=

	template <typename Type, std::enable_if_t<!is_string<Type>::value>* = nullptr>
	bool operator>=(Type value) const { return *reinterpret_cast<Type*>(pValue) >= value; }

	// <=

	template <typename Type, std::enable_if_t<!is_string<Type>::value>* = nullptr>
	bool operator<=(Type value) const { return *reinterpret_cast<Type*>(pValue) <= value; }

	// >

	template <typename Type, std::enable_if_t<!is_string<Type>::value>* = nullptr>
	bool operator>(Type value) const { return *reinterpret_cast<Type*>(pValue) > value; }

	// <

	template <typename Type, std::enable_if_t<!is_string<Type>::value>* = nullptr>
	bool operator<(Type value) const { return *reinterpret_cast<Type*>(pValue) < value; }

	// ==

	template <typename Type, std::enable_if_t<!is_string<Type>::value>* = nullptr>
	bool operator==(Type value) const { return *reinterpret_cast<Type*>(pValue) == value; }

	// !=

	template <typename Type, std::enable_if_t<!is_string<Type>::value>* = nullptr>
	bool operator!=(Type value) const { return *reinterpret_cast<Type*>(pValue) != value; }

	//------------------------------- PROPERTIES

	bool IsNull(void) const { return (pValue == nullptr); }

	template <typename Type>
	bool is_type(const btype_info<Type>& ti) const { return (type_name == ti.name()); }

	std::string get_type(void) const { return type_name; }

	void clear(void)
	{
		if (pValue && allocated)
			MatchType_CallFunc<void, DeleteType_Params>(DeleteType_Params());

		pValue = nullptr;
		allocated = false;
		type_name = "";
	}
};