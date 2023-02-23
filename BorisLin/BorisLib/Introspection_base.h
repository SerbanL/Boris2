#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <type_traits>
#include <utility>

////////////////////////////////////////////////////////////////////////////////////////////////// GENERALLY USEFUL CONSTRUCTS
//
//

//Function pointer to member function call helper
#define CALLFP(objectPointer, function) ((objectPointer)->*function)

//--------------------------------------------------
//Use for tag dispatching when std::true_type, std::false_type is not enough

template <int value>
struct int_tag {};

//use this in conjuction with a tag_selector built similar to e.g.:

/*
template <typename Type>
struct int_tag_select {

private:

//List all tests here, labelled as t1, t2, ... The tests should be either orthogonal or else increasing in specialisation.
static const bool t1 = is_streamable_in<std::stringstream, Type>::value;
static const bool t2 = is_pointer<Type>::value;
//...

static constexpr int parse_tests(void)
{
//... and here.
return get_value(0, 1, t1, t2);
}

//Do not change this
template <typename TType>
static constexpr int get_value(int tag, int test_number, TType test)
{
return (test == true ? test_number : tag);
}

template <typename TType, typename... PType>
static constexpr int get_value(int tag, int test_number, TType test, PType ... tests)
{
tag = get_value(tag, test_number, test);
return get_value(tag, ++test_number, tests...);
}

public:

//the tag value for tag dispatching
static const int value = parse_tests();
};

//e.g. call overloaded methods with an int_tag<number> tag as : int_tag< tag_selector<MyType>::value >()
*/

//--------------------------------------------------

template <typename Type>
struct VarInfo {

	//reference to variable
	Type& value;
	//stringified name of variable, used for loading values from a previously saved file
	std::string name;
	//if Type is a pointer the default behaviour is to remake it (delete then new). If keep_ptr == true then keep it as it is, just load anything in pointed-to-object if a complex object.
	bool keep_ptr = false;

	VarInfo(Type& value_, const std::string& name_, bool keep_ptr_ = false)
		: value(value_), keep_ptr(keep_ptr_)
	{
		//When used with ProgramState it's possible the variable scope is specified using "this->" or "basename<...>::"
		//In this case we only want the actual variable name, not the scope name as well
		//This happens if ProgramState is used on a templated class derived from a templated base, which requires specifying the templated base scope.
		size_t pos = name_.find("::");
		if (pos != std::string::npos) name = name_.substr(pos + 2);
		else {

			pos = name_.find("->");
			if (pos != std::string::npos) name = name_.substr(pos + 2);
			else name = name_;
		}
	}
};

template <typename Type>
struct ImplInfo {

	//keep implementation type
	Type* value = nullptr;
	//stringified name of implementation, used for loading from a previously saved file
	std::string name;

	ImplInfo(const std::string& name_)
		: name(name_)
	{}
};

//helper macro for constructing a VarInfo object
#define VINFO(variable) VarInfo<decltype(variable)>(variable, #variable)

#define VINFO_KEEPPTR(variable) VarInfo<decltype(variable)>(variable, #variable, true)

//#define IINFO(Implementation) VarInfo<Implementation>(Implementation(), #Implementation)
#define IINFO(Implementation) ImplInfo<Implementation>(#Implementation)

//convert anything to a std::string (name of variable, function, type)
#define STRINGIFY(name) #name

#define CLASS_STR(name) #name

//--------------------------------------------------

namespace Introspection {

	//The UnnamedType is a lambda which takes a parameter (the object for which we need to check if it has a given method), and return the type returned by the given method
	template <typename LambdaClosure>
	class Introspection_Helper {

	private:

		//in first field of decltype expression we try to do the following : instantiate the lambda and call it with an instance of the class we need to check if it has a given method.
		//the lambda return type is the type returned by the given method, thus if the method exists the first field in the decltype expression is valid : return std::true_type as per the last decltype field
		template <typename Type>
		constexpr auto testValidity(bool) const -> decltype(std::declval<LambdaClosure>()(std::declval<Type>()), std::true_type())
		{
			return std::true_type();
		}

		//SFINAE sink-hole : false if check above failed
		template <typename Type>
		constexpr auto testValidity(...) const -> std::false_type
		{
			return std::false_type();
		}

	public:

		//The check is done here
		template <typename Type>
		constexpr auto operator()(Type) const
		{
			return testValidity<Type>(true);
		}
	};

	//this method receives a lambda
	template <typename LambdaClosure>
	constexpr auto is_valid(LambdaClosure)
	{
		return Introspection_Helper<LambdaClosure>();
	}
};

#define HAS_METHOD(Type, method) Introspection::is_valid([](auto x) -> decltype(x.method()) {})(Type())
#define HAS_MEMBER(Type, member) Introspection::is_valid([](auto x) -> decltype(x.member) {})(Type())

#define POINTER_HAS_METHOD(pointer, method) Introspection::is_valid([](auto x) -> decltype(x->method()) {})(pointer)

//--------------------------------------------------
//check if type is a std::string or a reference to a std::string (irrespective of cv-qualifiers)

template <typename Type>
struct is_string :
	std::integral_constant
	<bool,
	std::is_same<typename std::remove_cv<Type>::type, std::string>::value ||
	std::is_same<Type, std::string&>::value ||
	std::is_same<Type, const std::string&>::value ||
	std::is_same<Type, volatile std::string&>::value
	>
	//std::remove_cv doesn't remove the const or volatile qualifiers from a const type& or a const type*
{};

//--------------------------------------------------
//check if type has >> or << operators for a given stream type

//helper
template <typename stream, typename Type>
struct __is_streamable_out {

private:

	template <typename stream_, typename Type_>
	static auto test(bool) -> decltype(std::declval<stream_&>() << std::declval<Type_&>(), std::true_type());

	template <typename, typename>
	static auto test(...)->std::false_type;

public:

	static const bool value = decltype(test<stream, Type>(true))::value;
};

//use this, e.g. is_streamable_out<std::string, double>() has base std::true_type(), etc.
template <typename stream, typename Type>
struct is_streamable_out : std::integral_constant<bool, __is_streamable_out<stream, Type>::value > {};

//---------

//helper
template <typename stream, typename Type>
struct __is_streamable_in {

private:

	template <typename stream_, typename Type_>
	static auto test(bool) -> decltype(std::declval<stream_&>() >> std::declval<Type_&>(), std::true_type());

	template <typename, typename>
	static auto test(...)->std::false_type;

public:

	static const bool value = decltype(test<stream, Type>(true))::value;
};

//use this
template <typename stream, typename Type>
struct is_streamable_in : std::integral_constant<bool, __is_streamable_in<stream, Type>::value > {};

//--------------------------------------------------
//check if type is indexable (has [] operator which can take an integer)

//helper
template <typename Type>
struct __is_indexable {

private:

	template <typename Type_>
	static auto test(bool) -> decltype(std::declval<Type_>()[std::declval<int>()], std::true_type());

	template <typename>
	static auto test(...)->std::false_type;

public:

	static const bool value = decltype(test<Type>(true))::value;
};

//use this
template <typename Type>
struct is_indexable : std::integral_constant<bool, __is_indexable<Type>::value > {};

//--------------------------------------------------
// check if type is a VarInfo type. Done as follows:
// 1. A VarInfo<Type> must have a member called value which has Type& type
// 2. Remove the reference
// 3. Construct a VarInfo with the type obtained above
// 4. This must have the same type as the template parameter. If any of these steps fail or final result is false, then template parameter it is not a VarInfo

//helper
template <typename Type>
struct __is_varinfo {

private:

	template <typename Type_>
	static auto test(bool)
		-> decltype(std::is_same<Type_, VarInfo< typename std::remove_reference< decltype(std::declval<Type_>().value) >::type >>());

	template <typename>
	static auto test(...)->std::false_type;

public:

	static const bool value = decltype(test<Type>(true))::value;
};

//use this
template <typename Type>
struct is_varinfo : std::integral_constant<bool, __is_varinfo<Type>::value > {};

//--------------------------------------------------
//get the type of stored elements in a mono container (here a container means anything that can be integer indexed : [index], i.e. is_indexable<VType>). Mono means all contained elements have the same type (e.g. std::string or std::vector etc.)
//CType is the container type. 
//contained_type<CType>::type is the contained element type
//if this is invoked with a type that cannot be indexed a compiler error is issued - good (so no need to check for is_indexable<CType>)

template <typename CType>
struct contained_type {

	//can index as [0] even for empty conatainers since in evaluating decltype(std::declval<CType>()[0]) no memory access actually occurs
	//need to std::remove_reference since the above returns a reference to the stored type in the container
	using type = typename std::remove_reference<decltype(std::declval<CType>()[0])>::type;
};

//--------------------------------------------------
//Use this to check if a concrete class type is public derived from an unspecified version of a template base class - single inheritance only, else ambiguous
//
// e.g. is_template_base_of<::Base, Derived>(); where Base is a templated class and Derived publicly inherits from a Base specialization
// :: is necessary when this expression is used within Derived or Base, otherwise Base template parameters are deduced from that context: the scope resolution operator stops that, and Base is then deduced from the type casting in the function call as intended (if valid)

//helper
template <template <typename...> class TemplateBase, typename Derived>
struct __is_template_base_of {

private:

	template <template <typename...> class __TemplateBase, typename... TemplateParams>
	static std::true_type test(const __TemplateBase<TemplateParams...>*);

	template <template <typename...> class __TemplateBase>
	static std::false_type test(...);

public:

	static const bool value = decltype(test<TemplateBase>(std::declval<Derived*>()))::value;
};

//use this
template <template <typename...> class TemplateBase, typename Derived>
struct is_template_base_of : std::integral_constant<bool, __is_template_base_of<TemplateBase, Derived>::value > {};

//--------------------------------------------------
//Check if Type has a constructor that takes the parameter types listed

//helper
template <typename Type, typename... PType>
struct __accepts_parameters {

private:

	template <typename Type_>
	static auto test(bool) -> decltype(Type_(std::declval<PType>()...), std::true_type());

	template <typename>
	static auto test(...) -> std::false_type;

public:

	static const bool value = decltype(test<Type>(true))::value;
};

//use this
template <typename Type, typename... PType>
struct accepts_parameters : std::integral_constant<bool, __accepts_parameters<Type, PType...>::value > {};

//--------------------------------------------------
//Check if Type has a constructor that takes N simple parameters, convertible to from an integer (so includes all fundamental types, std::strings)

//helper
template <typename Type, int N>
struct __accepts_N_parameters {

private:

	template <typename Type_, int... I>
	static auto test(std::integer_sequence<int, I...>) -> decltype(Type_(I...), std::true_type());

	template <typename>
	static auto test(...) -> std::false_type;

public:

	static const bool value = decltype(test<Type>(std::make_integer_sequence<int, N>{}))::value;
};

//use this
template <typename Type, int N>
struct accepts_N_parameters : std::integral_constant<bool, __accepts_N_parameters<Type, N>::value > {};

//--------------------------------------------------
//Get maximum number of parameters a Type constructor can take, where each parameter can be converted to from an integer

template <typename Type>
struct constructor_parameters {

private:

	static constexpr int find_max(void)
	{
		int max = 0;

		//allow maximum of 9 parameters checking - tried using index sequences to check each value but compiler kept complaining (something to be improved if needed)
		if (accepts_N_parameters<Type, 1>::value) max = 1;
		if (accepts_N_parameters<Type, 2>::value) max = 2;
		if (accepts_N_parameters<Type, 3>::value) max = 3;
		if (accepts_N_parameters<Type, 4>::value) max = 4;
		if (accepts_N_parameters<Type, 5>::value) max = 5;
		if (accepts_N_parameters<Type, 6>::value) max = 6;
		if (accepts_N_parameters<Type, 7>::value) max = 7;
		if (accepts_N_parameters<Type, 8>::value) max = 8;
		if (accepts_N_parameters<Type, 9>::value) max = 9;

		return max;
	}

public:

	static const int maximum = find_max();

	static void get_constructors(std::vector<int>& vec)
	{
		if (accepts_N_parameters<Type, 1>::value) vec.push_back(1);
		if (accepts_N_parameters<Type, 2>::value) vec.push_back(2);
		if (accepts_N_parameters<Type, 3>::value) vec.push_back(3);
		if (accepts_N_parameters<Type, 4>::value) vec.push_back(4);
		if (accepts_N_parameters<Type, 5>::value) vec.push_back(5);
		if (accepts_N_parameters<Type, 6>::value) vec.push_back(6);
		if (accepts_N_parameters<Type, 7>::value) vec.push_back(7);
		if (accepts_N_parameters<Type, 8>::value) vec.push_back(8);
		if (accepts_N_parameters<Type, 9>::value) vec.push_back(9);
	}
};

//--------------------------------------------------
//

//helper
template <typename BaseType, typename DerivedType, typename Param>
struct __can_instantiate_with_param {

private:

	template <typename BaseType_, typename DerivedType_, typename Param_>
	static auto test(bool) -> decltype(dynamic_cast<BaseType_*>(new DerivedType_(std::declval<Param_>())), std::true_type());

	template <typename BaseType_, typename DerivedType_, typename Param_>
	static auto test(...)->std::false_type;

public:

	static const bool value = decltype(test<BaseType, DerivedType, Param>(true))::value;
};

//use this
template <typename BaseType, typename DerivedType, typename Param>
struct can_instantiate_with_param : std::integral_constant<bool, __can_instantiate_with_param<BaseType, DerivedType, Param>::value > {};

template <typename BaseType, typename DerivedType>
struct __can_instantiate {

private:

	template <typename BaseType_, typename DerivedType_>
	static auto test(bool) -> decltype(dynamic_cast<BaseType_*>(new DerivedType_()), std::true_type());

	template <typename BaseType_, typename DerivedType_>
	static auto test(...)->std::false_type;

public:

	static const bool value = decltype(test<BaseType, DerivedType>(true))::value;
};

//use this
template <typename BaseType, typename DerivedType>
struct can_instantiate : std::integral_constant<bool, __can_instantiate<BaseType, DerivedType>::value > {};