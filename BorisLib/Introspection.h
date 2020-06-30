#pragma once

#include "Introspection_base.h"
#include "Types_VAL.h"

//--------------------------------------------------
//check if type is Id indexable (has [] operator which can take an INT2 as the index)

//helper
template <typename Type>
struct __is_indexable_byId {

private:

	template <typename Type_>
	static auto test(bool) -> decltype(std::declval<Type_>()[std::declval<INT2>()], std::true_type());

	template <typename>
	static auto test(...)->std::false_type;

public:

	static const bool value = decltype(test<Type>(true))::value;
};

//use this
template <typename Type>
struct is_indexable_byId : std::integral_constant<bool, __is_indexable_byId<Type>::value > {};

//--------------------------------------------------
//check if type is key indexable (has [] operator which can a std::string as parameter)

//helper
template <typename Type>
struct __is_indexable_bykey {

private:

	template <typename Type_>
	static auto test(bool) -> decltype(std::declval<Type_>()[std::declval<std::string>()], std::true_type());

	template <typename>
	static auto test(...)->std::false_type;

public:

	static const bool value = decltype(test<Type>(true))::value;
};

//use this
template <typename Type>
struct is_indexable_bykey : std::integral_constant<bool, __is_indexable_bykey<Type>::value > {};

//--------------------------------------------------
//check if type is a std::vector-like type : duck-typing test to allow for custom std::vector-like types. For a type to be considered a vector it must have the following:
//
// 1. be integer indexable : [index]
// 2. have these methods with no parameters : clear(), begin(), end(), size()
// 3. have the resize method which takes the same parameter type as that returned by size()
// 4. is not a std::string

//helper
template <typename Type>
struct __is_vector {

private:

	template <typename Type_>
	static auto test(bool) -> decltype(
		std::declval<Type_>().begin(),
		std::declval<Type_>().end(),
		std::declval<Type_>().clear(),
		std::declval<Type_>().size(),
		std::declval<Type_>().resize(std::declval<Type_>().size()),
		std::true_type()
		);

	template <typename>
	static auto test(...)->std::false_type;

public:

	static const bool value =
		is_indexable<Type>::value &&
		decltype(test<Type>(true))::value &&
		!is_string<Type>::value;
};

//use this
template <typename Type>
struct is_vector : std::integral_constant<bool, __is_vector<Type>::value > {};

//--------------------------------------------------
//check if type is a std::vector-like type with Id indexing: is_vector and is_indexable_byId

template <typename Type>
struct is_vector_withId : std::integral_constant<bool, is_vector<Type>::value && is_indexable_byId<Type>::value > {};

//--------------------------------------------------
//check if type is a std::vector-like type with key indexing: is_vector and is_indexable_bykey

template <typename Type>
struct is_vector_withkey : std::integral_constant<bool, is_vector<Type>::value && is_indexable_bykey<Type>::value > {};

