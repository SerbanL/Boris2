#pragma once

#include <string>

//This is a portable version of std::type_info
//Instead of using typeid(...).name(), make btype_info<...> object, then call name() on it. e.g. btype_info<double>().name()
//
//Since the names are hard-coded this is platform-independent.
//This is mainly intended for use with Types_Any.h, which will result in different ProgramState output save files on different OS 
//When you introduce a new type which Types_Any can handle, add a new specialization for btype_info etc. 
//It doesn't really matter what std::string you hard-code below, but to keep it consistent enter the std::string produced on Windows (naming is more sensible than on Linux anyway, although more verbose, but these are not used in performance-critical situations).

template <typename Type>
class btype_info {

private:

	std::string type_name = "";

public:

	btype_info() {}

	std::string name(void) const { return type_name; }
};

template <>
inline btype_info<void>::btype_info()
{
	type_name = "void";
}

template <>
inline btype_info<bool>::btype_info()
{
	type_name = "bool";
}

template <>
inline btype_info<char>::btype_info()
{
	type_name = "char";
}

template <>
inline btype_info<int>::btype_info()
{
	type_name = "int";
}

template <>
inline btype_info<float>::btype_info()
{
	type_name = "float";
}

template <>
inline btype_info<double>::btype_info()
{
	type_name = "double";
}

template <>
inline btype_info<class std::basic_string<char, struct std::char_traits<char>, class std::allocator<char> >>::btype_info()
{
	type_name = "class std::basic_string<char,struct std::char_traits<char>,class std::allocator<char> >";
}

template <>
inline btype_info<struct VAL2<int>>::btype_info()
{
	type_name = "struct VAL2<int>";
}

template <>
inline btype_info<struct VAL2<float>>::btype_info()
{
	type_name = "struct VAL2<float>";
}

template <>
inline btype_info<struct VAL2<double>>::btype_info()
{
	type_name = "struct VAL2<double>";
}

template <>
inline btype_info<struct VAL3<int>>::btype_info()
{
	type_name = "struct VAL3<int>";
}

template <>
inline btype_info<struct VAL3<float>>::btype_info()
{
	type_name = "struct VAL3<float>";
}

template <>
inline btype_info<struct VAL3<double>>::btype_info()
{
	type_name = "struct VAL3<double>";
}

template <>
inline btype_info<struct VAL4<int>>::btype_info()
{
	type_name = "struct VAL4<int>";
}

template <>
inline btype_info<struct VAL4<float>>::btype_info()
{
	type_name = "struct VAL4<float>";
}

template <>
inline btype_info<struct VAL4<double>>::btype_info()
{
	type_name = "struct VAL4<double>";
}

template <>
inline btype_info<struct __Box<void>>::btype_info()
{
	type_name = "struct __Box<void>";
}

template <>
inline btype_info<struct __Rect<void>>::btype_info()
{
	type_name = "struct __Rect<void>";
}

template <>
inline btype_info<class Sequence<double>>::btype_info()
{
	type_name = "class Sequence<double>";
}

template <>
inline btype_info<class Sequence<struct VAL3<double> >>::btype_info()
{
	type_name = "class Sequence<struct VAL3<double> >";
}

template <>
inline btype_info<class SequencePolar<void>>::btype_info()
{
	type_name = "class SequencePolar<void>";
}

template <>
inline btype_info<class CosSequence<double>>::btype_info()
{
	type_name = "class CosSequence<double>";
}

template <>
inline btype_info<class CosSequence<struct VAL3<double> >>::btype_info()
{
	type_name = "class CosSequence<struct VAL3<double> >";
}

template <>
inline btype_info<class SinOscillation<double>>::btype_info()
{
	type_name = "class SinOscillation<double>";
}

template <>
inline btype_info<class CosOscillation<double>>::btype_info()
{
	type_name = "class CosOscillation<double>";
}

template <>
inline btype_info<class SinOscillation<struct VAL3<double> >>::btype_info()
{
	type_name = "class SinOscillation<struct VAL3<double> >";
}

template <>
inline btype_info<class CosOscillation<struct VAL3<double> >>::btype_info()
{
	type_name = "class CosOscillation<struct VAL3<double> >";
}

template <>
inline btype_info<class StringSequence>::btype_info()
{
	type_name = "class StringSequence";
}

template <>
inline btype_info<class FileSequence<double>>::btype_info()
{
	type_name = "class FileSequence<double>";
}

template <>
inline btype_info<class FileSequence<struct VAL3<double> >>::btype_info()
{
	type_name = "class FileSequence<struct VAL3<double> >";
}