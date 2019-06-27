#pragma once

#include <tuple>

//TO BE REVISED IF NEEDED
/*
//Convenient object for handling tuples. Use as follows:

//First make a tuple_set object from your tuple (call it tup) as:

//				MAKE_TUPLE_SET(tup, tup_set);

//The size of the tuple can be obtained as:

//				tup_set.size()
//
//Set values as (also works in loops etc, just as you would expect):
//
//				tup_set[index] = value;
//
//To read values at a specific index, declare a variable then use the instruction below (type must match with stored type). 
//This is only useful really if all types stored in the tuple are the same, then reading can be done in a loop. For different types probably better to stick to auto value = get<index>(tup);
//
//				GET_TUPLE_ENTRY(tup_set, index, variable);
//
//Alternatively a void* can be returned using:
//
//				void* point_to_entry = tup_set(index);
//
//There's a specialisation for using VarInfo types in a tuple.
//
//The following gets the name in VarInfo entries (VarInfo::name string member)
//
//				tup_set.get_tuple_vinfo(index);
//
//Can also index VarInfo entries with a name, so a value can be set in the VarInfo::value member
//
//				tup_set[name] = value;
//
//Additionally, can search if a VarInfo with given name exists as:
//
//				int index = get_tuple_vinfo_index(name); //if index >= 0 the entry exists and is located at index
//
//When setting values in VarInfo entries, if a string is used, this is automatically converted to the VarInfo::value member if the conversion works (i.e. if ToNum works) :
//
//				tup_set[index] = value_as_a_string; 
//				tup_set[name] = value_as_a_string;

template <size_t N, typename Tuple>
class tuple_set {

private:

	Tuple & tup;
	int access_index = 0;

private:

	//------------------------------------------------------- Return void* to tuple entry at access_index. Use GET_TUPLE_ENTRY macro.

	template <std::size_t... I>
	void* get_tuple_entry(index_sequence<I...>)
	{
		return get_parse_tuple(0, get<I>(tup)...);
	}

	template <typename Type>
	void* get_parse_tuple(int parse_index, Type& entry)
	{
		if (parse_index == access_index) return &entry;
		else return nullptr;
	}

	template <typename Type, typename ... PType>
	void* get_parse_tuple(int parse_index, Type& entry, PType& ... further_entries)
	{
		if (parse_index == access_index) {

			return &entry;
		}
		else {

			parse_index++;
			return get_parse_tuple(parse_index, further_entries...);
		}
	}

	//------------------------------------------------------- Return name of variable stored in tuple entry at access_index, which is expected to be stored as VarInfo (return empty string if not a VarInfo). Use get_tuple_vinfo method

	template <std::size_t... I>
	string get_tuple_vinfo_entry(index_sequence<I...>)
	{
		return get_parse_vinfo_tuple(0, get<I>(tup)...);
	}

	template <typename Type>
	string get_parse_vinfo_tuple(int parse_index, Type& entry)
	{
		if (parse_index == access_index) return get_parse_vinfo_name(entry, is_varinfo<Type>());
		else return "";
	}

	template <typename Type, typename ... PType>
	string get_parse_vinfo_tuple(int parse_index, Type& entry, PType& ... further_entries)
	{
		if (parse_index == access_index) {
			
			return get_parse_vinfo_name(entry, is_varinfo<Type>());
		}
		else {

			parse_index++;
			return get_parse_vinfo_tuple(parse_index, further_entries...);
		}
	}

	//the entry is a VarInfo
	template <typename Type>
	string get_parse_vinfo_name(VarInfo<Type>& vi, true_type)
	{
		return vi.name;
	}

	//the entry is not a VarInfo
	template <typename Type>
	string get_parse_vinfo_name(Type& vi, false_type)
	{
		return "";
	}

	//------------------------------------------------------- Search for name in all VarInfo entries. Return index of first found name (-1 if not found). 
	//														  Index using a string : ["some name"], which should be followed by an assignment (if "some name" not found the assignment does nothing)

	template <std::size_t... I>
	int search_tuple_vinfo(const string& name, index_sequence<I...>)
	{
		return search_tuple_vinfo_parse(name, 0, get<I>(tup)...);
	}

	template <typename Type>
	int search_tuple_vinfo_parse(const string& name, int parse_index, Type& entry)
	{
		if (check_tuple_vinfo_name(name, entry, is_varinfo<Type>())) {

			return parse_index;
		}
		else return -1;	//not found
	}

	template <typename Type, typename ... PType>
	int search_tuple_vinfo_parse(const string& name, int parse_index, Type& entry, PType& ... further_entries)
	{
		if (check_tuple_vinfo_name(name, entry, is_varinfo<Type>())) {
			
			return parse_index; 
		}
		else {

			parse_index++;
			return search_tuple_vinfo_parse(name, parse_index, further_entries...);
		}
	}
	
	//the entry is a VarInfo
	template <typename Type>
	bool check_tuple_vinfo_name(const string& name, VarInfo<Type>& vi, true_type)
	{
		return (vi.name == name);
	}

	//the entry is not a VarInfo
	template <typename Type>
	bool check_tuple_vinfo_name(const string& name, Type& vi, false_type)
	{
		return false;
	}

	//------------------------------------------------------- Set value in tuple entry at access_index. Use as tup_set[index] = value;

	template <typename VType, std::size_t... I>
	void access_tuple_entry(const VType& setThis, index_sequence<I...>)
	{
		parse_tuple(setThis, 0, get<I>(tup)...);
	}

	template <typename VType, typename Type>
	void parse_tuple(const VType& setThis, int parse_index, Type& entry)
	{
		if (parse_index == access_index) set_tuple_entry_varinfo(setThis, entry, is_varinfo<Type>());
	}

	template <typename VType, typename Type, typename ... PType>
	void parse_tuple(const VType& setThis, int parse_index, Type& entry, PType& ... further_entries)
	{
		if (parse_index == access_index) {

			set_tuple_entry_varinfo(setThis, entry, is_varinfo<Type>());
		}
		else {

			parse_index++;
			parse_tuple(setThis, parse_index, further_entries...);
		}
	}

	//-----

	template <typename VType, typename Type>
	struct tag_selector {
		static const int value =
			1 * is_string<VType>::value +
			2 * is_string<Type>::value +
			//anything is convertible to an Any since it defines a templated conversion operator. Want to exclude this case here.
			4 * (is_convertible<VType, Type>::value && !is_same<Any, Type>::value) +
			8 * is_streamable_in<stringstream, Type>::value;
	};

	//a VarInfo
	template <typename VType, typename Type>
	void set_tuple_entry_varinfo(const VType& setThis, VarInfo<Type>& entry, true_type)
	{
		set_tuple_entry(setThis, entry.value, int_tag< tag_selector<VType, Type>::value >());
	}

	//Not a VarInfo
	template <typename VType, typename Type>
	void set_tuple_entry_varinfo(const VType& setThis, Type& entry, false_type)
	{
		set_tuple_entry(setThis, entry, int_tag< tag_selector<VType, Type>::value >());
	}

	//----- Possible cases:

	//string to string
	template <typename VType, typename Type>
	void set_tuple_entry(const VType& setThis, Type& value, int_tag<15>) { value = setThis; }

	//string to a non-string which can be streamed to
	template <typename VType, typename Type>
	void set_tuple_entry(const VType& setThis, Type& value, int_tag<9>)
	{
		stringstream ss(setThis);
		ss >> value;				//do not use ToNum as Type could be Any. >> operator will work with it also.
	}

	//convertible types where neither are strings (but may or may not be able to stream to Type)
	template <typename VType, typename Type>
	void set_tuple_entry(const VType& setThis, Type& value, int_tag<4>) { value = setThis; }	//implicit conversion (dont use c-style cast, or any cast, as Type may be any Any - in this case it will know what to do).
	template <typename VType, typename Type>
	void set_tuple_entry(const VType& setThis, Type& value, int_tag<12>) { value = setThis; }	//implicit conversion (dont use c-style cast, or any cast, as Type may be any Any - in this case it will know what to do).

	//sink-hole for the other tag values (some of which cannot occur)
	template <typename VType, typename Type>
	void set_tuple_entry(const VType& setThis, Type& value, ...) {}

public:

	tuple_set(Tuple& tup_) : tup(tup_), access_index(0)
	{}

	//indexing for writing value to tuple
	tuple_set& operator[](const int& access_index_)
	{
		access_index = access_index_;
		return *this;
	}

	//indexing for writing value to tuple - index using a string, specifically used for VarInfo entries, where the name will be searched against the name stored in VarInfo entries
	tuple_set& operator[](const string& name)
	{
		access_index = search_tuple_vinfo(name, make_index_sequence<N>{});
		
		return *this;
	}

	//write value to tuple (should index beforehand)
	template <typename VType>
	void operator=(const VType& setThis)
	{
		if(access_index >= 0 && access_index < N)
			access_tuple_entry(setThis, make_index_sequence<N>{});
	}

	//indexing for reading value from tuple (use GET_TUPLE_ENTRY macro instead of this, unless you just want the void*)
	void* operator()(const int& access_index_)
	{
		access_index = access_index_;
		return get_tuple_entry(make_index_sequence<N>{});
	}

	//indexing when the entry at access_index is a VarInfo : return the name of the variable stored in VarInfo
	string get_tuple_vinfo(const int& access_index_)
	{
		access_index = access_index_;
		return get_tuple_vinfo_entry(make_index_sequence<N>{});
	}

	int get_tuple_vinfo_index(const string& name)
	{
		access_index = search_tuple_vinfo(name, make_index_sequence<N>{});
		return access_index;
	}

	int size(void) { return N; }
};

#define MAKE_TUPLE_SET(tup, tup_set) tuple_set<tuple_size<decltype(tup)>::value, decltype(tup)> tup_set(tup)
#define GET_TUPLE_ENTRY(tup_set, index, to_variable) to_variable = *reinterpret_cast<decltype(to_variable)*>(tup_set(index))
*/