#pragma once

//
// Save/load "program state" to/from file
//
// What can be saved:
//
// 1. simple types (see is_simple_type trait checker below). (Save name then value. Load by checking name then value.)
// 2. complex types. (Save name then call SaveObjectState method for further action. Load by checking name then LoadObjectState for further action.)
// 3. pointers: (Save name for non-std::vector entries only. Load by checking name for non-vector entries only. For save/load do one of the following things below).
//  -> nullptr (Save "nullptr". Load by setting nullptr)
//  -> pointer to a complex type (save "!nullptr". Load by making instance with appropriate constructor - see int_tag_select_ctor below. Call Save/LoadObjectState method on that instance)
//  -> pointer to a non-complex type (save blank line. Load by making instance with appropriate constructor.) 
//  -> pointer to the non-complex base of a derived type, where one of the implementations is specified in the implementations std::tuple. (Save implementation name. Load by making instance of implementation with appropriate constructor. If implementation is a complex type then call Save/LoadObjectState method on that instance.)
// 4. vectors containing any of the types above. (Save name, size. Load by checking name then resize using size. Depending on entry type do one of the following things above).
//

#include <string>
#include <vector>
#include <tuple>
#include <cassert>

#include "Types.h"
#include "Types_Any.h"
#include "Introspection.h"
#include "Funcs_Vectors.h"

#define FILEROWCHARS	5000	//maximum number of characters per input file row

//need this so ProgramState can be specialized with 2 parameter packs ... also need it for the is_complex_type struct below.
template <typename ... >
class ProgramState {};

//---------------------------------------------------

//a simple type is defined as anything that is std::stringstream streamable in and out, is not a vector-like type and is not a pointer type.
//e.g. all fundamental types, strings and custom non-vector types with >> and << defined for std::stringstream
template <typename Type>
struct is_simple_type : std::bool_constant<
	is_streamable_in<std::stringstream, Type>::value &&
	is_streamable_out<std::stringstream, Type>::value &&
	!is_vector<Type>::value &&
	!std::is_pointer<Type>::value> {};

//a complex type is any class that inherits from ProgramState
template <typename Type>
struct is_complex_type : std::bool_constant< is_template_base_of<::ProgramState, Type>::value > {};

//---------------------------------------------------

template <typename OType, typename ... PType, typename ... IPType>
class ProgramState<OType, std::tuple<PType...>, std::tuple<IPType...> > {

private: //--------------------------------------------------- DATA

	//signal complex type end
	const std::string endType = "end type";

	//signal binary data block. This will be followed by number of BYTEs in the block (need this to keep compatibility with older save files - if this string encountered when not expected then jump over the binary block)
	const std::string binaryData = "bin data";

	//variables are saved using their names and variable references, and this info is contained in VarInfo. Need a std::tuple to expand the parameter pack.
	std::tuple< VarInfo<PType>... > objects;

	//list all possible implementations for abstract base classes here, available in this scope
	std::tuple< ImplInfo<IPType>... > implementations;

	//pointer to the owner object (this is the object which contains all the variables in the objects std::tuple. OType is the class which inherits this ProgramState)
	OType* pOwner;

	//test some traits of EntryType
	template <typename EntryType>
	struct int_tag_select
	{
	private:

		using EntryType_ = typename std::remove_reference<EntryType>::type;

		static const bool t1 = is_simple_type<EntryType_>::value;
		static const bool t2 = std::is_same<EntryType_, Any>::value;			//an Any is a simple type but needs special treatment (must save/load type_name)
		static const bool t3 = is_vector<EntryType_>::value;
		static const bool t4 = is_complex_type<EntryType_>::value;
		static const bool t5 = std::is_pointer<EntryType_>::value;

		static constexpr int parse_tests(void)
		{
			//... and here.
			return get_value(0, 1, t1, t2, t3, t4, t5);
		}

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

		//the tag value for tag dispatching. Highest test true value is taken (meaning highest by test number t1, t2, etc.)
		static const int value = parse_tests();
	};

	//test for particular constructors in CType - used to instantiate new objects when pointers are found in the objects std::tuple
	template <typename CType>
	struct int_tag_select_ctor
	{
	private:

		static const bool t1 = accepts_parameters<OType*>::value;			//has constructor taking a pointer of the owner type (see above OType* pOwner)
																			//also works if OType* is a pointer to a derived type with CType ctor defined for base pointer only
		static const bool t2 = accepts_parameters<const OType*>::value;
		static const bool t3 = accepts_parameters<CType>::value;			//has void constructor	- try this first if found

		static constexpr int parse_tests(void)
		{
			//... and here.
			return get_value(0, 1, t1, t2, t3);
		}

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

public:

	//use this so you don't have to type the template parameters more than once
	using ProgramStateNames = ProgramState<OType, std::tuple<PType...>, std::tuple<IPType...> >;

private: //--------------------------------------------------- METHODS

	//--------------------------------------------------- OBJECT FACTORY

	//OType* pOwner
	template <typename BaseType, typename DerivedType>
	void make_instance(BaseType** ppObj, int_tag<1> tag)
	{
		__make_instance<BaseType, DerivedType>(ppObj, tag, can_instantiate_with_param<BaseType, DerivedType, OType*>());
	}

	template <typename BaseType, typename DerivedType>
	void __make_instance(BaseType** ppObj, int_tag<1>, std::true_type)
	{
		*ppObj = dynamic_cast<BaseType*>(new DerivedType(pOwner));
	}

	//---

	//const OType* pOwner
	template <typename BaseType, typename DerivedType>
	void make_instance(BaseType** ppObj, int_tag<2> tag)
	{
		__make_instance<BaseType, DerivedType>(ppObj, tag, can_instantiate_with_param<BaseType, DerivedType, OType*>());
	}

	//const OType* pOwner
	template <typename BaseType, typename DerivedType>
	void __make_instance(BaseType** ppObj, int_tag<2>, std::true_type)
	{
		*ppObj = dynamic_cast<BaseType*>(new DerivedType(pOwner));
	}

	//---

	//void
	template <typename BaseType, typename DerivedType>
	void make_instance(BaseType** ppObj, int_tag<3> tag)
	{
		__make_instance<BaseType, DerivedType>(ppObj, tag, can_instantiate<BaseType, DerivedType>());
	}

	//void
	template <typename BaseType, typename DerivedType>
	void __make_instance(BaseType** ppObj, int_tag<3>, std::true_type)
	{
		*ppObj = dynamic_cast<BaseType*>(new DerivedType());
	}

	template <typename BaseType, typename DerivedType>
	void make_instance(...) {}

	template <typename BaseType, typename DerivedType>
	void __make_instance(...) {}

	//--------------------------------------------------- TUPLE INTERFACE PARSING TO FIND IMPLEMENTATION OF AN ABSTRACT BASE CLASS - SAVE

	//Start
	template <typename PointerType, typename Tuple, int... I>
	bool parse_implementations_save(std::ofstream &bdout, PointerType& pBase, Tuple& tup, std::index_sequence<I...> is)
	{
		return parse_implementations_save(bdout, pBase, tup, is, std::bool_constant< (bool)std::tuple_size<Tuple>::value >());
	}

	template <typename PointerType, typename Tuple, int... I>
	bool parse_implementations_save(std::ofstream &bdout, PointerType& pBase, Tuple& tup, std::index_sequence<I...>, std::true_type)
	{
		if (!parse_implementations_save(bdout, pBase, std::get<I>(tup)...)) {

			//implementation name not found
			bdout << endl;
			//complex or non-complex type?
			using CType = typename std::remove_pointer<typename std::remove_reference<PointerType>::type>::type;
			return parse_implementations_save(bdout, pBase, is_complex_type<CType>());
		}
		return true;
	}

	template <typename PointerType, typename Tuple, int... I>
	bool parse_implementations_save(std::ofstream &bdout, PointerType& pBase, Tuple& tup, std::index_sequence<I...>, std::false_type)
	{
		//zero length index sequence (empty std::tuple) - nothing to parse
		//complex or non-complex type?
		using CType = typename std::remove_pointer<typename std::remove_reference<PointerType>::type>::type;
		return parse_implementations_save(bdout, pBase, is_complex_type<CType>());
	}

	//Methods used to iterate over elements in std::tuple
	template <typename PointerType, typename Type>
	bool parse_implementations_save(std::ofstream &bdout, PointerType& pBase, Type& entry)
	{
		using IType = typename std::remove_pointer<typename std::remove_reference<decltype(entry.value)>::type>::type;	//just need the type of the interface
		auto pDownCast = dynamic_cast<IType*>(pBase);															//if dynamic_cast cannot match the two objects it will set pDownCast = nullptr

		if (pDownCast) {

			//Successful downcast of base class : found implementation. Save implementation name so we can remake this pointer when loading (it can change implementation or become nullptr)
			bdout << entry.name << endl;

			//complex or non-complex type?
			return parse_implementations_save(bdout, pDownCast, is_complex_type<IType>());
		}

		return false;
	}

	template <typename PointerType, typename Type, typename ... PType>
	bool parse_implementations_save(std::ofstream &bdout, PointerType& pBase, Type& entry, PType& ... further_entries)
	{
		if (parse_implementations_save(bdout, pBase, entry)) return true;
		else return parse_implementations_save(bdout, pBase, further_entries...);
	}

	//---

	//pointer to a complex type
	template <typename PointerType>
	bool parse_implementations_save(std::ofstream &bdout, PointerType& pointer, std::true_type)
	{
		pointer->SaveObjectState(bdout);
		return true;
	}

	//pointer but not to a complex type
	template <typename PointerType>
	bool parse_implementations_save(std::ofstream &bdout, PointerType& pointer, std::false_type)
	{
		return true;
	}

	//--------------------------------------------------- TUPLE INTERFACE PARSING TO FIND IMPLEMENTATION OF AN ABSTRACT BASE CLASS - LOAD

	//Start
	template <typename PointerType, typename Tuple, int... I>
	bool parse_implementations_load(std::ifstream &bdin, PointerType& pointer, std::string implementation_name, Tuple& tup, std::index_sequence<I...> is)
	{
		return parse_implementations_load(bdin, pointer, implementation_name, tup, is, std::bool_constant< (bool)std::tuple_size<Tuple>::value >());
	}

	template <typename PointerType, typename Tuple, int... I>
	bool parse_implementations_load(std::ifstream &bdin, PointerType& pointer, std::string implementation_name, Tuple& tup, std::index_sequence<I...>, std::true_type)
	{
		if (!parse_implementations_load(bdin, pointer, implementation_name, std::get<I>(tup)...)) {

			//implementation not found : pointer cannot be a pointer to an abstract base, and is also not a pointer to a complex type, go ahead and make an instance then return
			using TypeClass = typename std::remove_pointer<typename std::remove_reference<PointerType>::type>::type;

			if (pointer) {

				delete pointer;
				pointer = nullptr;
			}

			make_instance<TypeClass, TypeClass>(&pointer, int_tag< int_tag_select_ctor<TypeClass>::value >());
		}

		return true;
	}

	template <typename PointerType, typename Tuple, int... I>
	bool parse_implementations_load(std::ifstream &bdin, PointerType& pointer, std::string implementation_name, Tuple& tup, std::index_sequence<I...>, std::false_type)
	{
		//zero length index sequence (empty std::tuple) - nothing to do
		//complex or non-complex type?
		using CType = typename std::remove_pointer<typename std::remove_reference<PointerType>::type>::type;
		return parse_implementations_load(bdin, pointer, is_complex_type<CType>());
	}

	//---

	//Methods used to iterate over elements in std::tuple
	template <typename PointerType, typename Type>
	bool parse_implementations_load(std::ifstream &bdin, PointerType& pBase, std::string implementation_name, Type& entry)
	{
		if (entry.name == implementation_name) {

			using TypeClass = typename std::remove_pointer<typename std::remove_reference<PointerType>::type>::type;
			using IType = typename std::remove_pointer<typename std::remove_reference<decltype(entry.value)>::type>::type;

			if (pBase) {

				delete pBase;
				pBase = nullptr;
			}

			make_instance<TypeClass, IType>(&pBase, int_tag<int_tag_select_ctor<IType>::value>());

			auto pUpCast = dynamic_cast<IType*>(pBase);
			return parse_implementations_load(bdin, pUpCast, is_complex_type<IType>());
		}
		else return false;		//implementation not found
	}

	template <typename PointerType, typename Type, typename ... PType>
	bool parse_implementations_load(std::ifstream &bdin, PointerType& pBase, std::string implementation_name, Type& entry, PType& ... further_entries)
	{
		if (entry.name == implementation_name) {

			return parse_implementations_load(bdin, pBase, implementation_name, entry);
		}

		return parse_implementations_load(bdin, pBase, implementation_name, further_entries...);
	}

	//---

	//upcast pointer after finding and instantiating implementation : implementation is a complex type so load it
	template <typename PointerType>
	bool parse_implementations_load(std::ifstream& bdin, PointerType& pUpCast, std::true_type)
	{
		return pUpCast->LoadObjectState(bdin);
	}

	//upcast pointer after finding and instantiating implementation : implementation is not a complex type so nothing else to do here
	template <typename PointerType>
	bool parse_implementations_load(std::ifstream& bdin, PointerType& pUpCast, std::false_type)
	{
		return true;
	}

	//--------------------------------------------------- TUPLE PARSING FOR SAVING, ELEMENT BY ELEMENT IN ORDER

	//--------------
	//Parse
	template <typename Tuple, int... I>
	void parse_save_tuple(std::ofstream &bdout, Tuple& tup, std::index_sequence<I...> is)
	{
		//it is possible the objects std::tuple is empty 
		//(e.g. you have a std::vector containing pointers to abstract base classes with particular implementations.
		//It could be there's nothing to save in some or all of the implementations, but you still want this std::vector structure to be remade, so all implementations need to inherit from ProgramState.
		//Those that have nothing to save simply inherit as : public ProgramState< std::tuple<>, std::tuple<> >
		parse_save_tuple(bdout, tup, is, std::bool_constant< (bool)std::tuple_size<Tuple>::value >());
	}

	template <typename Tuple, int... I>
	void parse_save_tuple(std::ofstream &bdout, Tuple& tup, std::index_sequence<I...>, std::true_type)
	{
		parse_save_tuple(bdout, std::get<I>(tup)...);
	}

	template <typename Tuple, int... I>
	void parse_save_tuple(std::ofstream &bdout, Tuple& tup, std::index_sequence<I...>, std::false_type)
	{
		//zero length index sequence (empty std::tuple) - nothing to do
	}

	//Methods used to iterate over elements in std::tuple
	template <typename Type>
	void parse_save_tuple(std::ofstream &bdout, Type& entry)
	{
		bdout << entry.name << endl;
		save_tuple_entry(bdout, entry.value, int_tag< int_tag_select<decltype(entry.value)>::value >());
	}

	template <typename Type, typename ... PType>
	void parse_save_tuple(std::ofstream &bdout, Type& entry, PType& ... further_entries)
	{
		parse_save_tuple(bdout, entry);
		parse_save_tuple(bdout, further_entries...);
	}

	//--------------
	//Save std::tuple entry

	//entry value is a simple type
	template <typename Type>
	void save_tuple_entry(std::ofstream& bdout, Type& value, int_tag<1>, bool vector_entry = false)
	{
		//only vector entries (strings and Any excepted) are saved using binary
		if (std::is_same<Type, std::string>::value || !vector_entry) {

			bdout << value << endl;
		}
		else {

			char* binary_data = reinterpret_cast<char*>(&value);
			bdout.write(binary_data, sizeof(Type));
		}
	}

	//entry value is a simple type, in particular an Any
	template <typename Type>
	void save_tuple_entry(std::ofstream& bdout, Type& value, int_tag<2>, bool vector_entry = false)
	{
		bdout << value << endl;
		bdout << value.get_type() << endl;
	}

	//a std::vector type
	template <typename Type>
	void save_tuple_entry(std::ofstream& bdout, Type& value, int_tag<3>, bool vector_entry = false)
	{
		typedef std::remove_reference<Type>::type VType;
		typedef contained_type<VType>::type SType;

		VType vec = value;

		//stored types that can be saved are simple types, complex types (have SaveObjectState), or pointers (to complex types or base of a complex derived type).
		//Nested vectors not currently implemented.

		//signal std::vector start. size() returns the dimensions of the std::vector (the type is not necesarily an int - could be an INT3. It is the same type as that taken by the resize method, so can use ToNum method to convert when loading the std::vector later)
		bdout << vec.size() << endl;
		
		//all vectors next contain 2 lines : first there's the binary data string, then the number of BYTEs (excluding any key and id lines)
		size_t numBYTEs = (vec.end() - vec.begin()) * sizeof(SType);
		bdout << binaryData << endl;
		bdout << ToString(numBYTEs) << endl;

		bool vec_with_key = is_vector_withkey<VType>::value;
		bool vec_with_Id = is_vector_withId<VType>::value;

		//use iterators to parse the std::vector, as this method doesn't depend on the std::vector dimensions type, but works for all std::vector types. Thus need begin() and end().
		int idx = 0;
		for (auto it = vec.begin(); it != vec.end(); ++it, idx++) {

			//if this is a std::vector with key and/or Id then save key and/or Id in a separate line each
			if (vec_with_key) {

				save_tuple_entry_vector_key(bdout, vec, idx, is_vector_withkey<VType>());
			}

			if (vec_with_Id) {

				save_tuple_entry_vector_Id(bdout, vec, idx, is_vector_withId<VType>());
			}

			//depending on stored type must do different things : save in binary only if this is a simple vector (no key or Id) with simple entries
			save_tuple_entry(bdout, *it, int_tag< int_tag_select<SType>::value >(), !vec_with_key && !vec_with_Id);
		}
	}

	//-----

	//std::vector has key indexing
	template <typename Type>
	void save_tuple_entry_vector_key(std::ofstream& bdout, Type& vec, int idx, std::true_type)
	{
		std::string key = vec.get_key_from_index(idx);
		bdout << key << endl;
	}

	template <typename Type>
	void save_tuple_entry_vector_key(std::ofstream& bdout, Type& vec, int idx, std::false_type) {} //need to keep compiler happy

																								   //std::vector has Id indexing
	template <typename Type>
	void save_tuple_entry_vector_Id(std::ofstream& bdout, Type& vec, int idx, std::true_type)
	{
		INT2 Id = vec.get_id_from_index(idx);
		bdout << Id << endl;
	}

	template <typename Type>
	void save_tuple_entry_vector_Id(std::ofstream& bdout, Type& vec, int idx, std::false_type) {} //need to keep compiler happy

	//---

	//Complex object to be broken down further
	template <typename Type>
	void save_tuple_entry(std::ofstream &bdout, Type& value, int_tag<4>, bool vector_entry = false)
	{
		value.SaveObjectState(bdout);
	}

	//a pointer
	template <typename Type>
	void save_tuple_entry(std::ofstream& bdout, Type& value, int_tag<5>, bool vector_entry = false)
	{
		//class type pointed to
		using CType = typename std::remove_pointer<typename std::remove_reference<Type>::type>::type;

		//first check if nullptr. when loading this pointer will also be set to nullptr (first deleted if not nullptr on loading)
		if (value == nullptr) {

			bdout << "nullptr" << endl;
			return;
		}

		//next check if pointer to a complex type
		save_tuple_entry_pointer(bdout, value, is_complex_type<CType>());
	}

	//---Do something with the pointer depending on type pointed to

	//a pointer : to a complex type
	template <typename Type>
	void save_tuple_entry_pointer(std::ofstream& bdout, Type& value, std::true_type pcomplexType)
	{
		bdout << "!nullptr" << endl;		//means pointer to a complex object (not null)
		value->SaveObjectState(bdout);
	}

	//a pointer : not to a complex type. Check if this is a pointer to a base where one of the implementations is found in the implementations std::tuple. If so then save name, else save a blank line.
	template <typename Type>
	void save_tuple_entry_pointer(std::ofstream& bdout, Type& value, std::false_type pcomplexType)
	{
		parse_implementations_save(bdout, value, implementations, std::make_index_sequence<sizeof ...(IPType)>{});
	}

	//---

	//sink-hole
	template <typename Type>
	void save_tuple_entry(std::ofstream& bdout, Type& value, ...) {}

	//--------------------------------------------------- TUPLE PARSING FOR LOADING : FIND VARIABLE NAME

	template <typename Tuple, int... I>
	bool parse_load_tuple(std::ifstream& bdin, const std::string& var_name, int& entry_count, Tuple& tup, std::index_sequence<I...> is)
	{
		//see comments on parse_save_tuple for why this check is needed
		return parse_load_tuple(bdin, var_name, entry_count, tup, is, std::bool_constant< (bool)std::tuple_size<Tuple>::value >());
	}

	template <typename Tuple, int... I>
	bool parse_load_tuple(std::ifstream& bdin, const std::string& var_name, int& entry_count, Tuple& tup, std::index_sequence<I...>, std::true_type)
	{
		return parse_load_tuple(bdin, var_name, entry_count, std::get<I>(tup)...);
	}

	template <typename Tuple, int... I>
	bool parse_load_tuple(std::ifstream& bdin, const std::string& var_name, int& entry_count, Tuple& tup, std::index_sequence<I...>, std::false_type)
	{
		//zero length index sequence (empty std::tuple) - nothing to do
		return true;
	}

	//Methods used to iterate over elements in std::tuple
	template <typename Type>
	bool parse_load_tuple(std::ifstream& bdin, const std::string& var_name, int& entry_count, Type& entry)
	{
		if (var_name == entry.name) {

			//found variable name : try to load it
			entry_count++;
			return load_tuple_entry(bdin, entry.value, int_tag< int_tag_select<decltype(entry.value)>::value >(), entry.keep_ptr);
		}
		else return true;	//not found variable but still return true to attempt further loading (must be an old save version with this extra variable no longer in use)
	}

	template <typename Type, typename ... PType>
	bool parse_load_tuple(std::ifstream& bdin, const std::string& var_name, int& entry_count, Type& entry, PType& ... further_entries)
	{
		if (var_name == entry.name)
			return parse_load_tuple(bdin, var_name, entry_count, entry);
		else
			return parse_load_tuple(bdin, var_name, entry_count, further_entries...);
	}

	//--------------
	//Load std::tuple entry

	//entry value is a simple type
	template <typename Type>
	bool load_tuple_entry(std::ifstream& bdin, Type& value, int_tag<1>, bool keep_ptr = false, bool vector_entry = false)
	{
		//only simple vector entries (strings and Any excepted) are loaded from binary data
		if (is_string<Type>::value) {

			char line[FILEROWCHARS];
			if (bdin.getline(line, FILEROWCHARS)) {

				*reinterpret_cast<std::string*>(&value) = std::string(line);
			}
			else return false;
		}
		else {

			if (vector_entry) {

				bdin.read(reinterpret_cast<char*>(&value), sizeof(Type));
			}
			else {

				char line[FILEROWCHARS];
				if (bdin.getline(line, FILEROWCHARS)) {
					
					value = (Type)ToNum(std::string(line));
				}
				else return false;
			}
		}

		return true;
	}

	//entry value is a simple type, in particular an Any
	template <typename Type>
	bool load_tuple_entry(std::ifstream& bdin, Type& value, int_tag<2>, bool keep_ptr = false, bool vector_entry = false)
	{
		char line[FILEROWCHARS];

		if (bdin.getline(line, FILEROWCHARS)) {

			std::string read_text = std::string(line);
			if (bdin.getline(line, FILEROWCHARS)) {

				std::string type_name = std::string(line);
				value.convert_string_set_type(read_text, type_name);
				return true;
			}
			else return false;
		}
		else return false;
	}

	//load std::vector
	template <typename Type>
	bool load_tuple_entry(std::ifstream& bdin, Type& value, int_tag<3>, bool keep_ptr = false, bool vector_entry = false)
	{
		//stored elements type
		using SType = typename contained_type<Type>::type;

		char line[FILEROWCHARS];

		if (bdin.getline(line, FILEROWCHARS)) {

			std::string vec_size_text = std::string(line);

			clear_vector(value);						//deep clear the std::vector : if it contains pointers then delete them before clearing
			value.resize(ToNum(vec_size_text));

			//all vectors next contain 2 lines : first there's the binary data string, then the number of BYTEs (excluding any key and id lines). Skip over these as not needed.
			if (!bdin.getline(line, FILEROWCHARS)) return false;
			if (!bdin.getline(line, FILEROWCHARS)) return false;

			bool vec_with_key = is_vector_withkey<Type>::value;
			bool vec_with_Id = is_vector_withId<Type>::value;

			int idx = 0;
			for (auto it = value.begin(); it != value.end(); ++it, idx++) {

				if (vec_with_key) {

					if (bdin.getline(line, FILEROWCHARS)) load_vector_key(std::string(line), value, idx, is_vector_withkey<Type>());
					else return false;
				}

				if (vec_with_Id) {

					if (bdin.getline(line, FILEROWCHARS)) load_vector_Id(std::string(line), value, idx, is_vector_withId<Type>());
					else return false;
				}

				if (!load_tuple_entry(bdin, *it, int_tag< int_tag_select<SType>::value >(), keep_ptr, !vec_with_key && !vec_with_Id)) return false;
			}

			return true;
		}
		else return false;
	}

	//----

	//set key in std::vector
	template <typename Type>
	void load_vector_key(const std::string& line, Type& vec, int index, std::true_type)
	{
		vec.set_key_at_index(line, index);
	}

	template <typename Type>
	void load_vector_key(const std::string& line, Type& vec, int index, std::false_type) {} //never gets here but must keep compiler happy

																							//set Id in std::vector
	template <typename Type>
	void load_vector_Id(const std::string& line, Type& vec, int index, std::true_type)
	{
		vec.set_id_at_index(ToNum(line), index);
	}

	template <typename Type>
	void load_vector_Id(const std::string& line, Type& vec, int index, std::false_type) {} //never gets here but must keep compiler happy

	//-----

	//Complex object to be broken down further
	template <typename Type>
	bool load_tuple_entry(std::ifstream& bdin, Type& value, int_tag<4>, bool keep_ptr = false, bool vector_entry = false)
	{
		return value.LoadObjectState(bdin);
	}

	//pointer
	template <typename Type>
	bool load_tuple_entry(std::ifstream& bdin, Type& value, int_tag<5>, bool keep_ptr = false, bool vector_entry = false)
	{
		char line[FILEROWCHARS];
		if (bdin.getline(line, FILEROWCHARS)) {

			std::string read_text = std::string(line);
			if (read_text == "nullptr") {

				if (value) delete value;
				value = nullptr;
				return true;
			}
			else {

				//if read_text is "!nullptr" then it must be a pointer to a complex type - check type directly
				//if not read_text is blank then it must be a pointer to a non-complex type - again, check type directly
				//if read_text has a name then it must be a pointer to a base where the implementation name is held in the implementations std::tuple
				return load_pointer(bdin, value, read_text, is_complex_type<typename std::remove_pointer<Type>::type>(), keep_ptr);
			}
		}
		else return false;
	}

	//----

	//pointer to a complex type - remake object (implementation_name should contain "!nullptr", but not needed)
	template <typename Type>
	bool load_pointer(std::ifstream& bdin, Type& pointer, std::string implementation_name, std::true_type pcomplexType, bool keep_ptr = false)
	{
		using TypeClass = typename std::remove_pointer<typename std::remove_reference<Type>::type>::type;

		if (!keep_ptr) {

			//delete...
			if (pointer) {

				delete pointer;
				pointer = nullptr;
			}

			//... and remake
			make_instance<TypeClass, TypeClass>(&pointer, int_tag< int_tag_select_ctor<TypeClass>::value >());
		}

		if (pointer)
			return pointer->LoadObjectState(bdin);
		else return true;
	}

	//does not point to a complex type - try to downcast it to a possible implementation (from the implementations std::tuple)
	template <typename Type>
	bool load_pointer(std::ifstream& bdin, Type& pointer, std::string implementation_name, std::false_type pcomplexType, bool keep_ptr = false)
	{
		return parse_implementations_load(bdin, pointer, implementation_name, implementations, std::make_index_sequence<sizeof ...(IPType)>{});
	}

	//sink-hole
	template <typename Type>
	bool load_tuple_entry(std::ifstream& bdin, Type& value, ...) { return true; }

protected:

	ProgramState(OType* pOwner_, std::tuple< VarInfo<PType>... > variables, std::tuple< ImplInfo<IPType>... > implementations_)
		: pOwner(pOwner_), objects(variables), implementations(implementations_)
	{}

public:

	void SaveObjectState(std::ofstream &bdout)
	{
		parse_save_tuple(bdout, objects, std::make_index_sequence<sizeof...(PType)>{});
		bdout << endType << endl;
	}

	bool LoadObjectState(std::ifstream &bdin)
	{
		char line[FILEROWCHARS];

		for (int tup_entries_count = 0; tup_entries_count < std::tuple_size<decltype(objects)>::value;) {

			if (bdin.getline(line, FILEROWCHARS)) {

				if (std::string(line) == binaryData) {

					//not expecting to see this here, must be an older version save file (older than program version) : jump over the binary block (next line gives the number of BYTEs in the block)
					if (bdin.getline(line, FILEROWCHARS)) {

						int numBYTEs = ToNum(std::string(line));
						bdin.ignore(numBYTEs);
					}
					else return false;
				}

				//check for end of type - if found, return true to attempt further loading
				if (std::string(line) == endType) {

					RepairObjectState();
					return true;
				}

				if (!parse_load_tuple(bdin, std::string(line), tup_entries_count, objects, std::make_index_sequence<sizeof...(PType)>{})) {

					RepairObjectState();
					return false;
				}
			}
			else {

				RepairObjectState();
				return false;	//reached end of file, nothing else can be loaded
			}
		}

		//all possible loading done for this std::tuple, as the tup_set size limit has been reached. 
		//If this std::tuple was part of a higher-up object, there must be an endType std::string before we can return, so keep checking for this:

		while (bdin.getline(line, FILEROWCHARS)) {

			if (std::string(line) == endType) {

				RepairObjectState();
				return true;
			}
		}

		RepairObjectState();
		//reached end of file, nothing else can be loaded
		return false;
	}

	//after loading the object, some dependent quantities may need to be recalculated (dependent quantities are those that haven't been saved and can be set entirely from the ones that have been saved). Force client to consider this - must implement.
	virtual void RepairObjectState(void) = 0;
};