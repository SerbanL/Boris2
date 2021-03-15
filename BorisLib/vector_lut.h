////////////////////////////////////////////////////////////////////////////////////////////////// VECTOR with LOOK-UP TABLE
//
// std::vector storage with look-up table : there are two ways of indexing this: 1) as a normal std::vector, 2) using an identifier parameter (see below)
//
// each element in the std::vector has a type integer identifier (these will be a value from some enum and the caller specifies this value for each new entry) - this is called the major id
// there can be multiple entries of same type - these are not specified by the caller, but are assigned here and returned to caller - this is called the minor id
// vector_lut allows usage as with the usual stl std::vector, but in addition elements can be recalled using (major id, minor id) indexing : this is done internally using a look-up table
//
// !!!NOTE!!!: no checking for indexes out of range done here when indexing the std::vector: the user must ensure this

#pragma once

#include <string>
#include <algorithm>
#include <utility>

#include "Types_VAL.h"

template <typename Type> 
class vector_lut {

protected:

	//the actual std::vector
	std::vector<Type> storedObjects;

	//this has the same size as storedObjects with 1-2-1 correspondence between the entries: this LUT gives the majorId, minorId of the entry
	std::vector<INT2> LUT_byIndex;

	//look-up-table for storedObjects : first index is the majorId, second index is the minorId
	std::vector<std::vector<int>> LUT_byId;

private:

	//for given majorId return a minorId value and make sure LUT_byId has the correct size and is ready to accept an assignment at [INT2(majorId, minorId)] entry
	int setup_LUT_byID(int majorId);
	//make space for full Id
	void setup_LUT_byId(INT2 Id);

	//if no elements with majorId exist (and majorId is a valid index in LUT_byId) then clear the vector at LUT_byId[majorId]
	void clean_LUT_byId(int majorId = -1);

public:

	//------------------------------------ CONSTRUCTORS

	vector_lut(void) {}
	virtual ~vector_lut() {}

	//copy constructor for lvalues
	vector_lut(const vector_lut<Type>& copyThis) :
		storedObjects(copyThis.storedObjects),
		LUT_byIndex(copyThis.LUT_byIndex),
		LUT_byId(copyThis.LUT_byId)
	{}

	//overloaded copy constructor for rvalues using move semantics
	vector_lut(vector_lut<Type>&& moveThis) :
		storedObjects(std::move(moveThis.storedObjects)),
		LUT_byIndex(std::move(moveThis.LUT_byIndex)),
		LUT_byId(std::move(moveThis.LUT_byId))
	{}

	//assignment operator for lvalues
	vector_lut<Type>& operator=(const vector_lut<Type>& copyThis) 
	{
		storedObjects = copyThis.storedObjects;
		LUT_byIndex = copyThis.LUT_byIndex;
		LUT_byId = copyThis.LUT_byId;

		return *this;
	}

	//overloaded assignment operator for rvalues using move semantics
	vector_lut<Type>& operator=(vector_lut<Type>&& moveThis) 
	{
		storedObjects = std::move(moveThis.storedObjects);
		LUT_byIndex = std::move(moveThis.LUT_byIndex);
		LUT_byId = std::move(moveThis.LUT_byId);

		return *this;
	}

	//------------------------------------ INDEXING

	//usual std::vector indexing
	Type& operator[](const int &index) { return storedObjects[index]; }
	const Type& operator[](const int &index) const { return storedObjects[index]; }

	//indexing via look-up table with INT2(majorId, minorId). Note, can also call it with INT2(majorId) : this means minorId = 0, and is used when majorId is known to be unique.
	Type& operator[](const INT2 &id) { return storedObjects[LUT_byId[id.i][id.j]]; }
	const Type& operator[](const INT2 &id) const { return storedObjects[LUT_byId[id.i][id.j]]; }

	//indexing via major id only, where minorId is assumed to be zero (this is used for vector_lut objects which are guaranteed to have unique majorId entries)
	Type& operator()(const int &majorId) { return storedObjects[LUT_byId[majorId][0]]; }
	const Type& operator()(const int &majorId) const { return storedObjects[LUT_byId[majorId][0]]; }

	Type& front(void) { return storedObjects[0]; }
	Type& back(void) { return storedObjects[last()]; }

	//------------------------------------ ITERATORS

	Type* begin(void) { return &storedObjects[0]; }
	Type* data(void) { return &storedObjects[0]; }
	Type* end(void) { return &storedObjects[size()]; }

	//------------------------------------ SIZE MODIFICATION METHODS

	//add new object at the end and return the minor id
	int push_back(Type newObject, int majorId = 0);

	void pop_back(void);

	//insert object at given index and return minor id
	int insert(int index, Type newObject, int majorId = 0);

	//set value at full id entry: overwrite any existing entry, else create a new one.
	void set(Type newObject, INT2 id);

	void erase(INT2 id) { erase(LUT_byId[id.i][id.j]); }
	void erase(int index);

	void clear(void);

	void shrink_to_fit(void);

	//set new size. Any new entries beyond the old size will be assigned a majorId of 0.
	void resize(int new_size, Type value = Type());

	//------------------------------------ CONTROL OF ENTRIES

	//move element from source to destination index (i.e. ordering in std::vector changes)
	void move(int srcIdx, int dstIdx = 0);

	//------------------------------------ GET PROPERTIES

	//dimension of internal std::vector
	int size(void) const { return (int)storedObjects.size(); }

	//index of last element
	int last(void) const { if (size()) return ((int)size() - 1); else return -1; }

	//------------------------------------ MULTIPLE ENTRIES WITH SAME MAJOR ID

	//get number of elements with given major ID
	int get_num_IDs(int majorID) const;

	//get a minor id for element with majorID; return -1 if none found
	int get_next_minor_id(int majorID) const;

	//make all minor ids for elements with given major ID start from 0 and be numbered in steps of 1
	//use this e.g. when you want to ensure there's always a minor if of 0 if major id exists
	//DO NOT use if you are keeping track of minor ids externally
	void regularise_minor_ids(int majorID);

	//------------------------------------ CHECKERS

	//is the given id currently set? (i.e. is there a stored element with this id?)
	bool is_id_set(INT2 id) const { return (id.i >= 0 && id.i < (int)LUT_byId.size() && id.j >= 0 && id.j < (int)LUT_byId[id.i].size() && LUT_byId[id.i][id.j] >= 0); }

	//is the major id set ? (assume minor id = 0)
	bool is_ID_set(int majorID) const { return (majorID >= 0 && majorID < (int)LUT_byId.size() && LUT_byId[majorID][0] >= 0); }

	//does this std::vector contain the given value?
	bool has_value(Type value) const { for (int idx = 0; idx < (int)storedObjects.size(); idx++) if (storedObjects[idx] == value) return true; return false; }

	//------------------------------------ ID <--> INDEX GETTERS

	//--- from value

	//return Id of found value ((-1,-1) if not found)
	INT2 get_id_from_value(Type value) const { for (int idx = 0; idx < (int)storedObjects.size(); idx++) if (storedObjects[idx] == value) return get_id_from_index(idx); return INT2(-1, -1); }

	//return major ID of found value (-1 if not found)
	int get_ID_from_value(Type value) const { for (int idx = 0; idx < (int)storedObjects.size(); idx++) if (storedObjects[idx] == value) return get_ID_from_index(idx); return -1; }

	//--- from Id

	//return index of element with given id - this assumes the id exists
	int get_index_from_id(INT2 id) const { if (id.i >= 0 && id.i < (int)LUT_byId.size() && id.j >= 0 && id.j < (int)LUT_byId[id.i].size()) return LUT_byId[id.i][id.j]; else return -1; }

	//return index of element with given major id (assume minor id = 0) - this further assumes the id exists
	int get_index_from_ID(int majorID) const { if (majorID >= 0 && majorID < (int)LUT_byId.size()) return LUT_byId[majorID][0]; else return -1; }

	//--- from index

	//return id of element with given index - this assumes the index is valid
	INT2 get_id_from_index(int idx) const { if (idx < (int)LUT_byIndex.size()) return LUT_byIndex[idx]; else return INT2(-1, -1); }

	//return major ID of element with given index - this assumes the index is valid
	int get_ID_from_index(int idx) const { if (idx < (int)LUT_byIndex.size()) return LUT_byIndex[idx].major; else return -1; }

	//------------------------------------ ID SETTERS

	//over-ride existing Id with new one. No checking is done for clashing with existing Ids.
	void set_id_at_index(INT2 Id, int index);
};

//------------------------------------ PRIVATE

template <typename Type>
int vector_lut<Type>::setup_LUT_byID(int majorId)
{
	//find winId minor for given winId major
	int minorId = -1;

	//first must make sure LUT has the right size.
	if (majorId >= (int)LUT_byId.size()) {

		//This entry type has not been entered before: make new entry in LUT

		//make a single entry to hold the idx -1 : i.e. this entry does not hold a valid storedObjects index, meaning this type is not stored
		std::vector<int> unassignedEntry;
		unassignedEntry.push_back(-1);

		//make space in LUT with unassigned entries
		LUT_byId.resize(majorId + 1, unassignedEntry);

		minorId = 0;
	}
	else {

		//find first minorId available for this majorId type
		for (int i = 0; i < (int)LUT_byId[majorId].size(); i++) {

			if (LUT_byId[majorId][i] < 0) { minorId = i; break; }
		}

		//if a minorId value not found then need more space in LUT
		if (minorId == -1) {

			LUT_byId[majorId].push_back(-1);
			minorId = (int)LUT_byId[majorId].size() - 1;
		}
	}

	return minorId;
}

template <typename Type>
void vector_lut<Type>::setup_LUT_byId(INT2 Id)
{
	if (Id.major >= (int)LUT_byId.size()) {

		std::vector<int> unassignedEntry;
		unassignedEntry.push_back(-1);

		//make space in LUT with unassigned entries
		LUT_byId.resize(Id.major + 1, unassignedEntry);
	}

	if (Id.minor >= (int)LUT_byId[Id.major].size()) {

		//make space in LUT with unassigned entries
		LUT_byId[Id.major].resize(Id.minor + 1, -1);
	}
}

//if no elements with majorId exist (and majorId is a valid index in LUT_byId) then clear the vector at LUT_byId[majorId]
template <typename Type>
void vector_lut<Type>::clean_LUT_byId(int majorId)
{
	if (majorId >= 0 && majorId < LUT_byId.size() && get_num_IDs(majorId) == 0) {

		LUT_byId[majorId].clear();
	}

	//clean all
	if (majorId < 0) {

		for (int major = 0; major < LUT_byId.size(); major++) {

			if (get_num_IDs(major) == 0) LUT_byId[major].clear();
		}
	}
}

//------------------------------------ SIZE MODIFICATION METHODS

template <typename Type>
void vector_lut<Type>::resize(int new_size, Type value)
{
	int old_size = size();

	//all entries from index new_size up will become unassigned in LUT_byId
	for (int idx = new_size; idx < old_size; idx++) {

		INT2 Id = LUT_byIndex[idx];
		LUT_byId[Id.major][Id.minor] = -1;
	}

	clean_LUT_byId();

	//resize storage
	storedObjects.resize(new_size, value);

	//resize LUT by index : if new size is larger then all new entries are assigned a majorId of 0 in the loop below
	LUT_byIndex.resize(new_size, -1);

	for (int idx = new_size; idx < old_size; idx++) {

		int minorId = setup_LUT_byID(0);
		LUT_byIndex[idx] = INT2(0, minorId);
		LUT_byId[0][minorId] = idx;
	}
}

template <typename Type>
int vector_lut<Type>::push_back(Type newObject, int majorId)
{
	//make sure LUT_byId is set-up correctly for new entry
	int minorId = setup_LUT_byID(majorId);

	//push object in std::vector and store its id
	storedObjects.push_back(newObject);
	LUT_byIndex.push_back(INT2(majorId, minorId));

	//store index of pushed object in LUT
	LUT_byId[majorId][minorId] = (int)storedObjects.size() - 1;

	return minorId;
}

template <typename Type>
void vector_lut<Type>::pop_back(void)
{
	INT2 popId = LUT_byIndex[(int)storedObjects.size() - 1];
	LUT_byId[popId.i][popId.j] = -1;

	storedObjects.pop_back();

	LUT_byIndex.pop_back();

	clean_LUT_byId();
}

template <typename Type>
int vector_lut<Type>::insert(int index, Type newObject, int majorId)
{
	//make sure LUT_byId is set-up correctly for new entry
	int minorId = setup_LUT_byID(majorId);

	//insert object in std::vector at given index and store its id
	storedObjects.insert(storedObjects.begin() + index, newObject);
	LUT_byIndex.insert(LUT_byIndex.begin() + index, INT2(majorId, minorId));

	//currently, all indexes values stored in LUT_byId greater or equal to index need to be incremented by 1 since a new element was inserted at index position
	for (int idxMajor = 0; idxMajor < (int)LUT_byId.size(); idxMajor++) {
		for (int idxMinor = 0; idxMinor < (int)LUT_byId[idxMajor].size(); idxMinor++) {

			if (LUT_byId[idxMajor][idxMinor] >= index) LUT_byId[idxMajor][idxMinor]++;
		}
	}

	//store index of pushed object in LUT
	LUT_byId[majorId][minorId] = index;

	return minorId;
}

template <typename Type>
void vector_lut<Type>::set(Type newObject, INT2 id)
{
	//id must be positive
	if (id.major < 0 || id.minor < 0) return;

	if (is_id_set(id)) {

		//overwrite existing entry
		storedObjects[LUT_byId[id.i][id.j]] = newObject;
		return;
	}

	//entry does not exist: create it

	int majorId = id.major;

	while (true) {

		//create new entry for this majorId. Since entry doesn't exist then minorId <= id.minor.
		int minorId = setup_LUT_byID(majorId);

		//work up to minorId = id.minor by adding empty elements
		if (minorId != id.minor) storedObjects.push_back(Type());
		else storedObjects.push_back(newObject);
		
		LUT_byIndex.push_back(INT2(majorId, minorId));

		//store index of pushed object in LUT
		LUT_byId[majorId][minorId] = (int)storedObjects.size() - 1;

		//finish after adding element with required id
		if (minorId == id.minor) break;
	}
}

template <typename Type>
void vector_lut<Type>::erase(int index)
{
	if (index < 0 || index >= (int)storedObjects.size()) return;

	//mark in LUT the entry to be erased as not available any longer
	INT2 eraseId = LUT_byIndex[index];
	LUT_byId[eraseId.i][eraseId.j] = -1;

	clean_LUT_byId(eraseId.i);

	//erase the entry
	storedObjects.erase(storedObjects.begin() + index);
	LUT_byIndex.erase(LUT_byIndex.begin() + index);

	//all entries from index up must now also be updated in LUT_byId : must decrement entry there by 1 as those indexes have decreased by 1 now
	for (int i = index; i < storedObjects.size(); i++) {

		INT2 id = LUT_byIndex[i];
		LUT_byId[id.i][id.j]--;
	}
}

template <typename Type>
void vector_lut<Type>::clear(void)
{
	storedObjects.clear();
	LUT_byIndex.clear();
	LUT_byId.clear();
}

template <typename Type>
void vector_lut<Type>::shrink_to_fit(void)
{
	storedObjects.shrink_to_fit();
	LUT_byIndex.shrink_to_fit();
	LUT_byId.shrink_to_fit();
}

//------------------------------------ CONTROL OF ENTRIES

//move element from source to destination index
template <typename Type>
void vector_lut<Type>::move(int srcIdx, int dstIdx)
{
	if (srcIdx < 0 || dstIdx < 0 || srcIdx >= (int)storedObjects.size() || dstIdx >= (int)storedObjects.size() || srcIdx == dstIdx) return;

	//id of object to move
	INT2 idSrc = LUT_byIndex[srcIdx];
	//update LUT_byId of this object to its new destination : dstIdx
	LUT_byId[idSrc.i][idSrc.j] = dstIdx;

	if (dstIdx < srcIdx) {

		storedObjects.insert(storedObjects.begin() + dstIdx, storedObjects[srcIdx]);
		LUT_byIndex.insert(LUT_byIndex.begin() + dstIdx, idSrc);

		storedObjects.erase(storedObjects.begin() + srcIdx + 1);
		LUT_byIndex.erase(LUT_byIndex.begin() + srcIdx + 1);

		for (int i = dstIdx + 1; i <= srcIdx; i++) {

			INT2 idShift = LUT_byIndex[i];
			LUT_byId[idShift.i][idShift.j]++;
		}
	}
	else {

		storedObjects.insert(storedObjects.begin() + dstIdx + 1, storedObjects[srcIdx]);
		LUT_byIndex.insert(LUT_byIndex.begin() + dstIdx + 1, idSrc);

		storedObjects.erase(storedObjects.begin() + srcIdx);
		LUT_byIndex.erase(LUT_byIndex.begin() + srcIdx);

		for (int i = dstIdx - 1; i >= srcIdx; i--) {

			INT2 idShift = LUT_byIndex[i];
			LUT_byId[idShift.i][idShift.j]--;
		}
	}
}

//------------------------------------ MULTIPLE ENTRIES WITH SAME MAJOR ID

//get number of elements with given major ID
template <typename Type>
int vector_lut<Type>::get_num_IDs(int majorID) const
{
	if (majorID < (int)LUT_byId.size()) {

		int count = 0;

		for (int idx = 0; idx < LUT_byId[majorID].size(); idx++) {

			if (LUT_byId[majorID][idx] >= 0) count++;
		}

		return count;
	}
	else return 0;
}

//get a minor id for element with majorID; return -1 if none found
template <typename Type>
int vector_lut<Type>::get_next_minor_id(int majorID) const
{
	if (get_num_IDs(majorID) > 0) {

		for (int idx = 0; idx < LUT_byId[majorID].size(); idx++) {

			if (LUT_byId[majorID][idx] >= 0) return idx;
		}
	}
	else return -1;
}
//make all minor ids for elements with given major ID start from 0 and be numbered in steps of 1
//use this e.g. when you want to ensure there's always a minor if of 0 if major id exists
//DO NOT use if you are keeping track of minor ids externally
template <typename Type>
void vector_lut<Type>::regularise_minor_ids(int majorID)
{
	if (get_num_IDs(majorID) > 0) {

		int idx = 0;

		while (idx < LUT_byId[majorID].size()) {

			if (LUT_byId[majorID][idx] < 0) {

				//minor id not used -> erase it
				LUT_byId[majorID].erase(LUT_byId[majorID].begin() + idx);

				//do not increment idx
				continue;
			}

			idx++;
		}

		//now also adjust LUT_byIndex : minor ids have changed for all elements with majorID, so adjust them in LUT_byIndex
		for (idx = 0; idx < LUT_byId[majorID].size(); idx++) {

			//index in stored objects, hence also in LUT_byIndex
			int index = LUT_byId[majorID][idx];
			
			if (index >= 0 && index < LUT_byIndex.size())
				LUT_byIndex[index] = INT2(majorID, idx);
		}
	}
}

//------------------------------------ ID SETTERS

//over-ride existing Id with new one. No checking is done for clashing with existing Ids.
template <typename Type>
void vector_lut<Type>::set_id_at_index(INT2 Id, int index)
{
	if (!GoodIdx(size() - 1, index)) return;

	INT2 currId = LUT_byIndex[index];	//get currently set Id

	LUT_byIndex[index] = Id;			//set new Id in LUT index table

										//unassign old Id
	if (currId >= INT2(0) && currId.major < LUT_byId.size() && currId.minor < LUT_byId[currId.major].size())
		LUT_byId[currId.i][currId.j] = -1;

	setup_LUT_byId(Id);
	LUT_byId[Id.i][Id.j] = index;		//assign new Id to given index
}

////////////////////////////////////////////////////////////////////////
//
// vector_key : in addition to normal std::vector indexing this allows indexing by a text handle - similar to a map, but find this nicer to use.

template <typename Type> 
class vector_key {

private:

	//the actual std::vector and keys - don't use a pair to store them in a single std::vector as it complicates using when iterators to read elements sequentially
	std::vector<Type> storedObjects;
	std::vector<std::string> storedObjects_keys;

public:

	//------------------------------------ CONSTRUCTORS

	vector_key(void) {}
	~vector_key() {}

	//copy constructor for lvalues
	vector_key(const vector_key<Type>& copyThis) :
		storedObjects(copyThis.storedObjects),
		storedObjects_keys(copyThis.storedObjects_keys)
	{}

	//overloaded copy constructor for rvalues using move semantics
	vector_key(vector_key<Type>&& moveThis) :
		storedObjects(std::move(moveThis.storedObjects)),
		storedObjects_keys(std::move(moveThis.storedObjects_keys))
	{}

	//assignment operator for lvalues
	vector_key<Type>& operator=(const vector_key<Type>& copyThis)
	{
		storedObjects = copyThis.storedObjects;
		storedObjects_keys = copyThis.storedObjects_keys;

		return *this;
	}

	//overloaded assignment operator for rvalues using move semantics
	vector_key<Type>& operator=(vector_key<Type>&& moveThis)
	{
		storedObjects = std::move(moveThis.storedObjects);
		storedObjects_keys = std::move(moveThis.storedObjects_keys);

		return *this;
	}

	//------------------------------------ INDEXING

	//usual std::vector indexing
	Type& operator[](const int &index) { return storedObjects[index]; }
	const Type& operator[](const int &index) const { return storedObjects[index]; }

	Type& operator[](const std::string &key);
	const Type& operator[](const std::string &key) const;

	Type& front(void) { return storedObjects[0]; }
	Type& back(void) { return storedObjects[last()]; }

	//------------------------------------ ITERATORS

	Type* begin(void) { return &storedObjects[0]; }
	Type* data(void) { return &storedObjects[0]; }
	Type* end(void) { return &storedObjects[size()]; }

	//------------------------------------ CONTROL OF ENTRIES

	//move element from source to destination index (i.e. ordering in std::vector changes)
	void move(int srcIdx, int dstIdx = 0)
	{
		iter_swap(storedObjects.begin() + srcIdx, storedObjects.begin() + dstIdx);
		iter_swap(storedObjects_keys.begin() + srcIdx, storedObjects_keys.begin() + dstIdx);
	}

	//------------------------------------ SIZE MODIFICATION METHODS

	void push_back(Type newObject, const std::string& key);

	void pop_back(void);

	//make new entry at given index (if it exceeds size then resize, otherwise replace entry at the index)
	void insert(int index, Type newObject, const std::string& key);

	void erase(int index);

	void erase(const std::string& key);

	void clear(void);

	void shrink_to_fit(void);

	void resize(int new_size, Type newObject = Type());

	//------------------------------------ GET PROPERTIES

	//dimension of internal std::vector
	int size(void) const { return (int)storedObjects.size(); }

	//index of last element
	int last(void) const { if (size()) return ((int)size() - 1); else return 0; }

	//------------------------------------ CHECKERS

	bool has_key(const std::string& key) const;

	//------------------------------------ KEY <--> INDEX GETTERS

	int index_from_key(const std::string& key) const;

	std::string get_key_from_index(int idx) const { if (idx >= 0 && idx < (int)storedObjects.size()) return storedObjects_keys[idx]; else return ""; }

	//------------------------------------ KEY MODIFICATION

	//change key of entry in std::vector to new key - return true if succesful (i.e. old key exists)
	bool change_key(const std::string& old_key, const std::string& new_key);

	void set_key_at_index(const std::string& new_key, int index) { storedObjects_keys[index] = new_key; }

	//------------------------------------ ADVANCED PROPERTIES

	//find matches for a partial key start : return std::vector of indexes with keys containing keystart at the beginning
	std::vector<int> find_keystart(const std::string& keystart) const;
};

//------------------------------------ INDEXING

template <typename Type>
Type& vector_key<Type>::operator[](const std::string &key)
{
	for (int idx = 0; idx < (int)storedObjects_keys.size(); idx++) {

		if (storedObjects_keys[idx] == key) return storedObjects[idx];
	}

	return storedObjects[0];
}

template <typename Type>
const Type& vector_key<Type>::operator[](const std::string &key) const
{
	for (int idx = 0; idx < (int)storedObjects_keys.size(); idx++) {

		if (storedObjects_keys[idx] == key) return storedObjects[idx];
	}

	return storedObjects[0];
}


//------------------------------------ SIZE MODIFICATION METHODS

template <typename Type>
void vector_key<Type>::push_back(Type newObject, const std::string& key)
{
	storedObjects.push_back(newObject);
	storedObjects_keys.push_back(key);
}

template <typename Type>
void vector_key<Type>::pop_back(void)
{
	storedObjects.pop_back();
	storedObjects_keys.pop_back();
}

//make new entry at given index (if it exceeds size then resize, otherwise replace entry at the index)
template <typename Type>
void vector_key<Type>::insert(int index, Type newObject, const std::string& key)
{
	//resize if needed
	if (index >= (int)storedObjects.size()) {

		storedObjects.resize(index + 1);
		storedObjects_keys.resize(index + 1);
	}

	storedObjects[index] = newObject;
	storedObjects_keys[index] = key;
}

template <typename Type>
void vector_key<Type>::erase(int index)
{
	storedObjects.erase(storedObjects.begin() + index);
	storedObjects_keys.erase(storedObjects_keys.begin() + index);
}

template <typename Type>
void vector_key<Type>::erase(const std::string& key)
{
	int index = index_from_key(key);
	if (index >= 0) erase(index);
}

template <typename Type>
void vector_key<Type>::clear(void)
{
	storedObjects.clear();
	storedObjects_keys.clear();
}

template <typename Type>
void vector_key<Type>::shrink_to_fit(void)
{
	storedObjects.shrink_to_fit();
	storedObjects_keys.shrink_to_fit();
}

template <typename Type>
void vector_key<Type>::resize(int new_size, Type newObject)
{
	storedObjects.resize(new_size, newObject);
	storedObjects_keys.resize(new_size, "");
}


//------------------------------------ CHECKERS

template <typename Type>
bool vector_key<Type>::has_key(const std::string& key) const
{
	for (int idx = 0; idx < (int)storedObjects_keys.size(); idx++) {

		if (storedObjects_keys[idx] == key) return true;
	}

	return false;
}

//------------------------------------ KEY <--> INDEX GETTERS

template <typename Type>
int vector_key<Type>::index_from_key(const std::string& key) const
{
	for (int idx = 0; idx < (int)storedObjects_keys.size(); idx++) {

		if (storedObjects_keys[idx] == key) return idx;
	}

	return -1;
}

//------------------------------------ KEY MODIFICATION

//change key of entry in std::vector to new key - return true if succesful (i.e. old key exists)
template <typename Type>
bool vector_key<Type>::change_key(const std::string& old_key, const std::string& new_key)
{
	for (int idx = 0; idx < (int)storedObjects_keys.size(); idx++) {

		if (storedObjects_keys[idx] == old_key) {

			storedObjects_keys[idx] = new_key;
			return true;
		}
	}

	return false;
}

//------------------------------------ ADVANCED PROPERTIES

//find matches for a partial key start : return std::vector of indexes with keys containing keystart at the beginning
template <typename Type>
std::vector<int> vector_key<Type>::find_keystart(const std::string& keystart) const
{
	std::vector<int> indexes;

	for (int idx = 0; idx < (int)storedObjects_keys.size(); idx++) {

		std::string key = storedObjects_keys[idx];

		if (keystart.length() <= key.length() && key.substr(0, keystart.length()) == keystart)
			indexes.push_back(idx);
	}

	return indexes;
}

///////////////////////////////////////////////////////////////////////
//
// vector_key_lut : vector_key_lut can be indexed as a vector_lut, and in addition can be indexed by the std::string key.

template <typename Type> 
class vector_key_lut : 
	public vector_lut<Type> 
	//most of the code needed is in vector_lut. Don't need to also inherit from public_key as that complicates handling of storedObjects.
{		

private:

	std::vector<std::string> storedObjects_keys;

public:

	//------------------------------------ CONSTRUCTORS

	vector_key_lut(void) : vector_lut<Type>() {}

	//copy constructor for lvalues
	vector_key_lut(const vector_key_lut<Type>& copyThis) 
		: vector_lut<Type>(copyThis), storedObjects_keys(copyThis.storedObjects_keys)
	{}

	//overloaded copy constructor for rvalues using move semantics
	vector_key_lut(vector_key_lut<Type>&& moveThis)
		: vector_lut<Type>(std::move(moveThis)), storedObjects_keys(std::move(moveThis.storedObjects_keys))
	{}

	//assignment operator for lvalues
	vector_key_lut<Type>& operator=(const vector_key_lut<Type>& copyThis) 
	{ 
		vector_lut<Type>::storedObjects = copyThis.storedObjects;
		vector_lut<Type>::LUT_byIndex = copyThis.LUT_byIndex;
		vector_lut<Type>::LUT_byId = copyThis.LUT_byId;
		storedObjects_keys = copyThis.storedObjects_keys;

		return *this;
	}

	//overloaded assignment operator for rvalues using move semantics
	vector_key_lut<Type>& operator=(vector_key_lut<Type>&& moveThis)
	{ 
		vector_lut<Type>::storedObjects = std::move(moveThis.storedObjects);
		vector_lut<Type>::LUT_byIndex = std::move(moveThis.LUT_byIndex);
		vector_lut<Type>::LUT_byId = std::move(moveThis.LUT_byId);
		storedObjects_keys = std::move(moveThis.storedObjects_keys);

		return *this; 
	}

	//------------------------------------ INDEXING (identical to vector_lut operators - the base class).
	//NOTE: since we have an [] operator in this derived class (the one with const std::string& parameter below), then no [] base operators are inherited, so must redefine them here.

	//usual std::vector indexing
	Type& operator[](const int &index) { return vector_lut<Type>::storedObjects[index]; }
	const Type& operator[](const int &index) const { return vector_lut<Type>::storedObjects[index]; }

	//indexing via look-up table with INT2(majorId, minorId). Note, can also call it with INT2(majorId) : this means minorId = 0, and is used when majorId is known to be unique.
	Type& operator[](const INT2 &id) { return vector_lut<Type>::storedObjects[vector_lut<Type>::LUT_byId[id.i][id.j]]; }
	const Type& operator[](const INT2 &id) const { return vector_lut<Type>::storedObjects[vector_lut<Type>::LUT_byId[id.i][id.j]]; }

	//------------------------------------ INDEXING with key

	//indexing via key std::string
	Type& operator[](const std::string &key);
	const Type& operator[](const std::string &key) const;

	//------------------------------------ CONTROL OF ENTRIES

	//move element from source to destination index (i.e. ordering in std::vector changes)
	void move(int srcIdx, int dstIdx = 0)
	{
		vector_lut<Type>::move(srcIdx, dstIdx);
		iter_swap(storedObjects_keys.begin() + srcIdx, storedObjects_keys.begin() + dstIdx);
	}

	//------------------------------------ CHECKERS
	
	bool has_key(const std::string& key) const;
	
	//------------------------------------ ID <--> KEY <--> INDEX GETTERS (all involving keys, others are inherited from vector_lut)

	//return index of first element with key std::string. -1 if not found
	int get_index_from_key(const std::string& key) const;

	//get major ID of first element with key std::string. -1 if not found
	int get_ID_from_key(const std::string& key) const;

	//get full id of first element with key std::string. -1 if not found
	INT2 get_id_from_key(const std::string& key) const;

	//get key from index in std::vector
	std::string get_key_from_index(int index) const;

	//get key from major ID
	std::string get_key_from_ID(int majorID) const;

	//------------------------------------ KEY MODIFICATION

	//change key of entry in std::vector to new key - return true if succesful (i.e. old key exists)
	bool change_key(const std::string& old_key, const std::string& new_key);

	void set_key_at_index(const std::string& new_key, int index) { storedObjects_keys[index] = new_key; }

	//------------------------------------ SIZE MODIFICATION METHODS (overloaded vector_lut methods)

	int push_back(const std::string&, Type newObject, int majorId = 0);

	void pop_back(void);

	int insert(const std::string&, int index, Type newObject, int majorId = 0);

	//set new size. Any new entries beyond the old size will be assigned a majorId of 0 and empty key.
	void resize(int new_size, Type value = Type());

	void erase(INT2 id);
	void erase(int index);

	void clear(void);

	void shrink_to_fit(void);
};

//------------------------------------ INDEXING with key

//indexing via key std::string
template <typename Type>
Type& vector_key_lut<Type>::operator[](const std::string &key)
{
	for (int idx = 0; idx < (int)storedObjects_keys.size(); idx++)
		if (storedObjects_keys[idx] == key) return vector_lut<Type>::storedObjects[idx];

	return vector_lut<Type>::storedObjects[0];
}

template <typename Type>
const Type& vector_key_lut<Type>::operator[](const std::string &key) const
{
	for (int idx = 0; idx < (int)storedObjects_keys.size(); idx++)
		if (storedObjects_keys[idx] == key) return vector_lut<Type>::storedObjects[idx];

	return vector_lut<Type>::storedObjects[0];
}

//------------------------------------ CHECKERS

template <typename Type>
bool vector_key_lut<Type>::has_key(const std::string& key) const
{
	for (int idx = 0; idx < (int)storedObjects_keys.size(); idx++)
		if (storedObjects_keys[idx] == key) return true;

	return false;
}

//------------------------------------ ID <--> KEY <--> INDEX GETTERS (all involving keys, others are inherited from vector_lut)

//return index of first element with key std::string. -1 if not found
template <typename Type>
int vector_key_lut<Type>::get_index_from_key(const std::string& key) const
{
	for (int idx = 0; idx < (int)storedObjects_keys.size(); idx++)
		if (storedObjects_keys[idx] == key) return idx;

	return -1;
}

//get major ID of first element with key std::string. -1 if not found
template <typename Type>
int vector_key_lut<Type>::get_ID_from_key(const std::string& key) const
{
	for (int idx = 0; idx < (int)storedObjects_keys.size(); idx++)
		if (storedObjects_keys[idx] == key) return vector_lut<Type>::get_ID_from_index(idx);

	return -1;
}

//get full id of first element with key std::string. -1 if not found
template <typename Type>
INT2 vector_key_lut<Type>::get_id_from_key(const std::string& key) const
{
	for (int idx = 0; idx < (int)storedObjects_keys.size(); idx++)
		if (storedObjects_keys[idx] == key) return vector_lut<Type>::get_id_from_index(idx);

	return INT2(-1);
}

//get key from index in std::vector
template <typename Type>
std::string vector_key_lut<Type>::get_key_from_index(int index) const
{
	if (GoodIdx(vector_lut<Type>::size(), index)) return storedObjects_keys[index];

	return "";
}

//get key from major ID
template <typename Type>
std::string vector_key_lut<Type>::get_key_from_ID(int majorID) const
{
	if (vector_lut<Type>::is_ID_set(majorID)) return storedObjects_keys[vector_lut<Type>::LUT_byId[majorID][0]];

	return "";
}

//------------------------------------ KEY MODIFICATION

//change key of entry in std::vector to new key - return true if succesful (i.e. old key exists)
template <typename Type>
bool vector_key_lut<Type>::change_key(const std::string& old_key, const std::string& new_key)
{
	for (int idx = 0; idx < (int)storedObjects_keys.size(); idx++) {

		if (storedObjects_keys[idx] == old_key) {

			storedObjects_keys[idx] = new_key;
			return true;
		}
	}

	return false;
}


//------------------------------------ SIZE MODIFICATION METHODS (overloaded vector_lut methods)

template <typename Type>
int vector_key_lut<Type>::push_back(const std::string& key, Type newObject, int majorId)
{
	//store the key (last element)
	storedObjects_keys.push_back(key);

	//use base push_back method to store the value (last element) - return the minorId as well
	return vector_lut<Type>::push_back(newObject, majorId);
}

template <typename Type>
void vector_key_lut<Type>::pop_back(void)
{
	storedObjects_keys.pop_back();
	vector_lut<Type>::pop_back();
}

template <typename Type>
int vector_key_lut<Type>::insert(const std::string& key, int index, Type newObject, int majorId)
{
	//store the key
	storedObjects_keys.insert(storedObjects_keys.begin() + index, key);
	//store the value
	return vector_lut<Type>::insert(index, newObject, majorId);
}

//set new size. Any new entries beyond the old size will be assigned a majorId of 0 and empty key.
template <typename Type>
void vector_key_lut<Type>::resize(int new_size, Type value)
{
	storedObjects_keys.resize(new_size, "");
	vector_lut<Type>::resize(new_size, value);
}

template <typename Type>
void vector_key_lut<Type>::erase(INT2 id)
{
	int index = vector_lut<Type>::LUT_byId[id.i][id.j];
	erase(index);
}

template <typename Type>
void vector_key_lut<Type>::erase(int index)
{
	storedObjects_keys.erase(storedObjects_keys.begin() + index);
	vector_lut<Type>::erase(index);
}

template <typename Type>
void vector_key_lut<Type>::clear(void)
{
	storedObjects_keys.clear();
	vector_lut<Type>::clear();
}

template <typename Type>
void vector_key_lut<Type>::shrink_to_fit(void)
{
	storedObjects_keys.shrink_to_fit();
	vector_lut<Type>::shrink_to_fit();
}
