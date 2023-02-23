#pragma once

#include "BorisLib.h"

#include "MeshParams.h"

#include "MeshDefs.h"



//Manages materials data bases, allowing users to save/load materials definitions

//mdb files contain the following entries:

//0. Name - this is the material name, and it is the handle used with addmaterial command
//1. Formula - symbol equation for material (information only)
//2. Type - the type of mesh this material applies to (e.g. ferromagnetic, metal, insulator, ...)
//3. Description - brief description of entry, e.g. could indicate if this is applicable to a particular type of interface, such as Pt/Co, Ta/Co, or just bulk Co etc., or other relevant information so users know how they should use this entry
//4. Contributor - name of contributor for this entry (leave blank for Anonymous)
//5. State - is this a local, or an entry available in the official online data base. SHARED / LOCAL.
//
//Now follow all material parameters values. 
//The columns are labelled using the handles in MeshParams::meshParams, i.e. same handles visible in the console when using the params command.
//If the parameter has a unit then it is specified in brackets.
//e.g. a column label would look like: "Ms (A/m)", or "damping", etc. - without the quotation marks
//
//6a. Parameter - the parameter values are entered numerically in the units specified by the column descriptor
//6b. DOI - parameter values should have a valid literature reference; this is specified using a DOI, e.g. https://doi.org/10.1038/nphys3304
//
//7a. Parameter - //
//7b. DOI - //
//
//etc. for all parameters

//Fields character limits
//Enforced when loading entries from file, thus all sync requests also carry this limit
//If you edit an entry externally and go over the character limit, it will be truncated
#define FIELD_CHAR_LIMIT	200

class MaterialsDB {

	std::string default_databaseName = "BorisMDB.txt";

	//name of the currently selected database : this is a file path together with file name
	//the default setting uses the user documents path + "Boris Data" + "BorisMDB.txt"
	std::string databaseName_withpath;

	//all entries in the currently selected database
	//first row has column labels
	//parameter values start at params_start_idx (see below). Two columns each, one for the value, the other for the DOI
	std::vector<std::vector<std::string>> entries;

	//latest parameters loaded from data base - configured in the constructor to indicate enabled mesh parameters
	MeshParams dbEntry;

	//see entries above in database. This is the index (counting from 0) of first parameter value
	int params_start_idx = 6;

	//mesh type for which the dbEntry mesh parameter values apply
	//if set to MESH_SUPERMESH then nothing loaded
	MESH_ dbEntryType;

	//for all entries in MESH_ enum associate a handle which will be used with the Type field in data base entries
	vector_lut<std::string> meshTypeHandles;

private:

	//make a fresh mdb file at dbname, which must include a path, erasing any existing file
	//also loads entries
	BError make_mdb_file(std::string newdatabaseName);

	//make sure the mdb file has the correct entries for parameters as obtained from dbEntry
	//also loads entries
	BError fix_mdb_file(void);

	//load entries from currently selected database
	BError load_entries(void);

	//store entries in currently selected database
	BError store_entries(void);

	//find row index in entries for given materialName
	int find_entry_idx(std::string materialName);

	//replace certain characters with an ascii signalling code before sending the std::string with POST / GET requests
	void encode_characters_ascii(std::string& message);

public:

	MaterialsDB(std::vector<PARAM_>& enabledParams);

	//handles the materialsdatabase command - switches the currently selected database
	BError SwitchDataBase(std::string newdatabaseName);

	//get the currently set data base name
	std::string GetDataBaseName(void) { return databaseName_withpath; }

	//attempt to load entry from current database for material with given name
	//set in *pmeshType the type of material loaded, e.g. ferromagnetic, metal, insulator, etc. - see enum MESH_ in MeshDefs.h for possible values
	//pmeshType is nullptr then ignore this parameter
	//If succesful return true so the caller can create a mesh of the required type
	//The actual mesh parameter values will then be available here in dbEntry
	BError LoadMaterial(std::string materialName, int* pmeshType = nullptr);

	//copy parameters from current dbEntry to copy_to_this - e.g. a newly created mesh with addmaterial
	//can also be used with the setmaterial command, which overwrites mesh parameter values in an already created mesh
	void copy_parameters(MeshParams& copy_to_this);

	//add new entry in currently selected database - use with addmdbentry command
	BError AddMDBEntry(std::string materialName, MeshParams& meshParamsValues, MESH_ meshType);

	//delete entry in currently selected database
	BError DelMDBEntry(std::string materialName);

	//reload the currently selected database
	BError RefreshMDB(void);

	//From currently selected data base, request the named material is added to the main online shared materials data base. Handles the requestmdbsync command.
	//If successfully sent, the entry will be checked by a moderator, and if correct will become available for all users
	//If successfully sent, returnMessage is empty, otherwise it contains a message detailing why the sync request couldn't be sent (e.g. incorrectly formatted entry / missing parameters etc. - too specific to set this as a BError entry)
	//BError contains general info, e.g. couldn't connect etc. (or no error code if everything fine)
	BError RequestMDBSync(std::string materialName, std::string domain_name, std::string mdb_entry_handler, std::string emailContact, std::string *presponseMessage);

	//check if there are any new online materials data base entries - should only be called once at the start of the program to limit the potential for spamming the server
	//the server should also have a built-in timer to reject requests from users if sent too often from same ip
	BError CheckMDBLastUpdate(std::string *presponseMessage);

	//Update BorisMDB.txt with any new entries available in the online version. Handles the updatemdb command. 
	BError UpdateMDB(std::string domain_name, std::string mdb_handler, std::string* presponseMessage);
};
