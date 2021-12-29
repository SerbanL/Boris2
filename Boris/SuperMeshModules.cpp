#include "stdafx.h"
#include "SuperMesh.h"

//----------------------------------- MODULES CONTROL

BError SuperMesh::AddModule(std::string meshName, MOD_ moduleId)
{
	BError error(__FUNCTION__);
	
	if (moduleId <= MOD_ERROR) return error(BERROR_INCORRECTNAME);

	//can either a module to a specific mesh, or to all applicable meshes if meshName is empty, or a supermesh module
	if (!contains(meshName) && meshName != superMeshHandle && meshName.length()) return error(BERROR_INCORRECTNAME);

	auto adjust_mesh_modules = [&](std::string meshName, BError& error) -> BError {
		
		//check any super-mesh module was created correctly - if not, delete it
		if (pSMod.size() && meshName == superMeshHandle) {

			error = pSMod.back()->Error_On_Create();
			if (error) {

				pSMod.pop_back();
				return error;
			}
		}
		
		//now check list of exclusive super-mesh modules : any that are exclusive to the module just added must be removed from the list of active modules
		for (int idx = 0; idx < (int)superMeshExclusiveModules[moduleId].size(); idx++) {

			MOD_ module = superMeshExclusiveModules[moduleId][idx];

			if (meshName == superMeshHandle) {

				for (int idxMesh = 0; idxMesh < (int)pMesh.size(); idxMesh++) {

					pMesh[idxMesh]->DelModule(module);
				}
			}
			else if (IsSuperMeshModuleSet(module)) { delete pSMod[pSMod.get_index_from_ID(module)]; pSMod.erase(INT2(module, 0)); }
		}
		
		//check companion modules : if a module was added on a normal mesh, make sure any companion supermesh modules are also set
		if (meshName != superMeshHandle && meshName.length()) {

			for (int idx = 0; idx < (int)superMeshCompanionModules[moduleId].size(); idx++) {

				MOD_ module = superMeshCompanionModules[moduleId][idx];
				if (!error) error = AddModule(superMeshHandle, module);
			}
		}
		
		return error;
	};

	//add module to given mesh
	if (meshName != superMeshHandle && meshName.length()) {

		error = pMesh[meshName]->AddModule(moduleId);
		if (!error) error = adjust_mesh_modules(meshName, error);

		error = UpdateConfiguration(UPDATECONFIG_MODULEADDED);
	}
	//add module to supermesh. If module is already set then don't add another one : only one module of each type allowed on the supermesh
	else if (meshName == superMeshHandle && !IsSuperMeshModuleSet(moduleId)) {

		//add module to super-mesh
		switch (moduleId) {

		case MODS_SDEMAG:
			pSMod.push_back(new SDemag(this), MODS_SDEMAG);
			error = adjust_mesh_modules(superMeshHandle, error);
			break;

		case MODS_STRAYFIELD:
			pSMod.push_back(new StrayField(this), MODS_STRAYFIELD);
			error = adjust_mesh_modules(superMeshHandle, error);
			break;

		case MODS_STRANSPORT:
			pSMod.push_back(new STransport(this), MODS_STRANSPORT);
			error = adjust_mesh_modules(superMeshHandle, error);
			break;

		case MODS_OERSTED:
			pSMod.push_back(new Oersted(this), MODS_OERSTED);
			error = adjust_mesh_modules(superMeshHandle, error);
			break;

		case MODS_SHEAT:
			pSMod.push_back(new SHeat(this), MODS_SHEAT);
			error = adjust_mesh_modules(superMeshHandle, error);
			break;
		}

		error = UpdateConfiguration(UPDATECONFIG_MODULEADDED);
	}
	else if (!meshName.length()) {

		//try to add module to all meshes
		for (int idx = 0; idx < pMesh.size(); idx++) {

			error = pMesh[idx]->AddModule(moduleId);
			if (!error) error = adjust_mesh_modules(pMesh.get_key_from_index(idx), error);
			else error.reset();
		}

		error = UpdateConfiguration(UPDATECONFIG_MODULEADDED);
	}

	return error;
}

BError SuperMesh::DelModule(std::string meshName, MOD_ moduleId)
{
	BError error(__FUNCTION__);

	if (moduleId <= MOD_ERROR) return error(BERROR_INCORRECTNAME);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	if (meshName == superMeshHandle) {

		//delete module on super-mesh
		if (IsSuperMeshModuleSet(moduleId)) {

			if(pSMod(moduleId)) delete pSMod(moduleId);
			pSMod.erase(INT2(moduleId, 0));
		}
		
		//delete any companion modules on normal meshes
		for (int idx = 0; idx < (int)superMeshCompanionModules[moduleId].size(); idx++) {

			MOD_ module = superMeshCompanionModules[moduleId][idx];

			for (int idxMesh = 0; idxMesh < (int)pMesh.size(); idxMesh++) {

				pMesh[idxMesh]->DelModule(module);
			}
		}
	}
	else {

		//delete normal module
		pMesh[meshName]->DelModule(moduleId);

		//was this the last module of its type ?
		bool last_module = true;
		for (int idx = 0; idx < pMesh.size(); idx++) {

			if(pMesh[idx]->IsModuleSet(moduleId)) {
			
				last_module = false;
				break;
			}
		}

		if (last_module) {

			//check super mesh companion modules : if this was the last module of its type then delete the supermesh companion module too
			for (int idx = 0; idx < (int)superMeshCompanionModules[moduleId].size(); idx++) {

				MOD_ module = superMeshCompanionModules[moduleId][idx];

				//delete module on super-mesh
				if (IsSuperMeshModuleSet(module)) {

					if (pSMod(module)) delete pSMod(module);
					pSMod.erase(INT2(module, 0));
				}
			}
		}
	}

	error = UpdateConfiguration(UPDATECONFIG_MODULEDELETED);

	return error;
}