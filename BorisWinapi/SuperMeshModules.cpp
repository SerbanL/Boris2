#include "stdafx.h"
#include "SuperMesh.h"

//----------------------------------- MODULES CONTROL

BError SuperMesh::AddModule(string meshName, MOD_ moduleId)
{
	BError error(__FUNCTION__);

	if (moduleId <= MOD_ERROR) return error(BERROR_INCORRECTNAME);

	if (!contains(meshName) && meshName != superMeshHandle) return error(BERROR_INCORRECTNAME);

	//add module to given mesh
	if (meshName != superMeshHandle) {

		error = pMesh[meshName]->AddModule(moduleId);
	}
	//add module to supermesh. If module is already set then don't add another one : only one module of each type allowed
	else if (!IsSuperMeshModuleSet(moduleId)) {

		//add module to super-mesh
		switch (moduleId) {

		case MODS_SDEMAG:
			pSMod.push_back(new SDemag(this), MODS_SDEMAG);
			break;

		case MODS_STRAYFIELD:
			pSMod.push_back(new StrayField(this), MODS_STRAYFIELD);
			break;

		case MODS_STRANSPORT:
			pSMod.push_back(new STransport(this), MODS_STRANSPORT);
			break;

		case MODS_OERSTED:
			pSMod.push_back(new Oersted(this), MODS_OERSTED);
			break;

		case MODS_SHEAT:
			pSMod.push_back(new SHeat(this), MODS_SHEAT);
			break;
		}

		//check the super-mesh module was created correctly - if not, delete it
		error = pSMod.back()->Error_On_Create();

		if (error) pSMod.pop_back();
	}
	else return error(BERROR_INCORRECTACTION_SILENT);

	if (!error) {

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
		for (int idx = 0; idx < (int)superMeshCompanionModules[moduleId].size(); idx++) {

			MOD_ module = superMeshCompanionModules[moduleId][idx];

			if (meshName != superMeshHandle) {

				if (!error) error = AddModule(superMeshHandle, module);
			}
		}

		error = UpdateConfiguration();
	}

	return error;
}

BError SuperMesh::DelModule(string meshName, MOD_ moduleId)
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

	error = UpdateConfiguration();

	return error;
}