#include "stdafx.h"
#include "Simulation.h"

void Simulation::NewDataBoxField(DatumConfig dConfig) 
{
	//Can the DATA_ entry be displayed in the data box?
	if (dataDescriptor(dConfig.datumId).Label.length()) {

		int minorId = dataBoxList.push_back(dConfig);

		//Data box entry, showing the label of a given entry in Simulation::dataBoxList : minorId is the minor id of elements in Simulation::dataBoxList (major id there is always 0), auxId is the number of the interactive object in the list (i.e. entry number as it appears in data box in order). textId is the mesh name (if associated with this data type)
		string Label = "<b>[tc0,0,0,1/tc][io" + ToString((int)IOI_DATABOXFIELDLABEL) + "," + ToString(minorId) + "," + ToString(dataBoxList.last());

		//add mesh name if needed
		if (dConfig.meshName.length() && !dataDescriptor(dConfig.datumId).meshless) {

			Label += "," + dConfig.meshName + "/io]" + "<" + dConfig.meshName + "> ";
		}
		else Label += ", /io]";

		Label += dataDescriptor(dConfig.datumId).Label;

		//data box entries do not allow use of boxes on data - use whole mesh for these data
		dConfig.rectangle = Rect();

		string formattedTextEntry = Label + "</io></b>" + GetDataValueString(dConfig);

		BD.NewDataBoxField(formattedTextEntry);
	}
}

void Simulation::UpdateDataBox(void) 
{
	for (int idx = 0; idx < dataBoxList.size(); idx++) {

		BD.UpdateDataBoxField(idx, GetDataValueString(dataBoxList[idx]));
	}
}

void Simulation::RebuildDataBox(void) 
{
	BD.ClearDataBox();

	for (int idx = 0; idx < dataBoxList.size(); idx++) {

		DatumConfig dConfig = dataBoxList[idx];

		//Can the DATA_ entry be displayed in the data box?
		if (dataDescriptor(dConfig.datumId).Label.length()) {

			INT2 id = dataBoxList.get_id_from_index(idx);
			int minorId = id.minor;

			//Data box entry, showing the label of a given entry in Simulation::dataBoxList : minorId is the minor id of elements in Simulation::dataBoxList (major id there is always 0), auxId is the number of the interactive object in the list (i.e. entry number as it appears in data box in order). textId is the mesh name (if associated with this data type)
			string Label = "<b>[tc0,0,0,1/tc][io" + ToString((int)IOI_DATABOXFIELDLABEL) + "," + ToString(minorId) + "," + ToString(idx);

			//add mesh name if needed
			if (dConfig.meshName.length() && !dataDescriptor(dConfig.datumId).meshless) {

				Label += "," + dConfig.meshName + "/io]" + "<" + dConfig.meshName + "> ";
			}
			else Label += ", /io]";

			Label += dataDescriptor(dConfig.datumId).Label;

			//data box entries do not allow use of boxes on data - use whole mesh for these data
			dConfig.rectangle = Rect();

			string formattedTextEntry = Label + "</io></b>" + GetDataValueString(dConfig);

			BD.NewDataBoxField(formattedTextEntry);
		}
	}
}

void Simulation::DeleteDataBoxFields(string meshName) 
{
	for (int idx = dataBoxList.last(); idx >= 0; idx--) {

		if (dataBoxList[idx].meshName == meshName) dataBoxList.erase(idx);
	}
}

void Simulation::ChangeDataBoxLabels(string oldMeshName, string newMeshName) 
{
	for (int idx = 0; idx < dataBoxList.size(); idx++) {

		if (dataBoxList[idx].meshName == oldMeshName) dataBoxList[idx].meshName = newMeshName;
	}
}

Any Simulation::GetDataValue(DatumConfig dConfig) 
{
	//all possible DATA_ entries must be covered here

	switch (dConfig.datumId) {

	case DATA_STAGESTEP:
	{
		return Any(Check_and_GetStageStep());
	}
	break;

	case DATA_TIME:
	{
		return Any(SMesh.GetTime());
	}
	break;

	case DATA_STAGETIME:
	{
		return Any(SMesh.GetStageTime());
	}
	break;

	case DATA_ITERATIONS:
	{
		return Any(SMesh.GetIteration());
	}
	break;

	case DATA_SITERATIONS:
	{
		return Any(SMesh.GetStageIteration());
	}
	break;

	case DATA_MXH:
	{
		return Any(SMesh.Get_mxh());
	}
	break;

	case DATA_DMDT:
	{
		return Any(SMesh.Get_dmdt());
	}
	break;

	case DATA_DT:
	{
		return Any(SMesh.GetTimeStep());
	}
	break;

	case DATA_AVM:
	{
		return Any(SMesh[dConfig.meshName]->GetAverageMagnetisation(dConfig.rectangle));
	}
	break;

	case DATA_HA:
	{
		return Any(SMesh[dConfig.meshName]->CallModuleMethod(&Zeeman::GetField));
	}
	break;

	case DATA_JC:
	{
		return Any(SMesh[dConfig.meshName]->GetAverageChargeCurrentDensity(dConfig.rectangle));
	}
	break;

	case DATA_JSX:
	{
		return Any(SMesh[dConfig.meshName]->GetAverageSpinCurrentX(dConfig.rectangle));
	}
	break;

	case DATA_JSY:
	{
		return Any(SMesh[dConfig.meshName]->GetAverageSpinCurrentY(dConfig.rectangle));
	}
	break;

	case DATA_JSZ:
	{
		return Any(SMesh[dConfig.meshName]->GetAverageSpinCurrentZ(dConfig.rectangle));
	}
	break;

	case DATA_V:
	{
		return Any(SMesh[dConfig.meshName]->GetAverageElectricalPotential(dConfig.rectangle));
	}
	break;

	case DATA_S:
	{
		return Any(SMesh[dConfig.meshName]->GetAverageSpinAccumulation(dConfig.rectangle));
	}
	break;

	case DATA_ELC:
	{
		return Any(SMesh[dConfig.meshName]->GetAverageElectricalConductivity(dConfig.rectangle));
	}
	break;

	case DATA_POTENTIAL:
	{
		return Any(SMesh.CallModuleMethod(&STransport::GetPotential));
	}
	break;

	case DATA_CURRENT:
	{
		return Any(SMesh.CallModuleMethod(&STransport::GetCurrent));
	}
	break;

	case DATA_RESISTANCE:
	{
		return Any(SMesh.CallModuleMethod(&STransport::GetResistance));
	}
	break;

	case DATA_E_DEMAG:
	{
		if (!SMesh.IsSuperMeshModuleSet(MODS_SDEMAG)) {

			//read demag energy from named mesh if supermesh demag module not set
			return Any(SMesh[dConfig.meshName]->GetEnergy(MOD_DEMAG));
		}
		else {

			//get energy value from set supermesh demag module
			return SMesh.CallModuleMethod<double, SDemag>(&SDemag::GetEnergy);
		}
	}
	break;

	case DATA_E_EXCH:
	{
		if (SMesh[dConfig.meshName]->IsModuleSet(MOD_EXCHANGE6NGBR)) return Any(SMesh[dConfig.meshName]->GetEnergy(MOD_EXCHANGE6NGBR));
		else if (SMesh[dConfig.meshName]->IsModuleSet(MOD_DMEXCHANGE)) return Any(SMesh[dConfig.meshName]->GetEnergy(MOD_DMEXCHANGE));
		else return Any(SMesh[dConfig.meshName]->GetEnergy(MOD_IDMEXCHANGE));
	}
	break;

	case DATA_E_SURFEXCH:
	{
		return Any(SMesh[dConfig.meshName]->GetEnergy(MOD_SURFEXCHANGE));
	}
	break;

	case DATA_E_ZEE:
	{
		return Any(SMesh[dConfig.meshName]->GetEnergy(MOD_ZEEMAN));
	}
	break;

	case DATA_E_ANIS:
	{
		return Any(SMesh[dConfig.meshName]->GetEnergy(MOD_ANIUNI));
	}
	break;

	case DATA_E_ROUGH:
	{
		return Any(SMesh[dConfig.meshName]->GetEnergy(MOD_ROUGHNESS));
	}
	break;

	case DATA_E_TOTAL:
	{
		return Any(SMesh.GetTotalEnergy());
	}
	break;

	case DATA_DWSHIFT:
	{
		return Any(SMesh.Get_dwshift());
	}
	break;

	case DATA_SKYSHIFT:
	{
		return Any(SMesh[dConfig.meshName]->Get_skyshift(dConfig.rectangle));
	}
	break;

	case DATA_SKYPOS:
	{
		return Any(SMesh[dConfig.meshName]->Get_skypos_diameters(dConfig.rectangle));
	}
	break;

	case DATA_TRANSPORT_ITERSTOCONV:
	{
		return Any(SMesh.CallModuleMethod(&STransport::GetItersToConv));
	}
	break;

	case DATA_TRANSPORT_SITERSTOCONV:
	{
		return Any(SMesh.CallModuleMethod(&STransport::GetSItersToConv));
	}
	break;

	case DATA_TRANSPORT_CONVERROR:
	{
		return Any(SMesh.CallModuleMethod<double, STransport>(&STransport::GetEnergy));
	}
	break;

	case DATA_TRANSPORT_ASOR:
	{
		return Any(SMesh.CallModuleMethod<DBL2, STransport>(&STransport::Get_aSOR_Damping));
	}
	break;

	case DATA_TEMP:
	{
		return Any(SMesh[dConfig.meshName]->GetAverageTemperature(dConfig.rectangle));
	}
	break;

	case DATA_HEATDT:
	{
		return Any(SMesh.CallModuleMethod(&SHeat::get_heat_dT));
	}
	break;
	}

	return Any(0);
}

string Simulation::GetDataValueString(DatumConfig dConfig, bool ignore_unit) 
{
	string unit;
	if (!ignore_unit) unit = dataDescriptor(dConfig.datumId).unit;

	return GetDataValue(dConfig).convert_to_string(unit);
}

void Simulation::NewSaveDataEntry(DATA_ dataId, string meshName, Rect dataRect) 
{
	//if not meshless, make sure meshName is valid
	if (!dataDescriptor(dataId).meshless && !SMesh.contains(meshName)) return;

	if (dataDescriptor(dataId).meshless) meshName = "";

	//if not boxless and box is null, this actually means use the entire mesh - make this change (last check should not be really needed as you cannot have an entry which is meshless but not boxless)
	if (!dataDescriptor(dataId).boxless && dataRect.IsNull() && SMesh.contains(meshName)) { dataRect.e = SMesh[meshName]->GetMeshDimensions(); }

	if (dataDescriptor(dataId).boxless) dataRect = Rect();

	saveDataList.push_back(DatumConfig(dataId, meshName, dataRect));
}

void Simulation::EditSaveDataEntry(int index, DATA_ dataId, string meshName, Rect dataRect) 
{
	if (!GoodIdx(saveDataList.last(), index)) return;

	//if not meshless, make sure meshName is valid
	if (!dataDescriptor(dataId).meshless && !SMesh.contains(meshName)) return;

	if (dataDescriptor(dataId).meshless) meshName = "";

	//if not boxless and box is null, this actually means use the entire mesh - make this change (last check should not be really needed as you cannot have an entry which is meshless but not boxless)
	if (!dataDescriptor(dataId).boxless && dataRect.IsNull() && SMesh.contains(meshName)) { dataRect.e = SMesh[meshName]->GetMeshDimensions(); }

	if (dataDescriptor(dataId).boxless) dataRect = Rect();

	saveDataList[index] = DatumConfig(dataId, meshName, dataRect);
}

void Simulation::UpdateSaveDataEntries(string oldMeshName, string newMeshName) 
{
	for (int idx = 0; idx < saveDataList.size(); idx++) {

		if (saveDataList[idx].meshName == oldMeshName) saveDataList[idx].meshName = newMeshName;
	}
}

void Simulation::UpdateSaveDataEntries(Rect rect_old, Rect rect_new, string meshName) 
{
	for (int idx = 0; idx < saveDataList.size(); idx++) {

		if (!saveDataList[idx].rectangle.IsNull() && saveDataList[idx].meshName == meshName)
			saveDataList[idx].rectangle.resize(rect_old.size(), rect_new.size());
	}
}

void Simulation::DeleteSaveDataEntries(string meshName) 
{
	int idx = 0;

	while (true) {

		if (idx >= saveDataList.size()) break;

		if (saveDataList[idx].meshName == meshName) saveDataList.erase(idx);
		else idx++;
	}
}

void Simulation::SaveData(void) 
{
	//Actual data saving:

	if (saveDataFlag) {

		//we need a file name to save to
		if (savedataFile.size()) {

			ofstream bdout;

			//append to file or make a new one ?
			if (appendToDataFile) bdout.open((directory + savedataFile).c_str(), ios::out | ios::app);
			else {

				//Create new file
				bdout.open((directory + savedataFile).c_str(), ios::out);
				appendToDataFile = true;

				//Append header
				time_t rawtime;
				time(&rawtime);
				bdout << string(ctime(&rawtime)) << "\n";

				//List meshes
				for (int idx = 0; idx < SMesh.size(); idx++) {

					bdout << "<" + SMesh.key_from_meshIdx(idx) + "> : ";
					bdout << "Rectangle : " << ToString(SMesh[idx]->GetMeshRect(), "m") << ". ";
					bdout << "Cells : " << ToString(SMesh[idx]->GetMeshSize()) << ". ";
					bdout << "Cellsize : " << ToString(SMesh[idx]->GetMeshCellsize(), "m");
					bdout << "\n";
				}

				//Supermesh settings
				bdout << "<" + SMesh.superMeshHandle + "> : ";
				bdout << "Rectangle : " << ToString(SMesh.GetFMSMeshRect(), "m") << ". ";
				bdout << "Cells : " << ToString(SMesh.GetFMSMeshsize()) << ". ";
				bdout << "Cellsize : " << ToString(SMesh.GetFMSMeshCellsize(), "m");
				bdout << "\n";

				//List saved data labels and units
				bdout << "\nSaved data (dataname (unit) <meshname> (cells_rectangle)) : \n\n";

				for (int idx = 0; idx < saveDataList.size(); idx++) {

					//data name
					bdout << dataDescriptor.get_key_from_ID(saveDataList[idx].datumId);

					//unit?
					if (dataDescriptor(saveDataList[idx].datumId).unit.length())
						bdout << " (" << dataDescriptor(saveDataList[idx].datumId).unit << ")";

					//mesh?
					if (!dataDescriptor(saveDataList[idx].datumId).meshless)
						bdout << " <" << saveDataList[idx].meshName << ">";

					//box?
					if (!dataDescriptor(saveDataList[idx].datumId).boxless)
						bdout << " (" << ToString(saveDataList[idx].rectangle, "m") << ")";

					for (int tabs = 0; tabs < dataDescriptor(saveDataList[idx].datumId).components; tabs++)
						bdout << '\t';
				}

				bdout << "\n\n";
			}

			//save actual values as configured in saveDataList
			for (int idx = 0; idx < saveDataList.size(); idx++) {

				string value_string = GetDataValueString(saveDataList[idx], true);
				replaceall(value_string, ", ", "\t");

				bdout << value_string;
				if (idx != saveDataList.size() - 1) bdout << "\t";
			}

			bdout << "\n";
			bdout.close();
		}
	}

	//Image saving:

	if (saveImageFlag) {

		//update and refresh mesh display before saving image (the image is captured from on-screen image)
		UpdateMeshDisplay();
		RefreshScreen();

		string imageFile = directory + imageSaveFileBase + ToString(SMesh.GetIteration()) + ".png";
		BD.SaveMeshImage(imageFile);
	}
}