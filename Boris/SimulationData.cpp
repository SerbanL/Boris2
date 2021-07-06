#include "stdafx.h"
#include "Simulation.h"

void Simulation::NewDataBoxField(DatumConfig dConfig) 
{
	//Can the DATA_ entry be displayed in the data box?
	if (dataDescriptor(dConfig.datumId).Label.length()) {

		int minorId = dataBoxList.push_back(dConfig);

		//Data box entry, showing the label of a given entry in Simulation::dataBoxList : minorId is the minor id of elements in Simulation::dataBoxList (major id there is always 0), auxId is the number of the interactive object in the list (i.e. entry number as it appears in data box in order). textId is the mesh name (if associated with this data type)
		std::string Label = "<b>[tc0,0,0,1/tc][io" + ToString((int)IOI_DATABOXFIELDLABEL) + "," + ToString(minorId) + "," + ToString(dataBoxList.last());

		//add mesh name if needed
		if (dConfig.meshName.length() && !dataDescriptor(dConfig.datumId).meshless) {

			Label += "," + dConfig.meshName + "/io]" + "<" + dConfig.meshName + "> ";
		}
		else Label += ", /io]";

		Label += dataDescriptor(dConfig.datumId).Label;

		std::string formattedTextEntry = Label + "</io></b>" + GetDataValueString(dConfig);

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
			std::string Label = "<b>[tc0,0,0,1/tc][io" + ToString((int)IOI_DATABOXFIELDLABEL) + "," + ToString(minorId) + "," + ToString(idx);

			//add mesh name if needed
			if (dConfig.meshName.length() && !dataDescriptor(dConfig.datumId).meshless) {

				Label += "," + dConfig.meshName + "/io]" + "<" + dConfig.meshName + "> ";
			}
			else Label += ", /io]";

			Label += dataDescriptor(dConfig.datumId).Label;

			//data box entries do not allow use of boxes on data - use whole mesh for these data
			dConfig.rectangle = Rect();

			std::string formattedTextEntry = Label + "</io></b>" + GetDataValueString(dConfig);

			BD.NewDataBoxField(formattedTextEntry);
		}
	}
}

void Simulation::DeleteDataBoxFields(std::string meshName) 
{
	for (int idx = dataBoxList.last(); idx >= 0; idx--) {

		if (dataBoxList[idx].meshName == meshName) dataBoxList.erase(idx);
	}
}

void Simulation::ChangeDataBoxLabels(std::string oldMeshName, std::string newMeshName) 
{
	for (int idx = 0; idx < dataBoxList.size(); idx++) {

		if (dataBoxList[idx].meshName == oldMeshName) dataBoxList[idx].meshName = newMeshName;
	}
}

//update rectangles in dataBoxList (rect_old and rect_new are the old and new rectangles for the given mesh. All dataBoxList rectangles are scaled accordingly.
void Simulation::UpdateDataBoxEntries(Rect rect_old, Rect rect_new, std::string meshName)
{
	for (int idx = 0; idx < dataBoxList.size(); idx++) {

		if (!dataBoxList[idx].rectangle.IsNull() && dataBoxList[idx].meshName == meshName)
			dataBoxList[idx].rectangle.resize(rect_old.size(), rect_new.size());
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
		if (!SMesh[dConfig.meshName]->is_atomistic()) {

			return Any(dynamic_cast<Mesh*>(SMesh[dConfig.meshName])->GetAverageMagnetization(dConfig.rectangle));
		}
		else {

			//return average magnetization in the atomistic cell also : average atomistic moment (averaged over all participating unit cells), divided by unit cell volume
			return Any(dynamic_cast<Atom_Mesh*>(SMesh[dConfig.meshName])->GetAverageMoment(dConfig.rectangle) * MUB / SMesh[dConfig.meshName]->h.dim());
		}
	}
	break;

	case DATA_THAVM:
	{
		return Any((SMesh[dConfig.meshName])->GetThermodynamicAverageMagnetization(dConfig.rectangle));
	}
	break;

	case DATA_AVM2:
	{
		if (!SMesh[dConfig.meshName]->is_atomistic()) {

			return Any(dynamic_cast<Mesh*>(SMesh[dConfig.meshName])->GetAveragemagnetization2(dConfig.rectangle));
		}
		else return Any(DBL3());
	}
	break;

	case DATA_MX_MINMAX:
	{
		if (!SMesh[dConfig.meshName]->is_atomistic()) {

			return Any(dynamic_cast<Mesh*>(SMesh[dConfig.meshName])->GetmagnetizationXMinMax(dConfig.rectangle));
		}
		else {

			return Any(dynamic_cast<Atom_Mesh*>(SMesh[dConfig.meshName])->GetMomentXMinMax(dConfig.rectangle) * MUB / SMesh[dConfig.meshName]->h.dim());
		}
	}
	break;

	case DATA_MY_MINMAX:
	{
		if (!SMesh[dConfig.meshName]->is_atomistic()) {

			return Any(dynamic_cast<Mesh*>(SMesh[dConfig.meshName])->GetmagnetizationYMinMax(dConfig.rectangle));
		}
		else {

			return Any(dynamic_cast<Atom_Mesh*>(SMesh[dConfig.meshName])->GetMomentYMinMax(dConfig.rectangle) * MUB / SMesh[dConfig.meshName]->h.dim());
		}
	}
	break;

	case DATA_MZ_MINMAX:
	{
		if (!SMesh[dConfig.meshName]->is_atomistic()) {

			return Any(dynamic_cast<Mesh*>(SMesh[dConfig.meshName])->GetmagnetizationZMinMax(dConfig.rectangle));
		}
		else {

			return Any(dynamic_cast<Atom_Mesh*>(SMesh[dConfig.meshName])->GetMomentZMinMax(dConfig.rectangle) * MUB / SMesh[dConfig.meshName]->h.dim());
		}
	}
	break;

	case DATA_M_MINMAX:
	{
		if (!SMesh[dConfig.meshName]->is_atomistic()) {

			return Any(dynamic_cast<Mesh*>(SMesh[dConfig.meshName])->GetmagnetizationMinMax(dConfig.rectangle));
		}
		else {

			return Any(dynamic_cast<Atom_Mesh*>(SMesh[dConfig.meshName])->GetMomentMinMax(dConfig.rectangle) * MUB / SMesh[dConfig.meshName]->h.dim());
		}
	}
	break;

	case DATA_MONTECARLOPARAMS:
	{
		return Any(SMesh[dConfig.meshName]->Get_MonteCarlo_Params());
	}
	break;

	case DATA_HA:
	{
		return Any(SMesh[dConfig.meshName]->CallModuleMethod(&ZeemanBase::GetField));
	}
	break;

	case DATA_JC:
	{
		if (SMesh[dConfig.meshName]->IsModuleSet(MOD_TRANSPORT)) {

			return Any(SMesh[dConfig.meshName]->CallModuleMethod(&Transport::GetAverageChargeCurrent, dConfig.rectangle, {}));
		}
	}
	break;

	case DATA_JSX:
	{
		if (SMesh[dConfig.meshName]->IsModuleSet(MOD_TRANSPORT)) {

			return Any(SMesh[dConfig.meshName]->CallModuleMethod(&Transport::GetAverageSpinCurrent, 0, dConfig.rectangle, {}));
		}
	}
	break;

	case DATA_JSY:
	{
		if (SMesh[dConfig.meshName]->IsModuleSet(MOD_TRANSPORT)) {

			return Any(SMesh[dConfig.meshName]->CallModuleMethod(&Transport::GetAverageSpinCurrent, 1, dConfig.rectangle, {}));
		}
	}
	break;

	case DATA_JSZ:
	{
		if (SMesh[dConfig.meshName]->IsModuleSet(MOD_TRANSPORT)) {

			return Any(SMesh[dConfig.meshName]->CallModuleMethod(&Transport::GetAverageSpinCurrent, 2, dConfig.rectangle, {}));
		}
	}
	break;

	case DATA_RESPUMP:
	{
		return Any(SMesh[dConfig.meshName]->Average_mxdmdt(dConfig.rectangle));
	}
	break;

	case DATA_RESPUMP2:
	{
		return Any(SMesh[dConfig.meshName]->Average_mxdmdt2(dConfig.rectangle));
	}
	break;

	case DATA_IMSPUMP:
	{
		return Any(SMesh[dConfig.meshName]->Average_dmdt(dConfig.rectangle));
	}
	break;

	case DATA_IMSPUMP2:
	{
		return Any(SMesh[dConfig.meshName]->Average_dmdt2(dConfig.rectangle));
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
		if (!SMesh.IsSuperMeshModuleSet(MODS_SDEMAG) || SMesh.CallModuleMethod<bool, SDemag>(&SDemag::Get_Multilayered_Convolution_Status)) {

			//read demag energy from named mesh if supermesh demag module not set (also works for MOD_ATOM_DIPOLEDIPOLE since it's exclusive with MOD_DEMAG)
			//or else if SDemag set, but using multilayered convolution (then get energy directly from SDemag_Demag module
			return Any(SMesh[dConfig.meshName]->GetEnergyDensity(MOD_DEMAG, dConfig.rectangle));
		}
		else {

			//get energy value from set supermesh demag module, where we are using supermesh convolution
			return Any(SMesh.GetSuperMeshModule(MODS_SDEMAG)->GetEnergyDensity(dConfig.rectangle));
		}
	}
	break;

	case DATA_E_EXCH:
	{
		return Any(SMesh[dConfig.meshName]->GetEnergyDensity(MOD_EXCHANGE, dConfig.rectangle));
	}
	break;

	case DATA_T_EXCH:
	{
		return Any(SMesh[dConfig.meshName]->GetTorque(MOD_EXCHANGE, dConfig.rectangle));
	}
	break;

	case DATA_E_SURFEXCH:
	{
		return Any(SMesh[dConfig.meshName]->GetEnergyDensity(MOD_SURFEXCHANGE, dConfig.rectangle));
	}
	break;

	case DATA_T_SURFEXCH:
	{
		return Any(SMesh[dConfig.meshName]->GetTorque(MOD_SURFEXCHANGE, dConfig.rectangle));
	}
	break;

	case DATA_E_ZEE:
	{
		return Any(SMesh[dConfig.meshName]->GetEnergyDensity(MOD_ZEEMAN, dConfig.rectangle));
	}
	break;

	case DATA_T_ZEE:
	{
		return Any(SMesh[dConfig.meshName]->GetTorque(MOD_ZEEMAN, dConfig.rectangle));
	}
	break;

	case DATA_E_MOPTICAL:
	{
		return Any(SMesh[dConfig.meshName]->GetEnergyDensity(MOD_MOPTICAL, dConfig.rectangle));
	}
	break;

	case DATA_E_MELASTIC:
	{
		return Any(SMesh[dConfig.meshName]->GetEnergyDensity(MOD_MELASTIC, dConfig.rectangle));
	}
	break;

	case DATA_E_ANIS:
	{
		return Any(SMesh[dConfig.meshName]->GetEnergyDensity(MOD_ANIUNI, dConfig.rectangle));
	}
	break;

	case DATA_T_ANIS:
	{
		return Any(SMesh[dConfig.meshName]->GetTorque(MOD_ANIUNI, dConfig.rectangle));
	}
	break;

	case DATA_E_ROUGH:
	{
		return Any(SMesh[dConfig.meshName]->GetEnergyDensity(MOD_ROUGHNESS, dConfig.rectangle));
	}
	break;

	case DATA_E_TOTAL:
	{
		return Any(SMesh.GetTotalEnergyDensity());
	}
	break;

	case DATA_Q_TOPO:
	{
		return Any(SMesh[dConfig.meshName]->GetTopologicalCharge(dConfig.rectangle));
	}
	break;

	case DATA_DWSHIFT:
	{
		return Any(SMesh.Get_dwshift());
	}
	break;

	case DATA_DWPOS_X:
	{
		return Any(SMesh[dConfig.meshName]->FitDomainWall_X(dConfig.rectangle));
	}
	break;

	case DATA_DWPOS_Y:
	{
		return Any(SMesh[dConfig.meshName]->FitDomainWall_Y(dConfig.rectangle));
	}
	break;

	case DATA_DWPOS_Z:
	{
		return Any(SMesh[dConfig.meshName]->FitDomainWall_Z(dConfig.rectangle));
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
		return Any(SMesh.CallModuleMethod<double, STransport>(&STransport::GetEnergyDensity));
	}
	break;

	case DATA_TEMP:
	{
		return Any(SMesh[dConfig.meshName]->GetAverageTemperature(dConfig.rectangle));
	}
	break;

	case DATA_TEMP_L:
	{
		return Any(SMesh[dConfig.meshName]->GetAverageLatticeTemperature(dConfig.rectangle));
	}
	break;

	case DATA_HEATDT:
	{
		return Any(SMesh.CallModuleMethod(&SHeat::get_heat_dT));
	}
	break;

	case DATA_COMMBUFFER:
	{
		//the main purpose of this is to execute the command buffer during a simulation data saving schedule, so e.g. showdata commbuf is not useful
		if (!is_thread_running(THREAD_HANDLEMESSAGE)) {

			//This could do with more thinking. The problem is simulationMutex must be unlocked before we call RunCommandBuffer
			//try_lock to test if locked shouldn't be used this way since it can fail spuriously according to documentation - OS dependent.
			//The solution for now is : if we are here, then the only remaining possibility this has been called during a simulation schedule.
			//In this case simulationMutex is locked, so unlock it then lock back again.
			//Note, launching RunCommandBuffer asynchronously doesn't work since we have to complete running the buffered commands before proceeding - i.e. must be synchronous.
			simulationMutex.unlock();
			RunCommandBuffer();
			simulationMutex.lock();
		}

		return Any(command_buffer.size());
	}
	break;
	}

	return Any(0);
}

std::string Simulation::GetDataValueString(DatumConfig dConfig, bool ignore_unit) 
{
	std::string unit;
	if (!ignore_unit) unit = dataDescriptor(dConfig.datumId).unit;

	return GetDataValue(dConfig).convert_to_string(unit);
}

void Simulation::NewSaveDataEntry(DATA_ dataId, std::string meshName, Rect dataRect) 
{
	//if not meshless, make sure meshName is valid
	if (!dataDescriptor(dataId).meshless && !SMesh.contains(meshName)) return;

	if (dataDescriptor(dataId).meshless) meshName = "";

	//if not boxless and box is null, this actually means use the entire mesh - make this change (last check should not be really needed as you cannot have an entry which is meshless but not boxless)
	if (!dataDescriptor(dataId).boxless && dataRect.IsNull() && SMesh.contains(meshName)) { dataRect.e = SMesh[meshName]->GetMeshDimensions(); }

	if (dataDescriptor(dataId).boxless) dataRect = Rect();

	saveDataList.push_back(DatumConfig(dataId, meshName, dataRect));
}

void Simulation::EditSaveDataEntry(int index, DATA_ dataId, std::string meshName, Rect dataRect) 
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

void Simulation::UpdateSaveDataEntries(std::string oldMeshName, std::string newMeshName) 
{
	for (int idx = 0; idx < saveDataList.size(); idx++) {

		if (saveDataList[idx].meshName == oldMeshName) saveDataList[idx].meshName = newMeshName;
	}
}

void Simulation::UpdateSaveDataEntries(Rect rect_old, Rect rect_new, std::string meshName) 
{
	for (int idx = 0; idx < saveDataList.size(); idx++) {

		if (!saveDataList[idx].rectangle.IsNull() && saveDataList[idx].meshName == meshName)
			saveDataList[idx].rectangle.resize(rect_old.size(), rect_new.size());
	}
}

void Simulation::DeleteSaveDataEntries(std::string meshName) 
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
	//First build text to write to data file as a single row
	std::string row_text;

	//save actual values as configured in saveDataList
	for (int idx = 0; idx < saveDataList.size(); idx++) {

		std::string value_string = GetDataValueString(saveDataList[idx], true);
		replaceall(value_string, ", ", "\t");

		row_text += value_string;
		if (idx != saveDataList.size() - 1) row_text += "\t";
	}

	if (!is_thread_running(THREAD_DISKACCESS)) {
		
		if (savedata_diskbuffer_position < savedata_diskbuffer_size) savedata_diskbuffer[savedata_diskbuffer_position++] = row_text;
		if (savedata_diskbuffer_position == savedata_diskbuffer_size) {

			//flush overflow buffer synchronously if we have to: this contains entries written before the regular buffer could be completely flushed asynchronously, so come before any new entries in the regular buffer
			if (savedata_diskoverflowbuffer_position) SaveData_DiskBufferFlush(&savedata_diskoverflowbuffer, &savedata_diskoverflowbuffer_position);

			//flush regular buffer asynchronously
			while (single_call_launch<std::vector<std::string>*, int*>(&Simulation::SaveData_DiskBufferFlush, &savedata_diskbuffer, &savedata_diskbuffer_position, THREAD_DISKACCESS) != THREAD_DISKACCESS);
		}
	}
	else {

		//regular buffer still being written: add entries in overflow buffer
		if (savedata_diskoverflowbuffer_position < savedata_diskbuffer_size) savedata_diskoverflowbuffer[savedata_diskoverflowbuffer_position++] = row_text;
		else {

			//overflow buffer also full : flush it synchronously when the previous write operation finishes
			//stop disk thread: this command will complete when the writing operation finishes
			stop_thread(THREAD_DISKACCESS);
			//now flush overflow buffer synchronously when previous operation finishes
			SaveData_DiskBufferFlush(&savedata_diskoverflowbuffer, &savedata_diskoverflowbuffer_position);
		}
	}

	//Image saving:
	if (saveImageFlag) {

		//update and refresh mesh display before saving image (the image is captured from on-screen image)
		UpdateMeshDisplay();

		std::string imageFile = directory + imageSaveFileBase + ToString(SMesh.GetIteration()) + ".png";
		BD.SaveMeshImage(imageFile, image_cropping);
	}
}

//This method does the actual writing to disk : can be launched asynchronously from SaveData method
void Simulation::SaveData_DiskBufferFlush(std::vector<std::string>* pdiskbuffer, int* pdiskbuffer_position)
{
	if (saveDataFlag) {

		//we need a file name to save to
		if (savedataFile.size()) {

			std::ofstream bdout;

			//append to file or make a new one ?
			if (appendToDataFile) bdout.open((directory + savedataFile).c_str(), std::ios::out | std::ios::app);
			else {

				//Create new file
				bdout.open((directory + savedataFile).c_str(), std::ios::out);
				appendToDataFile = true;

				//Append header
				time_t rawtime;
				time(&rawtime);
				bdout << std::string(ctime(&rawtime)) << std::endl;

				//List meshes
				for (int idx = 0; idx < SMesh.size(); idx++) {

					bdout << "<" + SMesh.key_from_meshIdx(idx) + "> : ";
					bdout << "Rectangle : " << ToString(SMesh[idx]->GetMeshRect(), "m") << ". ";
					bdout << "Cells : " << ToString(SMesh[idx]->GetMeshSize()) << ". ";
					bdout << "Cellsize : " << ToString(SMesh[idx]->GetMeshCellsize(), "m");
					bdout << std::endl;
				}

				//Supermesh settings
				bdout << "<" + SMesh.superMeshHandle + "> : ";
				bdout << "Rectangle : " << ToString(SMesh.GetFMSMeshRect(), "m") << ". ";
				bdout << "Cells : " << ToString(SMesh.GetFMSMeshsize()) << ". ";
				bdout << "Cellsize : " << ToString(SMesh.GetFMSMeshCellsize(), "m");
				bdout << std::endl;

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

				bdout << std::endl << std::endl;
			}

			//write buffer entries
			for (int idx = 0; idx < *pdiskbuffer_position; idx++) {

				bdout << (*pdiskbuffer)[idx] << std::endl;
			}

			//buffer flushed now
			*pdiskbuffer_position = 0;

			bdout.close();
		}
	}
}