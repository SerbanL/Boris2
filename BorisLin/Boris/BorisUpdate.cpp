#include "stdafx.h"
#include "Simulation.h"
#include "Boris.h"

//-------------------------------------Program update checker

//check latest version available
int Simulation::CheckLatestVersion(void)
{
	HTTP httpComms(domain_name);

	std::string response_message;

	if (httpComms.is_available() && httpComms.http_post(version_checker, std::string("version=") + ToString(Program_Version), response_message)) {

		return ToNum(response_message);
	}
	else return 0;
}

//Check with "www.boris-spintronics.uk" if program version is up to date
void Simulation::CheckUpdate(void)
{
	std::string response_message;

	std::string console_update_message = "[tc0,0.5,0,1/tc]Program version update status : " + MakeIO(IOI_PROGRAMUPDATESTATUS);

	console_update_message += "</c>[tc0,0.5,0,1/tc] Check for updates on startup : " + MakeIO(IOI_UPDATESTATUSCHECKSTARTUP, start_check_updates);

	BD.DisplayFormattedConsoleMessage(console_update_message);

	//trying to connect
	version_checking = -1;

	HTTP httpComms(domain_name);

	if (httpComms.is_available() && httpComms.http_post(version_checker, std::string("version=") + ToString(Program_Version), response_message)) {

		int latest_version = ToNum(response_message);

		if (latest_version == Program_Version) {

			//program up to date
			version_checking = 1;
		}
		else {

			//update available
			version_checking = 2;
		}

		if (httpComms.is_available() && httpComms.http_get(mdb_lastupdate, response_message)) {

			BD.DisplayConsoleMessage("Online materials database last updated on (PT) : " + response_message);
		}
	}
	else {

		//connection failure
		version_checking = 0;
	}

	RefreshScreen();
}

void Simulation::OpenDownloadPage(void)
{
	open_web_page(domain_name + download_page);
}

#if GRAPHICS == 1 && OPERATING_SYSTEM == OS_WIN

//-------------------------------------Incremental Update : Windows only

#include <Urlmon.h>
#pragma comment(lib, "Urlmon.lib")

//for current program version check available incremental update versions and return them, including current version as first element.
//if incremental update not possible (e.g. because program is up to date) then returned vector will have only one element as the current version
//if error occured returned vector is empty
std::vector<int> Simulation::CheckDeltaUpdate(void)
{
	std::vector<int> update_versions;

	BD.DisplayConsoleMessage("Checking for incremental updates.");

	std::string from_version = ToString(Program_Version);

	HTTP httpComms(domain_name);
	std::string response_message;
	if (httpComms.is_available() && httpComms.http_post(delta_update_checker, std::string("from_version=") + from_version, response_message)) {

		std::vector<std::string> fields = split(response_message);

		int latest_version = Program_Version, next_version = Program_Version;
		
		if (fields.size() >= 1) {

			latest_version = ToNum(fields.back());
			if (fields.size() >= 2) next_version = ToNum(fields[1]);
		}

		if (latest_version == Program_Version && fields.size()) {

			BD.DisplayConsoleMessage("Program is up to date.");
			update_versions.push_back(Program_Version);
		}
		else if (next_version > Program_Version) {

			BD.DisplayConsoleMessage("Found incremental update versions : " + combine(subvec(fields, 1), ", "));
			
			for (int idx = 0; idx < fields.size(); idx++) update_versions.push_back(ToNum(fields[idx]));
		}
		else BD.DisplayConsoleError("Version database error.");
	}
	else BD.DisplayConsoleError("Couldn't connect.");

	return update_versions;
}

//check and download incremental update files; if target_version specified download only incremental update files up to target_version
std::vector<std::string> Simulation::DownloadDeltaUpdate(int target_version)
{
	std::vector<std::string> update_fileNames;

	//check available updates
	std::vector<int> update_versions = CheckDeltaUpdate();

	std::string program_directory = GetExeDirectory();
	std::vector<std::string> files_in_dir = GetFilesInDirectory(program_directory, "Boris", "");
	for (int idx = 0; idx < files_in_dir.size(); idx++) ExtractFilenameDirectory(files_in_dir[idx]);

	//Patch files are of the form:
	//BorisWinapi_Lite_patch_340_to_341
	//BorisWinapi_CUDA5_SP_patch_340_to_341
	//BorisWinapi_CUDA5_DP_patch_340_to_341
	//etc.
	auto download_patch_file = [&](int next_version, int previous_version) -> std::string {

		std::string update_file = "BorisWinapi";

#if COMPILECUDA == 1
		//Cuda version
		update_file += std::string("_CUDA") + ToString((double)__CUDA_ARCH__ / 100);

#if SINGLEPRECISION == 1
		//Single precision
		update_file += std::string("_SP");
#else
		//Double precision
		update_file += std::string("_DP");
#endif

#else
		//Lite version
		update_file += std::string("_Lite");
#endif

		update_file += std::string("_patch_") + ToString(previous_version) + std::string("_to_") + ToString(next_version);

		if (vector_contains(files_in_dir, update_file)) {

			BD.DisplayConsoleMessage("Patch file already exists in directory; not downloaded. " + update_file);
			return update_file;
		}

		BD.DisplayConsoleMessage("Downloading incremental update : " + update_file);

		std::string url = "http://" + domain_name + installers_dir + "/" + ToString((double)next_version / 100) + "/" + update_file;

		HRESULT res = URLDownloadToFile(nullptr, StringtoWCHARPointer(url), StringtoWCHARPointer(program_directory + update_file), 0, nullptr);

		if (res == S_OK) BD.DisplayConsoleMessage("Patch file downloaded : " + update_file);
		else if (res == INET_E_DOWNLOAD_FAILURE) { BD.DisplayConsoleError("Could not download file from : " + url); return ""; }
		else { BD.DisplayConsoleError("Download available, but error saving file (make sure you have privileges for writing to program folder): " + url); return ""; }

		return update_file;
	};

	for (int idx = 1; idx < update_versions.size(); idx++) {

		int next_version = update_versions[idx], previous_version = update_versions[idx - 1];

		if (!target_version || target_version >= next_version) {

			//download all relevant patch files

			std::string fileName = download_patch_file(next_version, previous_version);
			if (fileName.length()) update_fileNames.push_back(fileName);
		}
	}

	return update_fileNames;
}

//check, download, and install latest update file
int Simulation::InstallDeltaUpdate(int target_version)
{
	std::vector<std::string> patch_files = DownloadDeltaUpdate(target_version);

	if (patch_files.size()) {

		std::string all_patchFiles_string = combine(patch_files, " ");

		//downloaded incremental update. now apply patch.

		BD.DisplayConsoleMessage("Attempting to apply patches.");

		std::string program_directory = GetExeDirectory();
		std::string exe_fileName = GetExeFilename();

		//Launch PatchApply.exe which should be found in same executable directory. This takes parameters:
		//Name of current program executable
		//Name of patch to apply
		//PatchApply will then:
		//1) Wait for current Boris executable to finish
		//2) Apply patch using delta compression via xdelta utility which should be found in same executable directory (Using xdelta from https://github.com/jmacd/xdelta)
		//3) Restart Boris program, which should now be updated
		open_file(program_directory + patch_utility, exe_fileName + " " + all_patchFiles_string);
		
		//Finish current program by cleanly issuing WM_CLOSE command; this must be sent asynchronously
		std::function<void(HWND hWnd)> end_program = [](HWND hWnd) -> void { SendMessage(hWnd, WM_CLOSE, NULL, NULL); };
		std::thread end_program_thread(end_program, hWnd);
		end_program_thread.detach();
	}
	
	return 1;
}

//Roll back to a backup version using the PatchApply.exe utility
void Simulation::RollbackUpdate(int target_version)
{
	//downloaded incremental update. now apply patch.

	BD.DisplayConsoleMessage("Attempting to rollback version.");

	std::string program_directory = GetExeDirectory();
	std::string exe_fileName = GetExeFilename();

	//Launch PatchApply.exe which should be found in same executable directory. This takes parameters:
	//Name of current program executable
	//Rollback option
	//Version number to rollback to
	open_file(program_directory + patch_utility, exe_fileName + " rollback " + ToString(target_version));

	//Finish current program by cleanly issuing WM_CLOSE command; this must be sent asynchronously
	std::function<void(HWND hWnd)> end_program = [](HWND hWnd) -> void { SendMessage(hWnd, WM_CLOSE, NULL, NULL); };
	std::thread end_program_thread(end_program, hWnd);
	end_program_thread.detach();
}

#else

//Incremental update not available for Linux

std::vector<int> Simulation::CheckDeltaUpdate(void) { std::vector<int> empty; return empty; }
std::vector<std::string> Simulation::DownloadDeltaUpdate(int target_version) { std::vector<std::string> empty; return empty; }
int Simulation::InstallDeltaUpdate(int target_version) { return 0; }
void Simulation::RollbackUpdate(int target_version) {}

#endif