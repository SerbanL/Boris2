#pragma once

#include "resource.h"

#include "Simulation.h"
#include "BorisLib.h"
#include "Commands.h"

////////////////////////////////////
//Command-line arguments

//default arguments
static std::string server_port = "";
static std::string server_pwd = "changeme";
static int cudaDevice = -1;

//show window front or back
static std::string window_startup_option = "front";

////////////////////////////////////
//Version

static int Program_Version = 290;

////////////////////////////////////
//Program master object

static Simulation *pSim = nullptr;

