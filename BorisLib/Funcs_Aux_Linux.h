#pragma once

#include "BorisLib_Config.h"
#if OPERATING_SYSTEM == OS_LIN

#include <unistd.h>
#include <sys/sysinfo.h>
#include <time.h>

#include "Funcs_Aux_base.h"
#include "Types_VAL.h"

///////////////////////////////////////////////////////////////////////////////
//MEMORY

//Return total free memory in MB
inline size_t MemGetFreeMB(void)
{
	struct sysinfo myinfo;
    unsigned long total_bytes;

    sysinfo(&myinfo);

    total_bytes = myinfo.freeram;

    return total_bytes / (1024 * 1024);
}

//Return total memory in MB
inline size_t MemGetTotalMB(void)
{
    struct sysinfo myinfo;
    unsigned long total_bytes;

    sysinfo(&myinfo);

    total_bytes = myinfo.totalram;

    return total_bytes / (1024 * 1024);
}

inline unsigned GetSystemTickCount() 
{ 
	struct timespec ts;
    unsigned ticks = 0U;
    clock_gettime( CLOCK_REALTIME, &ts );
    ticks  = ts.tv_nsec / 1000000;
    ticks += ts.tv_sec * 1000;
    return ticks;
}

#endif