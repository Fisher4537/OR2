#include "pch.h"
#include "chrono.h"

double myWallTime()
{ 
#ifdef __APPLE__
	static double timeConvert = 0.0;
	if ( timeConvert == 0.0 )
	{
		mach_timebase_info_data_t timeBase;
		mach_timebase_info(&timeBase);
		timeConvert = (double)timeBase.numer / (double)timeBase.denom / 1000000000.0;
	}
	return mach_absolute_time() * timeConvert;
#elif __linux__
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double)ts.tv_sec + 1.0e-9*((double)ts.tv_nsec);
#elif _WIN32
	SYSTEMTIME sm;
	GetSystemTime(&sm);
	return (double)sm.wSecond + 1.0e-6 * ((double)sm.wMilliseconds);		// accurato ai millisecondi (non perfetto ma abbastanza corretto)
#endif																		// provare ad usare QueryPerformanceCounter
	return 0.0;
}

double second(){
#if 1
	double t = myWallTime();
	return(t);
#else
	return ((double)clock()/(double)CLK_TCK);
#endif
}