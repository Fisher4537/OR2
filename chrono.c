#include "pch.h"
#include <time.h>
#include <stdlib.h>
//#include <sys/resource.h>

#if defined(__MACH__) && defined(__APPLE__)
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

int clock_gettime_monotonic(struct timespec* tv);

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
#else
	struct timespec ts;
	clock_gettime_monotonic(&ts);
	return (double)ts.tv_sec + 1.0e-9*((double)ts.tv_nsec);
#endif // __APPLE__
	return 0.0;
}

#define MS_PER_SEC      1000ULL     // MS = milliseconds
#define US_PER_MS       1000ULL     // US = microseconds
#define HNS_PER_US      10ULL       // HNS = hundred-nanoseconds (e.g., 1 hns = 100 ns)
#define NS_PER_US       1000ULL

#define HNS_PER_SEC     (MS_PER_SEC * US_PER_MS * HNS_PER_US)
#define NS_PER_HNS      (100ULL)    // NS = nanoseconds
#define NS_PER_SEC      (MS_PER_SEC * US_PER_MS * NS_PER_US)

int clock_gettime_monotonic(struct timespec* tv)
{
	static LARGE_INTEGER ticksPerSec;
	LARGE_INTEGER ticks;
	double seconds;

	if (!ticksPerSec.QuadPart) {
		QueryPerformanceFrequency(&ticksPerSec);
		if (!ticksPerSec.QuadPart) {
			errno = ENOTSUP;
			return -1;
		}
	}

	QueryPerformanceCounter(&ticks);

	seconds = (double)ticks.QuadPart / (double)ticksPerSec.QuadPart;
	tv->tv_sec = (time_t)seconds;
	tv->tv_nsec = (long)((ULONGLONG)(seconds * NS_PER_SEC) % NS_PER_SEC);

	return 0;
}

double second(){
#if 1
	double t = myWallTime();
	return(t);
#else
	return ((double)clock()/(double)CLK_TCK);
#endif
}