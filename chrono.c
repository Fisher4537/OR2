#include "chrono.h"

double second(){

	#ifdef __APPLE__
		static double timeConvert = 0.0;
		if (timeConvert == 0.0){
			mach_timebase_info_data_t timeBase;
			mach_timebase_info(&timeBase);
			timeConvert = (double)timeBase.numer / (double)timeBase.denom / 1000000000.0;
		}
		return mach_absolute_time() * timeConvert;

	#elif __linux__
		struct timespec ts;
		clock_gettime(CLOCK_MONOTONIC, &ts);
		return (double)ts.tv_sec + 1.0e-9 * ((double)ts.tv_nsec);

	#elif _WIN32
		/* Metodo 1 -> QueryPerformanceCounter*/
		// granularity about 50 microsecs on my machine
		static LARGE_INTEGER freq, start;
		LARGE_INTEGER count;
		if (!QueryPerformanceCounter(&count))
			return 0.0;
		if (!freq.QuadPart) { // one time initialization
			if (!QueryPerformanceFrequency(&freq))
				return 0.0;
			start = count;
		}
		return (double)(count.QuadPart - start.QuadPart) / freq.QuadPart;

		/* Metodo 2
		SYSTEMTIME sm;
		GetSystemTime(&sm);
		return (double)sm.wSecond + 1.0e-6 * ((double)sm.wMilliseconds);
		*/

		/* Metodo 2	-> CPXgettime in tsp.c
		struct timespec ts, ts2;
		CPXgettime(env, &ts);
		mip_optimization(env, lp, inst, &status);
		CPXgettime(env, &ts2);
		inst->opt_time = (double)ts2.tv_nsec + 1.0e-9 * ((double)ts2.tv_nsec) - (double)ts.tv_nsec + 1.0e-9 * ((double)ts.tv_nsec);
		*/

	#endif
	
		return 0.0;

}
