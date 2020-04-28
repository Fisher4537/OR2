#include <time.h>
#include <stdlib.h>
//#include <sys/resource.h>

#ifdef _WIN32
#include <windows.h>
#include <sysinfoapi.h>
#endif
#if defined(__MACH__) && defined(__APPLE__)
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

double myWallTime();
double second();