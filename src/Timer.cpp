#include <fcl/SuperOvoid_global.h>

namespace fcl
{

	Timer::Timer()
	{
#ifdef _WIN32
		QueryPerformanceFrequency(&frequency);
		startCount.QuadPart = 0;
		endCount.QuadPart = 0;
#else
		startCount.tv_sec = startCount.tv_nsec = 0;
		endCount.tv_sec = endCount.tv_nsec = 0;
#endif

		stopped = 0;
		startTimeInMicroSec = 0;
		endTimeInMicroSec = 0;
	}


	Timer::~Timer()
	{
	}


	void Timer::start()
	{
		stopped = 0; // reset stop flag
#ifdef _WIN32
		QueryPerformanceCounter(&startCount);
#else
		// gettimeofday(&startCount, NULL);
		clock_gettime(CLOCK_THREAD_CPUTIME_ID, &startCount);
#endif
	}


	void Timer::stop()
	{
		stopped = 1; // set timer stopped flag

#ifdef _WIN32
		QueryPerformanceCounter(&endCount);
#else
		//gettimeofday(&endCount, NULL);
		clock_gettime(CLOCK_THREAD_CPUTIME_ID, &endCount);
#endif
	}


	double Timer::getElapsedTimeInMicroSec()
	{
#ifdef _WIN32
		if (!stopped)
			QueryPerformanceCounter(&endCount);

		startTimeInMicroSec = startCount.QuadPart * (1000000.0 / frequency.QuadPart);
		endTimeInMicroSec = endCount.QuadPart * (1000000.0 / frequency.QuadPart);
#else
		if (!stopped)
			clock_gettime(CLOCK_THREAD_CPUTIME_ID, &endCount);

		startTimeInMicroSec = (startCount.tv_sec * 1000000.0) + startCount.tv_nsec / 1000.0;
		endTimeInMicroSec = (endCount.tv_sec * 1000000.0) + endCount.tv_nsec / 1000.0;
#endif

		return endTimeInMicroSec - startTimeInMicroSec;
	}


	double Timer::getElapsedTimeInMilliSec()
	{
		return this->getElapsedTimeInMicroSec() * 0.001;
	}


	double Timer::getElapsedTimeInSec()
	{
		return this->getElapsedTimeInMicroSec() * 0.000001;
	}


	double Timer::getElapsedTime()
	{
		return this->getElapsedTimeInMilliSec();
	}
}