#ifndef FCL_SUPEROVOID_GLOBAL_H
#define FCL_SUPEROVOID_GLOBAL_H

#include <iostream>
#include <fcl/data_types.h>
#include <fcl/math/vec_3f.h>

#ifdef _WIN32
#define NOMINMAX  // required to avoid compilation errors with Visual Studio 2010
#include <windows.h>
#else
#include <time.h>
#endif

namespace fcl
{
	class Timer
	{
	public:
		Timer();
		~Timer();

		void start();                               ///< start timer
		void stop();                                ///< stop the timer
		double getElapsedTime();                    ///< get elapsed time in milli-second
		double getElapsedTimeInSec();               ///< get elapsed time in second (same as getElapsedTime)
		double getElapsedTimeInMilliSec();          ///< get elapsed time in milli-second
		double getElapsedTimeInMicroSec();          ///< get elapsed time in micro-second

	private:
		double startTimeInMicroSec;                 ///< starting time in micro-second
		double endTimeInMicroSec;                   ///< ending time in micro-second
		int stopped;                                ///< stop flag
#ifdef _WIN32
		LARGE_INTEGER frequency;                    ///< ticks per second
		LARGE_INTEGER startCount;
		LARGE_INTEGER endCount;
#else
		struct timespec startCount;
		struct timespec endCount;
#endif
	};


	struct NewtonRaphsonStats
	{
		// In
		double precision;    // Tolerance for NewtonRaphson algorithm.
		int maxIterations;   // Alg stops after this number of iterations even if 'precision' was not achieved.
		int guessQuality;    // Either the number of cells in OcTree, or Parametric QuadTree.

		bool returnBest;     // Alg returns the best result it found overall, if the last found iteration was worse.
		bool analytical;     // True to use analytical jacobian matrices. False uses numerical w/ finite differences.
		bool parametric;     // Parametric version of algorithm. If false, implicit version is used.
		bool forceImplicitNormals; // Uses implicit formula for normals even when parametric = true.
		bool superellipsoid; // Allow superovoids with taperingX = taperingY = 0 to have simpler math expressions (faster)

		enum INITIAL_GUESS { AVG_SPHERE, OCTREE, OBB, MESH, PARAMETRIC_QUADTREE };
		NewtonRaphsonStats::INITIAL_GUESS guessType;

		// Out
		// Time fields in microseconds
		double iterations;					// Number of executed iterations. If 0, broad-phase determined superovoids did not collide.
		long long setupTime;				// Not used.
		long long initialGuessTime;			// Time spent on estimation of initial iteration for NewtonRaphson.
		long long iterationsTime;           // Time spent on NewtonRaphson iterations.
		long long numericalJacobianTime;    // Time spent on computing Jacobian with finite differences.
		long long totalTime;                // Time spent overall, excluding broad-phase.

		int guessEstimations;				// Number of times a new guess was estimated
		int nonMinimumDistance;				// Number of times NewtonRaphson converged to a non-minimum distance solution
		int discardedGuess;					// Number of times NewtonRaphson diverged with cached guess and repeated with a new guess
		int didNotConverge;					// Number of times the NewtonRaphson algorithm diverged

		// Initial guess calculated inside the algorithm
		FCL_REAL guessA[3];
		FCL_REAL guessB[3];

		// Final result in Azimuth,Zenith parameters, if (parametric==1)
		FCL_REAL outParamsA[2];
		FCL_REAL outParamsB[2];

		void resetToDefault();
		void clearOutputs();
		void add(const NewtonRaphsonStats& other);
		void averageDivide(int n);

        Vec3f getGuessA() const
        {
            return Vec3f(guessA[0], guessA[1], guessA[2]);
        }

        Vec3f getGuessB() const
        {
            return Vec3f(guessB[0], guessB[1], guessB[2]);
        }
		
        void print(std::ostream &stream)
        {
            stream << "Iterations: " << iterations << std::endl
                << "Guess quality: " << guessQuality << std::endl
                << "Guess estimations: " << guessEstimations << std::endl
                << "NonMinimumDistance: " << nonMinimumDistance << std::endl
                << "DiscardedGuess: " << discardedGuess << std::endl
                << "DidNotConverge: " << didNotConverge << std::endl
                << "Guess A: " << getGuessA() << std::endl
                << "Guess B: " << getGuessB() << std::endl;
        }

		static void printHeader(std::ostream& out);

		friend std::ostream& operator<< (std::ostream& out, NewtonRaphsonStats& s);
		friend std::ostream& operator<< (std::ostream& out, const NewtonRaphsonStats& s);
	};

	// Global storage for the last collision detection using SuperOvoids
	extern bool g_lastStatsValid;
	extern NewtonRaphsonStats g_lastStats;
    extern int g_superovoidQueries;
}

#endif