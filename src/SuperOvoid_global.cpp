#include <fcl/SuperOvoid_global.h>
#include <cstring>

namespace fcl
{
	bool g_lastStatsValid = false;
	NewtonRaphsonStats g_lastStats;

	void NewtonRaphsonStats::resetToDefault()
	{
		// Clear everything
		memset(this, 0, sizeof(NewtonRaphsonStats));

		// Set some default values
		returnBest = true;
		precision = 1e-6;
		maxIterations = 30;
		guessQuality = 6;
		analytical = false;
		parametric = false;
		superellipsoid = true;
		guessType = NewtonRaphsonStats::INITIAL_GUESS::PARAMETRIC_QUADTREE;
	}

	void NewtonRaphsonStats::clearOutputs()
	{
		iterations = 0;
		setupTime = 0;
		initialGuessTime = 0;
		iterationsTime = 0;
		numericalJacobianTime = 0;
		totalTime = 0;
		guessEstimations = 0;
		nonMinimumDistance = 0;
		discardedGuess = 0;
		didNotConverge = 0;

		for (int i = 0; i < 3; i++)
		{
			guessA[i] = 0;
			guessB[i] = 0;
		}

		for (int i = 0; i < 2; i++)
		{
			outParamsA[i] = 0;
			outParamsB[i] = 0;
		}
	}

	void NewtonRaphsonStats::add(const NewtonRaphsonStats& other)
	{
		this->iterations += other.iterations;
		this->setupTime += other.setupTime;
		this->initialGuessTime += other.initialGuessTime;
		this->iterationsTime += other.iterationsTime;
		this->numericalJacobianTime += other.numericalJacobianTime;
		this->totalTime += other.totalTime;
		this->guessEstimations += other.guessEstimations;
		this->nonMinimumDistance += other.nonMinimumDistance;
		this->didNotConverge += other.didNotConverge;
		this->discardedGuess += other.discardedGuess;
	}

	void NewtonRaphsonStats::averageDivide(int n)
	{
		iterations = (iterations / n);
		setupTime = (long long)(1.0 * setupTime / n);
		initialGuessTime = (long long)(1.0 * initialGuessTime / n);
		iterationsTime = (long long)(1.0 * iterationsTime / n);
		numericalJacobianTime = (long long)(1.0 * numericalJacobianTime / n);
		totalTime = (long long)(1.0 * totalTime / n);
		guessEstimations = (int)(1.0 * guessEstimations / n);
		nonMinimumDistance = (int)(1.0 * nonMinimumDistance / n);
		didNotConverge = (int)(1.0 * didNotConverge / n);
		discardedGuess = (int)(1.0 * discardedGuess / n);
	}

	// Writes a header for a NewtonRaphsonStats in CSV format
	void NewtonRaphsonStats::printHeader(std::ostream& out)
	{
		out <<
			// In =======================================
			"precision," <<
			"maxIterations," <<
			"guessQuality," <<
			"returnBest," <<
			"analytical," <<
			"parametric," <<
			"superellipsoid," <<
			// Out ======================================
			"guessType," <<
			"iterations," <<
			"setupTime," <<
			"initialGuessTime," <<
			"iterationsTime," <<
			"numericalJacobianTime," <<
			"totalTime," <<
			"guessEstimations," <<
			"nonMinimumDistance," <<
			"discardedGuess," <<
			"didNotConverge," <<
			"guessA," <<
			"guessB";
	}

	// Writes NewtonRaphsonStats data in CSV format
	std::ostream& operator<< (std::ostream& out, NewtonRaphsonStats& s)
	{
		out <<
			// In =======================================
			s.precision << "," <<
			s.maxIterations << "," <<
			s.guessQuality << "," <<
			s.returnBest << "," <<
			s.analytical << "," <<
			s.parametric << "," <<
			s.superellipsoid << "," <<
			// Out ======================================
			s.guessType << "," <<
			s.iterations << "," <<
			s.setupTime << "," <<
			s.initialGuessTime << "," <<
			s.iterationsTime << "," <<
			s.numericalJacobianTime << "," <<
			s.totalTime << "," <<
			s.guessEstimations << "," <<
			s.nonMinimumDistance << "," <<
			s.discardedGuess << "," <<
			s.didNotConverge << "," <<
			s.getGuessA() << "," <<
			s.getGuessB();

		return out;
	}
}