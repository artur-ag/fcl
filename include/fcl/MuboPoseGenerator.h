#include <fcl/math/transform.h>
#include <fcl/shape/SuperOvoid.h>
#include <fcl/SuperOvoid_global.h>

namespace fcl
{
	class MuboPoseGenerator
	{
	private:
		/** State for the simple RNG */
		int randState;
		
	public:
		/** Predictable pseudo-random function to always generate the same output, regardless of running platform. */
		int next();
		
		void setSeed(int seed);
		float getRandomRange(FCL_REAL min, FCL_REAL max);
		float getNormalSample();
		void getRandomSuperovoid(FCL_REAL* dest);
		void randomizeSuperovoid(SuperOvoid& sov);
		void getRandomQuat(FCL_REAL* params);
		void setRandomRotation(CollisionObject* obj);
	};
}
