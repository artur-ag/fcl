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
		int next()
		{
			int const a = 1103515245;
			int const c = 12345;
			randState = a * randState + c;
			return (randState >> 16) & 0x7FFF;
		}

		void setSeed(int seed)
		{
			randState = seed;
		}

		float getRandomRange(FCL_REAL min, FCL_REAL max)
		{
			return min + (next() % 1024 / 1024.0) * (max - min);
		}

		float getNormalSample()
		{
			const int levels = 10240;
			const int samples = 16;

			FCL_REAL sum = 0;
			for (int i = 0; i < samples; i++)
				sum += (next() % levels) * 1.0 / levels;

			return sum / samples - 0.5;
		}

		void getRandomSuperovoid(FCL_REAL* dest)
		{
			// Debug: use spheres

			dest[0] = 1.0;
			dest[1] = 1.0;
			dest[2] = 1.0;

			dest[3] = 1;
			dest[4] = 1;

			dest[5] = 0;
			dest[6] = 0;

			//sov.a1 = getRandomRange(0.5, 1.5);
			//sov.a2 = getRandomRange(0.5, 1.5);
			//sov.a3 = getRandomRange(0.5, 1.5);

			//sov.epsilon1 = getRandomRange(0.3, 1.1);
			//sov.epsilon2 = getRandomRange(0.3, 1.1);

			//sov.taperingX = getRandomRange(-0.4, 0.4);
			//sov.taperingY = getRandomRange(-0.4, 0.4);
		}

		void randomizeSuperovoid(SuperOvoid& sov)
		{
			FCL_REAL params[7];
			getRandomSuperovoid(params);

			sov.a1 = params[0];
			sov.a2 = params[1];
			sov.a3 = params[2];

			sov.epsilon1 = params[3];
			sov.epsilon2 = params[4];

			sov.taperingX = params[5];
			sov.taperingY = params[6];
		}

		void getRandomQuat(FCL_REAL* params)
		{
			for (int i = 0; i < 4; i++)
				params[i] = getNormalSample();

			FCL_REAL length = 0;
			for (int i = 0; i < 4; i++)
				length += params[i] * params[i];
			length = std::sqrt(length);

			for (int i = 0; i < 4; i++)
				params[i] /= length;
		}

		void setRandomRotation(CollisionObject* obj)
		{
			FCL_REAL quat[4];
			getRandomQuat(quat);
			Quaternion3f rotation(quat[0], quat[1], quat[2], quat[3]);
			obj->setQuatRotation(rotation);
		}
	};
}
