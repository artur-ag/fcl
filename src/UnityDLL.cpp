#include "fcl/distance.h"
#include "fcl/collision.h"
#include "fcl/BVH/BVH_model.h"
#include "fcl/BV/BV.h"
#include "fcl/BV/OBBRSS.h"
#include "fcl/shape/geometric_shapes.h"
#include "fcl/narrowphase/narrowphase.h"
#include "fcl/shape/SuperOvoidDetails.h"
#include "fcl/MuboPoseGenerator.h"

#include "fcl/SuperOvoid_global.h"

using namespace fcl;

typedef BVHModel<OBBRSS> Model;

SuperOvoid* superOvoid1;
SuperOvoid* superOvoid2;

Model* mesh1;
Model* mesh2;

Model* getSuperOvoidMesh(SuperOvoid& s, int azSlices, int zeSlices, bool adaptiveMeshes)
{
	std::vector<Vec3f>* vertices = new std::vector<Vec3f>();
	s.storePoints(vertices, azSlices, zeSlices, !adaptiveMeshes);
	std::vector<Triangle>* triangles = new std::vector<Triangle>();
	s.storeTriangles(triangles, azSlices, zeSlices);

	Model* model = new Model();
	model->beginModel();
	model->addSubModel(*vertices, *triangles);
	model->endModel();

	return model;
}

Transform3f initTransform(double* values)
{
	return Transform3f(
		Quaternion3f(values[3], values[4], values[5], values[6]),
		Vec3f(values[0], values[1], values[2]));
}

extern "C"
{
	/// Helper function. Exposes the Superovoid-Superovoid minimum distance algorithm.
	/// s1 and s2:     Superovoid [a1 a2 a3 epsilon1 epsilon2 taperingX, taperingY]
	/// t1 and t2:     Transform (translation and rotation quaternion) [x y z][x y z w]
	/// initialGuess:  Initial guess in global coordinates, P and Q [x y z][x y z]. May be NULL.
	/// out:           Output. P and Q in global coordinates [x y z][x y z]
	__declspec(dllexport) bool doThing(double* s1, double* t1,  double* s2, double* t2,
		double* initialGuess, int guessQuality, double* out, NewtonRaphsonStats* stats)
	{
		if (guessQuality < 1)
			guessQuality = 1;

		Timer setupTimer = Timer();
		setupTimer.start();

		SuperOvoid* sov1 = new SuperOvoid(
			s1[0], s1[1], s1[2],
			s1[3], s1[4],
			s1[5], s1[6],
			guessQuality);
		Transform3f s1_transform(Quaternion3f(t1[3], t1[4], t1[5], t1[6]), Vec3f(t1[0], t1[1], t1[2]));

		SuperOvoid* sov2 = new SuperOvoid(
			s2[0], s2[1], s2[2],
			s2[3], s2[4],
			s2[5], s2[6],
			guessQuality);
		Transform3f s2_transform(Quaternion3f(t2[3], t2[4], t2[5], t2[6]), Vec3f(t2[0], t2[1], t2[2]));

		if (stats != NULL && stats->guessType == NewtonRaphsonStats::OCTREE)
		{
			sov1->computeLocalOctree(guessQuality);
			sov2->computeLocalOctree(guessQuality);
		}

		setupTimer.stop();
		if (stats != NULL)
			stats->setupTime = setupTimer.getElapsedTimeInMicroSec();

		if (initialGuess != NULL)
		{
			sov1->setCachedPoint(sov2, Vec3f(initialGuess[0], initialGuess[1], initialGuess[2]));
			sov2->setCachedPoint(sov1, Vec3f(initialGuess[3], initialGuess[4], initialGuess[5]));
		}

		FCL_REAL distance;
		Vec3f p, q;

		Timer collisionTimer = Timer();
		collisionTimer.start();

		bool separated = superOvoidSuperOvoidDistance(
			*sov1, s1_transform,
			*sov2, s2_transform,
			&distance, &p, &q, false, stats);

		collisionTimer.stop();
		if (stats != NULL)
			stats->totalTime = collisionTimer.getElapsedTimeInMicroSec();

		// Write results
		out[0] = p[0];
		out[1] = p[1];
		out[2] = p[2];

		out[3] = q[0];
		out[4] = q[1];
		out[5] = q[2];

		delete sov1;
		delete sov2;

        return separated;
	}



	__declspec(dllexport) void initSuperovoids(double* s1, double* s2)
	{
		delete superOvoid1;
		delete superOvoid2;

		superOvoid1 = new SuperOvoid(s1[0], s1[1], s1[2], s1[3], s1[4], s1[5], s1[6]);
		superOvoid2 = new SuperOvoid(s2[0], s2[1], s2[2], s2[3], s2[4], s2[5], s2[6]);
	}

	__declspec(dllexport) void initMeshSuperovoids(double* s1, double* s2, int* meshSettings)
	{
		delete mesh1;
		delete mesh2;

		SuperOvoid sov1 = SuperOvoid(s1[0], s1[1], s1[2], s1[3], s1[4], s1[5], s1[6]);
		SuperOvoid sov2 = SuperOvoid(s2[0], s2[1], s2[2], s2[3], s2[4], s2[5], s2[6]);

		mesh1 = getSuperOvoidMesh(sov1, meshSettings[0], meshSettings[1], meshSettings[2] == 1);
		mesh2 = getSuperOvoidMesh(sov2, meshSettings[3], meshSettings[4], meshSettings[5] == 1);
	}

	__declspec(dllexport) bool collideSuperovoids(double* t1, double* t2, double* point)
	{
		Transform3f transform1 = initTransform(t1);
		Transform3f transform2 = initTransform(t2);

		CollisionRequest request;
		request.num_max_contacts = 1;
		request.enable_contact = true;
		CollisionResult result;

		bool collided = collide(superOvoid1, transform1, superOvoid2, transform2, request, result);
		if (point != NULL && collided)
		{
			for (int i = 0; i < 3; i++)
			{
				point[i] = result.getContact(0).pos[i];
			}
		}

		return collided;
	}

	__declspec(dllexport) bool distanceSuperovoids(double* t1, double* t2, double* points)
	{
		Transform3f transform1 = initTransform(t1);
		Transform3f transform2 = initTransform(t2);

		DistanceRequest request;
		request.enable_nearest_points = true;
		DistanceResult result;

		bool distanced = distance(superOvoid1, transform1, superOvoid2, transform2, request, result);

		if (points != NULL)
		{
			for (int i = 0; i < 3; i++)
			{
				points[i+0] = result.nearest_points[0][i];
				points[i+3] = result.nearest_points[1][i];
			}
		}

		return distanced;
	}

	__declspec(dllexport) bool collideMeshes(double* t1, double* t2, double* point)
	{
		Transform3f transform1 = initTransform(t1);
		Transform3f transform2 = initTransform(t2);

		CollisionRequest request;
		request.num_max_contacts = 1;
		request.enable_contact = true;
		CollisionResult result;

		bool collided = collide(mesh1, transform1, mesh2, transform2, request, result);
		if (point != NULL && collided)
		{
			for (int i = 0; i < 3; i++)
			{
				point[i] = result.getContact(0).pos[i];
			}
		}

		return collided;
	}

	__declspec(dllexport) bool distanceMeshes(double* t1, double* t2, double* points)
	{
		Transform3f transform1 = initTransform(t1);
		Transform3f transform2 = initTransform(t2);

		DistanceRequest request;
		request.enable_nearest_points = true;
		DistanceResult result;

		bool distanced = distance(mesh1, transform1, mesh2, transform2, request, result);

		if (points != NULL)
		{
			for (int i = 0; i < 3; i++)
			{
				points[i + 0] = result.nearest_points[0][i];
				points[i + 3] = result.nearest_points[1][i];
			}
		}

		return distanced;
	}

	__declspec(dllexport) bool getStats(NewtonRaphsonStats* stats)
	{
		*stats = g_lastStats;
		return g_lastStatsValid;
	}


	enum VecType
	{
		TYPE_POINT,

		TYPE_NORMAL,
		TYPE_TANGENT,
		TYPE_BINORMAL,

		TYPE_D1_NORMAL,
		TYPE_D1_TANGENT,
		TYPE_D1_BINORMAL,

		TYPE_D2_NORMAL,
		TYPE_D2_TANGENT,
		TYPE_D2_BINORMAL
	};

	__declspec(dllexport) void getVector(double* superovoid, double* point, bool parametric, VecType type, double* out)
	{
		SuperOvoid s = SuperOvoid(
			superovoid[0],
			superovoid[1],
			superovoid[2],
			superovoid[3],
			superovoid[4],
			superovoid[5],
			superovoid[6]);

		Vec3f outVector;

		if (parametric)
		{
			switch (type)
			{
			case TYPE_POINT:
				outVector = s.getPoint(point[0], point[1]); break;

			case TYPE_NORMAL:
				outVector = s.getNormal(point[0], point[1]); break;
			case TYPE_TANGENT:
				outVector = s.getAzimuthTangent(point[0], point[1]); break;
			case TYPE_BINORMAL:
				outVector = s.getZenithTangent(point[0], point[1]); break;
			
			case TYPE_D1_NORMAL:
				outVector = s.getNormalDerivativePhi1(point[0], point[1]); break;
			case TYPE_D1_TANGENT:
				outVector = s.getAzimuthTangentDerivativePhi1(point[0], point[1]); break;
			case TYPE_D1_BINORMAL:
				outVector = s.getZenithTangentDerivativePhi1(point[0], point[1]); break;
			
			default:
				// TODO implement phi2 derivatives
				break;
			}
		}

		// Write output
		if (out != NULL)
			for (int i = 0; i < 3; i++)
				out[i] = outVector[i];
	}

	__declspec(dllexport) void getRandomMuboPose(double* sov1, double* quat1, double* sov2, double* quat2, int seed, int index)
	{
		if (index < 0)
			return;

		MuboPoseGenerator rand = MuboPoseGenerator();
		rand.setSeed(seed);

		// Force the RNG generator to advance to the correct state (index)
		for (int i = -1; i < index; i++)
		{
			rand.getRandomSuperovoid(sov1);
			
			rand.getRandomSuperovoid(sov2);
		}
	}
}