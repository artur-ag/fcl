#define BOOST_TEST_MODULE "FCL_SUPEROVOID"
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <fstream>
#include <cstdio>

#include "fcl/collision.h"
#include "fcl/shape/geometric_shapes.h"
#include "fcl/shape/SuperOvoid.h"
#include "fcl/shape/SuperOvoidDetails.h"
#include "fcl/SuperOvoid_global.h"

using namespace fcl;


void checkSuperEllipsoidCase(SuperOvoid ovoid, SuperOvoid ellipsoid);

BOOST_AUTO_TEST_CASE(superellipsoid)
{
    double
        a1 = 1,
        a2 = 1.4,
        a3 = 0.8;

    for (double e1 = 0.1; e1 < 1.9; e1 += 0.27)
    {
        for (double e2 = 0.11; e2 < 1.9; e2 += 0.24)
        {
            // Superovoid, never changed to superellipsoid
            SuperOvoid ovoid(a1, a2, a3, e1, e2, 0, 0);
            BOOST_CHECK(ovoid.isOvoid() == true);

            // Superellipsoid
            SuperOvoid ellipsoid(a1, a2, a3, e1, e2, 0, 0);
            ellipsoid.toSuperellipsoid();
            BOOST_CHECK(ellipsoid.isOvoid() == false);

            checkSuperEllipsoidCase(ovoid, ellipsoid);
        }
    }

    return;
}

void checkSuperEllipsoidCase(SuperOvoid ovoid, SuperOvoid ellipsoid)
{
    // Confirm that both SuperOvoid and SuperEllipsoid implementations return same values
    for (int a = -10; a < 10; a++)
    {
        double azimuth = a * 3.141592 / 10;

        for (int z = -10; z < 10; z++)
        {
            double zenith = z * 3.141592 / 2 / 10;

            //// Using azimuth and zenith directly, for normals
            //Vec3f ovoidNormal = ovoid.getNormal(azimuth, zenith, false);
            //Vec3f ellipsoidNormal = ellipsoid.getNormal(azimuth, zenith, false);
            //BOOST_CHECK(ovoidNormal.equal(ellipsoidNormal, 1e-9));

            // Using a 3D point, for normals
            Vec3f point = ovoid.getPoint(azimuth, zenith, false);
			Vec3f ovoidNormal = ovoid.getNormal(point);
			Vec3f ellipsoidNormal = ellipsoid.getNormal(point);
            BOOST_CHECK(ovoidNormal.equal(ellipsoidNormal, 1e-9));

            // Confirm that both Ovoid and SuperEllipsoid implementations return same implicit function
            double ovoidImplicit = ovoid.implicitFunction(point);
            double ellipsoidImplicit = ellipsoid.implicitFunction(point);
            BOOST_CHECK(std::abs(ovoidImplicit - ellipsoidImplicit) < 1e-9);
            // Because the points were generated on the surface, implicit should be equal to 0
            BOOST_CHECK(ovoidImplicit < 1e-9);
        }
    }

    BOOST_CHECK(ovoid.isOvoid() == true);
    BOOST_CHECK(ellipsoid.isOvoid() == false);
}

int checkParametricNormal(SuperOvoid ovoid)
{
    int fails = 0;
    FCL_REAL tolerance = 1e-3;

    // Confirm that both Implicit and Parametric implementations return same normals
    for (int a = -9; a < 9; a++)
    {
        double azimuth = a * 3.141592 / 10;

        for (int z = -9; z < 9; z++)
        {
            double zenith = z * 3.141592 / 2 / 10;

            Vec3f implicitNormal = ovoid.getNormal(azimuth, zenith, false);

            Vec3f paramTangent1 = ovoid.getAzimuthTangent(azimuth, zenith);
            Vec3f paramTangent2 = ovoid.getZenithTangent(azimuth, zenith);
            Vec3f paramNormal = paramTangent1.cross(paramTangent2).normalize();

            //std::cout << implicitNormal << std::endl;
            //std::cout << paramNormal << std::endl;
            //std::cout << implicitNormal.norm() << std::endl;
            //std::cout << paramNormal.norm() << std::endl;

            //BOOST_CHECK((implicitNormal - paramNormal).norm() < tolerance);

            if ((implicitNormal - paramNormal).norm() >= tolerance)
            {
                //std::cout << implicitNormal << std::endl;
                //std::cout << paramNormal << std::endl;
                //std::cout << implicitNormal.norm() << std::endl;
                //std::cout << paramNormal.norm() << std::endl;
                fails++;
            }
        }
    }

    return fails;
}

BOOST_AUTO_TEST_CASE(check_parametric_normal)
{
    double
        a1 = 1,
        a2 = 1.4,
        a3 = 0.8;

    int fails = 0;

    for (double e1 = 0.1; e1 < 1.9; e1 += 0.27)
    {
        for (double e2 = 0.11; e2 < 1.9; e2 += 0.24)
        {
            SuperOvoid ovoid(a1, a2, a3, e1, e2, 0.2, 0.2);
            fails += checkParametricNormal(ovoid);
        }
    }

    //std::cout << "Parametric normal fails: " << fails << std::endl;
	BOOST_CHECK(fails == 0);

    return;
}


BOOST_AUTO_TEST_CASE(selly_speed_implicitAndNormal)
{
	// Compare the general SuperOvoid case against the specific SuperEllipsoid case,
	// in which taperingX and taperingY are zero
	// (this simplifies maths behind getNormal and implicitFunction).

	// Approximate results on VS2013 x32 debug build, Intel Core i7 MQ-4700 2.40 GHz:
	// SuperOvoid:     45 ms
	// SuperEllipsoid: 30 ms

	{
		Timer timer;

		SuperOvoid ovoid(1.1, 1.2, 0.8, 0.6, 1.7, 0, 0);
		ovoid.toSuperellipsoid();

		timer.start();

		for (int a = 0; a < 100; a++)
		{
			double azimuth = a * 3.141592 / 100;

			for (int z = 0; z < 100; z++)
			{
				double zenith = z * 3.141592 / 2 / 100;

				// Using azimuth and zenith directly, for normals
				Vec3f ovoidNormal = ovoid.getNormal(azimuth, zenith, false);

				Vec3f point = ovoid.getPoint(azimuth, zenith, false);
				ovoidNormal = ovoid.getNormal(point);
				double ovoidImplicit = ovoid.implicitFunction(point);
			}
		}

		timer.stop();

		std::cout << "SuperEllipsoid specific: " << timer.getElapsedTimeInMilliSec() << " ms" << std::endl;
	}

	{
		Timer timer;

		SuperOvoid ovoid(1.1, 1.2, 0.8, 0.6, 1.7, 0, 0);

		timer.start();

		// Confirm that both Ovoid and SuperEllipsoid implementations return same values
		for (int a = 0; a < 100; a++)
		{
			double azimuth = a * 3.141592 / 100;

			for (int z = 0; z < 100; z++)
			{
				double zenith = z * 3.141592 / 2 / 100;

				// Using azimuth and zenith directly, for normals
				Vec3f ovoidNormal = ovoid.getNormal(azimuth, zenith, false);

				Vec3f point = ovoid.getPoint(azimuth, zenith, false);
				ovoidNormal = ovoid.getNormal(point);
				double ovoidImplicit = ovoid.implicitFunction(point);
			}
		}

		timer.stop();

		std::cout << "SuperOvoid general case: " << timer.getElapsedTimeInMilliSec() << " ms" << std::endl;
	}
}

void testDistance(SuperOvoid& o1, SuperOvoid& o2, NewtonRaphsonStats* options = NULL)
{
	Transform3f t1(Quaternion3f(), Vec3f(0, 0, 0));
	Transform3f t2(Quaternion3f(), Vec3f(0, 1.5, -3));

	o1.clearCachedPoint(&o2);
	o2.clearCachedPoint(&o1);

	const int N = 200;
	for (int i = 0; i < N; i++)
	{
		t2.setTranslation(t2.getTranslation() + Vec3f(0, 0, 4.0 / N));

		details::superOvoidSuperOvoidDistance(o1, t1, o2, t2, NULL, NULL, NULL, false, options);
	}
}

BOOST_AUTO_TEST_CASE(selly_speed_distance)
{
	SuperOvoid o1(1, 1, 1, 0.8, 0.8, 0, 0);
	SuperOvoid o2(1, 1, 1, 0.8, 0.8, 0, 0);

	Timer timer;
	double ovoidTime = 0, ellipsoidTime = 0;

	// SuperEllipsoid specific case
	g_lastStats.resetToDefault();
	g_lastStats.superellipsoid = true;
	g_lastStatsValid = true;

	timer.start();
	testDistance(o1, o2);
	timer.stop();
	ellipsoidTime = timer.getElapsedTimeInMicroSec();

	BOOST_CHECK(o1.isOvoid() == false);
	BOOST_CHECK(o2.isOvoid() == false);

	std::cout << "Distance with SuperEllipsoids: " << ellipsoidTime << " ms" << std::endl;

	// SuperOvoid general case
	g_lastStats.resetToDefault();
	g_lastStats.superellipsoid = false;
	g_lastStatsValid = true;

	timer.start();
	testDistance(o1, o2);
	timer.stop();
	ovoidTime = timer.getElapsedTimeInMicroSec();

	BOOST_CHECK(o1.isOvoid() == true);
	BOOST_CHECK(o2.isOvoid() == true);

	std::cout << "Distance with SuperOvoids:     " << ovoidTime << " ms" << std::endl;

	printf("Selly VS SOV: %f (%s)\n", ellipsoidTime / ovoidTime, (ellipsoidTime / ovoidTime < 1) ? "ok" : "not very good...");

}

BOOST_AUTO_TEST_CASE(implicit_vs_parametric)
{
	SuperOvoid o1(1, 1, 1, 0.6, 0.8, 0.2, 0.1);
	SuperOvoid o2(1, 1, 1, 0.8, 1.1, 0.4, -0.3);

	o1.computeLocalOctree(3);
	o2.computeLocalOctree(3);

	Timer timer;

	NewtonRaphsonStats options;
	options.precision = 1e-6;
	options.maxIterations = 100;
	options.guessQuality = 6;

	options.parametric = false;
	options.guessType = NewtonRaphsonStats::INITIAL_GUESS::OCTREE;

	timer.start();
	testDistance(o1, o2, &options);
	timer.stop();

	std::cout << "Distance with Implicit:   " << timer.getElapsedTimeInMilliSec() << " ms" << std::endl;

	options.parametric = true;
	options.guessType = NewtonRaphsonStats::INITIAL_GUESS::PARAMETRIC_QUADTREE;

	timer.start();
	testDistance(o1, o2, &options);
	timer.stop();

	std::cout << "Distance with Parametric: " << timer.getElapsedTimeInMilliSec() << " ms" << std::endl;

}

BOOST_AUTO_TEST_CASE(global_stats)
{
	std::memset(&g_lastStats, 0xBB, sizeof(NewtonRaphsonStats));
	g_lastStatsValid = false;

	SuperOvoid o1(1, 1, 1, 0.6, 0.8, 0.2, 0.1);
	SuperOvoid o2(1, 1, 1, 0.8, 1.1, 0.4, -0.3);

	{
		// Objects are close, and will collide
		Transform3f t1(Vec3f(0, 0, 0));
		Transform3f t2(Vec3f(0, 1.5, 0));

		BOOST_CHECK(g_lastStatsValid == false);

		bool collide = !details::superOvoidSuperOvoidDistance(o1, t1, o2, t2, NULL, NULL, NULL, true, NULL);

		BOOST_CHECK(collide == true);
		BOOST_CHECK(g_lastStatsValid == true);
		BOOST_CHECK(g_lastStats.iterations > 0);
		BOOST_CHECK(g_lastStats.iterations <= g_lastStats.maxIterations);
		BOOST_CHECK(g_lastStats.iterationsTime >= 0);
		BOOST_CHECK(g_lastStats.numericalJacobianTime >= 0);
	}

	{
		// Objects are far from each other, broad-phase will stop the algorithm
		// before any NewtonRaphson iteration takes place
		Transform3f t1(Vec3f(0, 0, 0));
		Transform3f t2(Vec3f(0, 5, 0));

		bool collide = !details::superOvoidSuperOvoidDistance(o1, t1, o2, t2, NULL, NULL, NULL, true, NULL);

		BOOST_CHECK(collide == false);
		BOOST_CHECK(g_lastStatsValid == true);
		BOOST_CHECK(g_lastStats.iterations < 1e-12);
		BOOST_CHECK(g_lastStats.iterationsTime < 1e-12);
		BOOST_CHECK(g_lastStats.numericalJacobianTime < 1e-12);
	}
}





#include "fcl/BVH/BVH_model.h"
#include "fcl/BV/BV.h"
#include "fcl/BV/OBBRSS.h"
#include "fcl/narrowphase/narrowphase.h"
#include "fcl/distance.h"

typedef BVHModel<OBBRSS> Model;
CollisionObject* getSuperOvoidMeshObject(SuperOvoid& s, Transform3f t, int azSlices, int zeSlices, bool adaptiveMeshes);

//Vec3f start1(0, 0, 0);
//Vec3f start2(2.36, 0.4, -0.3);
//Vec3f goal2( 2.36, 0.4, 0.9);

//// These give absurt 20+ ratios!
//Vec3f start1(0, 0, 0);
//Vec3f start2(2.22, 0.19, -0.74);
//Vec3f goal2(2.22, -0.19, 1.03);

// These give 1~6 ratio results
Vec3f start1(0, 0, 0);
Vec3f start2(2.22, 0.19, -0.3);
Vec3f goal2(0.6, -0.19, 1.8);

float testCollision(CollisionObject* a, CollisionObject* b, int N)
{
	int wasColliding = 0;
	Transform3f t1(Quaternion3f(), start1);
	Transform3f t2(Quaternion3f(), start2);
	a->setTransform(t1);

	//g_lastStatsValid = true;
	//g_lastStats.resetToDefault();
	//g_lastStats.guessType = NewtonRaphsonStats::OCTREE;

	//NewtonRaphsonStats stats;
	//stats = g_lastStats;
	
	for (int i = 0; i < N; i++)
	{
		float t = 1.0f * i / N;
		t2.setTranslation(start2 * (1 - t) + goal2 * t);
		b->setTransform(t2);

		CollisionRequest request;
		request.enable_contact = true;
		CollisionResult result;

		if (collide(a, b, request, result))
			wasColliding++;

		//if (g_lastStatsValid)
		//	stats.add(g_lastStats);
	}

	//stats.averageDivide(N);
	//g_lastStats = stats;

	return 1.0f * wasColliding / N;
}

float testDistance(CollisionObject* a, CollisionObject* b, int N)
{
	int collided = 0;
	Transform3f t1(Quaternion3f(), start1);
	Transform3f t2(Quaternion3f(), start2);
	a->setTransform(t1);

	if (!g_lastStatsValid)
	{
		g_lastStats.resetToDefault();
		g_lastStatsValid = true;
	}

	NewtonRaphsonStats accum;
	accum = g_lastStats;

	for (int i = 0; i < N; i++)
	{
		float t = 1.0f * i / N;
		t2.setTranslation(start2 * (1 - t) + goal2 * t);
		b->setTransform(t2);

		DistanceRequest request;
		request.enable_nearest_points = true;
		DistanceResult result;

		if (distance(a, b, request, result))
			collided++;

		if (g_lastStatsValid)
			accum.add(g_lastStats);
	}

	accum.averageDivide(N);
	g_lastStats = accum;

	return 1.0f * collided / N;
}

CollisionObject* getSuperOvoidMeshObject(SuperOvoid& s, int azSlices, int zeSlices, bool adaptiveMeshes)
{
	std::vector<Vec3f>* vertices = new std::vector<Vec3f>();
	s.storePoints(vertices, azSlices, zeSlices, !adaptiveMeshes);
	std::vector<Triangle>* triangles = new std::vector<Triangle>();
	s.storeTriangles(triangles, azSlices, zeSlices);

	Model* model = new Model();
	model->beginModel();
	model->addSubModel(*vertices, *triangles);
	model->endModel();

	CollisionObject* obj = new CollisionObject(boost::shared_ptr<CollisionGeometry>(model), Transform3f());
	return obj;
}

CollisionObject* getSuperOvoidObject(SuperOvoid* s)
{
	CollisionObject* obj = new CollisionObject(boost::shared_ptr<CollisionGeometry>(s), Transform3f());
	return obj;
}

BOOST_AUTO_TEST_CASE(vs_meshes)
{
	// Superovoid stuff
	SuperOvoid sov1(1, 1, 1, 0.6, 0.8, 0.2, 0.1);
	SuperOvoid sov2(1, 1, 1, 0.8, 1.1, 0.4, -0.3);
	CollisionObject* s1 = getSuperOvoidObject(&sov1);
	CollisionObject* s2 = getSuperOvoidObject(&sov2);

	// Mesh stuff
	const int azSlices = 8;
	const int zeSlices = 9;
	CollisionObject* o1 = getSuperOvoidMeshObject(sov1, azSlices, zeSlices, true);
	CollisionObject* o2 = getSuperOvoidMeshObject(sov2, azSlices, zeSlices, true);

	const int N = 100;
	Timer timer;
	double collisionSov, collisionMesh, distanceSov, distanceMesh;
	double collided = 0;

	// Time superovoid stuff
	timer.start();
	collided += testCollision(s1, s2, N);
	timer.stop();
	collisionSov = timer.getElapsedTimeInMicroSec();

	// Time mesh stuff
	timer.start();
	collided += testCollision(o1, o2, N);
	timer.stop();
	collisionMesh = timer.getElapsedTimeInMicroSec();

	// Time superovoid stuff
	timer.start();
	testDistance(s1, s2, N);
	timer.stop();
	distanceSov = timer.getElapsedTimeInMicroSec();

	// Time mesh stuff
	timer.start();
	testDistance(o1, o2, N);
	timer.stop();
	distanceMesh = timer.getElapsedTimeInMicroSec();

	printf("\n");
	printf("====    Collision and Distance times in ms\n");
	printf("====    Averaged with %d queries (%.0f%% collisions)   =====\n", N, collided / 2.0 * 100.0);
	printf("Col Sov  | %.4f\n", collisionSov / 1000.0);
	printf("Col Mesh | %.4f\n", collisionMesh / 1000.0);
	printf("Dis Sov  | %.4f\n", distanceSov / 1000.0);
	printf("Dis Mesh | %.4f\n", distanceMesh / 1000.0);
	printf("\n");
	printf("Collision speedup: %.2f\n", collisionMesh / collisionSov);
	printf("Distance speedup:  %.2f\n", distanceMesh / distanceSov);
	
}

void compareMeshVsSuperOvoid(SuperOvoid sov1, SuperOvoid sov2,
	float(*function)(CollisionObject*, CollisionObject*, int),
	const char* fileName, int meshIterations, int timeSteps, bool adaptiveMeshes)
{
	// Mesh stuff
	int azSlices = 2;
	int zeSlices = 5;

	// Open file to write results
	std::ofstream file;
	file.precision(6);
	file.open(fileName);
	if (!file.is_open())
	{
		std::cerr << "Could not open file." << std::endl;
		return;
	}

	file << "Vertices" << "," << "Triangles" << "," << "Superovoids" << "," << "Meshes" << std::endl;

	// Create superovoid objects, which are always the same
	CollisionObject* s1 = getSuperOvoidObject(&sov1);
	CollisionObject* s2 = getSuperOvoidObject(&sov2);

	for (int i = 0; i < meshIterations; i++)
	{
		sov1.clearCachedPoints();
		sov2.clearCachedPoints();

		// Create mesh objects
		CollisionObject* o1 = getSuperOvoidMeshObject(sov1, azSlices, zeSlices, adaptiveMeshes);
		CollisionObject* o2 = getSuperOvoidMeshObject(sov2, azSlices, zeSlices, adaptiveMeshes);

		const Model* model = (const Model*)o1->collisionGeometry().get();
		//printf("Mesh verts = %d    tris = %d\n", model->num_vertices, model->num_tris);

		double sovTime, meshTime;
		float wasColliding = 0;

		// Time superovoid stuff
		{
			Timer sovTimer = Timer();
			sovTimer.start();
			wasColliding = function(s1, s2, timeSteps);
			sovTimer.stop();
			sovTime = sovTimer.getElapsedTimeInMicroSec() / timeSteps;

			//printf("Collision with Sov : %04.4f us (%03.1f%% collisions)\n", sovTime, wasColliding * 100);
		}

		// Time mesh stuff
		{
			Timer meshTimer = Timer();
			meshTimer.start();
			wasColliding = function(o1, o2, timeSteps);
			meshTimer.stop();
			meshTime = meshTimer.getElapsedTimeInMicroSec() / timeSteps;

			//printf("Collision with Mesh: %04.4f us (%03.1f%% collisions)\n", meshTime, wasColliding * 100);
		}

		//printf("\t\t\tRatio: %.3f\n", meshTime / sovTime);

		file << model->num_vertices << "," << model->num_tris << "," << sovTime << "," << meshTime << std::endl;

		//std::cout << std::endl;
		azSlices += 2;
		zeSlices += 2;
	}

	file.close();
}

BOOST_AUTO_TEST_CASE(mesh_precision_vs_time)
{
	// Superovoid stuff
	SuperOvoid sov1(1.4, 0.96, 1.2, 0.7, 0.5, 0.2, 0.3);
	SuperOvoid sov2(1.1, 1.2, 1.0, 0.4, 1.1, 0.253, 0.243);

	const int timeSteps = 100;
	const int meshIterations = 15;

	compareMeshVsSuperOvoid(sov1, sov2, testCollision, "mesh_collide.csv", meshIterations, timeSteps, false);
	compareMeshVsSuperOvoid(sov1, sov2, testDistance, "mesh_distance.csv", meshIterations, timeSteps, false);

	printf("\n");

	compareMeshVsSuperOvoid(sov1, sov2, testCollision, "mesh_collide_adaptive.csv", meshIterations, timeSteps, true);
	compareMeshVsSuperOvoid(sov1, sov2, testDistance, "mesh_distance_adaptive.csv", meshIterations, timeSteps, true);
}

BOOST_AUTO_TEST_CASE(superovoid_vs_mesh_time)
{
	// Superovoid stuff
	SuperOvoid sov1(1.4, 0.96, 1.2, 0.7, 0.5, 0.2, 0.3);
	SuperOvoid sov2(1.1, 1.2, 1.0, 0.4, 1.1, 0.253, 0.243);

	const int timeSteps = 500;

	// Create superovoid objects, which are always the same
	CollisionObject* s1 = getSuperOvoidObject(&sov1);
	CollisionObject* s2 = getSuperOvoidObject(&sov2);

	// Open file to write results
	std::ofstream file;
	file.precision(6);
	file.open("superovoid_precision.csv");
	if (!file.is_open())
	{
		std::cerr << "Could not open file." << std::endl;
		return;
	}
	file << "Precision" << "," << "Iterations" << "," << "Time" << std::endl;

	double precision[7] = { 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2 };

	for (int i = 0; i < 7; i++)
	{
		sov1.clearCachedPoints();
		sov2.clearCachedPoints();

		double sovTime;

		g_lastStats.resetToDefault();
		g_lastStats.precision = precision[i];
		g_lastStatsValid = true;

		// Time superovoid stuff
		{
			Timer sovTimer = Timer();
			sovTimer.start();
			testDistance(s1, s2, timeSteps);
			sovTimer.stop();
			sovTime = sovTimer.getElapsedTimeInMicroSec() / timeSteps;

			file << g_lastStats.precision << "," << 1.0 * g_lastStats.iterations << "," << sovTime << std::endl;

			int iterations = (int)g_lastStats.iterations;

			printf("Distance with Sov : %04.4f us (iterations = %.2f, precision = %.10f)\n", sovTime, 1.0 * g_lastStats.iterations, g_lastStats.precision);
		}
	}

	file.close();
}

void compareMeshVsSuperOvoidPrecision(SuperOvoid sov1, SuperOvoid sov2, Transform3f tf1, Transform3f tf2, int meshIterations, bool adaptiveMeshes)
{
	// Mesh stuff
	int azSlices = 2;
	int zeSlices = 5;

	// Repetitions for average times
	int repetitions = 30;

	// Open precisionFile to write results
	std::ofstream precisionFile;
	precisionFile.precision(8);
	precisionFile.open(adaptiveMeshes ? "geom_precision_adaptive.csv" : "geom_precision_equallyspaced.csv");
	if (!precisionFile.is_open())
	{
		std::cerr << "Could not open precisionFile." << std::endl;
		return;
	}

	precisionFile << "Superovoid VS Mesh: geometric precision test on fcl::distance" << std::endl;
	precisionFile << "Vertices" << "," << "Triangles" << "," << "Superovoid A" << "," << "Superovoid B" << "," << "Mesh A" << "," << "Mesh B" << "," << "Error A" << "," << "Error B" << "," << "Time SOV" << "," << "Time Mesh" << std::endl;

	// Create superovoid objects, which are always the same
	CollisionObject* s1 = getSuperOvoidObject(&sov1);
	CollisionObject* s2 = getSuperOvoidObject(&sov2);
	s1->setTransform(tf1);
	s2->setTransform(tf2);

	for (int i = 0; i < meshIterations; i++)
	{
		sov1.clearCachedPoints();
		sov2.clearCachedPoints();

		// Create mesh objects
		CollisionObject* mesh1 = getSuperOvoidMeshObject(sov1, azSlices, zeSlices, adaptiveMeshes);
		CollisionObject* mesh2 = getSuperOvoidMeshObject(sov2, azSlices, zeSlices, adaptiveMeshes);
		mesh1->setTransform(tf1);
		mesh2->setTransform(tf2);

		// Get vertex and triangle count (equal for both meshes)
		const Model* model = (const Model*)mesh1->collisionGeometry().get();
		printf("Mesh verts = %d    tris = %d\n", model->num_vertices, model->num_tris);

		Vec3f sovA, sovB, meshA, meshB;
		double sovTime, meshTime;

		// Test superovoid stuff
		{
			Timer timer; timer.start();
			DistanceRequest request;
			request.enable_nearest_points = true;
			DistanceResult result;
			for (int i = 0; i < repetitions; i++)
				BOOST_CHECK(distance(s1, s2, request, result));
			timer.stop();
			
			sovA = result.nearest_points[0];
			sovB = result.nearest_points[1];
			sovTime = timer.getElapsedTimeInMicroSec() / repetitions;
			
			//std::cout << "Superovoid nearest points: " << sovA << "," << sovB << std::endl;
		}

		// Test mesh stuff
		{
			Timer timer; timer.start();
			DistanceRequest request;
			request.enable_nearest_points = true;
			DistanceResult result;
			for (int i = 0; i < repetitions; i++)
				BOOST_CHECK(distance(mesh1, mesh2, request, result));
			timer.stop();

			meshA = result.nearest_points[0];
			meshB = result.nearest_points[1];
			meshTime = timer.getElapsedTimeInMicroSec() / repetitions;

			//std::cout << "Superovoid nearest points: " << meshA << "," << meshB << std::endl;
		}

		double errorA = (sovA - meshA).length();
		double errorB = (sovB - meshB).length();
		//printf("Diff: %f - %f\n", errorA, errorB);

		precisionFile << model->num_vertices << "," << model->num_tris << "," << sovA << "," << sovB << "," << meshA << "," << meshB << "," << errorA << "," << errorB << "," << sovTime << "," << meshTime << std::endl;

		// Test with a more detailed mesh the next iteration
		azSlices += 2;
		zeSlices += 2;
	}

	precisionFile.close();
}

BOOST_AUTO_TEST_CASE(superovoid_vs_mesh_precision)
{
	// Superovoid stuff
	SuperOvoid sov1(1.4, 0.96, 1.2, 0.7, 0.5, 0.2, 0.3);
	SuperOvoid sov2(1.1, 1.2, 1.0, 0.4, 1.1, 0.253, 0.243);

	Transform3f tf1(Quaternion3f(), Vec3f(0, 0, 0));
	Transform3f tf2(Quaternion3f(), Vec3f(2.59, -0.3, -0.4));

	std::cout.precision(3);
	compareMeshVsSuperOvoidPrecision(sov1, sov2, tf1, tf2, 10, true);
	compareMeshVsSuperOvoidPrecision(sov1, sov2, tf1, tf2, 10, false);
}

BOOST_AUTO_TEST_CASE(octree_guess)
{
    // iCub thumb tip superellipsoids
    SuperOvoid s1(0.0064423, 0.016188, 0.0064423,   1, 0.64,   0, 0);
    SuperOvoid s2(0.0064423, 0.016188, 0.0064423,   1, 0.64,   0, 0);

    Transform3f t1(Quaternion3f(-0.19759836, -0.59332478, 0.75545207, 0.19548083), Vec3f(-0.23568268, -0.022810882, -0.029985857));
    Transform3f t2(Quaternion3f(0.19547792, 0.75545089, -0.59332602, -0.19760205), Vec3f(-0.22636475, -0.011973577, -0.029985381));

    DistanceRequest request;
    request.enable_nearest_points = true;
    DistanceResult result;

    g_lastStatsValid = true;
    g_lastStats.resetToDefault();
    g_lastStats.guessType = fcl::NewtonRaphsonStats::OCTREE;
    g_lastStats.guessQuality = 8;
    

    distance(&s1, t1, &s2, t2, request, result);

    BOOST_CHECK(g_lastStatsValid);
    g_lastStats.print(std::cout);

    BOOST_CHECK(g_lastStats.discardedGuess == 0);
    BOOST_CHECK(g_lastStats.didNotConverge == 0);

    BOOST_CHECK(s1.isCachedPointValid(&s2));
    BOOST_CHECK(s2.isCachedPointValid(&s1));
    
    distance(&s1, t1, &s2, t2, request, result);
    
    BOOST_CHECK(g_lastStats.discardedGuess == 0);
    BOOST_CHECK(g_lastStats.didNotConverge == 0);
    g_lastStats.print(std::cout);
    BOOST_CHECK(g_lastStats.iterations == 1);
}

BOOST_AUTO_TEST_CASE(mesh_small_memtest)
{
	// Superovoid stuff
	SuperOvoid sov1(1.4, 0.96, 1.2, 0.7, 0.5, 0.2, 0.3);
	SuperOvoid sov2(1.1, 1.2, 1.0, 0.4, 1.1, 0.253, 0.243);

	int azSlices = 4;
	int zeSlices = 7;

	// Create superovoid objects, which are always the same
	CollisionObject* s1 = getSuperOvoidMeshObject(sov1, azSlices, zeSlices, true);
	CollisionObject* s2 = getSuperOvoidMeshObject(sov2, azSlices, zeSlices, true);

	testDistance(s1, s2, 500);
}

BOOST_AUTO_TEST_CASE(mesh_large_memtest)
{
	// Superovoid stuff
	SuperOvoid sov1(1.4, 0.96, 1.2, 0.7, 0.5, 0.2, 0.3);
	SuperOvoid sov2(1.1, 1.2, 1.0, 0.4, 1.1, 0.253, 0.243);

	int azSlices = 16;
	int zeSlices = 21;

	// Create superovoid objects, which are always the same
	CollisionObject* s1 = getSuperOvoidMeshObject(sov1, azSlices, zeSlices, true);
	CollisionObject* s2 = getSuperOvoidMeshObject(sov2, azSlices, zeSlices, true);

	testDistance(s1, s2, 500);
}

std::vector<Vec3f>* getTestMesh()
{
	// Generated in Unity implementation
	// SuperOvoid (1, 1, 1, 0.64, 0.64, -0.28, -0.28);
	// Azimuth slicess = 9
	// Zenith slices = 21
	// Adaptive (not equally spaced)

	std::vector<Vec3f>* vertices = new std::vector<Vec3f>();

	vertices->push_back(Vec3f(0.0000, 0.0000, 1.0000));
	vertices->push_back(Vec3f(0.0000, 0.3437, 0.9684));
	vertices->push_back(Vec3f(0.1730, 0.3303, 0.9684));
	vertices->push_back(Vec3f(0.2591, 0.2898, 0.9684));
	vertices->push_back(Vec3f(0.3135, 0.2206, 0.9684));
	vertices->push_back(Vec3f(0.3404, 0.1121, 0.9684));
	vertices->push_back(Vec3f(0.3404, -0.1121, 0.9684));
	vertices->push_back(Vec3f(0.3135, -0.2206, 0.9684));
	vertices->push_back(Vec3f(0.2591, -0.2898, 0.9684));
	vertices->push_back(Vec3f(0.1730, -0.3303, 0.9684));
	vertices->push_back(Vec3f(0.0000, -0.3437, 0.9684));
	vertices->push_back(Vec3f(-0.1730, -0.3303, 0.9684));
	vertices->push_back(Vec3f(-0.2591, -0.2898, 0.9684));
	vertices->push_back(Vec3f(-0.3135, -0.2206, 0.9684));
	vertices->push_back(Vec3f(-0.3404, -0.1121, 0.9684));
	vertices->push_back(Vec3f(-0.3404, 0.1121, 0.9684));
	vertices->push_back(Vec3f(-0.3135, 0.2206, 0.9684));
	vertices->push_back(Vec3f(-0.2591, 0.2898, 0.9684));
	vertices->push_back(Vec3f(-0.1730, 0.3303, 0.9684));
	vertices->push_back(Vec3f(0.0000, 0.5377, 0.8732));
	vertices->push_back(Vec3f(0.2706, 0.5167, 0.8732));
	vertices->push_back(Vec3f(0.4052, 0.4534, 0.8732));
	vertices->push_back(Vec3f(0.4904, 0.3451, 0.8732));
	vertices->push_back(Vec3f(0.5325, 0.1754, 0.8732));
	vertices->push_back(Vec3f(0.5325, -0.1754, 0.8732));
	vertices->push_back(Vec3f(0.4904, -0.3451, 0.8732));
	vertices->push_back(Vec3f(0.4052, -0.4534, 0.8732));
	vertices->push_back(Vec3f(0.2706, -0.5167, 0.8732));
	vertices->push_back(Vec3f(0.0000, -0.5377, 0.8732));
	vertices->push_back(Vec3f(-0.2706, -0.5167, 0.8732));
	vertices->push_back(Vec3f(-0.4052, -0.4534, 0.8732));
	vertices->push_back(Vec3f(-0.4904, -0.3451, 0.8732));
	vertices->push_back(Vec3f(-0.5325, -0.1754, 0.8732));
	vertices->push_back(Vec3f(-0.5325, 0.1754, 0.8732));
	vertices->push_back(Vec3f(-0.4904, 0.3451, 0.8732));
	vertices->push_back(Vec3f(-0.4052, 0.4534, 0.8732));
	vertices->push_back(Vec3f(-0.2706, 0.5167, 0.8732));
	vertices->push_back(Vec3f(0.0000, 0.6992, 0.7117));
	vertices->push_back(Vec3f(0.3519, 0.6719, 0.7117));
	vertices->push_back(Vec3f(0.5269, 0.5895, 0.7117));
	vertices->push_back(Vec3f(0.6377, 0.4487, 0.7117));
	vertices->push_back(Vec3f(0.6923, 0.2280, 0.7117));
	vertices->push_back(Vec3f(0.6923, -0.2280, 0.7117));
	vertices->push_back(Vec3f(0.6377, -0.4487, 0.7117));
	vertices->push_back(Vec3f(0.5269, -0.5895, 0.7117));
	vertices->push_back(Vec3f(0.3519, -0.6719, 0.7117));
	vertices->push_back(Vec3f(0.0000, -0.6992, 0.7117));
	vertices->push_back(Vec3f(-0.3519, -0.6719, 0.7117));
	vertices->push_back(Vec3f(-0.5269, -0.5895, 0.7117));
	vertices->push_back(Vec3f(-0.6377, -0.4487, 0.7117));
	vertices->push_back(Vec3f(-0.6923, -0.2280, 0.7117));
	vertices->push_back(Vec3f(-0.6923, 0.2280, 0.7117));
	vertices->push_back(Vec3f(-0.6377, 0.4487, 0.7117));
	vertices->push_back(Vec3f(-0.5269, 0.5895, 0.7117));
	vertices->push_back(Vec3f(-0.3519, 0.6719, 0.7117));
	vertices->push_back(Vec3f(0.0000, 0.8405, 0.4716));
	vertices->push_back(Vec3f(0.4230, 0.8077, 0.4716));
	vertices->push_back(Vec3f(0.6334, 0.7087, 0.4716));
	vertices->push_back(Vec3f(0.7666, 0.5394, 0.4716));
	vertices->push_back(Vec3f(0.8323, 0.2741, 0.4716));
	vertices->push_back(Vec3f(0.8323, -0.2741, 0.4716));
	vertices->push_back(Vec3f(0.7666, -0.5394, 0.4716));
	vertices->push_back(Vec3f(0.6334, -0.7087, 0.4716));
	vertices->push_back(Vec3f(0.4230, -0.8077, 0.4716));
	vertices->push_back(Vec3f(0.0000, -0.8405, 0.4716));
	vertices->push_back(Vec3f(-0.4230, -0.8077, 0.4716));
	vertices->push_back(Vec3f(-0.6334, -0.7087, 0.4716));
	vertices->push_back(Vec3f(-0.7666, -0.5394, 0.4716));
	vertices->push_back(Vec3f(-0.8323, -0.2741, 0.4716));
	vertices->push_back(Vec3f(-0.8323, 0.2741, 0.4716));
	vertices->push_back(Vec3f(-0.7666, 0.5394, 0.4716));
	vertices->push_back(Vec3f(-0.6334, 0.7087, 0.4716));
	vertices->push_back(Vec3f(-0.4230, 0.8077, 0.4716));
	vertices->push_back(Vec3f(0.0000, 1.0000, 0.0000));
	vertices->push_back(Vec3f(0.5033, 0.9610, 0.0000));
	vertices->push_back(Vec3f(0.7536, 0.8432, 0.0000));
	vertices->push_back(Vec3f(0.9121, 0.6417, 0.0000));
	vertices->push_back(Vec3f(0.9903, 0.3261, 0.0000));
	vertices->push_back(Vec3f(0.9903, -0.3261, 0.0000));
	vertices->push_back(Vec3f(0.9121, -0.6417, 0.0000));
	vertices->push_back(Vec3f(0.7536, -0.8432, 0.0000));
	vertices->push_back(Vec3f(0.5033, -0.9610, 0.0000));
	vertices->push_back(Vec3f(0.0000, -1.0000, 0.0000));
	vertices->push_back(Vec3f(-0.5033, -0.9610, 0.0000));
	vertices->push_back(Vec3f(-0.7536, -0.8432, 0.0000));
	vertices->push_back(Vec3f(-0.9121, -0.6417, 0.0000));
	vertices->push_back(Vec3f(-0.9903, -0.3261, 0.0000));
	vertices->push_back(Vec3f(-0.9903, 0.3261, 0.0000));
	vertices->push_back(Vec3f(-0.9121, 0.6417, 0.0000));
	vertices->push_back(Vec3f(-0.7536, 0.8432, 0.0000));
	vertices->push_back(Vec3f(-0.5033, 0.9610, 0.0000));
	vertices->push_back(Vec3f(0.0000, 1.0963, -0.4716));
	vertices->push_back(Vec3f(0.5517, 1.0535, -0.4716));
	vertices->push_back(Vec3f(0.8262, 0.9244, -0.4716));
	vertices->push_back(Vec3f(0.9999, 0.7035, -0.4716));
	vertices->push_back(Vec3f(1.0856, 0.3575, -0.4716));
	vertices->push_back(Vec3f(1.0856, -0.3575, -0.4716));
	vertices->push_back(Vec3f(0.9999, -0.7035, -0.4716));
	vertices->push_back(Vec3f(0.8262, -0.9244, -0.4716));
	vertices->push_back(Vec3f(0.5517, -1.0535, -0.4716));
	vertices->push_back(Vec3f(0.0000, -1.0963, -0.4716));
	vertices->push_back(Vec3f(-0.5517, -1.0535, -0.4716));
	vertices->push_back(Vec3f(-0.8262, -0.9244, -0.4716));
	vertices->push_back(Vec3f(-0.9999, -0.7035, -0.4716));
	vertices->push_back(Vec3f(-1.0856, -0.3575, -0.4716));
	vertices->push_back(Vec3f(-1.0856, 0.3575, -0.4716));
	vertices->push_back(Vec3f(-0.9999, 0.7035, -0.4716));
	vertices->push_back(Vec3f(-0.8262, 0.9244, -0.4716));
	vertices->push_back(Vec3f(-0.5517, 1.0535, -0.4716));
	vertices->push_back(Vec3f(0.0000, 1.0472, -0.7117));
	vertices->push_back(Vec3f(0.5270, 1.0063, -0.7117));
	vertices->push_back(Vec3f(0.7892, 0.8829, -0.7117));
	vertices->push_back(Vec3f(0.9551, 0.6720, -0.7117));
	vertices->push_back(Vec3f(1.0369, 0.3415, -0.7117));
	vertices->push_back(Vec3f(1.0369, -0.3415, -0.7117));
	vertices->push_back(Vec3f(0.9551, -0.6720, -0.7117));
	vertices->push_back(Vec3f(0.7892, -0.8829, -0.7117));
	vertices->push_back(Vec3f(0.5270, -1.0063, -0.7117));
	vertices->push_back(Vec3f(0.0000, -1.0472, -0.7117));
	vertices->push_back(Vec3f(-0.5270, -1.0063, -0.7117));
	vertices->push_back(Vec3f(-0.7892, -0.8829, -0.7117));
	vertices->push_back(Vec3f(-0.9551, -0.6720, -0.7117));
	vertices->push_back(Vec3f(-1.0369, -0.3415, -0.7117));
	vertices->push_back(Vec3f(-1.0369, 0.3415, -0.7117));
	vertices->push_back(Vec3f(-0.9551, 0.6720, -0.7117));
	vertices->push_back(Vec3f(-0.7892, 0.8829, -0.7117));
	vertices->push_back(Vec3f(-0.5270, 1.0063, -0.7117));
	vertices->push_back(Vec3f(0.0000, 0.8857, -0.8732));
	vertices->push_back(Vec3f(0.4457, 0.8511, -0.8732));
	vertices->push_back(Vec3f(0.6675, 0.7468, -0.8732));
	vertices->push_back(Vec3f(0.8078, 0.5684, -0.8732));
	vertices->push_back(Vec3f(0.8771, 0.2889, -0.8732));
	vertices->push_back(Vec3f(0.8771, -0.2889, -0.8732));
	vertices->push_back(Vec3f(0.8078, -0.5684, -0.8732));
	vertices->push_back(Vec3f(0.6675, -0.7468, -0.8732));
	vertices->push_back(Vec3f(0.4457, -0.8511, -0.8732));
	vertices->push_back(Vec3f(0.0000, -0.8857, -0.8732));
	vertices->push_back(Vec3f(-0.4457, -0.8511, -0.8732));
	vertices->push_back(Vec3f(-0.6675, -0.7468, -0.8732));
	vertices->push_back(Vec3f(-0.8078, -0.5684, -0.8732));
	vertices->push_back(Vec3f(-0.8771, -0.2889, -0.8732));
	vertices->push_back(Vec3f(-0.8771, 0.2889, -0.8732));
	vertices->push_back(Vec3f(-0.8078, 0.5684, -0.8732));
	vertices->push_back(Vec3f(-0.6675, 0.7468, -0.8732));
	vertices->push_back(Vec3f(-0.4457, 0.8511, -0.8732));
	vertices->push_back(Vec3f(0.0000, 0.5995, -0.9684));
	vertices->push_back(Vec3f(0.3017, 0.5761, -0.9684));
	vertices->push_back(Vec3f(0.4518, 0.5055, -0.9684));
	vertices->push_back(Vec3f(0.5468, 0.3847, -0.9684));
	vertices->push_back(Vec3f(0.5937, 0.1955, -0.9684));
	vertices->push_back(Vec3f(0.5937, -0.1955, -0.9684));
	vertices->push_back(Vec3f(0.5468, -0.3847, -0.9684));
	vertices->push_back(Vec3f(0.4518, -0.5055, -0.9684));
	vertices->push_back(Vec3f(0.3017, -0.5761, -0.9684));
	vertices->push_back(Vec3f(0.0000, -0.5995, -0.9684));
	vertices->push_back(Vec3f(-0.3017, -0.5761, -0.9684));
	vertices->push_back(Vec3f(-0.4518, -0.5055, -0.9684));
	vertices->push_back(Vec3f(-0.5468, -0.3847, -0.9684));
	vertices->push_back(Vec3f(-0.5937, -0.1955, -0.9684));
	vertices->push_back(Vec3f(-0.5937, 0.1955, -0.9684));
	vertices->push_back(Vec3f(-0.5468, 0.3847, -0.9684));
	vertices->push_back(Vec3f(-0.4518, 0.5055, -0.9684));
	vertices->push_back(Vec3f(-0.3017, 0.5761, -0.9684));
	vertices->push_back(Vec3f(0.0000, 0.0000, -1.0000));

	return vertices;
}

BOOST_AUTO_TEST_CASE(mesh_vs_unity)
{
	SuperOvoid sov(1, 1, 1, 0.64, 0.64, -0.28, -0.28);
	std::vector<Vec3f> fclVertices;
	sov.storePoints(&fclVertices, 9, 21, false);

	std::vector<Vec3f>* unityVertices = getTestMesh();

	// Must have same number of vertices
	BOOST_CHECK(fclVertices.size() == unityVertices->size());

	if (fclVertices.size() == unityVertices->size())
	{
		// Vertices must have the same values
		int wrongVertices = 0;
		for (unsigned int i = 0; i < fclVertices.size(); i++)
		{
			Vec3f unity = (*unityVertices)[i];
			Vec3f fcl = fclVertices[i];

			Vec3f diff = unity - fcl;
			if (diff.length() > 0.0002)
				wrongVertices++;
		}

		BOOST_CHECK(wrongVertices == 0);
	}

	delete unityVertices;
}







/** State for the simple RNG */
int randState;

/** Predictable pseudo-random function to always generate the same output, regardless of running platform. */
int myRand()
{
	int const a = 1103515245;
	int const c = 12345;
	randState = a * randState + c;
	return (randState >> 16) & 0x7FFF;
}

void mySeedRand(int seed)
{
	randState = seed;
}

float getRandomRange(FCL_REAL min, FCL_REAL max)
{
	return min + (myRand() % 1024 / 1024.0) * (max - min);
}

float getNormalSample()
{
	const int levels = 10240;
	const int samples = 16;

	FCL_REAL sum = 0;
	for (int i = 0; i < samples; i++)
		sum += (myRand() % levels) * 1.0 / levels;

	return sum / samples - 0.5;
}

void setRandomRotation(CollisionObject* obj)
{
	FCL_REAL quat[4];
	for (int i = 0; i < 4; i++)
		quat[i] = getNormalSample();

	FCL_REAL length = 0;
	for (int i = 0; i < 4; i++)
		length += quat[i] * quat[i];
	length = std::sqrt(length);

	for (int i = 0; i < 4; i++)
		quat[i] /= length;

	Quaternion3f rotation(quat[0], quat[1], quat[2], quat[3]);
	obj->setQuatRotation(rotation);
}


void randomizeSuperovoid(SuperOvoid& sov)
{
	sov.a1 = 1.0;
	sov.a2 = 1.0;
	sov.a3 = 1.0;

	sov.epsilon1 = 1;
	sov.epsilon2 = 1;

	sov.taperingX = 0;
	sov.taperingY = 0;

	//sov.a1 = getRandomRange(0.5, 1.5);
	//sov.a2 = getRandomRange(0.5, 1.5);
	//sov.a3 = getRandomRange(0.5, 1.5);

	//sov.epsilon1 = getRandomRange(0.3, 1.1);
	//sov.epsilon2 = getRandomRange(0.3, 1.1);

	//sov.taperingX = getRandomRange(-0.4, 0.4);
	//sov.taperingY = getRandomRange(-0.4, 0.4);
}


bool isGuessEqual(NewtonRaphsonStats stats1, NewtonRaphsonStats stats2)
{
	return 
		stats1.getGuessA().equal(stats2.getGuessA(), 1e-12)
		&& stats1.getGuessB().equal(stats2.getGuessB(), 1e-12);
}

// Testing the implicit and parametric versions,
// starting from the same initial guess, obtained
// with the ParametricQuadtree method
BOOST_AUTO_TEST_CASE(mubo_implicit_vs_parametric_same_guess)
{
	std::ofstream fileParametric;
	fileParametric.open("mubo_parametric.csv");
	std::ofstream fileImplicit;
	fileImplicit.open("mubo_implicit.csv");
	
	// Write header
	fileParametric << "iteration,";
	fileImplicit << "iteration,";
	NewtonRaphsonStats::printHeader(fileParametric);
	NewtonRaphsonStats::printHeader(fileImplicit);
	fileParametric << std::endl;
	fileImplicit << std::endl;

	// Superovoid stuff
	SuperOvoid sov1(1.4, 0.96, 1.2, 0.7, 0.5, 0.2, 0.3);
	SuperOvoid sov2(1.1, 1.2, 1.0, 0.4, 1.1, 0.253, 0.243);

	// Create superovoid objects, which are always the same
	CollisionObject* s1 = getSuperOvoidObject(&sov1);
	CollisionObject* s2 = getSuperOvoidObject(&sov2);

	s1->setTranslation(Vec3f());
	s2->setTranslation(Vec3f(4, 0, 0));
	
	mySeedRand(1337);

	for (int i = 0; i < 5; i++)
	{
		randomizeSuperovoid(sov1);
		randomizeSuperovoid(sov2);
		setRandomRotation(s1);
		setRandomRotation(s2);

		g_lastStats.resetToDefault();
		g_lastStatsValid = true;
		
		// Use the same initial guess for both experiments
		g_lastStats.guessType = NewtonRaphsonStats::PARAMETRIC_QUADTREE;

		DistanceRequest request;
		DistanceResult parametricResult, implicitResult;
		NewtonRaphsonStats parametricStats, implicitStats;

		// Experiment 1: parametric
		sov1.clearCachedPoints();
		sov2.clearCachedPoints();
		request = DistanceRequest(true);
		g_lastStats.parametric = true;
		distance(s1, s2, request, parametricResult);
		parametricStats = g_lastStats;

		// Experiment 2: implicit
		sov1.clearCachedPoints();
		sov2.clearCachedPoints();
		request = DistanceRequest(true);
		g_lastStats.parametric = false;
		distance(s1, s2, request, implicitResult);
		implicitStats = g_lastStats;

		// Write results
		fileImplicit << i << "," << implicitStats << std::endl;
		fileParametric << i << "," << parametricStats << std::endl;

		//BOOST_CHECK(isGuessEqual(parametricStats, implicitStats));
		//BOOST_CHECK(std::abs(parametricResult.min_distance - implicitResult.min_distance) < 1e-4);
	}

	fileParametric.close();
	fileImplicit.close();
}