#define BOOST_TEST_MODULE "FCL_SUPEROVOID_CORRECTNESS"
#include <boost/test/unit_test.hpp>

#include "fcl/collision.h"
#include "fcl/shape/geometric_shapes.h"
#include "fcl/shape/SuperOvoid.h"
#include "fcl/shape/SuperOvoidDetails.h"

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

			// Using azimuth and zenith directly, for normals
			Vec3f ovoidNormal = ovoid.getNormal(azimuth, zenith, false);
			Vec3f ellipsoidNormal = ellipsoid.getNormal(azimuth, zenith, false);
			BOOST_CHECK(ovoidNormal.equal(ellipsoidNormal, 1e-9));

			// Using a 3D point, for normals
			Vec3f point = ovoid.getPoint(azimuth, zenith, false);
			ovoidNormal = ovoid.getNormal(point);
			ellipsoidNormal = ellipsoid.getNormal(point);
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

