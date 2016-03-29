#ifndef FCL_SUPEROVOID_DETAILS_H
#define FCL_SUPEROVOID_DETAILS_H

#include "fcl/shape/geometric_shapes.h"
#include "fcl/SuperOvoid_global.h"
#include "fcl/math/vec_3f.h"

namespace fcl
{
	namespace details
	{
		bool superOvoidSuperOvoidDistance(const SuperOvoid& s1, const Transform3f& tf1,
			const SuperOvoid& s2, const Transform3f& tf2,
			FCL_REAL* dist, Vec3f* p1, Vec3f* p2,
			bool collisionQuery,
			NewtonRaphsonStats* stats);

		bool isNaN(FCL_REAL* vector, int size);

		bool getNumericalJacobian(
			int size,
			const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2,
			void(*function)(const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, FCL_REAL* qk, FCL_REAL phi[6]),
			FCL_REAL* qk, FCL_REAL* jacobian);

		bool getAnalyticalParametricJacobian(
			const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2,
			FCL_REAL* qk, FCL_REAL* jacobian);
	}
}

#endif