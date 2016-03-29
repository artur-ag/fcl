/**
 * Implements analytical Jacobian matrices for the numerical
 * Newton-Raphson method.
 *
 * Authors: Artur Goncalves, Daniel Simoes Lopes
 */

#include "fcl/narrowphase/narrowphase.h"
#include "fcl/distance.h"
#include "fcl/collision.h"
#include "fcl/BVH/BVH_model.h"
#include "fcl/BV/BV.h"
#include "fcl/BV/OBBRSS.h"
#include "fcl/shape/SuperOvoid.h"
#include "fcl/shape/SuperOvoidDetails.h"

#include "fcl/SuperOvoid_global.h" // For NewtonRaphsonStats and timers

#include "lapacke.h"

namespace fcl
{
    namespace details
    {
		bool getNumericalJacobian(
			int size,
			const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2,
			void(*function)(const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, FCL_REAL* qk, FCL_REAL phi[6]),
			FCL_REAL* qk, FCL_REAL* jacobian)
		{
			assert(size <= 6);

			FCL_REAL phi[6];
			FCL_REAL qkPerturb[6];
			FCL_REAL phiPerturb[6];

			// Evaluate function at qk, store result in phi
			function(s1, tf1, s2, tf2, qk, phi);

			// Clear array to store Jacobian
			for (int i = 0; i < size * size; i++)
				jacobian[i] = 0;

			FCL_REAL perturb = 1e-6;
			for (int col = 0; col < size; col++)
			{
				for (int i = 0; i < size; i++)
					qkPerturb[i] = qk[i];
				qkPerturb[col] = qkPerturb[col] + perturb;
				function(s1, tf1, s2, tf2, qkPerturb, phiPerturb);

				for (int i = 0; i < size; i++)
				{
					jacobian[col * size + i] = (phiPerturb[i] - phi[i]) / perturb;
				}
			}

			// Return false if any entry in the Jaconian is NaN
			return !isNaN(jacobian, size * size);
		}

		bool getAnalyticalParametricJacobian(
			const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2,
			FCL_REAL* qk, FCL_REAL* jacobian)
		{
			// qk is of the form [p1i, p2i, p1j, p2j]
			FCL_REAL p1i = qk[0], p2i = qk[1], p1j = qk[2], p2j = qk[3];

			// Rotations of each surface
			Quaternion3f tf1Rotation = tf1.getQuatRotation();
			Quaternion3f tf2Rotation = tf2.getQuatRotation();

			// Normals, tangents, binormals, and respective derivatives, in global coordinates
			Vec3f n_OP =       tf1Rotation.transform(s1.getNormal(p1i, p2i));
			Vec3f dN_OP_dp1i = tf1Rotation.transform(s1.getNormalDerivativePhi1(p1i, p2i));
			Vec3f dN_OP_dp2i = tf1Rotation.transform(s1.getNormalDerivativePhi2(p1i, p2i));

			Vec3f t_OQ =       tf2Rotation.transform(s2.getAzimuthTangent(p1j, p2j));
			Vec3f dt_OQ_dp1j = tf2Rotation.transform(s2.getAzimuthTangentDerivativePhi1(p1j, p2j));
			Vec3f dt_OQ_dp2j = tf2Rotation.transform(s2.getAzimuthTangentDerivativePhi2(p1j, p2j));
			
			Vec3f b_OQ =       tf2Rotation.transform(s2.getZenithTangent(p1j, p2j));
			Vec3f db_OQ_dp1j = tf2Rotation.transform(s2.getZenithTangentDerivativePhi1(p1j, p2j));
			Vec3f db_OQ_dp2j = tf2Rotation.transform(s2.getZenithTangentDerivativePhi2(p1j, p2j));

			Vec3f d_PQ =        tf2.transform(s2.getPoint(p1j, p2j)) - tf1.transform(s1.getPoint(p1i, p2i));
			Vec3f d_d_PQ_dp1i = tf1Rotation.transform(-s1.getAzimuthTangent(p1i, p2i));
			Vec3f d_d_PQ_dp2i = tf1Rotation.transform(-s1.getZenithTangent(p1i, p2i));
			Vec3f d_d_PQ_dp1j = tf2Rotation.transform(s2.getAzimuthTangent(p1j, p2j));
			Vec3f d_d_PQ_dp2j = tf2Rotation.transform(s2.getZenithTangent(p1j, p2j));

			// 16 entries of jacobian matrix
			FCL_REAL dColumn1_dp1i = dN_OP_dp1i.dot(t_OQ);							// d(n_OP . t_OQ) / d p1i
			FCL_REAL dColumn2_dp1i = dN_OP_dp1i.dot(b_OQ);							// d(n_OP . b_OQ) / d p1i
			FCL_REAL dColumn3_dp1i = d_d_PQ_dp1i.dot(t_OQ);							// d(d_PQ . t_OQ) / d p1i
			FCL_REAL dColumn4_dp1i = d_d_PQ_dp1i.dot(b_OQ);							// d(d_PQ . b_OQ) / d p1i

			FCL_REAL dColumn1_dp2i = dN_OP_dp2i.dot(t_OQ);							// d(n_OP . t_OQ) / d p2i
			FCL_REAL dColumn2_dp2i = dN_OP_dp2i.dot(b_OQ);							// d(n_OP . b_OQ) / d p2i
			FCL_REAL dColumn3_dp2i = d_d_PQ_dp2i.dot(t_OQ);							// d(d_PQ . t_OQ) / d p2i
			FCL_REAL dColumn4_dp2i = d_d_PQ_dp2i.dot(b_OQ);							// d(d_PQ . b_OQ) / d p2i

			FCL_REAL dColumn1_dp1j = n_OP.dot(dt_OQ_dp1j);							// d(n_OP . t_OQ) / d p1j
			FCL_REAL dColumn2_dp1j = n_OP.dot(db_OQ_dp1j);							// d(n_OP . b_OQ) / d p1j
			FCL_REAL dColumn3_dp1j = d_d_PQ_dp1j.dot(t_OQ) + d_PQ.dot(dt_OQ_dp1j);	// d(d_PQ . t_OQ) / d p1j
			FCL_REAL dColumn4_dp1j = d_d_PQ_dp1j.dot(b_OQ) + d_PQ.dot(db_OQ_dp1j);	// d(d_PQ . b_OQ) / d p1j

			FCL_REAL dColumn1_dp2j = n_OP.dot(dt_OQ_dp2j);							// d(n_OP . t_OQ) / d p2j
			FCL_REAL dColumn2_dp2j = n_OP.dot(db_OQ_dp2j);							// d(n_OP . b_OQ) / d p2j
			FCL_REAL dColumn3_dp2j = d_d_PQ_dp2j.dot(t_OQ) + d_PQ.dot(dt_OQ_dp2j);	// d(d_PQ . t_OQ) / d p2j
			FCL_REAL dColumn4_dp2j = d_d_PQ_dp2j.dot(b_OQ) + d_PQ.dot(db_OQ_dp2j);	// d(d_PQ . b_OQ) / d p2j

			int i = 0;
			jacobian[i++] = dColumn1_dp1i;
			jacobian[i++] = dColumn2_dp1i;
			jacobian[i++] = dColumn3_dp1i;
			jacobian[i++] = dColumn4_dp1i;

			jacobian[i++] = dColumn1_dp2i;
			jacobian[i++] = dColumn2_dp2i;
			jacobian[i++] = dColumn3_dp2i;
			jacobian[i++] = dColumn4_dp2i;
			
			jacobian[i++] = dColumn1_dp1j;
			jacobian[i++] = dColumn2_dp1j;
			jacobian[i++] = dColumn3_dp1j;
			jacobian[i++] = dColumn4_dp1j;
			
			jacobian[i++] = dColumn1_dp2j;
			jacobian[i++] = dColumn2_dp2j;
			jacobian[i++] = dColumn3_dp2j;
			jacobian[i++] = dColumn4_dp2j;

			// Return false if any entry in the Jaconian is NaN
			return !isNaN(jacobian, 4 * 4);
		}

		// @brief Power a to b, b being a rational number, and always get a real result
		inline FCL_REAL spow(FCL_REAL a, FCL_REAL b)
		{
			return a * std::pow(std::pow(a, 2), (b - 1) / 2.0);
		}


    }
}
