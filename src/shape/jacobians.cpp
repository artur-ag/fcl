/**
 * Implements the Newton-Raphson numerical method to calculate
 * the minimum distance between two superovoid surfaces.
 *
 * Author: Artur Goncalves
 *
 * Based on MDC-ELLIPSOIDs (Daniel Simoes Lopes)
 * http://web.ist.utl.pt/ist151462/mdc-ellipsoids.html
 */

#include "fcl/shape/SuperOvoid.h"
#include "fcl/shape/SuperOvoidDetails.h"

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
	}	
}