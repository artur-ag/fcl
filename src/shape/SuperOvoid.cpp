/**
 * Implements the Newton-Raphson numerical method to calculate
 * the minimum distance between two superovoid surfaces.
 *
 * Author: Artur Goncalves
 *
 * Based on MDC-ELLIPSOIDs (Daniel Simoes Lopes)
 * http://web.ist.utl.pt/ist151462/mdc-ellipsoids.html
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

// #define FCL_SUPEROVOID_DEBUG_LOG 0
/* (additive)
0	Nothing
1	number of iterations
2	qk, phi, jacobian, delta_qk of each iteration
3	intermediate results of geometrical constraints
*/

/* Auxiliary routine: printing a matrix */
void print_matrix(char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda) {
    lapack_int i, j;
    printf("\n%s\n", desc);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
            printf("%.10f\n", a[j*lda + i]);
    }
    printf("\n");
}

/* Auxiliary routine: printing a vector of integers */
void print_vector(char* desc, lapack_int n, double* a) {
    lapack_int j;
    printf("\n%s\n", desc);
    for (j = 0; j < n; j++)
        printf("%.10f\n", a[j]);
    printf("\n");
}



namespace fcl
{
    namespace details
    {
        // Forward declarations
        void getInitialGuess          (const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, FCL_REAL qk[6], NewtonRaphsonStats* stats, bool parametricOutput);
        void getAvgSphereInitialGuess (const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, FCL_REAL qk[6], NewtonRaphsonStats* stats);
        void getParametricInitialGuess(const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, FCL_REAL qk[6], NewtonRaphsonStats* stats, bool parametricOutput);
        void getOctreeInitialGuess    (const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, FCL_REAL qk[6], NewtonRaphsonStats* stats);
        void getOBBInitialGuess       (const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, FCL_REAL qk[6], NewtonRaphsonStats* stats);
        void getMeshInitialGuess      (const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, FCL_REAL qk[6], NewtonRaphsonStats* stats);

        void computeTangentsWithHouseholder(const Vec3f normal, Vec3f* out_householder, Vec3f* out_tangent, Vec3f* out_binormal);
        void evaluateGeometricConstraintsImplicit(const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, FCL_REAL qk[6], FCL_REAL phi[6]);
        void evaluateGeometricConstraintsParametric(const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, FCL_REAL qk[4], FCL_REAL phi[4]);
        void evaluateGeometricConstraintsDistanceNormalThing(const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, FCL_REAL* qk, FCL_REAL phi[6]);
        bool solveNewtonRaphson(int size, FCL_REAL tolerance, int maxIterations, const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, void(*function)(const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, FCL_REAL* qk, FCL_REAL phi[6]), FCL_REAL* guess, NewtonRaphsonStats* stats);
        bool isMinimumDistance(const SuperOvoid& s1, const Transform3f& t1, const SuperOvoid& s2, const Transform3f& t2, const double* qk, bool parametric);

        /// @brief Returns true if any of the cells in a vector is NaN.
        bool isNaN(FCL_REAL* vector, int size)
        {
            for (int i = 0; i < size; i++)
                if (std::isnan(vector[i]))
                    return true;

            return false;
        }



        /// @brief Computes the distance between two Superovoids
        /// Returns true if the superovoids are separated. Otherwise, returns false.
        /// Either way, dist will be set to the signed distance between the two objects, and p1 and p2 the minimum distance points.
        bool superOvoidSuperOvoidDistance(
            const SuperOvoid& s1, const Transform3f& tf1,
            const SuperOvoid& s2, const Transform3f& tf2,
            FCL_REAL* dist, Vec3f* p1, Vec3f* p2,
            bool collisionQuery,
            NewtonRaphsonStats* stats)
        {
            // Try to use the global statistics for output
            if (stats == NULL)
            {
                stats = &g_lastStats;

                // Try to use as input as well; if they're invalid, reset to defaults
                if (!g_lastStatsValid)
                    stats->resetToDefault();
            }

            if (stats != NULL)
            {
                stats->clearOutputs();
            }

            // If the requested query was for collision,
            // we can run some broad-phase algs to avoid useless work
            if (collisionQuery)
            {
                // Simple sphere-sphere intersection
                if (s1.aabb_radius * s1.aabb_radius + s2.aabb_radius * s2.aabb_radius
                    < (tf2.getTranslation() - tf1.getTranslation()).sqrLength())
                {
                    // Save stats as a global variable, for benchmarks
                    if (stats != NULL)
                    {
                        g_lastStatsValid = true;
                        g_lastStats = *stats;
                    }
                    return true;
                }

                // Less-simple OBB intersection
                AABB aabb1 = s1.aabb_local;
                AABB aabb2 = s2.aabb_local;
                OBB obb1; convertBV(aabb1, tf1, obb1);
                OBB obb2; convertBV(aabb2, tf2, obb2);

                if (!obb1.overlap(obb2))
                {
                    // Save stats as a global variable, for benchmarks
                    if (stats != NULL)
                    {
                        g_lastStatsValid = true;
                        g_lastStats = *stats;
                    }
                    return true;
                }
            }

            Timer totalTimer = Timer();
            totalTimer.start();

            // Attempt to simplify superovoids into superellipsoids if their tapering values are zero
            if (stats != NULL)
            {
                if (stats->superellipsoid)
                {
                    const_cast<SuperOvoid*>(&s1)->toSuperellipsoid();
                    const_cast<SuperOvoid*>(&s2)->toSuperellipsoid();
                }
                else
                {
                    const_cast<SuperOvoid*>(&s1)->toSuperovoid();
                    const_cast<SuperOvoid*>(&s2)->toSuperovoid();
                }
            }

            bool useParametric = false;
            if (stats != NULL)
                useParametric = stats->parametric;

            // With the implicit version of the problem, the unknowns are
            // the X,Y,Z coordinates of the two points expressing the minimum
            // distance between the superovoids (6 unknowns).
            // In the parametric version, the unknowns are the Azimuth,Zenith
            // parameters for the same points (4 unknowns).
            int size = (useParametric ? 4 : 6);

            // Initial guess, in local coords
            FCL_REAL qk[6];
            bool useCachedGuesses = (s1.isCachedPointValid(&s2) && s2.isCachedPointValid(&s1));

            if (!useParametric)
            {
                if (!useCachedGuesses)
                {
                    getInitialGuess(s1, tf1, s2, tf2, qk, stats, useParametric);
                }
                else // if (useCachedGuesses)
                {
                    // Use the given initial guess: convert global to local coords
                    Vec3f localP = Transform3f(tf1).inverse().transform(s1.getCachedPoint(&s2));
                    Vec3f localQ = Transform3f(tf2).inverse().transform(s2.getCachedPoint(&s1));
                    qk[0] = localP[0]; qk[1] = localP[1]; qk[2] = localP[2];
                    qk[3] = localQ[0]; qk[4] = localQ[1]; qk[5] = localQ[2];
                }
            }
            else // if (useParametric)
            {
                if (!useCachedGuesses)
                {
					getInitialGuess(s1, tf1, s2, tf2, qk, stats, useParametric);
                    //getParametricInitialGuess(s1, tf1, s2, tf2, qk, stats, useParametric);
                }
                else // if (useCachedGuesses)
                {
                    qk[0] = s1.getCachedPoint(&s2)[0];
                    qk[1] = s1.getCachedPoint(&s2)[1];
                    qk[2] = s2.getCachedPoint(&s1)[0];
                    qk[3] = s2.getCachedPoint(&s1)[1];;
                }
            }

            FCL_REAL tolerance = 1e-6;
            int max_iterations = 30;
            if (stats != NULL)
            {
                tolerance = stats->precision;
                max_iterations = stats->maxIterations;
            }

            // The roots of this function satisfy the common normal condition for minimum distance, for either implicit or parametric
            void(*function)(const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, FCL_REAL* qk, FCL_REAL phi[6])
                = useParametric ? evaluateGeometricConstraintsParametric : evaluateGeometricConstraintsImplicit;

            // Solve the problem iteratively
            bool converged = solveNewtonRaphson(size, tolerance, max_iterations,
                s1, tf1, s2, tf2, function, qk, stats);

            bool isNonMinimum = !isMinimumDistance(s1, tf1, s2, tf2, qk, useParametric);
            if (stats != NULL && isNonMinimum)
                stats->nonMinimumDistance++;

            if (useCachedGuesses && (!converged || isNonMinimum))
            {
                if (stats != NULL)
                    stats->discardedGuess++;

                // Invalidate the cached guess and try again with a new guess
                const_cast<SuperOvoid*>(&s1)->clearCachedPoint(&s2);
                const_cast<SuperOvoid*>(&s2)->clearCachedPoint(&s1);

                getInitialGuess(s1, tf1, s2, tf2, qk, stats, useParametric);

                converged = solveNewtonRaphson(size, tolerance, max_iterations,
                    s1, tf1, s2, tf2, function, qk, stats);

                isNonMinimum = isMinimumDistance(s1, tf1, s2, tf2, qk, useParametric);
                if (stats != NULL && isNonMinimum)
                    stats->nonMinimumDistance++;

            }

            // Convert solution to global coordinates
            Vec3f localPi, localPj;
            Vec3f globalPi, globalPj;
            if (useParametric)
            {
                // Convert to cartesian first, then to global
                s1.getPoint(&localPi, qk[0], qk[1], false);
                s2.getPoint(&localPj, qk[2], qk[3], false);
                globalPi = tf1.transform(localPi);
                globalPj = tf2.transform(localPj);

                // Debug: copy result UV coords to stats
                if (stats != NULL)
                {
                    stats->outParamsA[0] = qk[0];
                    stats->outParamsA[1] = qk[1];
                    stats->outParamsB[0] = qk[2];
                    stats->outParamsB[1] = qk[3];
                }
            }
            else
            {
                localPi = Vec3f(qk[0], qk[1], qk[2]);
                localPj = Vec3f(qk[3], qk[4], qk[5]);

                // Convert to global
                globalPi = tf1.transform(localPi);
                globalPj = tf2.transform(localPj);
            }

            // Save solution as initial guess for future queries
            if (!useParametric)
            {
                const_cast<SuperOvoid*>(&s1)->setCachedPoint(&s2, globalPi);
                const_cast<SuperOvoid*>(&s2)->setCachedPoint(&s1, globalPj);
            }
            else // if (useParametric)
            {
                const_cast<SuperOvoid*>(&s1)->setCachedPoint(&s2, Vec3f(qk[0], qk[1], 0));
                const_cast<SuperOvoid*>(&s2)->setCachedPoint(&s1, Vec3f(qk[2], qk[3], 0));
            }

            // Save stats as a global variable, for benchmarks
            if (stats != NULL)
            {
                g_lastStatsValid = true;
                g_lastStats = *stats;
            }

            // Output values
            bool intersecting = (s1.implicitFunction(Transform3f(tf1).inverse().transform(globalPj))) < 0;

            FCL_REAL distance = (globalPj - globalPi).norm();
            if (dist)
                *dist = distance * (intersecting ? -1 : 1);

            if (p1)
                *p1 = globalPi;
            if (p2)
                *p2 = globalPj;

            totalTimer.stop();
            if (stats != NULL)
                stats->totalTime = totalTimer.getElapsedTimeInMicroSec();

            return !intersecting;
        }

        /// @brief Returns true if the normal at S1 has the same direction as d_pq, and the normal at S2 is opposite to d_pq
        bool isMinimumDistance(const SuperOvoid& s1, const Transform3f& t1, const SuperOvoid& s2, const Transform3f& t2, const double* qk, bool parametric)
        {
            Vec3f p(qk[0], qk[1], qk[2]);
            Vec3f q(qk[3], qk[4], qk[5]);

            if (parametric)
            {
                p = Vec3f(s1.getPoint(qk[0], qk[1], false));
                q = Vec3f(s2.getPoint(qk[2], qk[3], false));
            }

            bool intersect = s1.implicitFunction(Transform3f(t1).inverse().transform(t2.transform(q))) < 0;

            Vec3f d_pq = (intersect ? -1 : 1) * (t2.transform(q) - t1.transform(p)).normalize();

            Vec3f n_p = t1.getQuatRotation().transform(s1.getNormal(p)).normalize();
            Vec3f n_q = t2.getQuatRotation().transform(s2.getNormal(q)).normalize();

            double a = d_pq.dot(n_p);
            double b = d_pq.dot(n_q);

            if (a < 0.99 || b > -0.99)
                return false;

            return true;
        }

        bool solveNewtonRaphson(
            int size, FCL_REAL tolerance, int maxIterations,
            const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2,
            void(*function)(const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, FCL_REAL* qk, FCL_REAL phi[6]),
            FCL_REAL* guess,
            NewtonRaphsonStats* stats)
        {
            bool converged = true;

            assert(size <= 6);
            // Size can be either 4 or 6. These are statically allocated for performance, so they must be size 6.
            FCL_REAL phi[6];
            FCL_REAL qkPerturb[6];
            FCL_REAL phiPerturb[6];
            FCL_REAL jacobian[6 * 6];
            FCL_REAL bestQk[6];
            int pivots[6];

            Timer timer = Timer();
            timer.start();

            // Save the best solution found at each iteration, so it can be recovered if NR does not converge
            for (int i = 0; i < size; i++)
                bestQk[i] = guess[i];

            FCL_REAL* qk = guess;

            FCL_REAL bestNormDeltaQk = 3 * tolerance;
            FCL_REAL normDeltaQk = 2 * tolerance;
            int iteration = 0;

            while (normDeltaQk > tolerance)
            {
                if (iteration >= maxIterations)
                {
                    converged = false;
                    break;
                }

#if FCL_SUPEROVOID_DEBUG_LOG > 1
                printf("\n============================== Iteration %i ==========\n", iteration);
                print_vector("qk", 6, qk);
#endif

                function(s1, tf1, s2, tf2, qk, phi);

				Timer jacobianTimer = Timer();
				jacobianTimer.start();

				bool validJacobian;
				if (size == 6 || stats->analytical == false)
					validJacobian = getNumericalJacobian(size, s1, tf1, s2, tf2, function, qk, jacobian);
				else
				{
					validJacobian = getAnalyticalParametricJacobian(s1, tf1, s2, tf2, qk, jacobian);

					// Debug: compare analytical and numerical
					FCL_REAL numericalJacobian[16];
					FCL_REAL jacobianError[16];
					getNumericalJacobian(size, s1, tf1, s2, tf2, function, qk, numericalJacobian);

					for (int i = 0; i < size * size; i++)
					{
						FCL_REAL n = numericalJacobian[i];
						FCL_REAL a = jacobian[i];
						FCL_REAL error = std::abs(n - a) / std::abs(n) * 100;

						jacobianError[i] = error;
					}

					FCL_REAL maxError = 0;
					for (int i = 0; i < size * size; i++)
					{
						if (jacobianError[i] > maxError)
							maxError = jacobianError[i];
					}
					maxError = maxError;
				}

				jacobianTimer.stop();
				if (stats != NULL)
					stats->numericalJacobianTime += jacobianTimer.getElapsedTimeInMicroSec();

				if (!validJacobian)
				{
					converged = false;
					break;
				}

                // Debug
#if FCL_SUPEROVOID_DEBUG_LOG > 1
                print_vector("phi", size, phi);
                print_matrix("Jacobian matrix", size, size, J, size);
#endif

                // Solve linear system to get actualization term
                int linSolveInfo = LAPACKE_dgesv(
                    LAPACK_COL_MAJOR,
                    size, 1, jacobian, size,
                    pivots, phi, size);

                // delta_qk = -phi;

                if (linSolveInfo > 0)
                {
                    //printf("The diagonal element of the triangular factor of A,\n");
                    //printf("U(%i,%i) is zero, so that A is singular;\n", linSolveInfo, linSolveInfo);
                    //printf("the solution could not be computed.\n");

                    converged = false;
                    break;
                }

#if FCL_SUPEROVOID_DEBUG_LOG > 2
                print_vector("-delta_qk", size, phi);
#endif

                // If result is better than previous, store its corresponding qk
                normDeltaQk = LAPACKE_dlange(LAPACK_COL_MAJOR, 'F', size, 1, phi, size);

                if ((normDeltaQk < bestNormDeltaQk) || (stats != NULL && stats->returnBest == false)
                    && !isNaN(phi, size))
                {
                    bestNormDeltaQk = normDeltaQk;

                    for (int i = 0; i < size; i++)
                        bestQk[i] = qk[i];
                }

                // Get new approximation for next iteration
                for (int i = 0; i < size; i++)
                    qk[i] = qk[i] - phi[i];

                iteration++;
            }

#if FCL_SUPEROVOID_DEBUG_LOG > 0
            printf("Newton-Raphson iterations: %i\n", iteration);
#endif

            // Keep the best guess as final result
            for (int i = 0; i < size; i++)
                qk[i] = bestQk[i];

            // Keep track of time and write it to stats struct, if it exists
            timer.stop();

            if (stats != NULL)
            {
                stats->iterationsTime = timer.getElapsedTimeInMicroSec();
                stats->iterations = iteration;

                if (converged == false)
                    stats->didNotConverge++;
            }

            return converged;
        }

        // Prepare the system of geometric constraints that express the common normal concept.
        void evaluateGeometricConstraintsImplicit(
            const SuperOvoid& s1, const Transform3f& tf1,
            const SuperOvoid& s2, const Transform3f& tf2,
            FCL_REAL* qk,
            FCL_REAL phi[6])
        {
            // In local coordinates relative to each superovoid
            Vec3f s_iP(qk[0], qk[1], qk[2]);
            Vec3f s_jQ(qk[3], qk[4], qk[5]);

            // Implicit function value for each superovoid
            FCL_REAL F_i = s1.implicitFunction(s_iP);
            FCL_REAL F_j = s2.implicitFunction(s_jQ);

            // Normal vector in local coordinates
            Vec3f n_iP = s1.getNormal(s_iP);
            Vec3f n_jQ = s2.getNormal(s_jQ);

            // Tangent vectors in local coordinates
            Vec3f v_jQ;
            Vec3f t_jQ;
            Vec3f b_jQ;
            computeTangentsWithHouseholder(n_jQ, &v_jQ, &t_jQ, &b_jQ);

            // Above vectors expressed in global coordinates
            Quaternion3f tf1Rotation = tf1.getQuatRotation();
            Quaternion3f tf2Rotation = tf2.getQuatRotation();

            Vec3f n_OP = tf1Rotation.transform(n_iP);
            // Vec3f n_OQ = tf2Rotation.transform(n_jQ);
            // Vec3f v_OQ = tf2Rotation.transform(v_jQ);
            Vec3f t_OQ = tf2Rotation.transform(t_jQ);
            Vec3f b_OQ = tf2Rotation.transform(b_jQ);

            Vec3f r_OP = tf1.transform(s_iP);
            Vec3f r_OQ = tf2.transform(s_jQ);
            Vec3f d_PQ = r_OQ - r_OP;

            // Non-linear vector column Phi = Phi(q), with q = [s_iP; s_jQ], containing the geometric constraints
            phi[0] = n_OP.dot(t_OQ);
            phi[1] = n_OP.dot(b_OQ);
            phi[2] = d_PQ.dot(t_OQ);
            phi[3] = d_PQ.dot(b_OQ);
            phi[4] = F_i;
            phi[5] = F_j;

#if FCL_SUPEROVOID_DEBUG_LOG > 3
            std::cout << "F_i " << F_i << "\n";
            std::cout << "F_j " << F_j << "\n";
            std::cout << "\n";

            std::cout << "n_iP " << n_iP << "\n";
            std::cout << "n_jQ " << n_jQ << "\n";
            std::cout << "\n";

            std::cout << "v_jQ " << v_jQ << "\n";
            std::cout << "t_jQ " << t_jQ << "\n";
            std::cout << "b_jQ " << b_jQ << "\n";
            std::cout << "\n";

            std::cout << "n_OP " << n_OP << "\n";
            std::cout << "n_OQ " << n_OQ << "\n";
            std::cout << "v_OQ " << v_OQ << "\n";
            std::cout << "t_OQ " << t_OQ << "\n";
            std::cout << "b_OQ " << b_OQ << "\n";
            std::cout << "\n";

            std::cout << "r_OP " << r_OP << "\n";
            std::cout << "r_OQ " << r_OQ << "\n";
            std::cout << "d_PQ " << d_PQ << "\n";
            std::cout << "\n";

            std::cout << "phi1 " << phi[0] << "\n";
            std::cout << "phi2 " << phi[1] << "\n";
            std::cout << "phi3 " << phi[2] << "\n";
            std::cout << "phi4 " << phi[3] << "\n";
            std::cout << "phi5 " << phi[4] << "\n";
            std::cout << "phi6 " << phi[5] << "\n";
            std::cout << "\n";
#endif
        }

        // Prepare the system of geometric constraints that express the common normal concept.
        // TODO test this
        void evaluateGeometricConstraintsDistanceNormalThing(
            const SuperOvoid& s1, const Transform3f& tf1,
            const SuperOvoid& s2, const Transform3f& tf2,
            FCL_REAL* qk,
            FCL_REAL phi[6])
        {
            // In local coordinates relative to each superovoid
            Vec3f s_iP(qk[0], qk[1], qk[2]);
            Vec3f s_jQ(qk[3], qk[4], qk[5]);

            // Implicit function value for each superovoid
            FCL_REAL F_i = s1.implicitFunction(s_iP);
            FCL_REAL F_j = s2.implicitFunction(s_jQ);

            // Normal vector in local coordinates
            Vec3f n_iP = s1.getNormal(s_iP);
            Vec3f n_jQ = s2.getNormal(s_jQ);

            // Tangent vectors in local coordinates
            Vec3f v_jQ;
            Vec3f t_jQ;
            Vec3f b_jQ;
            computeTangentsWithHouseholder(n_jQ, &v_jQ, &t_jQ, &b_jQ);

            // Above vectors expressed in global coordinates
            Quaternion3f tf1Rotation = tf1.getQuatRotation();
            Quaternion3f tf2Rotation = tf2.getQuatRotation();

            Vec3f n_OP = tf1Rotation.transform(n_iP);
            Vec3f n_OQ = tf2Rotation.transform(n_jQ);
            // Vec3f v_OQ = tf2Rotation.transform(v_jQ);
            Vec3f t_OQ = tf2Rotation.transform(t_jQ);
            Vec3f b_OQ = tf2Rotation.transform(b_jQ);

            Vec3f r_OP = tf1.transform(s_iP);
            Vec3f r_OQ = tf2.transform(s_jQ);
            Vec3f d_PQ = (r_OQ - r_OP).normalize();

            // Non-linear vector column Phi = Phi(q), with q = [s_iP; s_jQ], containing the geometric constraints
            phi[0] = n_OP.dot(t_OQ);
            phi[1] = n_OP.dot(b_OQ);
            phi[2] = d_PQ.dot(n_OP) - 1;
            phi[3] = d_PQ.dot(n_OQ) + 1;
            phi[4] = F_i;
            phi[5] = F_j;

#if FCL_SUPEROVOID_DEBUG_LOG > 3
            std::cout << "F_i " << F_i << "\n";
            std::cout << "F_j " << F_j << "\n";
            std::cout << "\n";

            std::cout << "n_iP " << n_iP << "\n";
            std::cout << "n_jQ " << n_jQ << "\n";
            std::cout << "\n";

            std::cout << "v_jQ " << v_jQ << "\n";
            std::cout << "t_jQ " << t_jQ << "\n";
            std::cout << "b_jQ " << b_jQ << "\n";
            std::cout << "\n";

            std::cout << "n_OP " << n_OP << "\n";
            std::cout << "n_OQ " << n_OQ << "\n";
            std::cout << "v_OQ " << v_OQ << "\n";
            std::cout << "t_OQ " << t_OQ << "\n";
            std::cout << "b_OQ " << b_OQ << "\n";
            std::cout << "\n";

            std::cout << "r_OP " << r_OP << "\n";
            std::cout << "r_OQ " << r_OQ << "\n";
            std::cout << "d_PQ " << d_PQ << "\n";
            std::cout << "\n";

            std::cout << "phi1 " << phi[0] << "\n";
            std::cout << "phi2 " << phi[1] << "\n";
            std::cout << "phi3 " << phi[2] << "\n";
            std::cout << "phi4 " << phi[3] << "\n";
            std::cout << "phi5 " << phi[4] << "\n";
            std::cout << "phi6 " << phi[5] << "\n";
            std::cout << "\n";
#endif
        }

        // Prepare the system of geometric constraints that express the common normal concept, for parametric case
        void evaluateGeometricConstraintsParametric(
            const SuperOvoid& s1, const Transform3f& tf1,
            const SuperOvoid& s2, const Transform3f& tf2,
            FCL_REAL* qk,
            FCL_REAL phi[4])
        {
            // Normal vector in local coordinates
            Vec3f n_iP = s1.getNormal(qk[0], qk[1], false);
            Vec3f n_jQ = s2.getNormal(qk[2], qk[3], false);

            // Tangent vectors in local coordinates
            Vec3f v_jQ;
            Vec3f t_jQ;
            Vec3f b_jQ;
            computeTangentsWithHouseholder(n_jQ, &v_jQ, &t_jQ, &b_jQ);

            // Above vectors expressed in global coordinates
            Quaternion3f tf1Rotation = tf1.getQuatRotation();
            Quaternion3f tf2Rotation = tf2.getQuatRotation();

            Vec3f n_OP = tf1Rotation.transform(n_iP);
            // Vec3f n_OQ = tf2Rotation.transform(n_jQ);
            // Vec3f v_OQ = tf2Rotation.transform(v_jQ);
            Vec3f t_OQ = tf2Rotation.transform(t_jQ);
            Vec3f b_OQ = tf2Rotation.transform(b_jQ);

            // Points in local coordinates converted to global
            Vec3f s_iP; s1.getPoint(&s_iP, qk[0], qk[1], false);
            Vec3f s_jQ; s2.getPoint(&s_jQ, qk[2], qk[3], false);
            Vec3f r_OP = tf1.transform(s_iP);
            Vec3f r_OQ = tf2.transform(s_jQ);
            Vec3f d_PQ = r_OQ - r_OP;

            // Non-linear vector column Phi = Phi(q)  (with q = [a1, z1, a2, z2]) containing the geometric constraints
            phi[0] = n_OP.dot(t_OQ);
            phi[1] = n_OP.dot(b_OQ);
            phi[2] = d_PQ.dot(t_OQ);
            phi[3] = d_PQ.dot(b_OQ);

#if FCL_SUPEROVOID_DEBUG_LOG > 3
            std::cout << "n_iP " << n_iP << "\n";
            std::cout << "n_jQ " << n_jQ << "\n";
            std::cout << "\n";

            std::cout << "v_jQ " << v_jQ << "\n";
            std::cout << "t_jQ " << t_jQ << "\n";
            std::cout << "b_jQ " << b_jQ << "\n";
            std::cout << "\n";

            std::cout << "n_OP " << n_OP << "\n";
            std::cout << "n_OQ " << n_OQ << "\n";
            std::cout << "v_OQ " << v_OQ << "\n";
            std::cout << "t_OQ " << t_OQ << "\n";
            std::cout << "b_OQ " << b_OQ << "\n";
            std::cout << "\n";

            std::cout << "r_OP " << r_OP << "\n";
            std::cout << "r_OQ " << r_OQ << "\n";
            std::cout << "d_PQ " << d_PQ << "\n";
            std::cout << "\n";

            std::cout << "phi1 " << phi[0] << "\n";
            std::cout << "phi2 " << phi[1] << "\n";
            std::cout << "phi3 " << phi[2] << "\n";
            std::cout << "phi4 " << phi[3] << "\n";
            std::cout << "\n";
#endif
        }

        // Compute the tangent and binormal vectors, given one vector, using the Householder transformation
        void computeTangentsWithHouseholder(const Vec3f normal, Vec3f* out_householder, Vec3f* out_tangent, Vec3f* out_binormal)
        {
            // Householder vector
            FCL_REAL nLength = normal.norm();
            FCL_REAL h1 = std::max(normal[0] - nLength, normal[0] + nLength);
            FCL_REAL h2 = normal[1];
            FCL_REAL h3 = normal[2];

            out_householder->setValue(h1, h2, h3);

            float hSqr = (h1 * h1) + (h2 * h2) + (h3 * h3);

            // Second and third columns of Householder matrix
            // correspond to vectors orthogonal to the first (normal)
            out_tangent->setValue(
                -2 * h1 * h2 / hSqr,
                1 - 2 * h2 * h2 / hSqr,
                -2 * h2 * h3 / hSqr);

            out_binormal->setValue(
                -2 * h1 * h3 / hSqr,
                -2 * h2 * h3 / hSqr,
                1 - 2 * h3 * h3 / hSqr);
        }






        void getInitialGuess(const SuperOvoid& s1, const Transform3f& tf1, const SuperOvoid& s2, const Transform3f& tf2, FCL_REAL qk[6], NewtonRaphsonStats* stats, bool parametricOutput)
        {
            // Find an initial guess
            NewtonRaphsonStats::INITIAL_GUESS guessType = NewtonRaphsonStats::INITIAL_GUESS::OCTREE;
            if (stats != NULL)
            {
                guessType = stats->guessType;
                stats->guessEstimations++;
            }

            switch (guessType)
            {
            case NewtonRaphsonStats::AVG_SPHERE:
                getAvgSphereInitialGuess(s1, tf1, s2, tf2, qk, stats);
                break;
            case NewtonRaphsonStats::OCTREE:
                getOctreeInitialGuess(s1, tf1, s2, tf2, qk, stats);
                break;
            case NewtonRaphsonStats::OBB:
                getOBBInitialGuess(s1, tf1, s2, tf2, qk, stats);
                break;
            case NewtonRaphsonStats::MESH:
                getMeshInitialGuess(s1, tf1, s2, tf2, qk, stats);
                break;
            case NewtonRaphsonStats::PARAMETRIC_QUADTREE:
                getParametricInitialGuess(s1, tf1, s2, tf2, qk, stats, parametricOutput);
                break;

            default:
                std::cerr << "Initial guess type unknow: " << guessType << std::endl;
                getParametricInitialGuess(s1, tf1, s2, tf2, qk, stats, parametricOutput);
                break;
            }
        }

        // Estimate two points for the minimum distance between two spheres,
        // as an initial guess for the superovoid iterative algorithm.
        // Spheres' radii is the average of the half-lengths of the ovoid
        // in each axis.
        void getAvgSphereInitialGuess(
            const SuperOvoid& s1, const Transform3f& tf1,
            const SuperOvoid& s2, const Transform3f& tf2,
            FCL_REAL qk[6],
            NewtonRaphsonStats* stats)
        {
            Timer timer = Timer();
            timer.start();

            Vec3f r_Oi = tf1.getTranslation();
            Vec3f r_Oj = tf2.getTranslation();
            Vec3f d_vector = (r_Oj - r_Oi).normalize();

            // In global coordinates
            Vec3f p1 = r_Oi + ((s1.a1 + s1.a2 + s1.a3) / 3) * d_vector;
            Vec3f p2 = r_Oj - ((s2.a1 + s2.a2 + s2.a3) / 3) * d_vector;

            // Convert to local coordinates
            Transform3f inv1 = Transform3f(tf1).inverse();
            Transform3f inv2 = Transform3f(tf2).inverse();
            Vec3f s_iP = inv1.transform(p1);
            Vec3f s_jQ = inv2.transform(p2);

            qk[0] = s_iP[0]; qk[1] = s_iP[1]; qk[2] = s_iP[2];
            qk[3] = s_jQ[0]; qk[4] = s_jQ[1]; qk[5] = s_jQ[2];

            // Keep track of time and write it to stats struct, if it exists
            timer.stop();

            if (stats != NULL)
            {
                for (int i = 0; i < 3; i++)
                {
                    stats->guessA[i] = p1[i];
                    stats->guessB[i] = p2[i];
                }

                stats->initialGuessTime = timer.getElapsedTimeInMicroSec();
            }
        }

        void getOctreeInitialGuess(
            const SuperOvoid& s1, const Transform3f& tf1,
            const SuperOvoid& s2, const Transform3f& tf2,
            FCL_REAL qk[6],
            NewtonRaphsonStats* stats)
        {
            Timer timer = Timer();
            timer.start();

            int size = 6;
            if (stats != NULL)
                size = stats->guessQuality;

            if (s1.localOctree == NULL)
                const_cast<SuperOvoid*>(&s1)->computeLocalOctree(size);
            if (s2.localOctree == NULL)
                const_cast<SuperOvoid*>(&s2)->computeLocalOctree(size);

            OcTree* tree1 = s1.localOctree;
            OcTree* tree2 = s2.localOctree;

            CollisionObject tree1_obj(boost::shared_ptr<CollisionGeometry>(tree1, null_deleter()));
            CollisionObject tree2_obj(boost::shared_ptr<CollisionGeometry>(tree2, null_deleter()));

            tree1_obj.setTransform(tf1);
            tree2_obj.setTransform(tf2);

            DistanceRequest request;
            request.enable_nearest_points = true;
            DistanceResult result;

            CollisionRequest collRequest;
            collRequest.enable_contact = true;
            CollisionResult collResult;

            // Check if octrees are colliding
            int contacts = collide(&tree1_obj, &tree2_obj, collRequest, collResult);

            // Otherwise, check distance between the two octrees
            if (contacts == 0)
                distance(&tree1_obj, &tree2_obj, request, result);

            // Convert nearest_points (global) to qk (local)
            Transform3f inv1 = Transform3f(tf1).inverse();
            Transform3f inv2 = Transform3f(tf2).inverse();

            Vec3f p1, p2;
            if (contacts > 0)
            {
                p1 = p2 = collResult.getContact(0).pos;
                /*p1 = collResult.getContact(0).pos + collResult.getContact(0).normal * collResult.getContact(0).penetration_depth;
                p2 = collResult.getContact(0).pos - collResult.getContact(0).normal * collResult.getContact(0).penetration_depth;*/
            }
            else
            {
                p1 = result.nearest_points[0];
                p2 = result.nearest_points[1];
            }

            // Write output
            Vec3f s_iP = inv1.transform(p1);
            Vec3f s_jQ = inv2.transform(p2);
            qk[0] = s_iP[0]; qk[1] = s_iP[1]; qk[2] = s_iP[2];
            qk[3] = s_jQ[0]; qk[4] = s_jQ[1]; qk[5] = s_jQ[2];

            // Keep track of time and write it to stats struct, if it exists
            timer.stop();

            if (stats != NULL)
            {
                for (int i = 0; i < 3; i++)
                {
                    stats->guessA[i] = p1[i];
                    stats->guessB[i] = p2[i];
                }

                stats->initialGuessTime = timer.getElapsedTimeInMicroSec();
            }
        }

        // Get an initial estimate of the closest points between two superovoids
        void getOBBInitialGuess(
            const SuperOvoid& s1, const Transform3f& tf1,
            const SuperOvoid& s2, const Transform3f& tf2,
            FCL_REAL qk[6],
            NewtonRaphsonStats* stats)
        {
            Timer timer = Timer();
            timer.start();

            AABB box1 = s1.aabb_local;
            AABB box2 = s2.aabb_local;
            box1.min_ += tf1.getTranslation();
            box1.max_ += tf1.getTranslation();
            box2.min_ += tf2.getTranslation();
            box2.max_ += tf2.getTranslation();

            Vec3f p1, p2;
            box1.distance(box2, &p1, &p2);

            Vec3f s_iP = Transform3f(tf1).inverse().transform(p1);
            Vec3f s_jQ = Transform3f(tf2).inverse().transform(p2);
            qk[0] = s_iP[0]; qk[1] = s_iP[1]; qk[2] = s_iP[2];
            qk[3] = s_jQ[0]; qk[4] = s_jQ[1]; qk[5] = s_jQ[2];

            // Keep track of time and write it to stats struct, if it exists
            timer.stop();

            if (stats != NULL)
            {
                for (int i = 0; i < 3; i++)
                {
                    stats->guessA[i] = p1[i];
                    stats->guessB[i] = p2[i];
                }

                stats->initialGuessTime = timer.getElapsedTimeInMicroSec();
            }
        }



        CollisionObject* makeSuperOvoidObject(const SuperOvoid& s, const Transform3f& t, int azSlices, int zeSlices, bool adaptiveMeshes)
        {
            std::vector<Vec3f>* vertices = new std::vector<Vec3f>();
            s.storePoints(vertices, azSlices, zeSlices, !adaptiveMeshes);
            std::vector<Triangle>* triangles = new std::vector<Triangle>();
            s.storeTriangles(triangles, azSlices, zeSlices);

            BVHModel<OBBRSS>* model = new BVHModel<OBBRSS>();
            model->beginModel();
            model->addSubModel(*vertices, *triangles);
            model->endModel();

            CollisionObject* obj = new CollisionObject(boost::shared_ptr<CollisionGeometry>(model), t);
            return obj;
        }

        // Get an initial estimate of the closest points between two superovoids, via closest point between two low-poly meshes
        void getMeshInitialGuess(
            const SuperOvoid& s1, const Transform3f& tf1,
            const SuperOvoid& s2, const Transform3f& tf2,
            FCL_REAL qk[6],
            NewtonRaphsonStats* stats)
        {
            const int azSlices = 8;
            const int zeSlices = 13;
            bool adaptiveMeshes = true;

            Timer timer = Timer();
            timer.start();

            GJKSolver_libccd solver;

            // Generate two low-poly meshes of the superovoids
            CollisionObject* o1 = makeSuperOvoidObject(s1, tf1, azSlices, zeSlices, adaptiveMeshes);
            CollisionObject* o2 = makeSuperOvoidObject(s2, tf2, azSlices, zeSlices, adaptiveMeshes);

            // Compute minimum distance between the two meshes
            DistanceRequest request;
            request.enable_nearest_points = true;
            DistanceResult result;
            distance(o1, o2, request, result);

            // Convert to local coordinates (for each superovoid) and write to output
            Vec3f s_iP = Transform3f(tf1).inverse().transform(result.nearest_points[0]);
            Vec3f s_jQ = Transform3f(tf2).inverse().transform(result.nearest_points[1]);
            qk[0] = s_iP[0]; qk[1] = s_iP[1]; qk[2] = s_iP[2];
            qk[3] = s_jQ[0]; qk[4] = s_jQ[1]; qk[5] = s_jQ[2];

            // Keep track of time and write it to stats struct, if it exists
            timer.stop();

            if (stats != NULL)
            {
                for (int i = 0; i < 3; i++)
                {
                    stats->guessA[i] = result.nearest_points[0][i];
                    stats->guessB[i] = result.nearest_points[1][i];
                }

                stats->initialGuessTime = timer.getElapsedTimeInMicroSec();
            }
        }





        Vec3f getQuadtreeCenter(Vec3f min, Vec3f max, int corner)
        {
            switch (corner)
            {
            case 0: return Vec3f(min[0] * 0.75 + max[0] * 0.25, min[1] * 0.75 + max[1] * 0.25, 0); break;
            case 1: return Vec3f(min[0] * 0.75 + max[0] * 0.25, min[1] * 0.25 + max[1] * 0.75, 0); break;
            case 2: return Vec3f(min[0] * 0.25 + max[0] * 0.75, min[1] * 0.75 + max[1] * 0.25, 0); break;
            case 3: return Vec3f(min[0] * 0.25 + max[0] * 0.75, min[1] * 0.25 + max[1] * 0.75, 0); break;
            default: return Vec3f();
            }
        }

        void getQuadtreeCorner(Vec3f min, Vec3f max, int corner, Vec3f* outMin, Vec3f* outMax)
        {
            FCL_REAL
                minX = min[0],
                maxX = max[0],
                minY = min[1],
                maxY = max[1];
            FCL_REAL
                medX = (minX + maxX) / 2,
                medY = (minY + maxY) / 2;

            switch (corner)
            {
            case 0:
                *outMin = Vec3f(minX, minY, 0);
                *outMax = Vec3f(medX, medY, 0);
                break;
            case 1:
                *outMin = Vec3f(minX, medY, 0);
                *outMax = Vec3f(medX, maxY, 0);
                break;
            case 2:
                *outMin = Vec3f(medX, minY, 0);
                *outMax = Vec3f(maxX, medY, 0);
                break;
            case 3:
                *outMin = Vec3f(medX, medY, 0);
                *outMax = Vec3f(maxX, maxY, 0);
                break;

            default:
                return;
            }
        }

        // Initial guess for the parametric version of Newton-Raphson; based on quadtrees
        void getParametricInitialGuess(
            const SuperOvoid& s1, const Transform3f& tf1,
            const SuperOvoid& s2, const Transform3f& tf2,
            FCL_REAL qk[6],
            NewtonRaphsonStats* stats,
            bool parametricOutput)
        {
            // The parameter space is a 2D space with
            // azimuth: [-pi, pi]       a full circle of latitude
            // zenith:  [-pi/2, pi/2]   south to north pole

            Timer timer = Timer();
            timer.start();

            Vec3f min(-pi, -pi / 2, 0);
            Vec3f max(pi, pi / 2, 0);

            Vec3f minA = min, minB = min, maxA = max, maxB = max;

            Vec3f bestA, bestB;
            int bestCornerA = 0, bestCornerB = 0;
            FCL_REAL best;
            //FCL_REAL myQk[4], phi[4];

            int iterations = 6;
            if (stats != NULL)
                iterations = stats->guessQuality;

            for (int i = 0; i < iterations; i++)
            {
                best = NAN;

                // Test every pair of zones
                for (int a = 0; a < 4; a++)
                {
                    for (int b = 0; b < 4; b++)
                    {
                        Vec3f paramsA = getQuadtreeCenter(minA, maxA, a);
                        Vec3f paramsB = getQuadtreeCenter(minB, maxB, b);

                        Vec3f pointA = s1.getPoint(paramsA[0], paramsA[1], false);
                        Vec3f pointB = s2.getPoint(paramsB[0], paramsB[1], false);

                        // Cartesian points in global coordinates
                        pointA = tf1.transform(pointA);
                        pointB = tf2.transform(pointB);

                        // Cost is the Euclidean distance
                        // Minimizing this works well for separated superovoids,
                        // but not so well for intersecting superovoids
                        // because it will yield points at the intersection
                        // (a closed curve) and not in the "minimum distance" direction
                        FCL_REAL distanceSqr = (pointB - pointA).sqrLength();
                        //FCL_REAL radialDistanceSqr = (tf1.getTranslation() - tf2.getTranslation()).sqrLength();

                        // But we can use it to quickly discard bad guesses
                        // because the minimum distance will always be shorter than
                        // the distance between the centers of the superovoids
                        //if (distanceSqr > radialDistanceSqr)
                        //	continue;

                        //// Cost is the norm of the geometric constraints result vector
                        //myQk[0] = paramsA[0];
                        //myQk[1] = paramsA[1];
                        //myQk[2] = paramsB[0];
                        //myQk[3] = paramsB[1];
                        //evaluateGeometricConstraintsParametric(s1, tf1, s2, tf2, myQk, phi);
                        //FCL_REAL constraintsSqr = (phi[0] * phi[0]) + (phi[1] * phi[1]) + (phi[2] * phi[2]) + (phi[3] * phi[3]);

                        //// Normal vectors in local, and then global coordinates
                        //Vec3f n_iP = s1.getNormal(qk[0], qk[1], false);
                        //Vec3f n_jQ = s2.getNormal(qk[2], qk[3], false);
                        //Vec3f n_OP = tf1.getQuatRotation().transform(n_iP);
                        //Vec3f n_OQ = tf2.getQuatRotation().transform(n_jQ);

                        //Vec3f distanceVec = (pointB - pointA);
                        //distanceVec.normalize();
                        //FCL_REAL minDistFactorP = -distanceVec.dot(n_OP) + 1;  // Closer to 0 is better
                        //FCL_REAL minDistFactorQ = -distanceVec.dot(-n_OQ) + 1; // Closer to 0 is better

                        FCL_REAL cost;
                        //if (i < 1)
                        //	cost = distanceSqr + constraintsSqr + minDistFactorP + minDistFactorQ;
                        //else

                        cost = distanceSqr; // + constraintsSqr + minDistFactorP + minDistFactorQ;

                        //cost = constraintsSqr + (minDistFactorP * minDistFactorP) + (minDistFactorQ * minDistFactorQ);

                        if (cost < best || isnan(best))
                        {
                            best = cost;
                            bestA = paramsA;
                            bestB = paramsB;
                            bestCornerA = a;
                            bestCornerB = b;
                        }
                    }
                }

                // Update the "tree" with the best zones found
                getQuadtreeCorner(minA, maxA, bestCornerA, &minA, &maxA);
                getQuadtreeCorner(minB, maxB, bestCornerB, &minB, &maxB);
            }

            if (parametricOutput)
            {
                qk[0] = bestA[0];
                qk[1] = bestA[1];
                qk[2] = bestB[0];
                qk[3] = bestB[1];

                qk[4] = 1337;
                qk[5] = 1337;
            }
            else
            {
                // Calculate points in local coordinates
                Vec3f pointA = s1.getPoint(bestA[0], bestA[1], false);
                Vec3f pointB = s2.getPoint(bestB[0], bestB[1], false);

                // qk accepts local coordinates
                for (int i = 0; i < 3; i++)
                {
                    qk[i] = pointA[i];
                    qk[i + 3] = pointB[i];
                }
            }

            // Keep track of time and write it to stats struct, if it exists
            timer.stop();

            if (stats != NULL)
            {
                for (int i = 0; i < 3; i++)
                {
                    // Parameters in radians
                    //stats->guessA[i] = bestA[i];
                    //stats->guessB[i] = bestB[i];

                    // Guess converted to Cartesian coordinates
                    stats->guessA[i] = tf1.transform(s1.getPoint(bestA[0], bestA[1], false))[i];
                    stats->guessB[i] = tf2.transform(s2.getPoint(bestB[0], bestB[1], false))[i];

                }

                stats->initialGuessTime = timer.getElapsedTimeInMicroSec();
            }
        }







        // UNUSED

        //// Superovoid Tree initial guess estimation
        //void getTreeInitialGuess(
        //	const SuperOvoid& s1, const Transform3f& tf1,
        //	const SuperOvoid& s2, const Transform3f& tf2,
        //	FCL_REAL qk[6],
        //	NewtonRaphsonStats* stats)
        //{
        //	Timer timer = Timer();
        //	timer.start();

        //          AABB box1 = s1.aabb_local;
        //          AABB box2 = s2.aabb_local;

        //	const int iterations = 5;
        //	for (int i = 0; i < iterations; i++)
        //	{

        //	}


        //	// Keep track of time and write it to stats struct, if it exists
        //	timer.stop();

        //	if (stats != NULL)
        //	{
        //		for (int i = 0; i < 3; i++)
        //		{
        //			//stats->guessA[i] = p1[i];
        //			//stats->guessB[i] = p2[i];
        //		}

        //		stats->initialGuessTime = timer.getElapsedTimeInMicroSec();
        //	}
        //}

        // enum CORNER { BLD, BLU, BRD, BRU, FLD, FLU, FRD, FRU };

        //AABB getOctreeChild(AABB parent, CORNER corner)
        //{
        //	// An octree child has three common planes with its parent (one common corner)
        //	Vec3f min = parent.min_;
        //	Vec3f max = parent.max_;

        //	// In each dimension, the child will either inherit the min or the max coordinate of the parent
        //	bool xMax = ((int)corner & 0x01) != 0;
        //	bool yMax = ((int)corner & 0x02) != 0;
        //	bool zMax = ((int)corner & 0x04) != 0;
        //	bool maxes[] = { xMax, yMax, zMax };

        //	// For each dimension X,Y,Z
        //	for (int i = 0; i < 3; i++)
        //	{
        //		// If it inherits the 'max' coordinate, the min must increase
        //		if (maxes[i])
        //			min[i] = (min[i] + max[i]) / 2;
        //		// If it inherits the 'min' coordinate, the max decreases instead
        //		else
        //			max[i] = (min[i] + max[i]) / 2;
        //	}

        //	return AABB(min, max);
        //}

        //// Gets the explicit coordinates of a corner of an AABB
        //Vec3f getAABBCorner(AABB box, CORNER corner)
        //{
        //	Vec3f cornerPoint;

        //	bool xMax = ((int)corner & 0x01) != 0;
        //	bool yMax = ((int)corner & 0x02) != 0;
        //	bool zMax = ((int)corner & 0x04) != 0;
        //	bool maxes[] = { xMax, yMax, zMax };

        //	for (int i = 0; i < 3; i++)
        //		cornerPoint[i] = (maxes[i] ? box.min_[i] : box.max_[i]);

        //	return cornerPoint;
        //}

        //// Gets the closest pair of points, between two AABBs' corners
        //void getClosestPoints(AABB boxA, AABB boxB, int* outIndex1, int* outIndex2)
        //{
        //	int p1index = 0, p2index = 0;
        //	FCL_REAL sqrMinDistance = (boxA.min_ - boxB.min_).sqrLength();

        //	for (int a = 0; a < 8; a++)
        //	{
        //		for (int b = 0; b < 8; b++)
        //		{
        //			Vec3f cornerA = getAABBCorner(boxA, (CORNER)a);
        //			Vec3f cornerB = getAABBCorner(boxB, (CORNER)b);

        //			FCL_REAL dist = (cornerA - cornerB).sqrLength();

        //			if (dist < sqrMinDistance)
        //			{
        //				sqrMinDistance = dist;
        //				p1index = a;
        //				p2index = b;
        //			}
        //		}
        //	}

        //	*outIndex1 = p1index;
        //	*outIndex2 = p2index;
        //}
    }




    // From SuperOvoid.h
    void SuperOvoid::setCachedPoint(const SuperOvoid* other, Vec3f point)
    {
        // Clear the previous cached point for this pair of SuperOvoids, and add the new one
        cachedGuesses.erase(other);
        cachedGuesses.emplace(other, point);
    }

    void SuperOvoid::clearCachedPoint(const SuperOvoid* other)
    {
        cachedGuesses.erase(other);
    }

    void SuperOvoid::clearCachedPoints()
    {
        cachedGuesses.clear();
    }

    Vec3f SuperOvoid::getCachedPoint(const SuperOvoid* other) const
    {
        boost::unordered_map<const SuperOvoid*, Vec3f>::const_iterator t = cachedGuesses.find(other);

        if (t != cachedGuesses.end())
            return t->second;
        else
            return Vec3f();
    }

    bool SuperOvoid::isCachedPointValid(const SuperOvoid* other) const
    {
        return (cachedGuesses.find(other) != cachedGuesses.end());
    }

    std::ostream& operator<< (std::ostream& out, SuperOvoid& s)
    {
        out << "(" <<
            s.a1 << ", " <<
            s.a2 << ", " <<
            s.a3 << ", " <<
            s.epsilon1 << ", " <<
            s.epsilon2 << ", " <<
            s.taperingX << ", " <<
            s.taperingY << ", " <<
            ")";

        return out;
    }
}
