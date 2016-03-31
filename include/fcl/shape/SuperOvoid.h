#ifndef FCL_SUPEROVOID_H
#define FCL_SUPEROVOID_H

#include "fcl/shape/geometric_shapes.h"
#include "fcl/octree.h"
#include <iostream>
#include <cmath>
#include <boost/unordered_map.hpp>

#define SIGN(x) ( ((x) > 0) ? 1 : ((x) < 0 ? -1 : 0) )

struct null_deleter
{
	void operator()(void const *) const
	{
	}
};

namespace fcl
{
	const FCL_REAL pi = boost::math::constants::pi<FCL_REAL>();

	/// @brief SuperOvoid tapered in the Z direction, centered at 0.
	class SuperOvoid : public ShapeBase
	{
	public:
		SuperOvoid(
			FCL_REAL a1_, FCL_REAL a2_, FCL_REAL a3_,
			FCL_REAL epsilon1_, FCL_REAL epsilon2_,
			FCL_REAL taperingX_, FCL_REAL taperingY_,
			FCL_INT32 octreeSize = 2) : ShapeBase(),
			a1(a1_), a2(a2_), a3(a3_),
			epsilon1(epsilon1_), epsilon2(epsilon2_),
			taperingX(taperingX_), taperingY(taperingY_)
		{
			_isOvoid = true;

            computeLocalAABB();
            
			localOctree = NULL;
			octomapOcTree = NULL;
			// computeLocalOctree(octreeSize);
		}

		~SuperOvoid()
		{
			delete localOctree;
			localOctree = NULL;

			delete octomapOcTree;
			octomapOcTree = NULL;
		}

	public:
		/// @brief Radius in the X axis
		FCL_REAL a1;
		/// @brief Radius in the Y axis
		FCL_REAL a2;
		/// @brief Radius in the Z axis
		FCL_REAL a3;

		/// @brief Squareness in the XY plane
		FCL_REAL epsilon1;
		/// @brief Squareness in the Z axis
		FCL_REAL epsilon2;

		/// @brief Egg factor in the X axis
		FCL_REAL taperingX;
		/// @brief Egg factor in the Y axis
		FCL_REAL taperingY;

	private:
		/// @brief If false, taperingX and taperingY are ignored and this becomes a SuperEllipsoid.
		bool _isOvoid;

	public:
		bool isOvoid()
		{
			return _isOvoid;
		}

	private:
		octomap::OcTree* octomapOcTree;
	public:
		/// @brief Bounding octree around the superovoid, in local coords
		OcTree* localOctree;

		/// @brief Compute AABB 
		void computeLocalAABB()
		{
            aabb_local = getSuperovoidAABB();

			// Actually a sphere, not an AABB?
			aabb_center = aabb_local.center();
			aabb_radius = aabb_local.radius();
		}

    private:
		AABB getSuperovoidAABB() const
		{
			FCL_REAL xSize = a1 * std::pow(1.0 + std::abs(taperingX), 1 - epsilon2 / 2);
			FCL_REAL ySize = a2 * std::pow(1.0 + std::abs(taperingY), 1 - epsilon2 / 2);
			FCL_REAL zSize = a3;

			AABB aabb;
			aabb.min_ = Vec3f(-xSize, -ySize, -zSize);
			aabb.max_ = Vec3f(xSize, ySize, zSize);

			return aabb;
		}


    private:
        boost::unordered_map<const SuperOvoid*, Vec3f> cachedGuesses;

    public:
        void setCachedPoint(const SuperOvoid* other, Vec3f point);
        void clearCachedPoint(const SuperOvoid* other);
		void clearCachedPoints();
        Vec3f getCachedPoint(const SuperOvoid* other) const;
        bool isCachedPointValid(const SuperOvoid* other) const;

		/// @brief Convert this superovoid into a simpler superellipsoid. Returns true if this can actually be converted, false if taperingX or taperingY are not zero.
		/// This should speed things up a little bit.
		bool toSuperellipsoid()
		{
			// If ovoid parameters are zero, this is a regular SuperEllipsoid
			if (std::abs(taperingX) < 1e-6 && std::abs(taperingY) < 1e-6)
				_isOvoid = false;

			return !_isOvoid;
		}

		bool toSuperovoid()
		{
			_isOvoid = true;
			return _isOvoid;
		}

		inline FCL_REAL implicitFunction(Vec3f point) const;
		FCL_REAL implicitFunction(FCL_REAL x, FCL_REAL y, FCL_REAL z) const;
		void SuperOvoid::getPoint(Vec3f* target, FCL_REAL azimuth, FCL_REAL zenith, bool equallySpaced) const;
		Vec3f getNormal(Vec3f point) const;
		Vec3f getNormal(FCL_REAL azimuth, FCL_REAL zenith, bool useEquallySpacedTransform) const;

		void computeLocalOctree(FCL_REAL size)
		{
            FCL_REAL minLength = std::min(std::min(a1, a2), a3);

			FCL_REAL resolution = 1.0 / size * minLength;
			octomapOcTree = new octomap::OcTree(resolution);

			std::vector<Vec3f> points = std::vector<Vec3f>();
			storePoints(&points, 8, 13, true);

			for (unsigned int i = 0; i < points.size(); i++)
			{
				octomapOcTree->updateNode(
					octomap::point3d(points[i][0], points[i][1], points[i][2]),
					true);
			}

			octomapOcTree->updateInnerOccupancy();

			localOctree = new OcTree(boost::shared_ptr<const octomap::OcTree>(octomapOcTree, null_deleter()));
		}



		/// @brief Get node type
		NODE_TYPE getNodeType() const { return GEOM_SUPEROVOID; }

		//Matrix3f computeMomentofInertia() const
		//{
		//	FCL_REAL I = 0.4 * radius * radius * computeVolume();
		//	return Matrix3f(I, 0, 0,
		//		0, I, 0,
		//		0, 0, I);
		//}

		//FCL_REAL computeVolume() const
		//{
		//	return 4 * boost::math::constants::pi<FCL_REAL>() * radius * radius / 3;
		//}

	private:
		static FCL_REAL equallySpacedTransform(FCL_REAL angle, FCL_REAL epsilon)
		{
			// Subtract small margin between the epsilon value and the boundaries 0 and 2
			// otherwise numerical errors will screw up everything
			FCL_REAL margin = 0.05;
			epsilon = margin + epsilon * (1 - margin);

			return std::atan(SIGN(std::tan(angle)) * std::pow(std::abs(std::tan(angle)), 1.0 / epsilon));
		}

		static FCL_REAL lerp(FCL_REAL  start, FCL_REAL end, FCL_REAL amount)
		{
			return start + (end - start) * amount;
		}


	public:


		Vec3f getPoint(FCL_REAL azimuth, FCL_REAL zenith, bool equallySpaced) const
		{
			Vec3f target;
			getPoint(&target, azimuth, zenith, equallySpaced);
			return target;
		}

		void storePoints(std::vector<Vec3f>* vertices, int azimuthSlices, int zenithSlices, bool equallySpaced) const
		{
			/* Example layout for vertices[]:
			*
			* azimuthSlices = AS = 3
			* zenithSlices  = 13 (always odd)
			*
			*      1 for north pole              0
			*      AS * 2 first row              1  2  3  4  5  6
			*      AS * 2 second row             7  8  9 10 11 12
			*      AS * 2 equator               13 14 15 16 17 18
			*      AS * 2 second-to-last row    19 20 21 22 23 24
			*      AS * 2 last                  25 26 27 28 29 30
			*      1 for south pole             31
			*
			* vertices: ((ZS+1) / 2 - 2) * AS * 2 + 2
			*           middle rows * points per row + 2 poles
			*
			* Points in the rows will be in counterclockwise order
			*/

			//Vector3[] vertices = new Vector3[((zenithSlices + 1) / 2 - 2) * azimuthSlices * 2 + 2];
			//std::vector<Vec3f>* finalVerts = new std::vector<Vec3f>(((zenithSlices + 1) / 2 - 2) * azimuthSlices * 2 + 2);
			//std::vector<Vec3f> vertices = *finalVerts;

			//std::vector<Vec3f> vertices = *out_vertices;
			(*vertices).resize(((zenithSlices + 1) / 2 - 2) * azimuthSlices * 2 + 2);

			(*vertices)[0].setValue(0, 0, 1);

			// For zenithSlices = 13, iz should go between [1...5], because 0 and 6 are the poles
			for (int iz = 1; iz < (zenithSlices + 1) / 2 - 1; iz++)
			{
				FCL_REAL zenithA = lerp(pi / 2, -pi / 2, 1.0 * iz / (zenithSlices / 2));

				for (int ia = 0; ia < azimuthSlices; ia++)
				{
					FCL_REAL azimuth = lerp(pi / 2, -pi / 2,
						1.0 * ia / azimuthSlices);

					// Half of the points are mathematically generated, the other half is a symmetrical copy
					int index = 1 + (iz - 1) * azimuthSlices * 2 + ia;
					getPoint(&((*vertices)[index]), azimuth, zenithA, equallySpaced);

					int index2 = 1 + (iz - 1) * azimuthSlices * 2 + azimuthSlices + ia;
					(*vertices)[index2].setValue(-(*vertices)[index][0], -(*vertices)[index][1], (*vertices)[index][2]);
				}
			}

			(*vertices)[((zenithSlices + 1) / 2 - 2) * azimuthSlices * 2 + 2 - 1] = Vec3f(0, 0, -1);
		}

		static void storeTriangles(std::vector<Triangle>* triangles, int azimuthSlices, int zenithSlices)
		{
			int rows = (zenithSlices + 1) / 2;
			int pointsPerRow = azimuthSlices * 2;

			(*triangles).resize((rows - 2) * pointsPerRow * 2);
			int currentIndex = 0;

			// Fill north pole, using triangles with the pole as a common vertex
			{
				int northStart = 1;
				int northEnd = pointsPerRow + 1;

				for (int f = 0; f < pointsPerRow; f++)
				{
					(*triangles)[currentIndex].set(
						0,
						wrap(f + 1, northStart, northEnd),
						wrap(f + 2, northStart, northEnd)
						);
					currentIndex++;
				}
			}

			// Fill every side face, except the poles, using two coplanar triangles (quads)
			for (int r = 0; r < rows - 3; r++)
			{
				// Each row of points has a start and end points.
				// The last face of each row connects to the first, in a circular fashion.
				int start = r * pointsPerRow + 1;
				int end = (r + 1) * pointsPerRow + 1;

				int nextStart = (r + 1) * pointsPerRow + 1;
				int nextEnd = (r + 2) * pointsPerRow + 1;

				// For each face...
				for (int f = 0; f < pointsPerRow; f++)
				{
					// Each face has 2 triangles, each with 3 vertices
					(*triangles)[currentIndex].set(
						wrap(r * (pointsPerRow)+f, start, end),
						wrap((r + 1) * (pointsPerRow)+f, nextStart, nextEnd),
						wrap((r + 1) * (pointsPerRow)+f + 1, nextStart, nextEnd)
						);
					currentIndex++;

					(*triangles)[currentIndex].set(
						wrap(r * (pointsPerRow)+f, start, end),
						wrap((r + 1) * (pointsPerRow)+f + 1, nextStart, nextEnd),
						wrap(r * (pointsPerRow)+f + 1, start, end)
						);
					currentIndex++;
				}
			}

			// Fill south pole
			{
				int southStart = (rows - 3) * pointsPerRow + 1;
				int southEnd = (rows - 2) * pointsPerRow + 1;

				for (int f = 0; f < pointsPerRow; f++)
				{
					(*triangles)[currentIndex].set(
						southEnd,
						wrap(southStart + f + 1, southStart, southEnd),
						wrap(southStart + f, southStart, southEnd)
						);
					currentIndex++;
				}
			}
		}

		void printCoeffs()
		{
			printf(/*"a1: */ "%f\n", a1);
			printf(/*"a2: */ "%f\n", a2);
			printf(/*"a3: */ "%f\n", a3);
			printf(/*"e1: */ "%f\n", epsilon1);
			printf(/*"e2: */ "%f\n", epsilon2);
			printf(/*"tx: */ "%f\n", taperingX);
			printf(/*"ty: */ "%f\n", taperingY);
		}

        friend std::ostream& operator<< (std::ostream& out, SuperOvoid& s);
        friend std::ostream& operator<< (std::ostream& out, const SuperOvoid& s);

	private:
		static int wrap(int value, int start, int end)
		{
			if (value < start)
				value += (end - start);

			return start + (value - start) % (end - start);
		}

		static inline FCL_REAL signpow(FCL_REAL a, FCL_REAL b)
		{
			return a * std::pow(std::pow(a, 2), (b - 1) / 2.0);
		}
	};
}

#endif