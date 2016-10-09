/**
 * Implements the Newton-Raphson numerical method to calculate
 * the minimum distance between two superovoid surfaces.
 *
 * Author: Artur Goncalves
 *
 * Based on MDC-ELLIPSOIDs (Daniel Simoes Lopes)
 * http://web.ist.utl.pt/ist151462/mdc-ellipsoids.html
 */

#include "fcl/shape/superovoid.h"

namespace fcl
{
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






	// Previously on SuperOvoid.h

	/// @brief Evaluate implicit function at a given point, in local coordinates.  This returns 0 for points on the surface, >0 for points outside, <0 for points inside.
	FCL_REAL SuperOvoid::implicitFunction(Vec3f point) const
	{
		return implicitFunction(point[0], point[1], point[2]);
	}

	/// @brief Evaluate implicit function at a given point, in local coordinates.  This returns 0 for points on the surface, >0 for points outside, <0 for points inside.
	FCL_REAL SuperOvoid::implicitFunction(FCL_REAL x, FCL_REAL y, FCL_REAL z) const
	{
		x = x / a1;
		y = y / a2;
		z = z / a3;

		//if (_isOvoid)
			return
			pow(
			(pow(std::abs(x / (taperingX * z + 1)), 2 / epsilon1)
			+ pow(std::abs(y / (taperingY * z + 1)), 2 / epsilon1)),
			epsilon1 / epsilon2)
			+ pow(std::abs(z), 2 / epsilon2)
			- 1;
		/*else
			return
			pow(
			(pow(std::abs(x), 2 / epsilon1)
			+ pow(std::abs(y), 2 / epsilon1)),
			epsilon1 / epsilon2)
			+ pow(std::abs(z), 2 / epsilon2)
			- 1;*/
	}

	/// @brief Get local point in Cartesian coordinates with the given parametric coordinates (in radians).
	void SuperOvoid::getPoint(Vec3f* target, FCL_REAL azimuth, FCL_REAL zenith, bool equallySpaced) const
	{
		// Adjust parameters so that vertices are almost equally spaced, it's prettier
		if (equallySpaced)
		{
			azimuth = equallySpacedTransform(azimuth, epsilon1);
			zenith = equallySpacedTransform(zenith, epsilon2);
		}
		if (zenith > (constants::pi / 2))
			zenith += constants::pi;
		else if (zenith < (-constants::pi / 2))
			zenith -= constants::pi;

		FCL_REAL sinA = std::sin(azimuth),
			sinZ = std::sin(zenith),
			cosA = std::cos(azimuth),
			cosZ = std::cos(zenith);

		// Calculate coordinates
		FCL_REAL z = signpow(sinZ, epsilon2) * a3;

		FCL_REAL taperingFactorX = taperingX * z / a3 + 1;
		FCL_REAL taperingFactorY = taperingY * z / a3 + 1;

		FCL_REAL x = signpow(cosA, epsilon1)
			* signpow(cosZ, epsilon2)
			* taperingFactorX
			* a1;

		FCL_REAL y = signpow(sinA, epsilon1)
			* signpow(cosZ, epsilon2)
			* taperingFactorY
			* a2;

		target->setValue(x, y, z);
	}

	/// @brief Get local normal direction at a given point on the surface of the superovoid. Normalized.
	Vec3f SuperOvoid::getNormal(Vec3f point) const
	{
		FCL_REAL x = point[0] / a1;
		FCL_REAL y = point[1] / a2;
		FCL_REAL z = point[2] / a3;
		FCL_REAL absX = std::abs(x);
		FCL_REAL absY = std::abs(y);
		FCL_REAL absZ = std::abs(z);

		//if (_isOvoid)
		{
			FCL_REAL tfx = std::pow(taperingX * z + 1, -2 / epsilon1);
			FCL_REAL tfy = std::pow(taperingY * z + 1, -2 / epsilon1);
			FCL_REAL aux = std::pow(
				tfx * std::pow(absX, 2 / epsilon1)
				+ tfy * std::pow(absY, 2 / epsilon1),
				epsilon1 / epsilon2 - 1);

			FCL_REAL nX = tfx * std::pow(absX, 2 / epsilon1 - 1) * aux * SIGN(x);
			FCL_REAL nY = tfy * std::pow(absY, 2 / epsilon1 - 1) * aux * SIGN(y);

			FCL_REAL dtfx = taperingX * std::pow(taperingX * z + 1, -2 / epsilon1 - 1) * std::pow(absX, 2 / epsilon1);
			FCL_REAL dtfy = taperingY * std::pow(taperingY * z + 1, -2 / epsilon1 - 1) * std::pow(absY, 2 / epsilon1);

			FCL_REAL nZ = (-dtfx - dtfy) * aux + std::pow(absZ, 2 / epsilon2 - 1) * SIGN(z);

			Vec3f normal = Vec3f(nX / a1, nY / a2, nZ / a3);
			normal.normalize();
			return normal;
		}
		/*else
		{
			FCL_REAL aux = std::pow(
				std::pow(absX, 2 / epsilon1)
				+ std::pow(absY, 2 / epsilon1),
				epsilon1 / epsilon2 - 1);

			FCL_REAL nX = std::pow(absX, 2 / epsilon1 - 1) * aux * SIGN(x);
			FCL_REAL nY = std::pow(absY, 2 / epsilon1 - 1) * aux * SIGN(y);
			FCL_REAL nZ = std::pow(absZ, 2 / epsilon2 - 1) * SIGN(z);

			Vec3f normal = Vec3f(nX / a1, nY / a2, nZ / a3);
			normal.normalize();
			return normal;
		}*/
	}

	/// @brief Get local normal direction at a given set of coordinates on the surface of the superovoid (in radians). Normalized.
	Vec3f SuperOvoid::getNormal(FCL_REAL azimuth, FCL_REAL zenith, bool useEquallySpacedTransform) const
	{
		// Adjust parameters so that vertices are almost equally spaced, it's prettier
		if (useEquallySpacedTransform)
		{
			azimuth = equallySpacedTransform(azimuth, epsilon1);
			zenith = equallySpacedTransform(zenith, epsilon2);
		}
		if (zenith > (constants::pi / 2))
			zenith += constants::pi;
		else if (zenith < (-constants::pi / 2))
			zenith -= constants::pi;

		//if (_isOvoid)
		{
			//if (g_lastStatsValid && g_lastStats.forceImplicitNormals)
			//{
			//	// This uses the implicit formula
			//	Vec3f point;
			//	getPoint(&point, azimuth, zenith, false);
			//	return getNormal(point);
			//}
			//else
			{
				// This is the purely-parametric formula
				return getAzimuthTangent(azimuth, zenith).cross(getZenithTangent(azimuth, zenith)).normalize();
			}
		}
		/*else
		{
			double cosA = std::cos(azimuth);
			double sinA = std::sin(azimuth);
			double cosZ = std::cos(zenith);
			double sinZ = std::sin(zenith);

			double cosZPowered = SIGN(cosZ) * std::pow(std::abs(cosZ), 2 - epsilon2);
			double x = SIGN(cosA) * std::pow(std::abs(cosA), 2 - epsilon1) * cosZPowered;
			double y = SIGN(sinA) * std::pow(std::abs(sinA), 2 - epsilon1) * cosZPowered;
			double z = SIGN(sinZ) * std::pow(std::abs(sinZ), 2 - epsilon2);

			Vec3f normal = Vec3f(x / a1, y / a2, z / a3);
			normal.normalize();
			return normal;
		}*/
	}

	// @brief Get tangent in azimuth direction (tangent)
	Vec3f SuperOvoid::getAzimuthTangent(FCL_REAL azimuth, FCL_REAL zenith) const
	{
		Vec3f point = getPoint(azimuth, zenith, false);
		FCL_REAL z = point[2];

		FCL_REAL sinA = std::sin(azimuth);
		FCL_REAL cosA = std::cos(azimuth);
		FCL_REAL sinZ = std::sin(zenith);
		FCL_REAL cosZ = std::cos(zenith);

		FCL_REAL dx = a1 * (taperingX * z / a3 + 1) * signpow(cosZ, epsilon2) * (-sinA) * epsilon1 * std::pow(std::abs(cosA), epsilon1 - 1);
		FCL_REAL dy = a2 * (taperingY * z / a3 + 1) * signpow(cosZ, epsilon2) * (cosA)* epsilon1 * std::pow(std::abs(sinA), epsilon1 - 1);
		FCL_REAL dz = 0;

		return Vec3f(dx, dy, dz).normalize();
	}

	// @brief Get tangent in zenith direction (binormal)
	Vec3f SuperOvoid::getZenithTangent(FCL_REAL azimuth, FCL_REAL zenith) const
	{
		Vec3f point = getPoint(azimuth, zenith, false);
		FCL_REAL z = point[2];

		FCL_REAL sinA = std::sin(azimuth);
		FCL_REAL cosA = std::cos(azimuth);
		FCL_REAL sinZ = std::sin(zenith);
		FCL_REAL cosZ = std::cos(zenith);

		double dz = a3 * (cosZ)* epsilon2 * std::pow(std::abs(sinZ), epsilon2 - 1);

		double dEpsilon2Term = (-sinZ) * epsilon2 * std::pow(std::abs(cosZ), epsilon2 - 1);

		double dx = a1 * (taperingX * z / a3 + 1) * signpow(cosA, epsilon1) * dEpsilon2Term + dz * a1 * (taperingX / a3) * signpow(cosA, epsilon1) * signpow(cosZ, epsilon2);
		double dy = a2 * (taperingY * z / a3 + 1) * signpow(sinA, epsilon1) * dEpsilon2Term + dz * a2 * (taperingY / a3) * signpow(sinA, epsilon1) * signpow(cosZ, epsilon2);

		return Vec3f(dx, dy, dz).normalize();
	}

	// @brief dt / dphi1
	Vec3f SuperOvoid::getAzimuthTangentDerivativePhi1(FCL_REAL phi1, FCL_REAL phi2) const
	{
		// DSL notes, page 2 parametric
		FCL_REAL e1 = epsilon1, e2 = epsilon2, tx = taperingX, ty = taperingY;

		FCL_REAL dtx_dphi1 = a1 * (tx * signpow(sin(phi2), e2) + 1)
			* signpow(cos(phi2), e2) * e1
			* (-cos(phi1) * signpow(cos(phi1), e1 - 1)
			+ (-sin(phi1) * (e1 - 1) * (-sin(phi1)) * signpow(cos(phi1), e1 - 2)));

		FCL_REAL dty_dphi1 = a2 * (ty * signpow(sin(phi2), e2) + 1)
			* signpow(cos(phi2), e2) * e1
			* (-sin(phi1) * signpow(sin(phi1), e1 - 1)
			+ (cos(phi1) * (e1 - 1) * cos(phi1) * signpow(sin(phi1), e1 - 2)));

		return Vec3f(dtx_dphi1, dty_dphi1, 0).normalize();
	}

	// @brief dt / dphi2
	Vec3f SuperOvoid::getAzimuthTangentDerivativePhi2(FCL_REAL phi1, FCL_REAL phi2) const
	{
		// DSL notes, page 2 parametric
		FCL_REAL e1 = epsilon1, e2 = epsilon2, tx = taperingX, ty = taperingY;

		FCL_REAL sp1 = sin(phi1),
			sp2 = sin(phi2),
			cp1 = cos(phi1),
			cp2 = cos(phi2);

		FCL_REAL dtx_dphi2 = a1 * e1 * (-sp1) * signpow(cp1, e1 - 1)
			* (tx * e2 * cp2 * signpow(sp2, e2 - 1) * signpow(cp2, e2)
			+ (tx * signpow(sp2, e2) + 1) * e2 * -sp2 * signpow(cp2, e2 - 1));

		FCL_REAL dty_dphi2 = a2 * e1 * cp1 * signpow(sp1, e1 - 1)
			* (ty * e2 * cp2 * signpow(sp2, e2 - 1) * signpow(cp2, e2)
			+ (ty * signpow(sp2, e2) + 1) * e2 * -sp2 * signpow(cp2, e2 - 1));

		return Vec3f(dtx_dphi2, dty_dphi2, 0).normalize();
	}

	// @brief db / dphi1
	Vec3f SuperOvoid::getZenithTangentDerivativePhi1(FCL_REAL phi1, FCL_REAL phi2) const
	{
		// DSL notes, page 3 parametric
		FCL_REAL e1 = epsilon1, e2 = epsilon2, tx = taperingX, ty = taperingY;

		FCL_REAL sp1 = sin(phi1),
			sp2 = sin(phi2),
			cp1 = cos(phi1),
			cp2 = cos(phi2);

		FCL_REAL dbx_dp1 = a1 * e1 * (-sp1) * signpow(cp1, e1 - 1) *
			(tx * e2 * cp2 * signpow(sp2, e2 - 1) * signpow(cp2, e2) +
			(tx * signpow(sp2, e2) + 1) * e2 * (-sp2) * signpow(cp2, e2 - 1));

		FCL_REAL dby_dp1 = a2 * e1 * cp1 * signpow(sp1, e1 - 1) *
			(ty * e2 * cp2 * signpow(sp2, e2 - 1) * signpow(cp2, e2) + (ty * signpow(sp2, e2) + 1) * e2 * (-sp2) * signpow(cp2, e2 - 1));

		return Vec3f(dbx_dp1, dby_dp1, 0).normalize();
	}

	// @brief db / dphi2
	Vec3f SuperOvoid::getZenithTangentDerivativePhi2(FCL_REAL phi1, FCL_REAL phi2) const
	{
		// DSL notes, page 3 parametric
		FCL_REAL e1 = epsilon1, e2 = epsilon2, tx = taperingX, ty = taperingY;

		FCL_REAL sp1 = sin(phi1),
			sp2 = sin(phi2),
			cp1 = cos(phi1),
			cp2 = cos(phi2);

		// Huge common part between dbx_dp2 and dby_dp2
		FCL_REAL bigAsterisk = tx * e2 * ((e2 - 1) * cp2 * signpow(sp2, e2 - 2) * signpow(cp2, e2 + 1) + signpow(sp2, e2 - 1) * (e2 + 1) * (-sp2 * signpow(cp2, e2)))
			+ e2 * (tx * e2 * cp2 * signpow(sp2, e2 - 1) * (-sp2) * signpow(cp2, e2 - 1) + (tx * signpow(sp2, e2) + 1) * (-cp2) * signpow(cp2, e2 - 1))
			+ (tx * signpow(sp2, e2) + 1) * (-sp2) * (e2 - 1) * (-sp2) * signpow(cp2, e2 - 2);

		FCL_REAL dbx_dp2 = a1 * signpow(cp1, e1) * bigAsterisk;

		FCL_REAL dby_dp2 = a2 * signpow(sp1, e1) * bigAsterisk;

		FCL_REAL dbz_dp2 = -a3 * e2 * ((-sp2) * signpow(sp2, e2 - 1) + cp2 * (e2 - 1) * cp2 * signpow(sp2, e2 - 2));

		return Vec3f(dbx_dp2, dby_dp2, dbz_dp2).normalize();
	}

	// @brief dn / dphi1
	Vec3f SuperOvoid::getNormalDerivativePhi1(FCL_REAL phi1, FCL_REAL phi2) const
	{
		Vec3f t = getAzimuthTangent(phi1, phi2);
		Vec3f dT = getAzimuthTangentDerivativePhi1(phi1, phi2);
		Vec3f b = getZenithTangent(phi1, phi2);
		Vec3f dB = getZenithTangentDerivativePhi1(phi1, phi2);

		Vec3f dTnumerical = getAzimuthTangentDerivativePhi1(phi1, phi2);
		Vec3f dBnumerical = getZenithTangentDerivativePhi1(phi1, phi2);

		return (dT.cross(b) + t.cross(dB)).normalize();
	}

	// @brief dn / dphi2
	Vec3f SuperOvoid::getNormalDerivativePhi2(FCL_REAL phi1, FCL_REAL phi2) const
	{
		Vec3f t = getAzimuthTangent(phi1, phi2);
		Vec3f dT = getAzimuthTangentDerivativePhi2(phi1, phi2);
		Vec3f b = getZenithTangent(phi1, phi2);
		Vec3f dB = getZenithTangentDerivativePhi2(phi1, phi2);

		Vec3f dTnumerical = getAzimuthTangentDerivativePhi2(phi1, phi2);
		Vec3f dBnumerical = getZenithTangentDerivativePhi2(phi1, phi2);

		return dT.cross(b) + t.cross(dB);
	}

    // Related to minimum distance query cached results
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
        std::unordered_map<const SuperOvoid*, Vec3f>::const_iterator t = cachedGuesses.find(other);

        if (t != cachedGuesses.end())
            return t->second;
        else
            return Vec3f();
    }

    bool SuperOvoid::isCachedPointValid(const SuperOvoid* other) const
    {
        return (cachedGuesses.find(other) != cachedGuesses.end());
    }
}
