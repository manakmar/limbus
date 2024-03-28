#pragma once

#ifndef LIMBUS_INTERSECTING_BALLS_H
#define LIMBUS_INTERSECTING_BALLS_H

#include <limbus/Config.h>
#include <limbus/Logger.h>

#include <limbus/Vector.h>

#include <limbus/Ball.h>
#include <cmath>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Point_3.h>

// simplified and faster computation
#define SIMPLIFIED

LIMBUS_BEGIN_NAMESPACE

/** Defines methods for computing volume of intersection among one, two, three and four balls. 
 *
 * The volume is computed analytically. Empty intersections are <strong>NOT</strong>
 * handled in this class. It is supposed that this kind of information is handled on the 
 * topological level of weighted alpha shapes.
 *
 * Analytic equations come from the paper "Measuring Space Filling Diagrams and 
 * Voids" written by Herbert Edelsbrunner and Ping Fu in 1994 (available at 
 * https://pub.ista.ac.at/~edels/Papers/1994-07-MeasuringSpaceFillingDiagrams.pdf, 
 * tested in March 2024).
 */
class IntersectingBalls
{
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Point_3<K> Point;

private:

	/** Gets the value of PI */
	static double pi() {
		return boost::math::constants::pi<double>();
	}

public:
	
	/** Computes the area of a ball <c>i</c>.
	@param i the ball 
	@return the area of the ball 
	*/
	static double computeArea(const Ball& i) {
		assert(i.radius() >= 0 && "Negative radius of ball i.");
		return 4 * pi() * (i.radius() * i.radius());
	}
	
	/** Computes the volume of a ball <c>i</c>
	@param i the ball 
	@return the volume of the ball 
	*/
	static double computeVolume(const Ball& i) {
		assert(i.radius() >= 0 && "Negative radius of ball i.");
#ifdef SIMPLIFIED
		const double r = i.radius();
		return (4.0 / 3.0) * pi() * (r * r * r);
#else
		return i.radius() * computeArea(i) / 3.0;
#endif
	}

	
	/** Computes the volume of the intersection between balls <c>i</c> and <c>j</c>.
	@param i the first ball 
	@param j the second ball 
	@return the volume of the intersection
	*/
	static double computeVolume(const Ball& i, const Ball& j) {
		assert(i.radius() >= 0 && "Negative radius of ball i.");
		assert(j.radius() >= 0 && "Negative radius of ball j.");

		// Two intersecting balls form two caps.
		const double 
			cap_volume_ij = computeCapVolume(i, j),
			cap_volume_ji = computeCapVolume(j, i);

		// The overal volume is the sum of these two caps.
		return cap_volume_ij + cap_volume_ji;
	}

	/** Computes the volume of the intersection among balls <c>i</c>, <c>j</c> and <c>k</c>.
	@param i the first ball 
	@param j the second ball 
	@param k the third ball 
	@return the volume of the intersection
	*/
	static double computeVolume(const Ball& i, const Ball& j, const Ball& k) {
		assert(i.radius() >= 0 && "Negative radius of ball i.");
		assert(j.radius() >= 0 && "Negative radius of ball j.");
		assert(k.radius() >= 0 && "Negative radius of ball k.");

		const double 
			cap_volume_ijk = computeCapVolume(i, j, k),
			cap_volume_jik = computeCapVolume(j, i, k),
			cap_volume_kij = computeCapVolume(k, i, j);

		return cap_volume_ijk + cap_volume_jik + cap_volume_kij;
	}

	/** Computes the volume of the intersection among balls <c>i</c>, <c>j</c>, <c>k</c> and <c>l</c>.
	@param i the first ball 
	@param j the second ball 
	@param k the third ball 
	@param l the fourth ball 
	@return the volume of the intersection
	*/
	static double computeVolume(const Ball& i, const Ball& j, const Ball& k, const Ball& l) {
		assert(i.radius() >= 0 && "Negative radius of ball i.");
		assert(j.radius() >= 0 && "Negative radius of ball j.");
		assert(k.radius() >= 0 && "Negative radius of ball k.");
		assert(l.radius() >= 0 && "Negative radius of ball l.");

		const double
			cap_volume_ijkl = computeCapVolume(i, j, k, l),
			cap_volume_jikl = computeCapVolume(j, i, k, l),
			cap_volume_kijl = computeCapVolume(k, i, j, l),
			cap_volume_lijk = computeCapVolume(l, i, j, k);

		return cap_volume_ijkl + cap_volume_jikl + cap_volume_kijl + cap_volume_lijk;
	}

private:

	/// Computes the volume of the cap of the ball i which is inside ball j.
	/// The slicing plane defining the cap is the plane containing the circle
	/// of intersection between i and j.
	static double computeCapVolume(const Ball& i, const Ball& j) {

#ifdef SIMPLIFIED
		/* A) direct computation using a simple formula */
		const double h = computeCapHeight(i, j);
		const double v = ( pi() * h * h * (3 * i.radius() - h) ) / 3.0;
		assert(v >= 0 && "Negative cap volume");
		return v;
#else
		/* B) more-complex computation, using several sub-formulas */
		const double 
			cap_area_ij = computeCapArea(i, j),
			cap_height_ij = computeCapHeight(i, j),
			disk_area_ij = computeDiskArea((i, j);	

		// S is the volume of a cone having its apex in the center of the sphere,
		// plus the volume of the cap. Computed by integrating infinitesimaly 
		// small pyramids, each one having its apex in the center of the sphere 
		// and a base on the surface of the sphere. The domain of integration is
		// the area of the cap.
		const double S = i.radius() * cap_area_ij / 3.0;
		assert(S >= 0 && "Negative volume (S)");

		// C is the volume of the cone with apex in the center of the sphere.
		const double C = (i.radius() - cap_height_ij) * disk_area_ij / 3.0;
		assert(C >= 0 && "Negative volume (C)");

		// subtracting these two volumes yields the volume of the cap.
		assert(S >= C && "Negative cap volume.");
		return S - C;
#endif
	}

	/// Computes the volume of the intersection of two caps in the context of 
	/// ball i. The first cap is the cap (i, j) from ball i, the second cap is 
	/// the cap (i, k) from ball i.
	static double computeCapVolume(const Ball& i, const Ball& j, const Ball& k) {
		// Intersecting the surface of ball i with balls j and k yields 
		// the area of cap(i, j, k). S2 is the integral volume of cones 
		// sharing the common apex - the center of ball i - and having the 
		// base on the surface of ball i. The integral domain is the surface 
		// of cap(i, j, k). From this volume, two signed volumes must be 
		// subtracted in order to get the cap volume. These volumes - Cj and 
		// Ck - are the volumes of cones having their apex in the center of 
		// ball i and the bases are parts of the discs from the intersection 
		// of balls i and j, or i and k, respectivelly. 
		const double S2 = (i.radius() * computeCapArea(i, j, k)) / 3;
		const double Cj = (i.radius() - computeCapHeight(i, j)) * computeSegmentArea(i, j, k) / 3;
		const double Ck = (i.radius() - computeCapHeight(i, k)) * computeSegmentArea(i, k, j) / 3;
		const double v = S2 - Cj - Ck;

		if (v < 0) {
			LOGE("Error: " << "Negative volume of the intersection of two caps.");
			LOGE(" - ball i " << i.x() << ", " << i.y() << ", " << i.z() << ", " << i.radius());
			LOGE(" - ball j " << j.x() << ", " << j.y() << ", " << j.z() << ", " << j.radius());
			LOGE(" - ball k " << k.x() << ", " << k.y() << ", " << k.z() << ", " << k.radius());
		}
		assert(v >= 0 && "Negative volume of the intersection of two caps.");
		return v;
	}

	/// Computes the volume of the intersection of three caps in the context of 
	/// ball i. The first cap is the cap (i, j) from ball i, the second cap is 
	/// the cap (i, k) from ball i, the third is the cap (i, l) from ball i.
	static double computeCapVolume(const Ball& i, const Ball& j, const Ball& k, const Ball& l) {
		const double S3 = (i.radius() * computeCapArea(i, j, k, l)) / 3;
		const double Cj = (i.radius() - computeCapHeight(i, j)) * computeSegmentArea(i, j, k, l) / 3;
		const double Ck = (i.radius() - computeCapHeight(i, k)) * computeSegmentArea(i, k, j, l) / 3;
		const double Cl = (i.radius() - computeCapHeight(i, l)) * computeSegmentArea(i, l, j, k) / 3;
		const double v = S3 - Cj - Ck - Cl;

		assert(v >= 0 && "Negative volume of the intersection of three caps.");
		return v;
	}
	
	/// Computes the area of the cap of the ball i inside ball j.
	/// The slicing plane defining the cap is the plane containing the circle
	/// of intersection between i and j.
	static double computeCapArea(const Ball& i, const Ball& j) {

		const double cap_height_ij = computeCapHeight(i, j);

		// Rotating a circle segment (parametrized by spherical coordinates) 
		// creates a surface. Integrating this surface yields the cap area.
		const double cap_area_ij = 2 * pi() * i.radius() * cap_height_ij;
		assert(cap_area_ij >= 0 && "Negative cap area.");

		return cap_area_ij;
	}

	
	
	static double computeCapArea(const Ball& i, const Ball& j, const Ball& k) {

		// p_jk and p_kj are two orthogonal centers - the angle between any 
		// two centers of two different balls at an orthogonal center is
		// the right angle.
		const Vector
			p_jk(computeTriangleDual(i, j, k)), 
			p_kj(computeTriangleDual(i, k, j));

		const double 			
			l_j = computeSegmentAngle(i, j, k), 
			l_k = computeSegmentAngle(i, k, j);

		const Vector 
			s(i.center()), 
			t(j.center()), 
			u(k.center());

		const double 
			phi_jk = 0.5 - computeDihedralAngle(s, p_jk, t, u),
			phi_kj = 0.5 - computeDihedralAngle(s, p_kj, u, t);

		const double
			A1 = 0.5 * computeArea(i) * (phi_jk + phi_kj),
			A2 = 2 * pi() * i.radius() * l_j * (i.radius() - computeCapHeight(i, j)),
			A3 = 2 * pi() * i.radius() * l_k * (i.radius() - computeCapHeight(i, k));

		const double area = A1 - A2 - A3;
		assert(area >= 0 && "Negative area of the intersection of two caps");
		
		return area;
	}

	static double computeCapArea(const Ball& i, const Ball& j, const Ball& k, const Ball& l) {
		if (!isCcw(center(i), center(j), center(k), center(l))) {
			return computeCapAreaCcw(i, j, l, k);
		} else {
			return computeCapAreaCcw(i, j, k, l);
		}
	}

	static double computeCapAreaCcw(const Ball& i, const Ball& j, const Ball& k, const Ball& l) {
		const Vector 
			p_kj(computeTriangleDual(i, k, j)),
			p_lk(computeTriangleDual(i, l, k)),
			p_jl(computeTriangleDual(i, j, l));

		const double
			l_j = computeSegmentAngle(i, j, k, l),
			l_k = computeSegmentAngle(i, k, l, j),
			l_l = computeSegmentAngle(i, l, j, k);

		const Vector 
			s(i.center()), 
			t(j.center()), 
			u(k.center()), 
			v(l.center());

		const double 
			phi_kj = 0.5 - computeDihedralAngle(s, p_kj, u, t),
			phi_lk = 0.5 - computeDihedralAngle(s, p_lk, v, u),
			phi_jl = 0.5 - computeDihedralAngle(s, p_jl, t, v);

		const double
			A1 = 0.5 * computeArea(i) * (phi_kj + phi_lk + phi_jl - 0.5),
			A2 = 2 * pi() * i.radius() * l_j * (i.radius() - computeCapHeight(i, j)),
			A3 = 2 * pi() * i.radius() * l_k * (i.radius() - computeCapHeight(i, k)),
			A4 = 2 * pi() * i.radius() * l_l * (i.radius() - computeCapHeight(i, l));

		const double area = A1 - A2 - A3 - A4;
		assert(area >= 0 && "Negative area of the intersection of three caps");

		return area;
	}

	
	/// Computes the height of a cap of the ball i. 
	/// The cap is the one being inside ball j, the slicing plane defining 
	/// the cap is the plane of the circle of intersection between i and j.
	static double computeCapHeight(const Ball& i, const Ball& j) {
		Vector	s(center(i)),		// the center of ball i
			Y(computeOrthogonalCenter(i, j));	// the orthogonal center between i and j

		// Y lies in the slicing plane (as it is the center of the circle 
		// between i and j).
		
		// Determine if the center of i is inside the cap (i.e. behind the
		// slicing plane). The cap height will be computed according it.
		bool isHidden = Vector::dot(center(j) - s, Y - s) < 0;

		const double dist_sY = distance(s, Y);

		double cap_height_ij;
		if (isHidden) {
			cap_height_ij = i.radius() + dist_sY;
		} else {
			cap_height_ij = i.radius() - dist_sY;
		}

		LOGD("Ball i: " << i.x() << ", " << i.y() << ", " << i.z() << ", " << i.radius());
		LOGD("Ball j: " << j.x() << ", " << j.y() << ", " << j.z() << ", " << j.radius());
		LOGD("Cap height: " << cap_height_ij << ", 2r: " << 2 * i.radius());
		assert(cap_height_ij >= 0 && cap_height_ij <= 2 * i.radius() && "Cap height is out of bounds.");
		return cap_height_ij;
	}

	
	static double computeDiskArea(const Ball& i, const Ball& j) {
		const double computeDiskRadius_ij = computeDiskRadius(i, j);
		const double disk_area_ij = pi() * computeDiskRadius_ij * computeDiskRadius_ij;
		assert(disk_area_ij >= 0 && "Negative area of a disk.");
		return disk_area_ij;
	}

	
	static double computeDiskLength(const Ball& i, const Ball& j) {
		return 2 * pi() * computeDiskRadius(i, j);
	}

	static double computeDiskRadius(const Ball& i, const Ball& j) {
		const double cap_height_ij = computeCapHeight(i, j);
		return sqrt( cap_height_ij * (2 * i.radius() - cap_height_ij) );
	}

	
	static double computeSegmentHeight(const Ball& i, const Ball& j, const Ball& k) {
		const Vector 
			Y3(computeOrthogonalCenter(i, j, k)), 
			Y2(computeOrthogonalCenter(i, j)),
			c_k(center(k));
		const double 
			r_ij = computeDiskRadius(i, j), // the radius of the disk defined by balls i and j
			Y_23 = distance(Y2, Y3);

		bool isHidden = Vector::dot(c_k - Y2, Y3 - Y2) < 0;

		const double h = r_ij + (isHidden ? Y_23 : -Y_23);
		assert(0 <= h && h <= 2 * r_ij && "Segment height is out of bounds.");

		return h;
	}


	static double computeSegmentArea(const Ball& i, const Ball& j, const Ball& k) {
		const double 
			S = 0.5 * computeDiskRadius(i, j) * computeSegmentLength(i, j, k);
		
		const Vector 
			p_jk(computeTriangleDual(i, j, k)),
			p_kj(computeTriangleDual(i, k, j));

		const double
			H = computeDiskRadius(i, j) - computeSegmentHeight(i, j, k),
			T = 0.5 * H * distance(p_jk, p_kj);

		const double area = S - T;
		assert(0 <= area && "Negative area of a segment.");

		return area;
	}

	static double computeSegmentArea(const Ball& i, const Ball& j, const Ball& k, const Ball& l) {
		if (!isCcw(center(i), center(j), center(k), center(l))) {
			return isSegmentAreaCcw(i, j, l, k);
		} else {
			return isSegmentAreaCcw(i, j, k, l);
		}
	}

	static double isSegmentAreaCcw(const Ball& i, const Ball& j, const Ball& k, const Ball& l) {
		const Vector 
			p_jk = computeTriangleDual(i, j, k),
			p_kj = computeTriangleDual(i, k, j),
			p_jl = computeTriangleDual(i, j, l),
			p_lj = computeTriangleDual(i, l, j),
			Y(computeOrthogonalCenter(i, j, k, l));

		const double
			h_k = computeSegmentHeight(i, j, k),
			h_l = computeSegmentHeight(i, j, l),
			r_ij = computeDiskRadius(i, j),
			S = 0.5 * r_ij * computeSegmentLength(i, j, k, l),
			T_k = 0.5 * (r_ij - h_k) * distance(p_kj, Y),
			T_l = 0.5 * (r_ij - h_l) * distance(p_jl, Y);

		const double area = S - T_k - T_l;

		assert(0 <= area && "Negative area of an intersection of two segments.");
		return area;
	}

	static double computeSegmentLength(const Ball& i, const Ball& j, const Ball& k) {
		return computeSegmentAngle(i, j, k) * computeDiskLength(i, j);
	}

	static double computeSegmentLength(const Ball& i, const Ball& j, const Ball& k, const Ball& l) {
		return computeSegmentAngle(i, j, k, l) * computeDiskLength(i, j);
	}
	
	static double computeSegmentAngle(const Ball& i, const Ball& j, const Ball& k) {
		const Vector
			p_jk(computeTriangleDual(i, j, k)),
			p_kj(computeTriangleDual(i, k, j)),
			s(center(i)), 
			t(center(j)), 
			u(center(k));
		const double 
			a1 = computeDihedralAngle(s, t, u, p_jk),
			a2 = computeDihedralAngle(s, t, u, p_kj);

		const double angle = a1 + a2;
		assert(0.0 <= angle && angle <= 1.0 && "Angle out of range.");

		return angle;
	}


	
	static double computeSegmentAngle(const Ball& i, const Ball& j, const Ball& k, const Ball& l) {
		const Vector 
			p_jl(computeTriangleDual(i, j, l)),
			p_kj(computeTriangleDual(i, k, j)),
			s(center(i)), 
			t(center(j)),
			u(center(k)),
			v(center(l));
		const double 
			a1 = computeDihedralAngle(s, t, u, p_kj),
			a2 = computeDihedralAngle(s, t, v, p_jl),
			a3 = computeDihedralAngle(s, t, u, v);

		const double angle = a1 + a2 - a3;
		assert(0.0 <= angle && angle <= 1.0 && "Angle out of range.");

		return angle;
	}
	
	static Vector computeTriangleDual(const Ball& i, const Ball& j, const Ball& k) {
		const Vector 
			Y(computeOrthogonalCenter(i, j, k)),
			s(center(i)),
			t(center(j)),
			u(center(k)),
			N(Vector::cross(t - s, u - s));
		const double
			S1 = Vector::dot(Y - s, N),
			S2 = Vector::dot(N, N),
			S3 = Vector::dot(Y - s, Y - s),
			r = i.radius(),
			xi = (-S1 + sqrt(S1*S1 - S3*S2 + r*r*S2)) / S2;
		assert(xi >= 0);
		return Y + xi * N;
	}

	static Vector computeOrthogonalCenter(const Ball& i, const Ball& j, const Ball& k, const Ball& l) {
		const Vector 
			ci(center(i)), 
			cj(center(j)), 
			ck(center(k)), 
			cl(center(l));

		const double 
			ri2 = i.radius() * i.radius(), 
			rj2 = j.radius() * j.radius(), 
			rk2 = k.radius() * k.radius(), 
			rl2 = l.radius() * l.radius();

		const double
			i0 = 0.5 * (ri2 - Vector::dot(ci, ci)),
			j0 = 0.5 * (rj2 - Vector::dot(cj, cj)),
			k0 = 0.5 * (rk2 - Vector::dot(ck, ck)),
			l0 = 0.5 * (rl2 - Vector::dot(cl, cl));

		const double
			D0 = CGAL::determinant(
				ci.x, ci.y, ci.z, 1.0,
				cj.x, cj.y, cj.z, 1.0,
				ck.x, ck.y, ck.z, 1.0,
				cl.x, cl.y, cl.z, 1.0),
			Dx = CGAL::determinant(
				-i0, ci.y, ci.z, 1.0,
				-j0, cj.y, cj.z, 1.0,
				-k0, ck.y, ck.z, 1.0,
				-l0, cl.y, cl.z, 1.0),
			Dy = CGAL::determinant(
				ci.x, -i0, ci.z, 1.0,
				cj.x, -j0, cj.z, 1.0,
				ck.x, -k0, ck.z, 1.0,
				cl.x, -l0, cl.z, 1.0),
			Dz = CGAL::determinant(
				ci.x, ci.y, -i0, 1.0,
				cj.x, cj.y, -j0, 1.0,
				ck.x, ck.y, -k0, 1.0,
				cl.x, cl.y, -l0, 1.0);

		const Vector Y(Dx/D0, Dy/D0, Dz/D0);
		return Y;
	}

	static Vector computeOrthogonalCenter(const Ball& i, const Ball& j, const Ball& k) {
		const Vector 
			ci(center(i)), 
			cj(center(j)), 
			ck(center(k));

		const double 
			ri2 = i.radius() * i.radius(), 
			rj2 = j.radius() * j.radius(), 
			rk2 = k.radius() * k.radius();

		const double
			i0 = 0.5 * (ri2 - Vector::dot(ci, ci)),
			j0 = 0.5 * (rj2 - Vector::dot(cj, cj)),
			k0 = 0.5 * (rk2 - Vector::dot(ck, ck));

		const double
			A1 = CGAL::determinant(
				ci.y, ci.z, 1.0,
				cj.y, cj.z, 1.0,
				ck.y, ck.z, 1.0),
			A2 = CGAL::determinant(
				ci.z, ci.x, 1.0,
				cj.z, cj.x, 1.0,
				ck.z, ck.x, 1.0),
			A3 = CGAL::determinant(
				ci.x, ci.y, 1.0,
				cj.x, cj.y, 1.0,
				ck.x, ck.y, 1.0),
			A4 = CGAL::determinant(
				ci.x, ci.y, ci.z,
				cj.x, cj.y, cj.z,
				ck.x, ck.y, ck.z),

			D0 = CGAL::determinant(
				ci.x, ci.y, ci.z, 1.0,
				cj.x, cj.y, cj.z, 1.0,
				ck.x, ck.y, ck.z, 1.0,
				  A1,   A2,   A3, 0.0),
			Dx = CGAL::determinant(
				-i0, ci.y, ci.z, 1.0,
				-j0, cj.y, cj.z, 1.0,
				-k0, ck.y, ck.z, 1.0,
				 A4,   A2,   A3, 0.0),
			Dy = CGAL::determinant(
				ci.x, -i0, ci.z, 1.0,
				cj.x, -j0, cj.z, 1.0,
				ck.x, -k0, ck.z, 1.0,
				  A1,  A4,   A3, 0.0),
			Dz = CGAL::determinant(
				ci.x, ci.y, -i0, 1.0,
				cj.x, cj.y, -j0, 1.0,
				ck.x, ck.y, -k0, 1.0,
				  A1,   A2,  A4, 0.0);

		const Vector Y(Dx/D0, Dy/D0, Dz/D0);
		return Y;
	}

	/// Computes the orthogonal center between balls i and j.
	static Vector computeOrthogonalCenter(const Ball& i, const Ball& j) {
		const Vector 
			c_i(center(i)),  // the center of ball i
			c_j(center(j)),  // the center of ball j
			c_ij(c_j - c_i); // the Vector from c_i to c_j

		// The orthogonal center Y is on the line from c_i to c_j:
		// Y = lambda * c_i + (1 - lambda) * c_j
		// Putting this into equations
		// |Y - i|^2 = r_i^2 - r_y^2
		// |Y - j|^2 = r_j^2 - r_y^2
		// leads to the linear equation in lambda:
		// lambda = ((r_j^2 - r_i^2) / |j - i|^2 + 1) / 2

		const double 
			ir2 = i.radius() * i.radius(), 
			jr2 = j.radius() * j.radius();

		const double
			aux = Vector::dot(c_ij, c_ij),
			lambda = (1 + (jr2 - ir2)/ aux) * 0.5;

		return lambda * c_i + (1 - lambda) * c_j;
	}

	/// Computes the Euclidean distance between two vectors u and v.
	static double distance(const Vector& u, const Vector& v) {
		return Vector::distance(u, v);
	}

	/// Gets the center of ball i.
	static Vector center(const Ball& i) {
		return i.center();
	}
	
	/// Computes the dihedral angle between two triangles (s, t, u) and (s, t, v). 
	/// The angle is normalized from the interval [0, 2pi) to the interval [0, 1).
	/// From the two possible angles, the minimal one is returned, so the returned 
	/// value is always from 0.0 to 0.5.
	static double computeDihedralAngle(const Vector& s, const Vector& t, const Vector& u, const Vector& v) {
		// Compute normals
		const Vector 
			Mu( Vector::cross(u - s, u - t) ), Nu( (1.0 / sqrt(Vector::dot(Mu, Mu))) * Mu ),
			Mv( Vector::cross(v - s, v - t) ), Nv( (1.0 / sqrt(Vector::dot(Mv, Mv))) * Mv );

		// Compute the cosine of the angle between these normals.
		const double cos_uv = Vector::dot(Nu, Nv);
		assert(cos_uv >= -1.0 && cos_uv <= 1.0 && "Cosine out of range [-1.0, 1.0]");

		// Get the angle [0.0, pi] and normalize it to [0.0, 0.5]
		const double a = acos( cos_uv ) / (2 * pi());

		assert(a >= 0.0 && a <= 0.5 && "Angle out of range [0.0, 0.5]");
		return a <= 0.5 ? a : 0.5; // clamping avoids numerical instabilities
	}

	

	static bool isCcw(const Vector& i, const Vector& j, const Vector& k, const Vector& l) {
		const Point 
			pi(i.x, i.y, i.z), 
			pj(j.x, j.y, j.z), 
			pk(k.x, k.y, k.z), 
			pl(l.x, l.y, l.z);
		CGAL::Orientation ori = CGAL::orientation(pj, pk, pl, pi);
		//CGAL::Orientation ori = CGAL::orientation(pi, pj, pk, pl);
		
		assert((ori == CGAL::COUNTERCLOCKWISE || ori == CGAL::CLOCKWISE) && "Undecideable orientation.");
		return ori == CGAL::COUNTERCLOCKWISE;
		// TODO: throw an exception when the orientation is undecideable
	}
};

LIMBUS_END_NAMESPACE

#endif
