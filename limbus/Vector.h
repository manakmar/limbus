#pragma once

#ifndef LIMBUS_VECTOR_H
#define LIMBUS_VECTOR_H

#include <limbus/Config.h>

#include <cmath>


LIMBUS_BEGIN_NAMESPACE

/** Defines a 3D vector. */
struct Vector
{
public:
	/** The x-coordinate. */
	double x;

	/** The y-coordinate. */
	double y;

	/** The z-coordinate. */
	double z;

	/** Creates a zero vector. */
	Vector(void) : x(0), y(0), z(0)
	{
	}

		
	/** Creates a vector by specifying its coordinates. 
	@param x the value of the x-coordinate
	@param y the value of the y-coordinate
	@param z the value of the z-coordinate
	*/
	Vector(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

		
	/** Copy constructor. */
	Vector(const Vector& other) {
		x = other.x;
		y = other.y;
		z = other.z;
	}
		
	/** Computes a cross-product of two vectors.
	@param u the first vector
	@param v the second vector
	@return a vector representing the cross product <c>u cross v</c>
	*/
	static Vector cross(const Vector& u, const Vector& v) {
		return Vector(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
	}

	/** Computes a dot-product of two vectors.
	@param u the first vector
	@param v the second vector
	@return a scalar representing the dot-product <c>u dot v</c>
	*/
	static double dot(const Vector& u, const Vector& v) {
		return u.x * v.x + u.y * v.y + u.z * v.z;
	}

		
	/** Computes the Euclidean distance between two vectors.
	@param u the first vector
	@param v the second vector
	@return the Euclidean distance between <c>u</c> and <c>v</c>
	*/
	static double distance(const Vector& u, const Vector& v) {
		const Vector w(u - v);
		return sqrt(dot(w, w));
	}

		
	/** Subtracts another vector from this vector and returns the result. 
	@param other the vector to be subtracted from this vector
	@return a vector representing the subtraction <c>this - other</c>
	*/
	Vector operator-(const Vector& other) const {
		return Vector(x - other.x, y - other.y, z - other.z);
	}

		
	/** Adds another vector to this vector and returns the result. 
	@param other the vector to be added to this vector
	@return a vector representing the addition <c>this + other</c>
	*/
	Vector operator+(const Vector& other) const {
		return Vector(x + other.x, y + other.y, z + other.z);
	}

		
	/** Multiplies this vector by a scalar value and returns the result. 
	@param scalar the scalar by which to multiply this vector
	@return a vector representing the multiplication <c>this * scalar</c>
	*/
	Vector operator*(double scalar) const {
		return Vector(x * scalar, y * scalar, z * scalar);
	}

	Vector operator-() const {
		return Vector(-x, -y, -z);
	}
};

/** Multiplies this vector by a scalar value from left and returns the result. 

This is equivalent to the multiplication by a scalar from the right.

@param scalar the scalar by which to multiply this vector
@return a vector representing the multiplication <c>scalar * this</c>
*/
inline Vector operator*(double scalar, const Vector& u) {
	return u * scalar;
}

LIMBUS_END_NAMESPACE

#endif