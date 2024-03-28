#pragma once

#ifndef LIMBUS_BALL_H
#define LIMBUS_BALL_H

#include <limbus/Config.h>

#include <limbus/Vector.h>

LIMBUS_BEGIN_NAMESPACE

/** Defines a spherical ball by its center and radius. */
class Ball {
private:
	/** The center of the ball. */
	double x_, y_, z_;

	/** The radius of the ball. */
	double r_;

public:

	/** Empty constructor - initializes all fields to zero. */
	Ball() : 
		x_(0.0), 
		y_(0.0), 
		z_(0.0), 
		r_(0.0)
	{
	}

	/** Creates a ball by specifying its center coordinates and the radius. 
	 * @param centerX the x-coordinate of the center of the ball.
	 * @param centerY the y-coordinate of the center of the ball. 
	 * @param centerZ the z-coordinate of the center of the ball.
	 * @param radius the radius of the ball.
	 */
	Ball(double centerX, double centerY, double centerZ, double radius): 
		x_(centerX), 
		y_(centerY), 
		z_(centerZ), 
		r_(radius)
	{
	}

	/** Creates a ball by specifying its center and the radius. 
	 * @param center the center of the ball.
	 * @param radius the radius of the ball.
	 */
	Ball(const Vector& center, double radius):
		x_(center.x), 
		y_(center.y),
		z_(center.z),
		r_(radius)
	{
	}

	/** Copy constructor. 
	 *
	 * @param ball the other ball from which to copy the center and radius.
	 */
	Ball(const Ball & ball):
		x_(ball.x_), y_(ball.y_), z_(ball.z_), r_(ball.r_)
	{
	}
		
	/** Gets the center. 
	 * 
	 * @return the center as a vector.
	 */
	Vector center() const {
		return Vector(x_, y_, z_);
	}

	/** Gets the radius. 
	 * 
	 * @return the radius value.
	 */
	 double radius() const {
		return r_;
	}

	/** Gets the x-coordinate of the center. 
	 * 
	 * @return the x-coordinate value of the center.
	 */
	const double x() const {
		return x_;
	}

	/** Gets the y-coordinate of the center. 
	 *
	 * @return the y-coordinate value of the center.
	 */
	const double y() const {
		return y_;
	}

	/** Gets the z-coordinate of the center. 
	 *
	 * @return the z-coordinate value of the center.
	 */
	const double z() const {
		return z_;
	}
};

LIMBUS_END_NAMESPACE

#endif
