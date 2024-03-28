#pragma once

#ifndef LIMBUS_RANDOM_VECTOR_GENERATOR_H
#define LIMBUS_RANDOM_VECTOR_GENERATOR_H

#include <limbus/Config.h>

#include <boost/random.hpp>

//#define DEFAULT_RANDOM_SEED 346785441

LIMBUS_BEGIN_NAMESPACE

/** Defines a generator of random 3D vectors with floating point coordinates ranging from 0 to 1. */
class RandomVectorGenerator	{
private:
	// Define the type of the random number generator.
	typedef boost::mt19937 base_generator_type;

	// Define the type of the distribution - produces "double"
	// values between 0 and 1 (0 inclusive, 1 exclusive).
	typedef boost::uniform_01<> distribution_type;

	base_generator_type generator;
	distribution_type distribution;

	// Define a random number generator and initialize it with a reproducible
	// seed. (The seed is unsigned, otherwise the wrong overload may be 
	// selected when using mt19937 as the base_generator_type.)
	// Note:: seed initialization in ctor.
	boost::variate_generator<base_generator_type&, distribution_type > random;
public:

	/** Creates a generator. */
	RandomVectorGenerator(void)
		: generator(), 
		random(generator, distribution)
	{
	}

	/** Destroys the generator. */
	~RandomVectorGenerator(void)
	{
	}

	/** Generates a random vector. 
	@return the random vector
	*/
	Vector operator()() {
		return Vector(random(), random(), random());
	}
};

LIMBUS_END_NAMESPACE

#endif
