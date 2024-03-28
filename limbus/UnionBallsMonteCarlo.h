#pragma once

#ifndef UNION_BALLS_MONTE_CARLO_H
#define UNION_BALLS_MONTE_CARLO_H

#include <limbus/Config.h>
#include <limbus/Logger.h>

#include <vector>
#include <iterator>
#include <algorithm>
#include <limits>
#include <iomanip>

#include <boost/functional/hash.hpp>
#include <boost/unordered/unordered_map.hpp>

#include <limbus/Ball.h>
#include <limbus/RandomVectorGenerator.h>

#define DEFAULT_CELL_SIZE 2.0
#define DEFAULT_SAMPLES_COUNT 1000
#define SHOW_VOLUME_EVERY_NSAMPLES 100000

using namespace std;

LIMBUS_BEGIN_NAMESPACE

/** Represents the union of a set of balls and offers methods 
to compute the volume of the union by a randomized Monte Carlo algorithm. 
*/
class UnionBallsMonteCarlo {
private:

	typedef std::vector<Ball> balls_container;
	typedef balls_container::iterator balls_iterator;
	typedef balls_container::const_iterator balls_const_iterator;

	balls_container balls;
		
	typedef boost::unordered_multimap<size_t, const Ball*> ball_lists;

	ball_lists ballLists;

	RandomVectorGenerator genRandomVector01;


public:
	/** The type for counting samples. */
	typedef long long sample_size_t;

	/**
	* Creates the union of balls from given parameters. 
	* The balls are copied into an internal buffer.
	*/
	UnionBallsMonteCarlo(const Ball* balls, int n);

	/** Destroys the object and releases the associated memory. */
	~UnionBallsMonteCarlo(void);

	/** Computes the volume of the union of the ball set. */
	double computeVolume();

private:
	/**
		* Computes the bounding box for the set of balls. The bounding box 
		* is given by its length and the offset of the point with minimal 
		* coordinates.
		*/
	void computeBoundingBox(Vector & offset, Vector & length) const;

	/**
		* Converts the coordinates of a given point to the grid coordinates.
		*/
	void toGridCoordinates(const Vector & p, int & gridX, int & gridY, int & gridZ) const;

	size_t toListIndex(int gridX, int gridY, int gridZ) const {
		std::size_t seed = 35719465; // chosen arbitrarily
		boost::hash_combine(seed, gridX);
		boost::hash_combine(seed, gridY);
		boost::hash_combine(seed, gridZ);
		return seed;
	}

		
	/**
		* Determines if the given sample hits some ball in the list at the specified grid cell.
		* Returns true if the sample is inside some ball of the list, false otherwise.
		*/
	bool testSampleHit(const Vector & sample, int gridX, int gridY, int gridZ) const;

	Vector getHalfRange() const {
		double r = getCellSize();
		return Vector(r, r, r);		// TODO: configuration?
	}

	int getSamplesCountPerVolumeUnit() const {
		return DEFAULT_SAMPLES_COUNT; // TODO: configuration
	}

	/** Gets the total number of balls in the set. */
	int getBallsCount() const {
		return balls.size();
	}
		
	/** Gets the ball at the specified index. */
	Ball getBall(int ballIndex) const {
		return balls[ballIndex];
	}
		
	/**
	* Gets the size of the grid cell.
	*/
	double getCellSize() const {
		return DEFAULT_CELL_SIZE; // TODO: configurable in a better way
	}


	Vector randomVector01() {
		return genRandomVector01();
	}

	/** Inserts the input set of balls to the grid. */
	void putBallsToGrid();

	sample_size_t countRandomHitsBruteForce(sample_size_t nSamples, const Vector & bBoxOffset, const Vector & bBoxLength);
	sample_size_t countRandomHitsUsingGrid(sample_size_t nSamples, const Vector & bBoxOffset, const Vector & bBoxLength);
		
	static bool isSampleInsideBall(const Vector& sample, const Ball & b) { 
		const Vector v(sample - b.center());
		return Vector::dot(v, v) <= b.radius() * b.radius();
	}
};

inline UnionBallsMonteCarlo::UnionBallsMonteCarlo(const Ball* balls, int n) {
	this->balls.reserve(n);
	for (int i = 0; i < n; ++i) {
		this->balls.push_back(balls[i]);
	}
}

inline UnionBallsMonteCarlo::~UnionBallsMonteCarlo(void) {
}

inline void UnionBallsMonteCarlo::computeBoundingBox(Vector & offset, Vector & length) const {

	const int n = getBallsCount();
	if (n == 0) {
		offset = Vector(0, 0, 0);
		length = Vector(0, 0, 0);
		return;
	}

	// create an initial bounding box from the first ball
	const Ball b0 = getBall(0);
	const Vector c0(b0.center());
	Vector 
		boxMin(c0.x - b0.radius(), c0.y - b0.radius(), c0.z - b0.radius()), 
		boxMax(c0.x + b0.radius(), c0.y + b0.radius(), c0.z + b0.radius());

	// update the bounding box by all remaining balls
	for (int i = 1; i < n; ++i) {
		const Ball b = getBall(i);
		Vector c(b.center());

		boxMin.x = min(boxMin.x, c.x - b.radius()); boxMax.x = max(boxMax.x, c.x + b.radius());
		boxMin.y = min(boxMin.y, c.y - b.radius()); boxMax.y = max(boxMax.y, c.y + b.radius());
		boxMin.z = min(boxMin.z, c.z - b.radius()); boxMax.z = max(boxMax.z, c.z + b.radius());	
	}

	// compute the offset and length
	offset = boxMin;
	length = boxMax - boxMin;

	assert(length.x >= 0 && length.y >= 0 && length.z >= 0 && "Negative bounding box.");
}

inline void UnionBallsMonteCarlo::toGridCoordinates(const Vector & p, int & gridX, int & gridY, int & gridZ) const {
	const double cellSize = getCellSize();
	gridX = (int)floor(p.x / cellSize);
	gridY = (int)floor(p.y / cellSize);
	gridZ = (int)floor(p.z / cellSize);
}

inline double UnionBallsMonteCarlo::computeVolume(void) {
	Vector bBoxOffset, bBoxLength;
	computeBoundingBox(bBoxOffset, bBoxLength);
		
	const double bBoxVolume = bBoxLength.x * bBoxLength.y * bBoxLength.z;

	sample_size_t nSamples = std::max<sample_size_t>(1, (sample_size_t)ceil(bBoxVolume * getSamplesCountPerVolumeUnit()));
	if (nSamples < 0)
		nSamples = std::numeric_limits<sample_size_t>::max();

	//sample_size_t nHits = countRandomHitsUsingGrid(nSamples, bBoxOffset, bBoxLength);
	sample_size_t nHits = countRandomHitsBruteForce(nSamples, bBoxOffset, bBoxLength);

	// estimate the volume from the relative number of hits and the bounding box volume
	const double estimatedVolume = (bBoxVolume * nHits) / nSamples;
	return estimatedVolume;
}

	
inline UnionBallsMonteCarlo::sample_size_t 
UnionBallsMonteCarlo::countRandomHitsBruteForce(sample_size_t nSamples, const Vector & bBoxOffset, const Vector & bBoxLength) {

	const double bBoxVolume = bBoxLength.x * bBoxLength.y * bBoxLength.z;

	// generate random samples
	sample_size_t nHits = 0;
	for (sample_size_t i = 1; i <= nSamples; ++i) {
		// generate a random sample in the bounding box
		const Vector rnd01(randomVector01());
		const Vector sample(
			bBoxOffset.x + rnd01.x * bBoxLength.x, 
			bBoxOffset.y + rnd01.y * bBoxLength.y, 
			bBoxOffset.z + rnd01.z * bBoxLength.z);

		// if the sample is inside at least one ball, the sample counts as a hit.
		for (balls_const_iterator it(balls.begin()), ite(balls.end()); it != ite; ++it) {
			if (isSampleInsideBall(sample, *it)) {
				++nHits;
				break;
			}
		}

#if SHOW_VOLUME_EVERY_NSAMPLES > 0
		if (i % SHOW_VOLUME_EVERY_NSAMPLES == 0) {
			LOG(i << ", " << ((bBoxVolume * nHits) / i) << ", progress: " << 
				(int)((100.0 * i) / nSamples) << "%");
		}
#endif
	}

	return nHits;
}

inline UnionBallsMonteCarlo::sample_size_t 
UnionBallsMonteCarlo::countRandomHitsUsingGrid(sample_size_t nSamples, const Vector & bBoxOffset, const Vector & bBoxLength) {
	const Vector halfRangeBox(1, 1, 1); // TODO: specify the half range

	// generate random samples
	sample_size_t nHits = 0;
	for (sample_size_t i = 1; i <= nSamples; ++i) {

		// generate a random sample in the bounding box
		const Vector rnd01(randomVector01());
		const Vector sample(
			bBoxOffset.x + rnd01.x * bBoxLength.x, 
			bBoxOffset.y + rnd01.y * bBoxLength.y, 
			bBoxOffset.z + rnd01.z * bBoxLength.z);
			
		// compute the 3D cube representing the space which needs to be tested
		const Vector sampleFrom(sample - halfRangeBox);
		const Vector sampleTo(sample + halfRangeBox);

		// get the 3D cube in cell coordinates
		int x1, y1, z1, x2, y2, z2;
		toGridCoordinates(sampleFrom, x1, y1, z1);
		toGridCoordinates(sampleTo, x2, y2, z2);

		// enumerate all cells which needs to be tested (in the inclusive range [sampleFrom, sampleTo])
		// and test if the sample hits any ball in any list.
		bool hit = false;
		for (int ix = x1; ix <= x2; ++ix) {
			for (int iy = y1; iy <= y2; ++iy) {
				for (int iz = z1; iz <= z2; ++iz) {
					// get the list of the cell (ix, iy, iz) and test the sample against it
					hit = testSampleHit(sample, ix, iy, iz);
					if (hit) { break; }
				}
				if (hit) { break; }
			}
			if (hit) { break; }
		}

		nHits += hit ? 1 : 0;
	}

	return nHits;
}

inline bool UnionBallsMonteCarlo::testSampleHit(const Vector & sample, int gridX, int gridY, int gridZ) const {
	size_t listIndex = toListIndex(gridX, gridY, gridZ);
		
	// Get the list. It might be empty (e.g. undefined).
	// Iterate all balls of the list:
	// - test each ball if the sample is inside
	// - if it was, then return hit=true
	// - else return hit=false

	std::pair<ball_lists::const_iterator, ball_lists::const_iterator> range = ballLists.equal_range(listIndex);
	for (ball_lists::const_iterator it = range.first, ite = range.second; it != it; ++it) {
		const Ball & b = *(it->second);
		if (isSampleInsideBall(sample, b)) {
			return true;
		}
	}

	return false;
}

inline void UnionBallsMonteCarlo::putBallsToGrid() {
	// For each ball:
	// - determine its bounding box

	// - translate the box to the grid coordinates
	// - for each grid cell:
	//   - put the ball to the list of this cell
	//   - optionally, skip the subsequent lists with the same index (i.e. add the ball only once)

	for (balls_const_iterator it = balls.begin(), ite = balls.end(); it != ite; ++it) {
		const Ball * pball = &(*it);
		const Vector c(pball->center());
		const Vector r(pball->radius(), pball->radius(), pball->radius());

		// get the bounding box of the ball in grid coordinates
		int x1, y1, z1, x2, y2, z2;
		toGridCoordinates(c - r, x1, y1, z1);
		toGridCoordinates(c + r, x2, y2, z2);

		// enumerate all cells covering the ball and put the ball to their lists
		for (int ix = x1; ix <= x2; ++ix) {
			for (int iy = y1; iy <= y2; ++iy) {
				for (int iz = z1; iz <= z2; ++iz) {
					size_t listIndex = toListIndex(ix, iy, iz);
					// TODO: avoid inserting a ball to a list more than once
					ballLists.insert(ball_lists::value_type(listIndex, pball));
				}
			}
		}						
	}
}

LIMBUS_END_NAMESPACE

#endif
