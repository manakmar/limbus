#pragma once

#ifndef LIMBUS_UNION_BALLS_H
#define LIMBUS_UNION_BALLS_H

#include <limbus/Config.h>

#include <limits>
#include <vector>
#include <functional>

#include <boost/math/constants/constants.hpp>

#pragma warning( push )
#pragma warning( disable : 4146 ) // <gmp.h> unary minus operator applied to unsigned type, result still unsigned
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/number_utils.h>
#pragma warning( pop )

#include <limbus/Logger.h>
#include <limbus/Ball.h>
#include <limbus/IntersectingBalls.h>
#include <limbus/VoidSpace.h>
#include <limbus/AlphaShapeTraits.h>

LIMBUS_BEGIN_NAMESPACE

/** Represents the union of a set of balls and offers methods 
 * to compute the volume of the union analytically via weighted alpha shapes. 
 */
class UnionBalls : private AlphaShapeTraits {
private:

	#pragma region Private Types and Classes
	/** A shortcut to the class for computing the intersection of several balls. */
	typedef IntersectingBalls Intersection;
		
	#pragma endregion

	/** The pointer to the alpha shape instance. */
	Alpha_shape_ptr mAlphaShape;

public:

	/** Creates the union of balls from given parameters. 
	 *
	 * @param balls the pointer to the first of <c>n</c> subsequent balls
	 * @param n the number of balls
	 */
	UnionBalls(const Ball* balls, size_t n):
	  mAlphaShape(createAlphaShape(balls, balls + n))
	{
	}

	/** Creates the union of balls from given parameters. 
	 *
	 * The balls in the range [begin, end) will be copied into an 
	 * internal buffer. The <c>value_type</c> of <c>InputIterator</c> 
	 * must be <c>Ball</c>.
	 *
	 * @param begin the iterator to the begin of the range of balls
	 * @param end the iterator to the end of the range of balls
	 */
	template < class InputIterator >
	UnionBalls(const InputIterator& begin, const InputIterator& end):
		mAlphaShape(createAlphaShape(begin, end))
	{
	}

	/** Destroys the object and releases the associated memory. */
	virtual ~UnionBalls(void) {}

	/** Computes the volume of the union of balls.
	 *
	 * @return the volume
	 */
	inline double computeVolume()
	{
          const Alpha_shape_ptr as = mAlphaShape;
          assert(&as != 0 && "Null alpha shape.");

          double v = 0.0;

          // volume for vertices
          VolumeFunction verticesVolumeAccumulator(as);
          for (Alpha_shape::Finite_vertices_iterator it(as->finite_vertices_begin()), ite(as->finite_vertices_end()); it != ite; ++it) {
	          verticesVolumeAccumulator(it);
          }
          v += verticesVolumeAccumulator.getVolume();
          LOG("sum vertices: +" << verticesVolumeAccumulator.getVolume());

          // volume of edges
          const double edgesVolume = std::for_each(as->finite_edges_begin(), as->finite_edges_end(), VolumeFunction(as)).getVolume();
          v -= edgesVolume;
          LOG("sum edges: -" << edgesVolume);

          // volume of facets
          const double facetsVolume = std::for_each(as->finite_facets_begin(), as->finite_facets_end(), VolumeFunction(as)).getVolume();
          v += facetsVolume;
          LOG("sum facets: +" << facetsVolume);

          // volume of cells
          VolumeFunction cellsVolumeAccumulator(as);
          for (Alpha_shape::Finite_cells_iterator it(as->finite_cells_begin()), ite(as->finite_cells_end()); it != ite; ++it) {
	          verticesVolumeAccumulator(it);
          }
          v -= cellsVolumeAccumulator.getVolume();
          LOG("sum cells: -" << cellsVolumeAccumulator.getVolume());

          return v;
        }

	/** Finds the configurations of balls forming internal voids.
	 *
	 * @param out output iterator for the class VoidSpace
	 */
	template <class OutputIterator>
	void findVoids(OutputIterator out) {
		VoidSpace::findVoidSpaces(mAlphaShape, out);
	}

private:

	/** Creates a weighted point from the given ball. 
	 *
	 * The point will have the same center as the ball but the weight will be 
	 * set to the square of the radius.
	 *
	 * @param b the ball
	 * @return the weighted point
	 */
	inline static Weighted_point createWeightedPoint(const Ball & b)
	{
	  const Vector & c(b.center());
	  return Weighted_point(Bare_point(c.x, c.y, c.z), CGAL::square(b.radius()));
        }

	/** Creates a weighted alpha shape from the given set of weighted points.
	 * 	
	 * Furthermore, each vertex of the alpha shape will get a ball.
	 *
	 * @return the weighted alpha shape
	 */
	/*
	static Alpha_shape_ptr createAlphaShape(const std::vector<Weighted_point> & points);
	*/
	template <class BallInputIterator>
	static Alpha_shape_ptr createAlphaShape(const BallInputIterator& begin, const BallInputIterator& end);

	/** Defines a functional for volume computation via basic alpha shape vertices, edges, facets and cells. 
	 *
	 * Each alpha shape element represents an intersection of one, two, three or four balls. The volume of
	 * the intersection is computed by this functional.
	 */
	struct VolumeFunction {

	private:
		/** The alpha shape. */
		const Alpha_shape_ptr mAlphaShape;

		/** The volume computed so far. */
		double mVolume;
	public:

		/** Creates a functional for volume computation. */
		VolumeFunction(const Alpha_shape_ptr alphaShape):
			mAlphaShape(alphaShape), mVolume(0.0) {
		}

		/** Computes the volume of the ball associted with the given vertex. 
		 * @param vh the handle to the vertex of an alpha shape
		 * @see computeVolume()
		 */
		void operator()(const Vertex_handle& vh) {
			const Alpha_shape::Classification_type cls = mAlphaShape->classify(vh);
			LOGD("Vertex class: " << cls);

			// skip facets that do not belong to the alpha-complex
			if (cls == Alpha_shape::EXTERIOR) {
				return;
			}

			const Ball& b_i(vh->info().getBall());
			const double v = Intersection::computeVolume(b_i);

			LOGD("union(b_i)");
			LOGD("  b_i: " << b_i.x() << ", " << b_i.y() << ", " << b_i.z() << "; " << b_i.radius());
			LOGD("  volume: " << v);


			mVolume += v;
		}

			
		/** Computes the volume of the intersection of two balls associated with 
		 * vertices of the given edge.
		 * @param e the edge of an alpha shape
		 * @see computeVolume()
		 */
		void operator()(const Edge& e) {

			const Alpha_shape::Classification_type cls = mAlphaShape->classify(e);
			LOGD("Edge class: " << cls);

			// skip facets that do not belong to the alpha-complex
			if (cls == Alpha_shape::EXTERIOR) {
				return;
			}

			const Cell_handle 
				cell = e.first;
			const int 
				i = e.second, 
				j = e.third;
			const Vertex_handle 
				vh_i = cell->vertex(i), 
				vh_j = cell->vertex(j);
			const Ball& 
				b_i(vh_i->info().getBall()), 
				b_j(vh_j->info().getBall());

			const double v = Intersection::computeVolume(b_i, b_j);
				
			LOGD("union(b_i, b_j)");
			LOGD("  b_i: " << b_i.x() << ", " << b_i.y() << ", " << b_i.z() << "; " << b_i.radius());
			LOGD("  b_j: " << b_j.x() << ", " << b_j.y() << ", " << b_j.z() << "; " << b_j.radius());
			LOGD("  volume: " << v);

			mVolume += v;
		}

			
		/** Computes the volume of the intersection of three balls associated with 
		 * vertices of the given facet.
		 * @param f the facet of an alpha shape
		 * @see computeVolume()
		 */
		void operator()(const Facet& f) {
			const Alpha_shape::Classification_type cls = mAlphaShape->classify(f);
			LOGD("Facet class: " << cls);

			// skip facets that do not belong to the alpha-complex
			if (cls == Alpha_shape::EXTERIOR) {
				return;
			}

			const Cell_handle cell = f.first;
			const int i = f.second;

			Ball b[3];
			for (int j = 0, k = 0; j < 4; ++j) {
				if (j != i) {
					const Vertex_handle vh_j = cell->vertex(j);
					b[k++] = vh_j->info().getBall();
				}
			}

			const double v = Intersection::computeVolume(b[0], b[1], b[2]);
				
			LOGD("union(b_0, b_1, b_2)");
			LOGD("  b_0: " << b[0].x() << ", " << b[0].y() << ", " << b[0].z() << "; " << b[0].radius());
			LOGD("  b_1: " << b[1].x() << ", " << b[1].y() << ", " << b[1].z() << "; " << b[1].radius());
			LOGD("  b_2: " << b[2].x() << ", " << b[2].y() << ", " << b[2].z() << "; " << b[2].radius());
			LOGD("  volume: " << v);

			mVolume += v;
		}

			
		/** Computes the volume of the intersection of four balls associated with 
		 * vertices of the given cell.
		 * @param c the handle to the cell of an alpha shape
		 * @see computeVolume()
		 */
		void operator()(const Cell_handle& c) {
			const Alpha_shape::Classification_type cls = mAlphaShape->classify(c);
			LOGD("Cell class: " << cls);

			// skip facets that do not belong to the alpha-complex
			if (cls == Alpha_shape::EXTERIOR) {
				return;
			}

			Ball b[4];
			for (int j = 0; j < 4; ++j) {
				const Vertex_handle vh_j = c->vertex(j);
				b[j] = vh_j->info().getBall();
			}

			const double v = Intersection::computeVolume(b[0], b[1], b[2], b[3]);
				
			LOGD("union(b_0, b_1, b_2, b_3)");
			LOGD("  b_0: " << b[0].x() << ", " << b[0].y() << ", " << b[0].z() << "; " << b[0].radius());
			LOGD("  b_1: " << b[1].x() << ", " << b[1].y() << ", " << b[1].z() << "; " << b[1].radius());
			LOGD("  b_2: " << b[2].x() << ", " << b[2].y() << ", " << b[2].z() << "; " << b[2].radius());
			LOGD("  b_3: " << b[3].x() << ", " << b[3].y() << ", " << b[3].z() << "; " << b[3].radius());
			LOGD("  volume: " << v);

			mVolume += v;
		}
			
		/** Gets the sum of all volumes computated so far. 
		 * @return the volume
		 */
		double getVolume() const { return mVolume; }
	};
	
	static const double getAlpha() { return 0.0; }

	static const Alpha_shape::Mode getMode() { return Alpha_shape::GENERAL; }
};

template <class BallInputIterator>
UnionBalls::Alpha_shape_ptr 
UnionBalls::
createAlphaShape(const BallInputIterator& begin, const BallInputIterator& end)
{
	// todo: spatial sort of the input points

	Triangulation t;
	
	// create a regular triangulation and initialize its vertices
	Vertex_handle v = NULL;
	size_t vertexId = 0;
	for (BallInputIterator it = begin; it != end; ++it, vertexId++) {
		const Ball& b(*it);
		Weighted_point p = createWeightedPoint(b);
		v = t.insert(p, v);
		if (v != Vertex_handle()) {
			v->info().setId(vertexId);
			v->info().setBall(b);
		} 
		else {
			LOGD("Could not insert a ball to the triangulation:");
			LOGD("  ball: " << b.x() << ", " << b.y() << ", " << b.z() << ", " << b.radius());
			LOGD("  id: " << vertexId);
		}
	}

	// initialize the cells of the triangulation. Each cell gets a number, even the infinite cells.
	size_t cellId(0);
	for (Triangulation::Cell_iterator it = t.cells_begin(), ite = t.cells_end(); it != ite; ++it) {
		it->info().setId(cellId++);
	}

	// set up the infinite vertex (it is a single vertex instance in the triangulation outside the vertex-container)
	double nan = std::numeric_limits<double>::quiet_NaN();
	t.infinite_vertex()->info().setId(-1);
	t.infinite_vertex()->info().setBall(Ball(nan, nan, nan, nan));

	// Create the alpha shape from the triangulation. These two triangulations will be 
	// swapped in constant time.
	Alpha_shape_ptr as(new Alpha_shape(t, getAlpha(), getMode()));
	return as;
}

LIMBUS_END_NAMESPACE

#endif
