#pragma once

#ifndef LIMBUS_VOID_SPACE_H
#define LIMBUS_VOID_SPACE_H

#include <limbus/Config.h>

#include <iterator>
#include <list>
#include <stack>
#include <stdexcept>
#include <limits>
#include <limits.h>
#include <unordered_set>
#include <cstddef>

#include <boost/iterator/function_output_iterator.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <limbus/AlphaShapeTraits.h>
#include <limbus/UnionBalls.h>

LIMBUS_BEGIN_NAMESPACE

/** Represents a void space in a system of balls.
 *
 * A void is represented as a set of tetrahedral cells. Each cell has four 
 * vertices corresponding to the balls defining the boundary of the void.
 * <p>
 * Implementation note: Alpha shapes are used for searching voids. All voids 
 * found in a single alpha shape share a common instance of the alpha shape.
 * Visiting cells is then performed as a graph search in the alpha shape.
 */
class VoidSpace : 
	public AlphaShapeTraits 
{
private:

	typedef std::stack<Cell_handle>     Cell_handle_stack;
	typedef std::unordered_set<Cell*> Cell_ptr_set;
	typedef boost::shared_ptr<Cell_ptr_set> Cell_ptr_set_shared_ptr;
	typedef boost::shared_ptr<Cell_handle_stack> Cell_handle_stack_shared_ptr;

	static const int CELL_NEIGHBORS_COUNT = 4;

private:
	Cell_handle     mStartCell;
	Alpha_shape_ptr mAlphaShape;
	int				mNumCells;
	bool			mIsInfinite;
	size_t			mMinVertexId;
	size_t			mMaxVertexId;

public:

	/** The classification of cells, edges and vertices of a void space. */
	enum ClassificationType 
	{
		EXTERIOR = Alpha_shape::EXTERIOR,
		SINGULAR = Alpha_shape::SINGULAR,
		REGULAR = Alpha_shape::REGULAR,
		INTERIOR = Alpha_shape::INTERIOR,
	};

	/** Defines an iterator for a connected component of EXTERIOR cells (a void space) in an alpha shape. */
	struct CellIterator 
	{
	  using iterator_category = std::input_iterator_tag;
          using difference_type   = std::ptrdiff_t;
          using value_type        = Cell_handle;
          using pointer           = Cell_handle*;  // or also value_type*
          using reference         = Cell_handle&;  // or also value_type&
	private:
		#pragma region Private member variables

		/** Pointer to an alpha shape object. The alpha shape is used for the classification of cells 
		  * and their triangular faces. This tells the iteration process where to go and where to not 
		  * go during a graph search of cells. 
		  */
		Alpha_shape_ptr mAlphaShape;

		/** Explicit pointer to the current cell. There can be many iterators sharing a search state. */
		Cell_handle mCurrentCell;
		
		/** Smart pointer to a set of marked cells which have been processed or which are on the search schedule (required for a graph search). */
		Cell_ptr_set_shared_ptr  mMarkedCells;
		
		/** Smart pointer to a stack of cells representing the current search schedule (required for a graph search). */
		Cell_handle_stack_shared_ptr mScheduledCells; 

		#pragma endregion
	
	public:

		/** Creates a NULL iterator. */
		explicit CellIterator() :
			mAlphaShape(boost::make_shared<Alpha_shape>()),  
			mMarkedCells(), 
			mScheduledCells() 
		{
		  mCurrentCell = NULL;
		}
	
		/** creates a new iterator and initializes its state.
		 *
		 * @param startCell a handle of a cell from which to start the iteration process
		 * @param alphaShape a pointer to an alpha for the classification of simplices during the iteration
		 */
		CellIterator(Cell_handle startCell, Alpha_shape_ptr alphaShape) :
			mAlphaShape(alphaShape), 
			mCurrentCell(startCell), 
			mMarkedCells(boost::make_shared<Cell_ptr_set>()), 
			mScheduledCells(boost::make_shared<Cell_handle_stack>()) 
		{
			schedule(startCell);
		}

		/** Creates a weak-copy iterator, sharing the search state with another iterator except 
		 *  the current result of the other iterator, which is always stored explicitly.
		 *
		 *  @param other another iterator
		 */
		CellIterator(const CellIterator& other) : 
			mAlphaShape(other.mAlphaShape), 
			mCurrentCell(other.mCurrentCell),
			mMarkedCells(other.mMarkedCells), 
			mScheduledCells(other.mScheduledCells) {
		}

		/** The prefix ++ operator (++it) advances the iterator to the next cell handle 
		  * in the iteration process and returns an iterator to the new current cell handle. 
		  */
		CellIterator& operator++() { increment(); return *this; }

		/** The postfix ++ operator (it++) advances the iterator to the next cell handle 
		  * in the iteration process and returns an iterator to the previously current cell handle. 
		  */
		CellIterator operator++(int) { CellIterator tmp(*this); increment(); return tmp;}

		/** Determines whether this iterator equals to another iterator by comparing their explicitly stored current results. */
		bool operator==(const CellIterator& rhs) { return mCurrentCell == rhs.mCurrentCell; }
		
		/** Determines whether this iterator differs from another iterator by comparing their explicitly stored current results. */
		bool operator!=(const CellIterator& rhs) { return mCurrentCell != rhs.mCurrentCell; }

		/** Returns the current iterator result. */
		reference operator*() { return mCurrentCell; }

		/** Returns a pointer to the current iterator result. */
		pointer operator->() { return &mCurrentCell; }

	private:
		#pragma region Private Methods

		/** Advances to the next cell handle in the iteration order and sets the current cell handle variable to the next handle.
		*/
		void increment() {

			if (isEmptySchedule())
				return;

			// get a cell from the schedule 
			Cell_handle pCell(get());

			// schedule the neighbors of the cell which are EXTERIOR, NOT MARKED and have a common EXTERIOR triangular facet
			for (int neighborId = 0; neighborId < CELL_NEIGHBORS_COUNT; ++neighborId) {

				// the neighbor shall not be scheduled if the triangular facet between the cell and the neighbor is not EXTERIOR
				const Alpha_shape::Classification_type ct1(mAlphaShape->classify(pCell, neighborId));
				if (ct1 != Alpha_shape::EXTERIOR)
					continue;

				Cell_handle neighborCell = pCell->neighbor(neighborId);

				// the neighbor shall not be scheduled if it is not EXTERIOR or if it is MARKED
				const Alpha_shape::Classification_type ct2(mAlphaShape->classify(neighborCell));
				if (ct2 != Alpha_shape::EXTERIOR || isMarked(neighborCell))
					continue;

				// schedule the neighbor as it is not MARKED, it is EXTERIOR and it is connected to the cell by an EXTERIOR triangular facet
				schedule(neighborCell);
			}

			mCurrentCell = isEmptySchedule() ? 0 : peek();
		}

		/** Determines whether there are some cell handles left in the schedule. */
		bool isEmptySchedule() const {
			return mScheduledCells->empty();
		}

		/** Determines whether the given cell handle is marked or not. */
		bool isMarked(Cell_handle pCell) {
			return mMarkedCells->end() != mMarkedCells->find(&(*pCell));
		}

		/** Marks the given cell handle. */
		void setMarked(Cell_handle pCell) {
			mMarkedCells->insert(&(*pCell));
		}

		/** Puts the given cell handle to the schedule and marks it. */
		void schedule(Cell_handle pCell) {
			setMarked(pCell);
			mScheduledCells->push(pCell);
		}

		/** Gets the next cell handle in the iteration order without removing it from the schedule. */
		Cell_handle peek() const {
			return mScheduledCells->top();
		}

		/** Gets the next cell handle in the iteration order and removes it from the schedule. */
		Cell_handle get() {
			Cell_handle tmp(peek());
			mScheduledCells->pop();
			return tmp;
		}
		#pragma endregion
	};

	/** Gets an iterator for discovering cell handles defining this void space.
	 *
	 * @return the iterator
	 */
	CellIterator 
	getCellsBegin() const 
	{
		return CellIterator(mStartCell, mAlphaShape);
	}

	/** Gets an iterator representing the end of the iteration sequence. 
	 *
	 * @return the iterator
	 */
	CellIterator 
	getCellsEnd() const 
	{
		return CellIterator();
	}

	/** Gets the number of cells of this void. 
	 *
	 * @return the number of cells.
	 */
	int 
	getCellsCount() const 
	{
		return mNumCells;
	}

	/** Determines if this void contains an infinite cell (a cell adjacent to an infinite vertex).
	 *
	 * @return true if the void contains an infinite cell, false if it does not contain an infinite cell.
	 */
	bool 
	isInfinite() const 
	{
		return mIsInfinite;
	}

	size_t
	getMinVertexId(void) const
	{
		return mMinVertexId;
	}

	size_t
	getMaxVertexId(void) const
	{
		return mMaxVertexId;
	}

	/** Determines whether the given cell is closed. */
	ClassificationType 
	classify(const Cell_handle& cell) const  
	{
		return static_cast<ClassificationType>(mAlphaShape->classify(cell));
	}

	/** Determines whether the given facet is closed. 
	 * 
	 * A cell without the specified vertex index constitutes a facet.
	 * 
	 * @param cell the cell of a facet
	 * @param i the index of a vertex in the cell (0 ... 4)
	 */
	ClassificationType 
	classify(const Cell_handle& cell, int i) const  
	{
		return static_cast<ClassificationType>(mAlphaShape->classify(cell, i));
	}

	/** Determines whether the given edge is closed.
	 *
	 * Two vertex indices in the context of a cell define an edge.
	 *
	 * @param cell the cell of an edge
	 * @param i the index of the first vertex in the cell (0 ... 4)
	 * @param j the index of the second vertex in the cell (0 ... 4)
	 */
	ClassificationType 
	classify(const Cell_handle& cell, int i, int j) const 
	{
		return static_cast<ClassificationType>(mAlphaShape->classify(cell, i, j));
	}

	/** Finds the vertices of this void.
	 * 
	 * @tparam ConstVertexPtrOutputIterator Output iterator for the type <c>const Vertex *</c>.
	 * @param out The output iterator for vertices found by the search.
	 * @see Vertex
	 */
	template <class ConstVertexPtrOutputIterator>
	void 
	findVertices(ConstVertexPtrOutputIterator out) const {
		std::unordered_set<const Vertex*> marked;

		const int N_SLOTS = 4;
		for (CellIterator it(getCellsBegin()), ite(getCellsEnd()); it != ite; ++it) {
			for (int i = 0; i < N_SLOTS; ++i) {
				const Vertex * v = &(*((*it)->vertex(i)));
				bool isVertexMarked = (marked.find(v) != marked.end());
				if (isVertexMarked)
					continue;

				marked.insert(v);

				*out = v;
				++out;
			}
		}
	}

private:

	/** Creates a void instance from an EXTERIOR cell and an alpha shape. 
	 *
	 * @param startCell the pointer to a start cell from which the void can be discovered. If the cell can not 
	 * be classified as EXTERIOR, the function throws an <c>invalid_argument</c> exception.
	 * @param alphaShape the pointer to the alpha shape object used for the classification of cells
	 */
	VoidSpace(Cell_handle startCell, Alpha_shape_ptr alphaShape):
		mStartCell(startCell),
		mAlphaShape(alphaShape),
		mNumCells(0),
		mIsInfinite(false)
	{
		if (mAlphaShape->classify(mStartCell) != Alpha_shape::EXTERIOR) {
			throw std::invalid_argument("The start cell was not classified as EXTERIOR.");
		}
	}

	/** Sets the number of cells of this void. 
	 * 
	 * @param nCells The number of cells.
	 */
	void 
	setCellsCount(int nCells) 
	{
		mNumCells = nCells;
	}

	/** Sets the flag which determines whether this void is finite or infinite.
	 * 
	 * @param isInfinite true if the void is infinite, false if it is finite.
	 */
	void 
	setInfinite(bool isInfinite) 
	{
		mIsInfinite = isInfinite;
	}

	void
	setMinVertexId(size_t minVertexId)
	{
		mMinVertexId = minVertexId;
	}

	void
	setMaxVertexId(size_t maxVertexId)
	{
		mMaxVertexId = maxVertexId;
	}

	/** Finds a new void reachable from a given cell handle and skips cells 
	 * that already defined some void. Hence, new voids are set to the
	 * output iterator and the iterator is incremented but all cells of any 
	 * old void are skipped.
	 */
	template <class OutputIterator>
	struct VoidFinder 
	{
	private:
		typedef VoidSpace::Cell_ptr_set Cell_ptr_set;
		typedef VoidSpace::Cell_ptr_set_shared_ptr Cell_ptr_set_shared_ptr;
		
		Alpha_shape_ptr mAlphaShape;
		OutputIterator mOut;
		Cell_ptr_set_shared_ptr mMarkedCells;

	public:

		VoidFinder(Alpha_shape_ptr alphaShape, OutputIterator out):
			mAlphaShape(alphaShape),
			mOut(out),
			mMarkedCells(boost::make_shared<Cell_ptr_set>())
		{
		}

		VoidFinder(const VoidFinder& other):
			mOut(other.mOut),
			mAlphaShape(other.mAlphaShape),
			mMarkedCells(other.mMarkedCells)
		{
		}

		/**
		 * Finds the void reachable from a start cell, skips old void cells.
		 */
		void operator()(const Cell_handle& startCell) 
		{	
			// skip the cells which already define some previously discovered void space
			if (isMarked(startCell))
				return;

			// create a void from the given start cell
	        VoidSpace v(startCell, mAlphaShape);

			bool isInfinite(false);
			int nCells(0);
			
			#pragma warning( "Better initialisation, this probably won't work." )
			size_t 
				minVertexId(UINT_MAX), 
				maxVertexId = 0;

			// find all the cells defining the void and mark them - these will be skipped in the future
			// NOTE: Idea: Maybe this could be done by some union-find approach without 
			//       visiting all cells of a void (by merging neighbors of a cell to a set).
			for (VoidSpace::CellIterator it = v.getCellsBegin(), ite = v.getCellsEnd(); it != ite; ++it) {
				Cell_handle cell = *it;
				setMarked(cell);

				// count the cells and determine if the void is connected to an infinite cell.
				++nCells;
				isInfinite |= mAlphaShape->is_infinite(cell);

				// find the minimum and maximum vertex index
				for (int j = 0; j < 4; ++j) {
					size_t id = cell->vertex(j)->info().getId();
					if (id < minVertexId)
						minVertexId = id;
					if (id > maxVertexId)
						maxVertexId = id;
				}
			}

			v.setCellsCount(nCells);
			v.setInfinite(isInfinite);
			v.setMinVertexId(minVertexId);
			v.setMaxVertexId(maxVertexId);

			*mOut = v;
			++mOut;
	    }

	private:
		bool isMarked(Cell_handle cell) const {
			return mMarkedCells->end() != mMarkedCells->find(&(*cell));
		}

		void setMarked(Cell_handle cell) const {
			mMarkedCells->insert(&(*cell));
		}
	};

public:

	/** Finds the voids in the given alpha shape and outputs them using an output iterator.
	 *
	 * @tparam OutputIterator output iterator class applicable to the class VoidSpace.
	 * @param[in] alphaShape the alpha shape for searching the voids.
	 * @param out the iterator used for the output of VoidSpace instances
	 */
	template <class OutputIterator>
	static void 
	findVoidSpaces(Alpha_shape_ptr alphaShape, OutputIterator out) 
	{	
		VoidFinder<OutputIterator> voidFinder(alphaShape, out);
		alphaShape->get_alpha_shape_cells(boost::make_function_output_iterator(voidFinder), Alpha_shape::EXTERIOR);
	}
};

LIMBUS_END_NAMESPACE

#endif
