#pragma once

#ifndef LIMBUS_ALPHA_SHAPE_TRAITS_H
#define LIMBUS_ALPHA_SHAPE_TRAITS_H

#include <limbus/Config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

#include <boost/shared_ptr.hpp>

#include <limbus/Ball.h>

LIMBUS_BEGIN_NAMESPACE

/** Contains types useful for working with CGAL alpha shapes. */
class AlphaShapeTraits
{
public:
	/** A shortcut to a CGAL kernel. */
	typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

	class VertexInfo 
	{
	private:
		size_t mId;
		Ball mBall;

	public:

		void setId(size_t id) {
			mId = id;
		}

		size_t getId() const { 
			return mId; 
		}

		void setBall(const Ball& ball) {
			mBall = ball;
		}

		const Ball& getBall() const {
			return mBall;
		}
	};

	class CellInfo 
	{
	private:
		size_t mId;

	public:

		size_t getId() { 
			return mId; 
		}

		void setId(size_t id) { 
			mId = id; 
		}
	};

	typedef CGAL::Regular_triangulation_vertex_base_3<Kernel> Vb0;
	typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, Kernel, Vb0> Vb1;
	/** Vertex base. */
	typedef CGAL::Alpha_shape_vertex_base_3<Kernel, Vb1> Vb;
	
	
	typedef CGAL::Regular_triangulation_cell_base_3<Kernel> Cb0;
	typedef CGAL::Triangulation_cell_base_with_info_3<CellInfo, Kernel, Cb0> Cb1; 
	/** Cell base. */
	typedef CGAL::Alpha_shape_cell_base_3<Kernel, Cb1>       Cb;

	/** Triangulation data structure. */
	typedef CGAL::Triangulation_data_structure_3<Vb,Cb> Tds;

	/** Triangulation (regular). */
	typedef CGAL::Regular_triangulation_3<Kernel,Tds> Triangulation;
		
	/** CGAL Alpha shape. */
	typedef CGAL::Alpha_shape_3<Triangulation>        Alpha_shape;

	typedef boost::shared_ptr<Alpha_shape>            Alpha_shape_ptr;
		
	/** The type of a cell handle of an alpha shape. */
	typedef Alpha_shape::Cell_handle          Cell_handle;
		
	/** The type of a vertex handle of an alpha shape. */
	typedef Alpha_shape::Vertex_handle        Vertex_handle;
		
	/** The type of a vertex of an alpha shape. */
	typedef Alpha_shape::Vertex               Vertex;
		
	/** The type of an edge of an alpha shape. */
	typedef Alpha_shape::Edge                 Edge;
		
	/** The type of a facet of an alpha shape. */
	typedef Alpha_shape::Facet                Facet;
		
	/** The type of a cell of an alpha shape. */
	typedef Alpha_shape::Cell                 Cell;

	/** The type of a weighted point from geometric traits. */
	typedef Triangulation::Weighted_point                Weighted_point;
		
	/** The type of a bare point from geometric traits. */
	typedef Triangulation::Bare_point                    Bare_point;

	/** Iterator for vertices of an alpha shape. */
	typedef Alpha_shape::Alpha_shape_vertices_iterator Vertices_iterator;
		
	/** Iterator for finite edges of an alpha shape. */
	typedef Alpha_shape::Finite_edges_iterator         Finite_edges_iterator;
		
	/** Iterator for finite facets of an alpha shape. */
	typedef Alpha_shape::Alpha_shape_facets_iterator   Facets_iterator;
		
	/** Pointer to the alpha shape. */
	typedef std::unique_ptr<Alpha_shape> Alpha_shape_auto_ptr;
};

LIMBUS_END_NAMESPACE

#endif
