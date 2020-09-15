#ifndef POLYGONRASTERIZE_H_
#define POLYGONRASTERIZE_H_

#include "SDOT/RegularGrid.h"
#include "SDOT/BoundingBox.h"

#include <vector>
#include <memory>
#include <deque>

// standard includes
#include <iostream>
#include <cassert>

// includes for defining the Voronoi diagram adaptor
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/number_utils.h>


namespace sdot {

  class PolygonRasterizeIter{
  public:

      // typedefs for defining the adaptor
      typedef CGAL::Exact_predicates_exact_constructions_kernel K;

      // typedef for the result type of the point location
      typedef K::Point_2                    Point_2;
      typedef K::Vector_2                   Vector_2;
      typedef CGAL::Polygon_2<K>            Polygon_2;
      typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;
      typedef std::pair<Polygon_2::Vertex_const_iterator,Polygon_2::Vertex_const_iterator> PolygonEdge;

      enum Direction{
        CW,
        CCW
      };

      PolygonRasterizeIter(std::shared_ptr<RegularGrid> const& gridIn,
                           std::shared_ptr<Polygon_2>   const& polyIn);

      PolygonRasterizeIter& Increment();

      /** Returns whether the indices are valid (i.e., in the polygon) or not.
          If called immediately after construction, this can be used to detect if
          any grid cells are inside the polygon.  Similarly, it can be used to
          detect when all indices inside the polygon have been detected.

          @code{.cpp}
          InteriorPolygonGrid iter(grid,poly);
          while(iter.IsValid()){
            iter.Increment();
          }
          @endcode
      */
      bool IsValid() const;

      std::pair<int,int> const& Indices() const;

      /** Returns whether or not the current grid cell intersects the polygon boundary. */
      bool IsBoundary() const{return overlapPoly!=nullptr;};

      /** Returns a polygon containing the current grid clipped to the polygon. */
      std::shared_ptr<Polygon_2> const& OverlapPoly() const{return overlapPoly;};

      std::shared_ptr<BoundingBox> const& Cell() const{return cellBox;};

      double cellBoundTime = 0;
      double edgeTime = 0;
      double xIndTime = 0;
      double interTime = 0;
      unsigned int numInters = 0;
      unsigned int numIncrs = 0;
      
    private:

      const double compTol = 1e-14;
      const double horizTol = 1e-10;

      bool isValid;

      // Boolean denoting if current grid cell is entirely within the polygon
      bool isBoundary;

      std::pair<int,int> indices;

      // The following doubles hold the bounds of the current grid cell
      std::shared_ptr<BoundingBox> cellBox;

      std::shared_ptr<RegularGrid> grid;
      std::shared_ptr<Polygon_2> poly;

      std::shared_ptr<Polygon_2> overlapPoly;

      /** Edges moving clockwise around the Laguerre cell that might possibly
         intersect the current grid cell.
      */
      std::deque<PolygonEdge> edgesCW;
      std::deque<double> slopesCW; // inverse slopes (dx/dy) for clockwise edges

      /** Edges moving counter-clockwise around the Laguerre cell that might possibly
         intersect the current grid cell.
      */
      std::deque<PolygonEdge> edgesCCW;
      std::deque<double> slopesCCW; // inverse slopes (dx/dy) for counter-clockwise edges


      // Stores the maximum x index that could potentially overlap with the left boundary (CW)
      unsigned int maxLeftIndBC;

      // Stores the minimum x index that could potentially overlap with the right boundary (CCW)
      unsigned int minRightIndBC;

      // Stores the last xInd for the current y
      unsigned int xInd_end;

      /** Returns the next edge (pair of vertex iterators) in either the CW or CCW
          direction starting from the currIter vertex.  If horizCheck is true, this
          will return the next edge that is NOT horizontal.
      */
      std::pair<Polygon_2::Vertex_const_iterator,Polygon_2::Vertex_const_iterator> NextEdge(Polygon_2::Vertex_const_iterator const& currIter,
                                                    Direction                               dir,
                                                    bool                                    horizCheck=false);

       Polygon_2::Vertex_const_iterator NextVertex(Polygon_2::Vertex_const_iterator const& currIter,
                                                   Direction                               dir);

       /** Returns bounds on where boundary cells might occur.  The output of this
          function is a tuple [leftMin, leftMax, rightMin, rightMax], containing
          ranges of the xindices that start a row [leftMin,leftMax], and the x
          indices that end a row [rightMin, rightMax]
       */
       std::tuple<unsigned int, unsigned int, unsigned int, unsigned int> GetXInds(double y) const;


       std::shared_ptr<Polygon_2> BoundaryIntersection();

       /** Updates the edge queues (if necessary) when the y index is updated.
       */
       void UpdateEdges(unsigned int yInd);

       /** Computes the points where a particular Laguerre cell edge crosses the
           boundary of the current grid cell.  The first result is where the Laguerre
           edge enters the cell.  The second result is where the Laguerre cells leaves
           the cell.   If only one crossing exists (i.e., the  edge starts  or ends
           in the grid cell), then a nullptr is used to indicate the missing crossing.
           The `includeEdges` input gives control over what points are considered "inside"
           the grid cell.  If includeEdges==true, then the grid cell edges are inclusive.  Otherwise
           the grid cell edges are not inclusive.  For example, if includeEdges==true
           and an Laguerre edge target point lies on the grid cell edge, then no
           exit edge crossing will be returned because the Laguerre edge doesn't actually
           leave the grid cell.  If includeEdges==false, then the Laguerre edge target
           point would be considered outside the grid cell and it would be returned
           as an exit point.
       */
       // std::pair<std::shared_ptr<Point_2>, std::shared_ptr<Point_2>> EdgeCrossings(PolygonEdge const& edge,
       //                                                                             bool               includeEdges);


       /** Given two intersection points on the boundary of the current grid cell, this function
           returns a vector of corner points between the intersection points if the
           grid cell boundary was traversed in the counter-clockwise direction
           between prevCross and nextCross.
       */
       void AddCorners(Point_2 const& nextPt, std::vector<Point_2> &polyPts) const;


       void UpdateCellBounds();

  };

} // namespace sdot


#endif // POLYGONRASTERIZE_H_
