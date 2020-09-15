#ifndef BOUNDINGBOX_H
#define BOUNDINGBOX_H

// includes for defining the Voronoi diagram adaptor
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>



namespace sdot {

  /** @class BoundingBox
      This class provides tools for working with a rectangular bounding box.
      Specifically, clipping line segments and polygons to the bounding box.
  */
  class BoundingBox{
  public:

    // typedefs for defining the adaptor
    typedef CGAL::Exact_predicates_exact_constructions_kernel K;

    // typedef for the result type of the point location
    typedef K::Point_2                    Point_2;
    typedef K::Vector_2                   Vector_2;
    typedef CGAL::Polygon_2<K>            Polygon_2;
    typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;

    enum class Side { Left, Right, Bottom, Top };

    BoundingBox(double xmin_, double xmax_,
                double ymin_, double ymax_);

    /** Given a line segment from a point src to a point tgt, this function clips
        the line segment to the bounding box.  Returns true if the segment intersects
        the bounding box and false if the line segment is completely outside.
    */
    bool ClipSegment(Point_2 &src, Point_2 &tgt) const;

    /** Computes the intersection of a polygon with the bounding box.  Uses
        an implementation of the Sutherland–Hodgman algorithm.
    */
    std::shared_ptr<Polygon_2> ClipPolygon(std::shared_ptr<Polygon_2> const& poly) const;

    /** Assuming the last point in polyPts and nextPt lie on the bounding box of
        the grid, this function adds any corner points needed to move around the
        bounding box from polyPts.back() to nextPt in a counter-clockwise direction.
    */
    void AddCorners(Point_2 const& nextPt, std::vector<Point_2> &polyPts) const;

    bool IsInside(Point_2 const& pt) const;

    bool OnEdge(Point_2 const& pt) const;

    const double xMax, xMin, yMin, yMax;

  private:

    std::vector<Point_2> ClipRight(std::vector<Point_2> const& poly) const;
    std::vector<Point_2> ClipLeft(std::shared_ptr<Polygon_2> const& poly) const;
    std::vector<Point_2> ClipTop(std::vector<Point_2> const& poly) const;
    std::vector<Point_2> ClipBottom(std::vector<Point_2> const& poly) const;

    Polygon_2 boxPoly;
    const double compTol = 1e-12;

  }; // class BoundingBox

};

#endif // #ifndef BOUNDINGBOX_H
