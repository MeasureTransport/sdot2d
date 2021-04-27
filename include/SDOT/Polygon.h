#ifndef POLYGON_H_
#define POLYGON_H_

#include <Eigen/Core>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <vector>

namespace sdot{

/** @class Polygon
    @brief Provides a wrapper, with some extra capabilities, around the CGAL POLYGON_2 class.
    @details This is useful for exposing the CGAL polygon object to python, as
             well as partitioning arbitrary polygons into convex components.
*/
class Polygon{

public:

  typedef CGAL::Exact_predicates_exact_constructions_kernel K;
  typedef K::Point_2                    Point_2;
  typedef CGAL::Polygon_2<K>            Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;

  /** Construct an SDOT::Polygon from polygon vertices.  If the vertices are
      in clockwise order, the order is reversed to make the orientation counter-clockwise.
  */
  Polygon(Eigen::Matrix2Xd const& pts);

  /** Construct an SDOT::Polygon from a shared pointer to a CGAL::Polygon_2 object. */
  Polygon(std::shared_ptr<Polygon_2> const& cgalPolyIn);

  /** Returns the vertices of the polygon. */
  Eigen::Matrix2Xd GetVertices() const;

  /** Uses CGAL to construct a convex partition of the polygon.  If this polygon
      is alredy convex, a pointer to *this is returned.
  */
  std::vector<std::shared_ptr<Polygon>> ConvexPartition() const;

  /** Returns the convex hull of the polygon.  If this polygon is already convex,
      a pointer to *this is returned.
  */
  std::shared_ptr<Polygon> ConvexHull() const;

  /** Returns true if the polygon is convex.  False otherwise. */
  bool IsConvex() const{return cgalPoly->is_convex();};

  /** Returns the CGAL Polygon_2 object. */
  std::shared_ptr<Polygon_2> ToCGAL(){return cgalPoly;};

  /** Returns the intersection polygon of this polygon with another. */
  Polygon Intersection(Polygon const& otherPoly) const;

  /** Returns true if this polygon intersects the other polygon and false otherwise. */
  bool Intersects(Polygon const& otherPoly) const;

  /** Returns the unweighted area of the polygon. */
  double Area() const{return area;};

  /** Returns the second moment of area around an axis perpendicular to 2d plane
      and passing through the centroid of the polygon.   Multiplying the output
      of this function by the mass of the particle gives the moment of inertia.
  */
  double SecondMoment() const;

  /** Returns the centroid of the polygon. */
  Eigen::Vector2d const& Centroid() const{return centroid;};

  /** Translate all vertices in the polygon. */
  Polygon& Translate(Eigen::Vector2d const& step);

  /** Rotate the polygon around its centroid. */
  Polygon& Rotate(double angle);

private:
  std::shared_ptr<Polygon_2> cgalPoly;

  double xmin, ymin, xmax, ymax; // The bounding box of the polygon
  Eigen::Vector2d centroid;
  double maxRadius=-1; // The radius of a circle centered at the centroid that contains all vertices in the polygon
  double area=-1; // The area of the polygon


  /** Computes the centroid and area  of the polygon. */
  void _ComputeCentroidArea();
  void _ComputeRadius();

}; // class Polygon

} // namespace sdot

#endif // #ifndef POLYGON_H_
