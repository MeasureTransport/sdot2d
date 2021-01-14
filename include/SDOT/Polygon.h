#ifndef POLYGON_H_
#define POLYGON_H_

#include <Eigen/Core>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

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
  typedef K::Point_2                   Point_2;
  typedef CGAL::Polygon_2<K>           Polygon_2;

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

private:
  std::shared_ptr<Polygon_2> cgalPoly;

}; // class Polygon

} // namespace sdot

#endif // #ifndef POLYGON_H_
