#include <Eigen/Core>

#include <vector>
#include <memory>

// includes for defining the Voronoi diagram adaptor
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Regular_triangulation_adaptation_traits_2.h>
#include <CGAL/Regular_triangulation_adaptation_policies_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>


namespace sdot{

/**
@class LaguerreDiagram
@brief Contains tools for constructing and querying a bounded Laguerre diagram in 2d.

*/
class LaguerreDiagram
{
  // typedefs for defining the adaptor
  typedef CGAL::Exact_predicates_exact_constructions_kernel           K;
  typedef CGAL::Regular_triangulation_2<K>                            DT;
  typedef CGAL::Regular_triangulation_adaptation_traits_2<DT>         AT;
  typedef CGAL::Regular_triangulation_degeneracy_removal_policy_2<DT> AP;
  typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                           PowerDiagram;

  // typedef for the result type of the point location
  typedef AT::Site_2                    Site_2;
  typedef AT::Point_2                   Point_2;
  typedef K::Vector_2                   Vector_2;
  typedef CGAL::Polygon_2<K>            Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;

  typedef PowerDiagram::Locate_result             Locate_result;
  typedef PowerDiagram::Vertex_handle             Vertex_handle;
  typedef PowerDiagram::Face_handle               Face_handle;
  typedef PowerDiagram::Halfedge_handle           Halfedge_handle;
  typedef PowerDiagram::Ccb_halfedge_circulator   Ccb_halfedge_circulator;

public:

  /** Constructs the bounded Laguerre diagram.
      @param[in] pts A 2xN matrix with each column defining a seed point in the diagram.
      @param[in] costs A vector of length N containing the costs (squared radius) associated with each point.
      @param[in] bndryPts A 2xM matrix defining a boundary polygon around the domain.  Each column represents a point, which we assume are ordered counter-clockwise.
  */
  LaguerreDiagram(Eigen::Matrix2Xd const& pts,
                  Eigen::VectorXd  const& costs,
                  Eigen::Matrix2Xd const& bndryPts);

  std::shared_ptr<Polygon_2> GetCell(int ind){return laguerreCells.at(ind);};

  int NumCells() const{return laguerreCells.size();};

private:

  void CreateBoundaryPolygon(Eigen::Matrix2Xd const& bndryPts);

  void CreateUnboundedDiagram(Eigen::Matrix2Xd const& pts,
                              Eigen::VectorXd  const& costs);

  void CreateBoundedCells(Eigen::Matrix2Xd const& pts);

  std::shared_ptr<Polygon_2> BoundOneCell(PowerDiagram::Face_handle const& face);

  /// The number of points used to construct the Laguere diagram.
  int numPts;

  /// A vector of bounded Laguerre cells.  These defines bounded diagram.
  std::vector<std::shared_ptr<Polygon_2>> laguerreCells;

  /// Stores a polygon defining the bounds of the domain of interest
  Polygon_2 boundPoly;

  /// Stores the unbounded polygon created by CGAL
  PowerDiagram unboundedDiagram;

  const double infDist = 1e5;

};

} // namespace sdot
