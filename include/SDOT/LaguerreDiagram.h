#ifndef LAGUERREDIAGRAM_H
#define LAGUERREDIAGRAM_H

#include <Eigen/Core>

#include <vector>
#include <memory>
#include <unordered_map>

// includes for defining the Voronoi diagram adaptor
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Regular_triangulation_adaptation_traits_2.h>
#include <CGAL/Regular_triangulation_adaptation_policies_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include "SDOT/BoundingBox.h"

namespace sdot{

/**
@class LaguerreDiagram
@brief Contains tools for constructing and querying a bounded Laguerre diagram in 2d.

*/
class LaguerreDiagram
{

public:

  // See https://doc.cgal.org/latest/Triangulation_2/Triangulation_2_2info_insert_with_pair_iterator_regular_2_8cpp-example.html#_a1
  // typedefs for defining the adaptor
  typedef CGAL::Exact_predicates_exact_constructions_kernel              K;
  typedef CGAL::Regular_triangulation_vertex_base_2<K>                   Vbase;
  typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K,Vbase> Vb;
  typedef CGAL::Regular_triangulation_face_base_2<K>                     Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                    Tds;
  typedef CGAL::Regular_triangulation_2<K, Tds>                          DT;

//  typedef CGAL::Regular_triangulation_2<K, Tds>                       DT;
  typedef CGAL::Regular_triangulation_adaptation_traits_2<DT>         AT;
  //typedef CGAL::Regular_triangulation_degeneracy_removal_policy_2<DT> AP;
  typedef CGAL::Voronoi_diagram_2<DT,AT>                           PowerDiagram;

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


  /** Constructs the bounded Laguerre diagram.
      @param[in] pts A 2xN matrix with each column defining a seed point in the diagram.
      @param[in] costs A vector of length N containing the
       costs (squared radius) associated with each point.
      @param[in] bndryPts A 2xM matrix defining a boundary polygon around the domain.  Each column represents a point, which we assume are ordered counter-clockwise.
  */
  LaguerreDiagram(double xBndLeftIn,   double xBndRightIn,
                  double yBndBottomIn, double yBndTopIn,
                  Eigen::Matrix2Xd  const& pts,
                  Eigen::VectorXd  const& costs);

  LaguerreDiagram(BoundingBox const& bbox,
                  Eigen::Matrix2Xd  const& pts,
                  Eigen::VectorXd  const& costs);

  /** Returns a CGAL Polygon_2 object representing one of the Laguerre cells. */
  std::shared_ptr<Polygon_2> const& GetCell(int ind) const{return laguerreCells.at(ind);};

  /**
  Returns the vertices of a single cell in the Laguerre diagram.
  @param[in] ind The index of the cell
  @return A matrix containing cell vertices.  Each column is a point in 2d.
  */
  Eigen::Matrix2Xd GetCellVertices(int ind) const;


  /**
    Internal edges in the Laguerre diagram are the intersection of two Laguerre
    cells.  The edge is therefore the dividing line between the two cells.  This
    function returns an internal edge given the indices of two cells in the
    diagram.
  */
  std::vector<std::tuple<unsigned int, Point_2, Point_2>> const& InternalEdges(unsigned int cellInd) const{return internalEdges.at(cellInd);};

  /** Returns the center of mass of the cell.  Using in Lloyd's algorithm for
      centroidal diagrams.  For a polygon covering a region $\Omega$, the center
      of mass is given by $\int_\Omega x dx / \int_\Omega dx$, where $x\in\mathbb{R}^2$ is the position.
      The GetCellCenter function computes this integral by breaking the convex
      polygon into triangles and computing summing the integrals over each
      triangle.
      @param[in] cellInd The index of the Laguerre cell of interest.
      @return The center of mass of the Laguerre cell.
  */
  Eigen::Vector2d CellCentroid(unsigned int cellInd) const;

  /** Repeatedly calls CellCentroid to compute the centers of mass for all
      Laguerre cells.  Column $i$ of the output matrix contains the centroid of
      the ith Laguerre cell.
  */
  Eigen::Matrix2Xd Centroids() const;

  /** Returns the seed points that were used to construct the Laguerre diagram. */
  Eigen::Matrix2Xd SeedPts() const;

  /** Returns the number of Laguerre cells. */
  int NumCells() const{return laguerreCells.size();};

  /** Returns the underlying CGAL PowerDiagram object. */
  PowerDiagram* BaseDiagram(){return &unboundedDiagram;};

  struct ConstructionException : public std::runtime_error
  {
    ConstructionException(std::string msg) : std::runtime_error(msg.c_str()){};
  };

  BoundingBox const& BoundBox() const{return bbox;};

private:

  const double compTol = 1e-14;

  void CreateBoundaryPolygon(Eigen::Matrix2Xd const& bndryPts);

  void CreateUnboundedDiagram(Eigen::Matrix2Xd const& pts,
                              Eigen::VectorXd  const& costs);

  void CreateBoundedCells(Eigen::Matrix2Xd const& pts);

  std::shared_ptr<Polygon_2> BoundOneCell(PowerDiagram::Face_handle const& face);

  void AddInternalEdges(PowerDiagram::Face_handle const& face);

  /** Returns true if the halfEdge has a source node inside the bounding box.
  */
  bool HasInternalSource(Ccb_halfedge_circulator halfEdge) const;

  /** Returns true if the halfEdge has a target node inside the bounding box.
  */
  bool HasInternalTarget(Ccb_halfedge_circulator halfEdge) const;

  /// The number of points used to construct the Laguere diagram.
  int numPts;

  /// A vector of bounded Laguerre cells.  These defines bounded diagram.
  std::vector<std::shared_ptr<Polygon_2>> laguerreCells;

  /// Stores a polygon defining the bounds of the domain of interest
  BoundingBox bbox;

  // double xBndLeft, xBndRight, yBndBottom, yBndTop;
//  Polygon_2 boundPoly;

  /// Stores the unbounded polygon created by CGAL
  PowerDiagram unboundedDiagram;

  /* Stores each internal edge of the Laguerre diagram internalEdges[i][j] returns
     the edge in the Laguerre diagram that divides cell i and cell j.
  */
  std::vector<std::vector<std::tuple<unsigned int, Point_2,Point_2>>> internalEdges;

  const double infDist;

};

} // namespace sdot

#endif
