#ifndef LAGUERREDIAGRAM_H
#define LAGUERREDIAGRAM_H

#include <Eigen/Core>

#include <vector>
#include <memory>
#include <unordered_map>
#include <chrono>

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
#include "SDOT/Distribution2d.h"
#include "SDOT/Polygon.h"
#include "SDOT/PolygonRasterize.h"

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


  /** Creates a centroidal Voronoi/power diagram using Lloyd's algorithm.
      @param[in] bbox A bounding box over the region of interest.
      @param[in] initialPts The points used to construct the initial diagram.
      @param[in] prices A vector of prices for each point.  If all the prices are
                        a Voronoi diagram is constructed.
      @param[in] maxIts The maximum number of iterations in Lloyd's algorithm
      @param[in] tol Convergence tolerance based on maximum absolute difference
                     between cell centroids and seed points.
  */
  static std::shared_ptr<LaguerreDiagram> BuildCentroidal(BoundingBox const& bbox,
                                                          Eigen::Matrix2Xd  const& initialPts,
                                                          Eigen::VectorXd  const& prices,
                                                          unsigned int maxIts=200,
                                                          double tol=1e-3,
                                                          std::shared_ptr<Distribution2d> const& dist=nullptr);

  /** Creates a centroidal Voronoi/Power diagram using Lloyd's algorithm.  Initializes
      the points with latin hypercube sampling.
      @param[in] bbox The bounding box around the domain of interest.
      @param[in] numPts The number of points in the diagram.
      @param[in] maxIts The maximum number of iterations in Lloyd's algorithm
      @param[in] tol Convergence tolerance based on maximum absolute difference
                     between cell centroids and seed points.
  */
  static std::shared_ptr<LaguerreDiagram> BuildCentroidal(BoundingBox const& bbox,
                                                          unsigned int       numPts,
                                                          unsigned int maxIts=200,
                                                          double tol=1e-3,
                                                          std::shared_ptr<Distribution2d> const& dist=nullptr);

  /** Generates numPts samples in the bounding box using lating hypercube sampling. */
  static Eigen::Matrix2Xd LatinHypercubeSample(BoundingBox const& bbox,
                                               unsigned int       numPts);

  /** Returns a CGAL Polygon_2 object representing one of the Laguerre cells. */
  std::shared_ptr<Polygon> GetCell(int ind) const{return std::make_shared<Polygon>(laguerreCells.at(ind));};

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
      The CellCentroid function computes this integral by breaking the convex
      polygon into triangles and computing summing the integrals over each
      triangle.
      @param[in] cellInd The index of the Laguerre cell of interest.
      @return The center of mass of the Laguerre cell.
  */
  Eigen::Vector2d CellCentroid(unsigned int cellInd) const;

  /** Computes the weighted center of mass of the cell, i.e.,
  \f[
  C = \frac{\int_\Omega x \rho(x) dx}{\int_\Omega \rho(x) dx}
  \f]

  @param[in] cellInd The index of the Laguerre cell of interest.
  @param[in] dist The density (discretized onto grid cells) defining \f$\rho(x)\f$.
  */
  Eigen::Vector2d CellCentroid(unsigned int cellInd, std::shared_ptr<Distribution2d> const& dist) const;

  /** Returns the area of one of the Laguerre cells, possibly weighted by the
      the two-dimensional distribution.
  */
  double CellArea(unsigned int cellInd,
                  std::shared_ptr<Distribution2d> const& dist=nullptr) const;

  /** Integrates a function f(x) defined on the Laguerre cell.   If the dist
      input is defined, then this functions computed the weighted integral
      $\int_A f(x)w(x) dx$ for the distribution $w(x)$.
  */

  template<typename TriType>
  auto IntegrateOverCell(unsigned int cellInd,
                           TriType triFunc) const;

  template<typename TriType, typename RectType>
  auto IntegrateOverCell(unsigned int cellInd,
                           TriType triFunc,
                           RectType rectFunc) const{return IntegrateOverCell(cellInd,triFunc);};

  template<typename TriType, typename RectType>
  auto IntegrateOverCell(unsigned int cellInd,
                           TriType triFunc,
                           RectType rectFunc,
                           std::shared_ptr<Distribution2d> const& dist) const;

  /** Repeatedly calls CellCentroid to compute the centers of mass for all
      Laguerre cells.  Column $i$ of the output matrix contains the centroid of
      the ith Laguerre cell.
      @param[in] dist (Optional) Discretized density.  If provided, the centroids will be centers of mass wrt to this density.
  */
  Eigen::Matrix2Xd Centroids(std::shared_ptr<Distribution2d> const& dist=nullptr) const;

  /** Repeatedly calls CellArea to compute the areas for all Laguerre cells.
      @param[in] dist (Optional) Discretized density.  If provided, the mass of
                 each cell obtained by integrating this density will be returned.
  */
  Eigen::VectorXd Areas(std::shared_ptr<Distribution2d> const& dist=nullptr) const;

  /** Returns the seed points that were used to construct the Laguerre diagram. */
  Eigen::Matrix2Xd SeedPts() const;

  /** Returns the prices associated with each Laguerre cell. */
  Eigen::VectorXd const& Prices() const;

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

  void CreateBoundedCells();

  /** This function clips a polygon to the halfspace defined by the Laguerre diagram face handle.
      This sdefines one part of the Sutherland-Hodgman algorthm used here to bound
      the unbounded Laguerre diagram.
      @param[in] poly The polygon we want to clip.
      @param[in] halfEdge The half edge describing one edge (possibly unbounded) in the Laguerre cell.
      @returns A shared pointer to the polygon.     Note that this might just be the same as the input if no clipping is necessary.
  */
  std::shared_ptr<LaguerreDiagram::Polygon_2> ClipToHalfspace(std::shared_ptr<LaguerreDiagram::Polygon_2> const& poly,
                                                              Ccb_halfedge_circulator halfEdge);



  std::shared_ptr<Polygon_2> BoundOneCell(PowerDiagram::Face_handle const& face);

  void AddInternalEdges(PowerDiagram::Face_handle const& face);

  /** Returns true if the halfEdge has a source node inside the bounding box.
  */
  bool HasInternalSource(Ccb_halfedge_circulator halfEdge) const;

  /** Returns true if the halfEdge has a target node inside the bounding box.
  */
  bool HasInternalTarget(Ccb_halfedge_circulator halfEdge) const;

  /** Returns the vertices on the ends of a halfEdge.  Fills in "large" values
      for edges whose source or target is at infinity.
  */
  std::pair<Point_2,Point_2> GetEdgeVerts(Ccb_halfedge_circulator& halfEdge);

  /// The number of points used to construct the Laguere diagram.
  int numPts;

  Eigen::Matrix2Xd seedPts;

  Eigen::VectorXd prices;

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



template<typename TriType, typename RectType>
auto LaguerreDiagram::IntegrateOverCell(unsigned int cellInd,
                                          TriType triFunc,
                                          RectType rectFunc,
                                          std::shared_ptr<Distribution2d> const& dist) const
{

  auto start = std::chrono::steady_clock::now();
  if(dist==nullptr)
    return IntegrateOverCell(cellInd, triFunc);

  auto& grid = dist->Grid();

  Eigen::Vector2d pt1(2), pt2(2), pt3(2);
  pt1 << 0.0, 0.0;

  // Should set result to 0
  auto zeroVal = rectFunc(pt1,pt1);

  // Loop over the grid cells in this Laguerre cell
  if(laguerreCells.at(cellInd)==nullptr)
    return zeroVal;

  if(laguerreCells.at(cellInd)->size()==0)
    return zeroVal;

  auto start2 = std::chrono::steady_clock::now();
  PolygonRasterizeIter gridIter(dist->Grid(), laguerreCells.at(cellInd));
  auto end2 = std::chrono::steady_clock::now();
  double setup_time = std::chrono::duration<double>(end2-start2).count();

  unsigned int xInd, yInd;
  double gridCellDens;
  auto result = zeroVal;


  double bndry_total = 0;
  int num_bndry = 0;
  double int_total = 0;
  double incr_total = 0;
  int num_int = 0;
  int num_incr = 0;

  while(gridIter.IsValid()){

    xInd = gridIter.Indices().first;
    yInd = gridIter.Indices().second;

    // The probability in this grid cell
    gridCellDens = dist->Density(xInd,yInd);

    if(gridCellDens>0.0){

      // The area of the intersection of this grid cell and the Laguerre cell
      auto interResult = zeroVal;

      if(gridIter.IsBoundary()){
        start2 = std::chrono::steady_clock::now();

        // Break the intersection polygon into triangles and add contributions from each triangle
        std::shared_ptr<PolygonRasterizeIter::Polygon_2> overlapPoly = gridIter.OverlapPoly();
        assert(overlapPoly);
        assert(overlapPoly->size()>2); // <- Makes sure there is at least 3 nodes in the polygon

        auto beginVert = overlapPoly->vertices_begin();
        auto vert1 = beginVert;
        vert1++;
        auto vert2 = vert1;
        vert2++;

        pt1[0] = CGAL::to_double( beginVert->x() );
        pt1[1] = CGAL::to_double( beginVert->y() );

        for(; vert2!=overlapPoly->vertices_end(); vert2++, vert1++)
        {
          pt2[0] = CGAL::to_double( vert1->x() );
          pt2[1] = CGAL::to_double( vert1->y() );
          pt3[0] = CGAL::to_double( vert2->x() );
          pt3[1] = CGAL::to_double( vert2->y() );

          interResult += triFunc(pt1,pt2,pt3);
        }

        end2 = std::chrono::steady_clock::now();
        bndry_total += std::chrono::duration<double>(end2-start2).count();
        num_bndry++;

      }else{
        start2 = std::chrono::steady_clock::now();

        pt1[0] = grid->xMin + grid->dx*xInd;
        pt1[1] = grid->yMin + grid->dy*yInd;
        pt2[0] = grid->xMin + grid->dx*(xInd+1);
        pt2[1] = grid->yMin + grid->dy*(yInd+1);

        interResult += rectFunc(pt1,pt2);

        end2 = std::chrono::steady_clock::now();
        int_total += std::chrono::duration<double>(end2-start2).count();
        num_int++;

      }

      result += interResult*gridCellDens;
    }

    start2 = std::chrono::steady_clock::now();

    gridIter.Increment();

    end2 = std::chrono::steady_clock::now();
    incr_total += std::chrono::duration<double>(end2-start2).count();
    num_incr++;
  }

  return result;
}

/* Implementation of templated function. */
template<typename TriType>
auto LaguerreDiagram::IntegrateOverCell(unsigned int cellInd, TriType triFunc) const
{
  Eigen::Vector2d pt1(2), pt2(2), pt3(2);
  pt1 << 0.0, 0.0;
  if(laguerreCells.at(cellInd)==nullptr){
    return triFunc(pt1,pt1,pt1);
  }

  auto beginVert = laguerreCells.at(cellInd)->vertices_begin();
  auto vert1 = beginVert;
  vert1++;
  auto vert2 = vert1;
  vert2++;

  pt1(0) = CGAL::to_double( beginVert->x() );
  pt1(1) = CGAL::to_double( beginVert->y() );

  //interArea = gridCellDens*CGAL::to_double( overlapPoly->area() );

  pt2(0) = CGAL::to_double( vert1->x() );
  pt2(1) = CGAL::to_double( vert1->y() );
  pt3(0) = CGAL::to_double( vert2->x() );
  pt3(1) = CGAL::to_double( vert2->y() );

  auto result = triFunc(pt1,pt2,pt3);
  vert2++;
  vert1++;

  for(; vert2!=laguerreCells.at(cellInd)->vertices_end(); vert2++, vert1++)
  {
    pt2(0) = CGAL::to_double( vert1->x() );
    pt2(1) = CGAL::to_double( vert1->y() );
    pt3(0) = CGAL::to_double( vert2->x() );
    pt3(1) = CGAL::to_double( vert2->y() );

    result += triFunc(pt1,pt2,pt3);
  }

  return result;
}


} // namespace sdot

#endif
