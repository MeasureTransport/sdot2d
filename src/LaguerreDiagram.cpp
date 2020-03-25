#include "SDOT/LaguerreDiagram.h"

// standard includes
#include <iostream>

// Additional CGAL includes
#include <CGAL/squared_distance_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Boolean_set_operations_2.h>

using namespace sdot;

void LaguerreDiagram::CreateBoundaryPolygon(Eigen::Matrix2Xd const& bndryPts)
{
  // Set up the domain boundary.  Points should go around the bounding polygon in counter-clockwise order
  for(int i=0; i<bndryPts.cols(); ++i)
    boundPoly.push_back(Point_2(bndryPts(0,i),bndryPts(1,i)));
}

void LaguerreDiagram::CreateUnboundedDiagram(Eigen::Matrix2Xd const& pts,
                                             Eigen::VectorXd  const& costs)
{
  assert(costs.size()==numPts);

  std::vector<Site_2> wpoints(numPts);
  for(int i=0; i<numPts; ++i)
    wpoints.at(i) = Site_2(Point_2(pts(0,i), pts(1,i)), costs(i));

  unboundedDiagram.insert(wpoints.begin(), wpoints.end());

  assert(unboundedDiagram.is_valid());
}


std::shared_ptr<LaguerreDiagram::Polygon_2> LaguerreDiagram::BoundOneCell(PowerDiagram::Face_handle const& face)
{

      // Get the point in the Regular triangulation at the center of this Laguerre cell
      auto vs = face->dual()->point();

      // Loop over the edges in the face
      Ccb_halfedge_circulator halfEdgeStart = face->ccb();
      Ccb_halfedge_circulator halfEdge = halfEdgeStart;

      std::vector<Point_2> polyPts;

      // Loop over all of the edges in the face
      do {

        // Get the node in the triangulation on the other side of this edge
        auto vt = halfEdge->twin()->face()->dual()->point();

        // Compute the direction of dividing line between this cell and the other cell
        Vector_2 dir(vt.y()-vs.y(), -(vt.x()-vs.x()));

        // Compute a vector that points away from the edge
        Vector_2 perpDir = vs.point() - vt.point();

        // If the edge has a source, add the source point to the polygon
        if(halfEdge->has_source()){
          polyPts.push_back(halfEdge->source()->point());

          if(!halfEdge->has_target()){
            /* If it doesn't have a target point, then it will be outside the boundary.
              To accomplish this, we add a point far far outside the bounding polygon
              and then later do an intersection with the bounding poly.
             */
            polyPts.push_back(halfEdge->source()->point() - infDist*dir);
            polyPts.push_back(halfEdge->source()->point() - infDist*dir + infDist*perpDir);
          }

        // This edge starts at infinity, but has a target
        }else if(halfEdge->has_target()){
          //polyPts.push_back(halfEdge->target()->point() + infDist*dir + infDist*perpDir);
          polyPts.push_back(halfEdge->target()->point() + infDist*dir);
          //polyPts.push_back(halfEdge->target()->point());

        // This edge has neither a source nor a target
        }else{

          // Find the weighted mid point between the triangulation sites
          double ts = CGAL::to_double( CGAL::squared_distance(vt.point(), vs.point()).exact() );
          double h = (0.5/ts)*CGAL::to_double(ts + vs.weight() - vt.weight());
          Point_2 midPt = vs.point() + h*(vt.point()-vs.point());

          polyPts.push_back(midPt - infDist*dir);
          polyPts.push_back(midPt - infDist*dir + infDist*perpDir);
          polyPts.push_back(midPt + infDist*dir + infDist*perpDir);
          polyPts.push_back(midPt + infDist*dir);

        }

      } while ( ++halfEdge != halfEdgeStart);

      // bigLag is a polygon containing the laguerre cell with large values replacing infinity
      Polygon_2 bigLag(polyPts.begin(), polyPts.end());
      std::list<Polygon_with_holes_2> temp;
      CGAL::intersection(bigLag, boundPoly, std::back_inserter(temp));

      // There should only be one polygon intersection because everything is simple and convex
      assert(std::distance(temp.begin(),temp.end())==1);

      return std::make_shared<Polygon_2>(temp.begin()->outer_boundary());
}

void LaguerreDiagram::CreateBoundedCells(Eigen::Matrix2Xd const& pts)
{

  // Figure out which faces correspond to the original points
  std::vector<Face_handle> faces;
  for(int i=0; i<numPts; ++i){
    Locate_result lr = unboundedDiagram.locate( Point_2(pts(0,i), pts(1,i)) );
    Face_handle* f = boost::get<Face_handle>(&lr);
    assert(f);
    faces.push_back(*f);
  }

  // Loop over the cells in the Laguerre diagram
  for(auto face : faces)
    laguerreCells.push_back( BoundOneCell(face) );

}

LaguerreDiagram::LaguerreDiagram(Eigen::Matrix2Xd const& pts,
                                 Eigen::VectorXd  const& costs,
                                 Eigen::Matrix2Xd const& bndryPts)
{
  numPts = pts.cols();
  assert(costs.size()==numPts);

  CreateBoundaryPolygon(bndryPts);

  CreateUnboundedDiagram(pts,costs);

  CreateBoundedCells(pts);
}
