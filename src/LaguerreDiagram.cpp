#include "SDOT/LaguerreDiagram.h"

// standard includes
#include <iostream>

// Additional CGAL includes
#include <CGAL/squared_distance_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Boolean_set_operations_2.h>

using namespace sdot;

LaguerreDiagram::LaguerreDiagram(double xBndLeftIn,   double xBndRightIn,
                                 double yBndBottomIn, double yBndTopIn,
                                 Eigen::Matrix2Xd const& pts,
                                 Eigen::VectorXd  const& costs) : xBndLeft(xBndLeftIn),
                                                                  xBndRight(xBndRightIn),
                                                                  yBndBottom(yBndBottomIn),
                                                                  yBndTop(yBndTopIn)
{
  numPts = pts.cols();
  assert(costs.size()==numPts);

  Eigen::Matrix2Xd bndryPts(2,4);
  bndryPts << xBndLeft,   xBndRight,  xBndRight, xBndLeft,
              yBndBottom, yBndBottom, yBndTop,   yBndTop;

  CreateBoundaryPolygon(bndryPts);

  CreateUnboundedDiagram(pts,costs);

  CreateBoundedCells(pts);
}

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

  //std::vector<Site_2> wpoints(numPts);
  for(unsigned int i=0; i<numPts; ++i)
    unboundedDiagram.insert( Site_2(Point_2(pts(0,i), pts(1,i)), costs(i)) )->dual()->info()= i;
    //wpoints.at(i) = ;

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
  internalEdges.clear();
  internalEdges.resize(numPts);

  // Figure out which faces correspond to the original points
  std::vector<Face_handle> faces(numPts);
  for(auto faceIter = unboundedDiagram.faces_begin(); faceIter!=unboundedDiagram.faces_end(); ++faceIter){
    faces.at(faceIter->dual()->info()) = *faceIter;
  }

  for(auto& face : faces){
    laguerreCells.push_back( BoundOneCell(face) );
    AddInternalEdges(face);
  }
}

void LaguerreDiagram::AddInternalEdges(PowerDiagram::Face_handle const& face)
{

  // Get the point in the Regular triangulation at the center of this Laguerre cell
  auto vs = face->dual()->point();

  // Loop over the edges in the face
  Ccb_halfedge_circulator halfEdgeStart = face->ccb();
  Ccb_halfedge_circulator halfEdge = halfEdgeStart;

  Point_2 srcPt, tgtPt;

  // Loop over all of the edges in the face
  do {

    // Get the node in the triangulation on the other side of this edge
    auto vt = halfEdge->twin()->face()->dual()->point();

    // Compute the direction of dividing line between this cell and the other cell
    Vector_2 dir(vt.y()-vs.y(), -(vt.x()-vs.x()));


    // If the edge has a source, add the source point to the polygon
    if(halfEdge->has_source()){
      srcPt = halfEdge->source()->point();

      if(halfEdge->has_target()){
        tgtPt = halfEdge->target()->point();
      }else{
        tgtPt = halfEdge->source()->point() - infDist*dir;
      }

    // This edge starts at infinity, but has a target
    }else if(halfEdge->has_target()){
      srcPt = halfEdge->target()->point() + infDist*dir;
      tgtPt = halfEdge->target()->point();

    // This edge has neither a source nor a target
    }else{

      // Find the weighted mid point between the triangulation sites
      double ts = CGAL::to_double( CGAL::squared_distance(vt.point(), vs.point()).exact() );
      double h = (0.5/ts)*CGAL::to_double(ts + vs.weight() - vt.weight());
      Point_2 midPt = vs.point() + h*(vt.point()-vs.point());

      srcPt = midPt - infDist*dir;
      tgtPt = midPt + infDist*dir;
    }

    bool hasIntersection = ClipToBoundary(srcPt, tgtPt);

    // If the edge does not intersect the domain at all, just skip it
    if(hasIntersection){
      // This edge has both a source and target, so it is an internal edge.  Save it!
      unsigned int currInd = face->dual()->info();
      unsigned int otherInd = halfEdge->twin()->face()->dual()->info();
      internalEdges.at(currInd).push_back(std::make_tuple(otherInd, srcPt, tgtPt));
    }

  } while ( ++halfEdge != halfEdgeStart);

}


bool LaguerreDiagram::ClipToBoundary(Point_2& srcPt, Point_2& tgtPt) const
{

  double sx = CGAL::to_double(srcPt.x());
  double sy = CGAL::to_double(srcPt.y());
  double tx = CGAL::to_double(tgtPt.x());
  double ty = CGAL::to_double(tgtPt.y());

  // Check if it's vertical
  if(std::abs(sx-tx)<compTol){

    // Check to see if it crosses the domain at all
    if((sx>xBndLeft-compTol)&&(sx<xBndRight+compTol)&&(std::max(sy,ty)>yBndBottom-compTol)&&(std::min(sy,ty)<yBndTop+compTol)){
      srcPt = Point_2(sx, std::min(yBndTop, std::max(yBndBottom, sy)));
      tgtPt = Point_2(tx, std::min(yBndTop, std::max(yBndBottom, ty)));
      return true;
    }else{
      return false;
    }

  // Check if it's horizontal
  }else if(std::abs(sy-ty)<compTol){

    if((sy>yBndBottom-compTol)&&(sy<yBndTop-compTol)&&(std::max(sx,tx)>xBndLeft-compTol)&&(std::min(sx,tx)<xBndRight+compTol)){
      srcPt = Point_2(std::min(xBndRight, std::max(xBndLeft, sx)), sy);
      tgtPt = Point_2(std::min(xBndRight, std::max(xBndLeft, tx)), ty);
      return true;
    }else{
      return false;
    }

  // Otherwise it has a non zero finite slope
  }else{


    // Parameterize the line as srcPt + t*(tgtPt-srcPt) with t in [0,1]
    Eigen::Matrix<double,4,1> qs(4), ps(4);

    ps(0) = sx-tx;
    ps(1) = tx-sx;
    ps(2) = sy-ty;
    ps(3) = ty-sy;

    qs(0) = (sx-xBndLeft);///(tx-sx);
    qs(1) = (xBndRight-sx);//(tx-sx);
    qs(2) = (sy-yBndBottom);///(ty-sy);
    qs(3) = (yBndTop-sy);///(ty-sy);


    double u1 = 0.0;
    for(unsigned int i=0; i<4; ++i){
      if(ps(i)<0.0)
        u1 = std::max(u1, qs(i)/ps(i));
    }

    double u2 = 1.0;
    for(unsigned int i=0; i<4; ++i){
      if(ps(i)>0.0)
        u2 = std::min(u2, qs(i)/ps(i));
    }

    if(u1>u2){
      return false;
    }else{
      srcPt = Point_2(sx + u1*ps(1), sy + u1*ps(3));
      tgtPt = Point_2(sx + u2*ps(1), sy + u2*ps(3));
      return true;
    }
  }

}
