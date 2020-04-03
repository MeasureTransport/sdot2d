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
                                 Eigen::VectorXd  const& costs) : bbox(xBndLeftIn, xBndRightIn, yBndBottomIn, yBndTopIn),
                                                                  infDist(2.0*std::sqrt(std::pow(xBndRightIn-xBndLeftIn,2.0) + std::pow(yBndTopIn-yBndBottomIn,2.0)))
{
  numPts = pts.cols();
  assert(costs.size()==numPts);

  CreateUnboundedDiagram(pts,costs);
  CreateBoundedCells(pts);
}

// void LaguerreDiagram::CreateBoundaryPolygon(Eigen::Matrix2Xd const& bndryPts)
// {
//   // Set up the domain boundary.  Points should go around the bounding polygon in counter-clockwise order
//   for(int i=0; i<bndryPts.cols(); ++i)
//     boundPoly.push_back(Point_2(bndryPts(0,i),bndryPts(1,i)));
// }

void LaguerreDiagram::CreateUnboundedDiagram(Eigen::Matrix2Xd const& pts,
                                             Eigen::VectorXd  const& costs)
{
  assert(costs.size()==numPts);

  //std::vector<Site_2> wpoints(numPts);
  for(unsigned int i=0; i<numPts; ++i){
    assert(costs(i)>=0);
    unboundedDiagram.insert( Site_2(Point_2(pts(0,i), pts(1,i)), std::sqrt(costs(i))) )->dual()->info()= i;
  }
    //wpoints.at(i) = ;

  assert(unboundedDiagram.is_valid());
}

bool LaguerreDiagram::HasInternalSource(Ccb_halfedge_circulator halfEdge) const
{
  if(!halfEdge->has_source()){
    return false;
  }else{
    double x = CGAL::to_double(halfEdge->source()->point().x());
    double y = CGAL::to_double(halfEdge->source()->point().y());
    return (x<bbox.xMax+compTol)&&(x>bbox.xMin-compTol)&&(y<bbox.yMax+compTol)&&(y>bbox.yMin-compTol);
  }
}

bool LaguerreDiagram::HasInternalTarget(Ccb_halfedge_circulator halfEdge) const
{
  if(!halfEdge->has_target()){
    return false;
  }else{
    double x = CGAL::to_double(halfEdge->target()->point().x());
    double y = CGAL::to_double(halfEdge->target()->point().y());
    return (x<bbox.xMax+compTol)&&(x>bbox.xMin-compTol)&&(y<bbox.yMax+compTol)&&(y>bbox.yMin-compTol);
  }
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
        // std::cout << "\nLaguerre edge: ";
        // if(!halfEdge->has_source()){
        //   std::cout << "(infty) -> ";
        // }else{
        //   std::cout << "(" << CGAL::to_double(halfEdge->source()->point().x()) << "," << CGAL::to_double(halfEdge->source()->point().y()) << ") -> ";
        // }
        // if(!halfEdge->has_target()){
        //   std::cout << "(infty)" << std::endl;
        // }else{
        //   std::cout << "(" << CGAL::to_double(halfEdge->target()->point().x()) << "," << CGAL::to_double(halfEdge->target()->point().y()) << ")" << std::endl;
        // }
        // std::cout << "HasInternalSource = "<< HasInternalSource(halfEdge) << std::endl;
        // std::cout << "HasInternalTarget = "<< HasInternalTarget(halfEdge) << std::endl;

        // Get the node in the triangulation on the other side of this edge
        auto vt = halfEdge->twin()->face()->dual()->point();

        // Compute the direction of dividing line between this cell and the other cell
        Vector_2 dir(vs.y()-vt.y(), -(vs.x()-vt.x()));
        dir /= std::sqrt(CGAL::to_double(dir.squared_length()));

        // If the edge doesn't have a source or a target...
        if( (!halfEdge->has_source()) && (!halfEdge->has_target()) ){

          // Find the weighted mid point between the triangulation sites
          double ts = CGAL::to_double( CGAL::squared_distance(vt.point(), vs.point()).exact() );
          double h = (0.5/ts)*CGAL::to_double(ts + vs.weight() - vt.weight());
          Point_2 midPt = vs.point() + h*(vt.point()-vs.point());

          Point_2 tgtPt = midPt + infDist*dir;
          Point_2 srcPt = midPt - infDist*dir;
          bbox.ClipSegment(srcPt,tgtPt);

          polyPts.push_back(tgtPt);
          bbox.AddCorners(srcPt, polyPts);
          polyPts.push_back(srcPt);

        }else{

          // If there is a source node, add it!
          bool hasInSrc = HasInternalSource(halfEdge);
          if(hasInSrc){
            polyPts.push_back(halfEdge->source()->point());
          }

          if((!HasInternalTarget(halfEdge)) && hasInSrc){
            /* If this halfedge doesn't have a target node, then the next half
               edge won't have a source node and we need to add in all boundary
               nodes (including corners) need to bound the Laguerre cell.
            */

            // Figure out where the current ray crosses the boundary
            Point_2 curr_src = halfEdge->source()->point();
            Point_2 curr_tgt = curr_src + infDist*dir;
            bbox.ClipSegment(curr_src,curr_tgt);
            polyPts.push_back(curr_tgt);

            // Figure out where the next ray crosses the boundary
            auto nextHalfEdge = halfEdge;
            nextHalfEdge++;
            while(!HasInternalTarget(nextHalfEdge)){
              nextHalfEdge++;
            }
            assert(!HasInternalSource(nextHalfEdge));

            auto next_vt = nextHalfEdge->twin()->face()->dual()->point();
            Vector_2 next_dir(vs.y()-next_vt.y(), -(vs.x()-next_vt.x()));
            next_dir /= std::sqrt(CGAL::to_double(next_dir.squared_length()));

            Point_2 next_tgt = nextHalfEdge->target()->point();
            Point_2 next_src = next_tgt - infDist*next_dir;
            bbox.ClipSegment(next_src, next_tgt);

            // Add corner points (if necessary), until we reach the next ray's source.
            bbox.AddCorners(next_src, polyPts);

            // Add the next edges source and target
            polyPts.push_back(next_src);
            polyPts.push_back(next_tgt);
          }
        }

      } while ( ++halfEdge != halfEdgeStart);

      auto outPoly = std::make_shared<Polygon_2>(polyPts.begin(), polyPts.end());//temp.begin()->outer_boundary());
      assert(outPoly);

      if(outPoly->size()<3){
        std::stringstream msg;
        msg << "Could not construct Laguerre diagram. Found empty Laguerre cell.";
        throw LaguerreDiagram::ConstructionException(msg.str());
      }

      // Remove near duplicate vertices
      auto currVert = outPoly->vertices_begin();
      auto nextVert = currVert;
      nextVert++;
      while(nextVert != outPoly->vertices_end()){

        if( CGAL::squared_distance(*nextVert, *currVert)< compTol){
          nextVert = outPoly->erase(nextVert);
        }else{
          currVert = nextVert;
          nextVert++;
        }
      }

      if(CGAL::squared_distance(*outPoly->vertices_begin(), *currVert)< compTol){
        outPoly->erase(currVert);
      }

      // Sanity check to make sure there are no duplicate vertices or edge crossings
      if(!outPoly->is_simple()){
        std::cout << "Created a polygon that is not simple:\n";
        for(auto& p : polyPts){
          std::cout << "(" << CGAL::to_double(p.x()) << "," << CGAL::to_double(p.y()) << ") -> ";
        }
        std::cout << std::endl;
        assert(outPoly->is_simple());
      }

      return outPoly;
}


void LaguerreDiagram::CreateBoundedCells(Eigen::Matrix2Xd const& pts)
{
  internalEdges.clear();
  internalEdges.resize(numPts);

  // Figure out which faces correspond to the original points
  std::vector<Face_handle> faces(numPts);
  unsigned int faceCount = 0;
  for(auto faceIter = unboundedDiagram.faces_begin(); faceIter!=unboundedDiagram.faces_end(); ++faceIter){
    faces.at(faceIter->dual()->info()) = *faceIter;
    faceCount++;
  }

  if(faceCount!=numPts){
    std::stringstream msg;
    msg << "Could not construct Laguerre diagram. ";
    msg << "There should be " << numPts << " faces, but only " << faceCount << " faces exist in the Laguerre diagram.";
    throw LaguerreDiagram::ConstructionException(msg.str());
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

    bool hasIntersection = bbox.ClipSegment(srcPt, tgtPt);

    // If the edge does not intersect the domain at all, just skip it
    if(hasIntersection){
      // This edge has both a source and target, so it is an internal edge.  Save it!
      unsigned int currInd = face->dual()->info();
      unsigned int otherInd = halfEdge->twin()->face()->dual()->info();
      internalEdges.at(currInd).push_back(std::make_tuple(otherInd, srcPt, tgtPt));
    }

  } while ( ++halfEdge != halfEdgeStart);

}


// bool LaguerreDiagram::ClipToBoundary(Point_2& srcPt, Point_2& tgtPt) const
// {
//
//   double sx = CGAL::to_double(srcPt.x());
//   double sy = CGAL::to_double(srcPt.y());
//   double tx = CGAL::to_double(tgtPt.x());
//   double ty = CGAL::to_double(tgtPt.y());
//
//   // Check if it's vertical
//   if(std::abs(sx-tx)<compTol){
//
//     // Check to see if it crosses the domain at all
//     if((sx>xBndLeft-compTol)&&(sx<xBndRight+compTol)&&(std::max(sy,ty)>yBndBottom-compTol)&&(std::min(sy,ty)<yBndTop+compTol)){
//       srcPt = Point_2(sx, std::min(yBndTop, std::max(yBndBottom, sy)));
//       tgtPt = Point_2(tx, std::min(yBndTop, std::max(yBndBottom, ty)));
//       return true;
//     }else{
//       return false;
//     }
//
//   // Check if it's horizontal
//   }else if(std::abs(sy-ty)<compTol){
//
//     if((sy>yBndBottom-compTol)&&(sy<yBndTop-compTol)&&(std::max(sx,tx)>xBndLeft-compTol)&&(std::min(sx,tx)<xBndRight+compTol)){
//       srcPt = Point_2(std::min(xBndRight, std::max(xBndLeft, sx)), sy);
//       tgtPt = Point_2(std::min(xBndRight, std::max(xBndLeft, tx)), ty);
//       return true;
//     }else{
//       return false;
//     }
//
//   // Otherwise it has a non zero finite slope
//   }else{
//
//
//     // Parameterize the line as srcPt + t*(tgtPt-srcPt) with t in [0,1]
//     Eigen::Matrix<double,4,1> qs(4), ps(4);
//
//     ps(0) = sx-tx;
//     ps(1) = tx-sx;
//     ps(2) = sy-ty;
//     ps(3) = ty-sy;
//
//     qs(0) = (sx-xBndLeft);///(tx-sx);
//     qs(1) = (xBndRight-sx);//(tx-sx);
//     qs(2) = (sy-yBndBottom);///(ty-sy);
//     qs(3) = (yBndTop-sy);///(ty-sy);
//
//
//     double u1 = 0.0;
//     for(unsigned int i=0; i<4; ++i){
//       if(ps(i)<0.0)
//         u1 = std::max(u1, qs(i)/ps(i));
//     }
//
//     double u2 = 1.0;
//     for(unsigned int i=0; i<4; ++i){
//       if(ps(i)>0.0)
//         u2 = std::min(u2, qs(i)/ps(i));
//     }
//
//     if(u1>u2){
//       return false;
//     }else{
//       srcPt = Point_2(sx + u1*ps(1), sy + u1*ps(3));
//       tgtPt = Point_2(sx + u2*ps(1), sy + u2*ps(3));
//       return true;
//     }
//   }
//
// }
