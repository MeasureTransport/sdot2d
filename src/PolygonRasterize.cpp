#include "SDOT/PolygonRasterize.h"

#include <CGAL/Boolean_set_operations_2.h>

using namespace sdot;


PolygonRasterizeIter::PolygonRasterizeIter(std::shared_ptr<RegularGrid> const& gridIn,
                                           std::shared_ptr<Polygon_2>   const& polyIn) : grid(gridIn),
                                                                                         poly(polyIn)
{

  if(poly->size()<=1){
    throw std::runtime_error("Invalid Polygon passed to \"PolygonRasterizeIter\".");
  }

  // Figure out the first y cell in the regular grid that lies within this polygon
  auto bottomVert = polyIn->bottom_vertex();
  int yInd = grid->BottomNode( CGAL::to_double(bottomVert->y()) );// std::ceil( CGAL::to_double(edgesCW.at(0).first->y() - grid->yMin)/grid->dy ) + 1;

  if(yInd < grid->NumCells(1)){
    isValid = true;
  }else{
    isValid = false;
  }

  // Set up the initial edges.  Do not include horizontal edges.
  // Recall that the  polygon orders vertices in counter-clockwise order
  auto temp = NextEdge(bottomVert, CW, true);

  edgesCW.push_back(  temp );

  slopesCW.push_back( CGAL::to_double( (edgesCW.front().second->x()-edgesCW.front().first->x()) / (edgesCW.front().second->y()-edgesCW.front().first->y())));


  edgesCCW.push_back( NextEdge(bottomVert, CCW, true) );
  slopesCCW.push_back( CGAL::to_double( (edgesCCW.front().second->x()-edgesCCW.front().first->x()) / (edgesCCW.front().second->y()-edgesCCW.front().first->y())));

  UpdateEdges(yInd);

  // The y values bounding the current row of grid cells
  double y = grid->yMin + yInd*grid->dy;

  unsigned int xInd_begin;
  std::tie(xInd_begin, maxLeftIndBC, minRightIndBC, xInd_end) = GetXInds(y);

  indices = std::make_pair(xInd_begin,yInd);
  UpdateCellBounds();
  overlapPoly = BoundaryIntersection();

}

std::tuple<unsigned int, unsigned int, unsigned int, unsigned int> PolygonRasterizeIter::GetXInds(double y) const
{
  // Compute the minimum x value for any CW edge obtained within this row of grid cells
  double leftMinX = std::numeric_limits<double>::infinity();
  double leftMaxX = -std::numeric_limits<double>::infinity();
  double rightMinX = std::numeric_limits<double>::infinity();
  double rightMaxX = -std::numeric_limits<double>::infinity();

  double leftMinY = CGAL::to_double( edgesCW.front().first->y() );
  double leftMaxY = CGAL::to_double( edgesCW.back().second->y() );

  double rightMinY = CGAL::to_double( edgesCCW.front().first->y() );
  double rightMaxY = CGAL::to_double( edgesCCW.back().second->y() );

  for(int i=0; i<edgesCW.size(); ++i){

    // The intersection of each Laguerre cell edge with a row of grid cells forms a line segment.
    double x1; // The x value at the start of the line segment.
    double x2; // The x value at the end of the line segment.

    // if the edge starts before this row of grid cells
    if(edgesCW.at(i).first->y() < y-compTol){
      x1 = CGAL::to_double( edgesCW.at(i).first->x() + slopesCW.at(i)*(y-edgesCW.at(i).first->y()) );

    // Otherwise the edge starts inside the grid cell
    }else{
      x1 = CGAL::to_double( edgesCW.at(i).first->x() );
    }

    // If the edge ends after this row of grid cells
    if(edgesCW.at(i).second->y() > y+grid->dy+compTol){
      x2 = CGAL::to_double( edgesCW.at(i).first->x() + slopesCW.at(i)*(y+grid->dy-edgesCW.at(i).first->y()) );

    // Otherwise the edge ends inside the grid cell
    }else{
      x2 = CGAL::to_double( edgesCW.at(i).second->x() );
    }

    leftMinX = std::min(leftMinX, std::min(x1,x2));
    leftMaxX = std::max(leftMaxX, std::max(x1,x2));
  }


  for(int i=0; i<edgesCCW.size(); ++i){

    // The intersection of each Laguerre cell edge with a row of grid cells forms a line segment.
    double x1; // The x value at the start of the line segment.
    double x2; // The x value at the end of the line segment.

    // if the edge starts before this row of grid cells
    if(edgesCCW.at(i).first->y() < y-compTol){
      x1 = CGAL::to_double( edgesCCW.at(i).first->x() + slopesCCW.at(i)*(y-edgesCCW.at(i).first->y()) );

    // Otherwise the edge starts inside the grid cell
    }else{
      x1 = CGAL::to_double( edgesCCW.at(i).first->x() );
    }

    // If the edge ends after this row of grid cells
    if(edgesCCW.at(i).second->y() > y+grid->dy+compTol){
      x2 = CGAL::to_double( edgesCCW.at(i).first->x() + slopesCCW.at(i)*(y+grid->dy-edgesCCW.at(i).first->y()) );

    // Otherwise the edge ends inside the grid cell
    }else{
      x2 = CGAL::to_double( edgesCCW.at(i).second->x() );
    }

    rightMaxX = std::max(rightMaxX, std::max(x1,x2));
    rightMinX = std::min(rightMinX, std::min(x1,x2));
  }

//  std::cout << "y = " << y << ", x bounds = " << leftMinX << ", " << rightMaxX << std::endl;

  // round the min and max to nearest grid cell boundaries
  unsigned int leftInd1 = grid->LeftNode(leftMinX);
  unsigned int leftInd2 = std::max(int(grid->RightNode(leftMaxX))-1,0);
  unsigned int rightInd1 = grid->LeftNode(rightMinX);
  unsigned int rightInd2 = std::max(int(grid->RightNode(rightMaxX))-1,0);
  leftInd2 = std::max(leftInd1, leftInd2);
  rightInd1 = std::min(rightInd1, rightInd2);

  // If either the CW edges or the CCW edges start or end inside this row
  if((leftMinY>y+compTol)||(rightMinY>y+compTol)||(leftMaxY<y+grid->dy-compTol)||(rightMaxY<y+grid->dy-compTol)){
    return std::make_tuple(leftInd1, rightInd2, leftInd2, rightInd2);
  }else{
    return std::make_tuple(leftInd1,leftInd2,rightInd1,rightInd2);
  }
}



void PolygonRasterizeIter::UpdateEdges(unsigned int yInd)
{
  double ymin = yInd*grid->dy + grid->yMin;
  double ymax = ymin + grid->dy;

  // Remove any edges that will no longer intersect grid cells
  while(!edgesCW.empty()){

    if(edgesCW.front().second->y()>ymin-compTol){
      break;
    }else{
      edgesCW.pop_front();
      slopesCW.pop_front();
    }
  }

  // Now add any new edges that might possibly intersect grid cells
  while(edgesCW.back().second->y() < ymax+compTol){
    auto nextEdge = NextEdge(edgesCW.back().second, CW, true);

    if(nextEdge.first->y() > nextEdge.second->y()){
      break;
    }else{
      edgesCW.push_back(nextEdge);
      slopesCW.push_back( CGAL::to_double( (nextEdge.second->x()-nextEdge.first->x()) / (nextEdge.second->y()-nextEdge.first->y())));
    }
  }


  // Remove any edges that will no longer intersect grid cells
  while(!edgesCCW.empty()){

    if(edgesCCW.front().second->y()>ymin-compTol){
      break;
    }else{
      edgesCCW.pop_front();
      slopesCCW.pop_front();
    }
  }

  // Now add any new edges that might possibly intersect grid cells
  while(edgesCCW.back().second->y() < ymax+compTol){
    auto nextEdge = NextEdge(edgesCCW.back().second, CCW, true);

    if(nextEdge.first->y() > nextEdge.second->y()){
      break;
    }else{
      edgesCCW.push_back(nextEdge);
      slopesCCW.push_back( CGAL::to_double( (nextEdge.second->x()-nextEdge.first->x()) / (nextEdge.second->y()-nextEdge.first->y())));
    }
  }
}

void PolygonRasterizeIter::UpdateCellBounds()
{
  cellBox = std::make_shared<BoundingBox>(indices.first * grid->dx + grid->xMin,
                                          (indices.first+1) * grid->dx + grid->xMin,
                                          indices.second * grid->dy + grid->yMin,
                                          (indices.second+1) * grid->dy + grid->yMin);
}

// std::pair<std::shared_ptr<PolygonRasterizeIter::Point_2>,
//           std::shared_ptr<PolygonRasterizeIter::Point_2>> PolygonRasterizeIter::EdgeCrossings(PolygonEdge const& edge,
//                                                                                               bool includeEdges)
// {
//   // locations of the source and target vertices of the edge
//   double tx = CGAL::to_double( edge.second->x() );
//   double ty = CGAL::to_double( edge.second->y() );
//   double sx = CGAL::to_double( edge.first->x() );
//   double sy = CGAL::to_double( edge.first->y() );
//
//   std::shared_ptr<Point_2> enterPt, exitPt;
//
//   double enterX, enterY;
//   double exitX, exitY;
//
//   // The edge is entirely to the left or entirely to the right of the grid cell
//   if(includeEdges){
//     if(((tx>xright+compTol)&&(sx>xright+compTol))||((tx<xleft-compTol)&&(sx<xleft-compTol))){
//       enterPt=nullptr;
//       exitPt=nullptr;
//       return std::make_pair(enterPt, exitPt);
//     }
//   }else{
//     if(((tx>xright-compTol)&&(sx>xright-compTol))||((tx<xleft+compTol)&&(sx<xleft+compTol))){
//       enterPt=nullptr;
//       exitPt=nullptr;
//       return std::make_pair(enterPt, exitPt);
//     }
//   }
//
//   // The edge is entirely above or below the grid cell
//   if(includeEdges){
//     if(((ty>ytop+compTol)&&(sy>ytop+compTol))||((ty<ybot-compTol)&&(sy<ybot-compTol))){
//       enterPt=nullptr;
//       exitPt=nullptr;
//       return std::make_pair(enterPt, exitPt);
//     }
//   }else{
//     if(((ty>ytop-compTol)&&(sy>ytop-compTol))||((ty<ybot+compTol)&&(sy<ybot+compTol))){
//       enterPt=nullptr;
//       exitPt=nullptr;
//       return std::make_pair(enterPt, exitPt);
//     }
//   }
//
//   // If it's vertical...
//   if(std::abs(tx-sx)<horizTol){
//
//     // if it overlaps the boundary exactly
//     if((tx<xleft+compTol)||(tx > xright-compTol)){
//       enterPt=nullptr;
//       exitPt=nullptr;
//       return std::make_pair(enterPt, exitPt);
//     }
//
//     enterX = sx;
//     exitX = sx;
//     if(sy>ty){
//       enterY = ytop;
//       exitY = ybot;
//     }else{
//       enterY = ybot;
//       exitY = ytop;
//     }
//
//   // If it's horizontal
//   }else if(std::abs(ty-sy)<horizTol){
//
//     // if it doesn't intersect the grid cell
//     if((ty<ybot+compTol)||(ty>ytop-compTol)){
//       enterPt=nullptr;
//       exitPt=nullptr;
//       return std::make_pair(enterPt, exitPt);
//     }
//
//     enterY = sy;
//     exitY = sy;
//     if(tx<sx){ // leave left side
//       enterX = xright;
//       exitX = xleft;
//
//     }else{ // leave right side
//       enterX = xleft;
//       exitX = xright;
//     }
//
//   // If it's not vertical or horizontal...
//   }else{
//
//     double slope = (ty-sy)/(tx-sx);
//
//     // If slope is positive, an intersection will only occur if the bottom right point is below the line segment.
//     if(slope>0){
//
//       if((ybot > sy + slope*(xright-sx)-compTol)||(ytop < sy + slope*(xleft-sx)+compTol)){
//         enterPt=nullptr;
//         exitPt=nullptr;
//         return std::make_pair(enterPt, exitPt);
//       }
//
//     // Similarly, if the slope is negative, an intersection will only occur if the bottom left point is below the line segment.
//     }else if(slope<0){
//       if((ybot > sy + slope*(xleft-sx)-compTol)||(ytop < sy + slope*(xright-sx)+compTol)){
//         enterPt=nullptr;
//         exitPt=nullptr;
//         return std::make_pair(enterPt, exitPt);
//       }
//
//     }else{
//       std::cerr << "I shouldn't be here...  The slop should always be nonzero at this stage." << std::endl;
//       assert(slope!=0);
//     }
//
//     // could enter from left and exit from right
//     if(tx>sx){
//       enterY = std::min(ytop, std::max(ybot, sy + slope*(xleft-sx)));
//       exitY = std::min(ytop, std::max(ybot, sy + slope*(xright-sx)));
//
//     // could enter from right and exit from left
//     }else{
//       enterY = std::min(ytop, std::max(ybot, sy + slope*(xright-sx)));
//       exitY = std::min(ytop, std::max(ybot, sy+slope*(xleft-sx)));
//     }
//
//     // Could enter from top and leave from bottom
//     if(ty<sy){
//       enterX = std::min(xright, std::max(xleft, sx + (ytop-sy)/slope));
//       exitX = std::min(xright, std::max(xleft, sx + (ybot-sy)/slope));
//
//     // Could enter from bottom and leave from top
//     }else{
//       enterX = std::min(xright, std::max(xleft, sx + (ybot-sy)/slope));
//       exitX = std::min(xright, std::max(xleft, sx + (ytop-sy)/slope));
//     }
//
//   }
//
//   bool srcIsInside, tgtIsInside;
//
//   if(includeEdges){
//     srcIsInside = (sx>xleft-compTol)&&(sx<xright+compTol)&&(sy>ybot-compTol)&&(sy<ytop+compTol);
//     tgtIsInside = (tx>xleft-compTol)&&(tx<xright+compTol)&&(ty>ybot-compTol)&&(ty<ytop+compTol);
//   }else{
//     srcIsInside = (sx>xleft+compTol)&&(sx<xright-compTol)&&(sy>ybot+compTol)&&(sy<ytop-compTol);
//     tgtIsInside = (tx>xleft+compTol)&&(tx<xright-compTol)&&(ty>ybot+compTol)&&(ty<ytop-compTol);
//   }
//
//   if(tgtIsInside){
//     exitPt = nullptr;
//   }else{
//     exitPt = std::make_shared<Point_2>(exitX,exitY);
//   }
//
//   if(srcIsInside){
//     enterPt = nullptr;
//   }else{
//     enterPt = std::make_shared<Point_2>(enterX,enterY);
//   }
//
//   return std::make_pair(enterPt, exitPt);
// }


// void PolygonRasterizeIter::AddCorners(Point_2 const& nextPt, std::vector<Point_2> &polyPts) const
// {
//   // const double xmin = indices.first * grid->dx + grid->xMin;
//   // const double xmax = xmin + grid->dx;
//   // const double ymin = indices.second * grid->dy + grid->yMin;
//   // const double ymax = ymin + grid->dy;
//
//   double backX = CGAL::to_double( polyPts.back().x() );
//   double backY = CGAL::to_double( polyPts.back().y() );
//
//   double enterX = CGAL::to_double( nextPt.x() );
//   double enterY = CGAL::to_double( nextPt.y() );
//
//   // Check to see if we need to add add any corners.
//   // If the xs and the ys are different OR if the points are far enough apart, then we need to add corner
//   bool sameXdiffY = (std::abs(backX-enterX)<compTol)&&(std::abs(backY-enterY)>compTol);
//   bool sameYdiffX = (std::abs(backY-enterY)<compTol)&&(std::abs(backX-enterX)>compTol);
//   bool sameEdge = sameXdiffY || sameYdiffX;
//
//   if(sameEdge)
//     return;
//
//   for(unsigned int i=0; i<4; ++i){
//   //while(!sameEdge){
//     // Add corner point
//
//     // Currently on the right
//     if((std::abs(backX-xright)<compTol)&&(backY<ytop-compTol)){
//       polyPts.push_back( Point_2(xright,ytop) );
//
//     // Currently on the top
//     }else if((std::abs(backY-ytop)<compTol)&&(backX>xleft+compTol)){
//       polyPts.push_back( Point_2(xleft, ytop) );
//
//     // Currently on the left
//     }else if((std::abs(backX-xleft)<compTol)&&(backY>ybot+compTol)){
//       polyPts.push_back( Point_2(xleft, ybot) );
//
//     // On the bottom
//     }else if((std::abs(backY-ybot)<compTol)&&(backX<xright-compTol)){
//       polyPts.push_back( Point_2(xright, ybot) );
//
//     }else{
//       std::cerr << "I shouldn't be here... " << std::endl;
//       for(auto& pt : polyPts){
//         std::cout << "  " << pt << std::endl;
//       }
//       assert(false);
//     }
//
//     backX = CGAL::to_double( polyPts.back().x() );
//     backY = CGAL::to_double( polyPts.back().y() );
//
//     sameXdiffY = (std::abs(backX-enterX)<compTol)&&(std::abs(backY-enterY)>compTol);
//     sameYdiffX = (std::abs(backY-enterY)<compTol)&&(std::abs(backX-enterX)>compTol);
//     sameEdge = sameXdiffY || sameYdiffX;
//
//     if(sameEdge){
//       break;
//     }
//   }
// }

std::shared_ptr<PolygonRasterizeIter::Polygon_2> PolygonRasterizeIter::BoundaryIntersection()
{

  return cellBox->ClipPolygon(poly);
  //
  // //////////////////////////////////////////////////////////////////////////////
  // // USING CGAL THIS WAY IS REALLY SLOW, BUT CAN BE A GOOD SANITY CHECK
  // //
  // Polygon_2 gridCellPoly;
  // gridCellPoly.push_back(Point_2(indices.first*grid->dx+grid->xMin, indices.second*grid->dy+grid->yMin));
  // gridCellPoly.push_back(Point_2((indices.first+1)*grid->dx+grid->xMin, indices.second*grid->dy+grid->yMin));
  // gridCellPoly.push_back(Point_2((indices.first+1)*grid->dx+grid->xMin, (indices.second+1)*grid->dy+grid->yMin));
  // gridCellPoly.push_back(Point_2(indices.first*grid->dx+grid->xMin, (indices.second+1)*grid->dy+grid->yMin));
  //
  // // The edgesCW and edgesCCW containers contain edges that could overlap with
  // // each grid cell
  // if(!poly->is_simple()){
  //   std::cout << "Found a polygon that is not simple: " << *poly << std::endl;
  //   assert(poly->is_simple());
  // }
  //
  // std::list<Polygon_with_holes_2> temp;
  // CGAL::intersection(*poly, gridCellPoly, std::back_inserter(temp));
  //
  // // There should only be one polygon intersection because everything is simple and convex
  // assert(std::distance(temp.begin(),temp.end())==1);
  //
  // auto outPoly = std::make_shared<Polygon_2>(temp.begin()->outer_boundary());
  //
  // return outPoly;
  // //
  // //////////////////////////////////////////////////////////////////////////////
  //
  //
  // std::vector<Point_2> polyPts;
  // std::shared_ptr<Point_2> enterPt, exitPt;
  //
  // std::shared_ptr<Point_2> nextCorner;
  //
  // // start with the ccw edges
  // for(auto& edge : edgesCCW){
  //   std::tie(enterPt, exitPt) = EdgeCrossings(edge, false);
  //
  //   // std::cout << "CCW Edge: " << *edge.first << " -> " << *edge.second << std::endl;
  //
  //   // Add the entrance and exit points
  //   if(enterPt){
  //     // std::cout << "  enter pt: " << *enterPt << std::endl;
  //     if(polyPts.size()>0){
  //
  //       // Check to make sure this is a new point
  //       if(polyPts.back()!=(*enterPt)){
  //         AddCorners(*enterPt, polyPts);
  //         polyPts.push_back(*enterPt);
  //       }
  //
  //     }else{
  //       polyPts.push_back(*enterPt);
  //     }
  //
  //   // If no entrance point, then add source node if it is inside the grid cell
  //   }else if((edge.first->x()<xright-compTol)&&(edge.first->x()>xleft+compTol)&&(edge.first->y()>ybot+compTol)&&(edge.first->y()<ytop+compTol)){
  //     if(polyPts.size()>0){
  //       if(polyPts.back()!=(*edge.first)){
  //         polyPts.push_back(*edge.first);
  //       }
  //     }else{
  //       polyPts.push_back(*edge.first);
  //     }
  //   }
  //
  //   if(exitPt){
  //     // std::cout << "  exit pt: " << *exitPt << std::endl;
  //     if(polyPts.size()>0){
  //       if(polyPts.back()!=(*exitPt)){
  //         //AddCorners(*exitPt, polyPts);
  //         polyPts.push_back(*exitPt);
  //       }
  //     }else{
  //       polyPts.push_back(*exitPt);
  //     }
  //
  //   // If no exit point, check to see if target is inside
  //   }else if((edge.second->x()<xright-compTol)&&(edge.second->x()>xleft+compTol)&&(edge.second->y()>ybot+compTol)&&(edge.second->y()<ytop+compTol)){
  //     if(polyPts.size()>0){
  //       if(polyPts.back()!=(*edge.second)){
  //         polyPts.push_back(*edge.second);
  //       }
  //     }else{
  //       polyPts.push_back(*edge.second);
  //     }
  //   }
  // }
  //
  // // std::cout << "PolyPts 0:" << std::endl;
  // // for(auto& p : polyPts)
  // //   std::cout << "  " << p << std::endl;
  //
  // /* If the last CW edge and last CCW edge end at the same y-value, then
  //    there is a horizontal edge that might cross this grid cell */
  //  double ccwy = CGAL::to_double( edgesCCW.back().second->y() );
  //  double cwy = CGAL::to_double( edgesCW.back().second->y() );
  //  double ccwx = CGAL::to_double( edgesCCW.back().second->x() );
  //  double cwx = CGAL::to_double( edgesCW.back().second->x() );
  //  if((std::abs(ccwy-cwy)<compTol)&&(ccwx>xleft+compTol)&&(cwx<xright-compTol)){
  //    if((ccwy>ybot+compTol)&&(ccwy<ytop-compTol)){
  //      Point_2 newPt(std::min(xright,ccwx), ccwy);
  //      if(polyPts.size()>0){
  //        if(newPt != polyPts.back())
  //          polyPts.push_back(newPt);
  //      }else{
  //        polyPts.push_back( newPt );
  //      }
  //
  //      polyPts.push_back( Point_2(std::max(xleft,cwx), cwy) );
  //    }
  //  }
  //
  //  // std::cout << "PolyPts 1:" << std::endl;
  //  // for(auto& p : polyPts)
  //  //   std::cout << "  " << p << std::endl;
  //
  // // Add the CW edges, but go in counter-clockwise order
  // for(auto edgeIter = edgesCW.rbegin(); edgeIter!=edgesCW.rend(); edgeIter++){
  //   // std::cout << "CW Edge: " << *edgeIter->second << " -> " << *edgeIter->first << std::endl;
  //   std::tie(enterPt, exitPt) = EdgeCrossings(std::make_pair(edgeIter->second, edgeIter->first), false);
  //
  //   if(enterPt){
  //     // std::cout << "  enter pt: " << *enterPt << std::endl;
  //
  //     if(polyPts.size()>0){
  //       if(polyPts.back()!=(*enterPt)){
  //         AddCorners(*enterPt, polyPts);
  //         polyPts.push_back(*enterPt);
  //       }
  //     }else{
  //       polyPts.push_back(*enterPt);
  //     }
  //
  //   // If no enter point, check edge
  //   }else if((edgeIter->second->x()<xright-compTol)&&(edgeIter->second->x()>xleft+compTol)&&(edgeIter->second->y()>ybot+compTol)&&(edgeIter->second->y()<ytop-compTol)){
  //     if(polyPts.size()>0){
  //       if(polyPts.back()!=(*edgeIter->second)){
  //         polyPts.push_back(*edgeIter->second);
  //       }
  //     }else{
  //       polyPts.push_back(*edgeIter->second);
  //     }
  //   }
  //
  //   if(exitPt){
  //     // std::cout << "  exit pt: " << *exitPt << std::endl;
  //
  //     if(polyPts.size()>0){
  //       if(polyPts.back()!=(*exitPt)){
  //         //AddCorners(*exitPt, polyPts);
  //         polyPts.push_back(*exitPt);
  //       }
  //     }else{
  //       polyPts.push_back(*exitPt);
  //     }
  //   }else if((edgeIter->first->x()<xright-compTol)&&(edgeIter->first->x()>xleft+compTol)&&(edgeIter->first->y()>ybot+compTol)&&(edgeIter->first->y()<ytop+compTol)){
  //     if(polyPts.size()>0){
  //       if(polyPts.back()!=(*edgeIter->first)){
  //         polyPts.push_back(*edgeIter->first);
  //       }
  //     }else{
  //       polyPts.push_back(*edgeIter->first);
  //     }
  //   }
  // }
  //
  // // std::cout << "PolyPts 2:" << std::endl;
  // // for(auto& p : polyPts)
  // //   std::cout << "  " << p << std::endl;
  //
  //
  // /* If the first CW edge and first CCW edge start at the same y-value, then
  //    there is a horizontal edge we might need to account for. Note that the edge
  //    must cross one or both of the grid cell sides */
  // ccwy = CGAL::to_double( edgesCCW.front().first->y() );
  // cwy = CGAL::to_double( edgesCW.front().first->y() );
  // ccwx = CGAL::to_double( edgesCCW.front().first->x() );
  // cwx = CGAL::to_double( edgesCW.front().first->x() );
  //
  // if((std::abs(ccwy-cwy)<compTol)&&(ccwx>xleft+compTol)&&(cwx<xright-compTol)){
  //   if(ccwy>ybot+compTol){
  //     if(polyPts.size()>0){
  //       Point_2 newPt(std::min(xright,ccwx), ccwy);
  //       if(polyPts.front() != newPt)
  //         polyPts.insert(polyPts.begin(), newPt);
  //       polyPts.insert(polyPts.begin(), Point_2(std::max(xleft,cwx), cwy));
  //     }else{
  //       polyPts.push_back(Point_2(std::max(xleft,cwx), cwy));
  //       polyPts.push_back(Point_2(std::min(xright,ccwx), ccwy));
  //     }
  //   }
  // }
  //
  // // std::cout << "PolyPts 3:" << std::endl;
  // // for(auto& p : polyPts)
  // //   std::cout << "  " << p << std::endl;
  //
  // // Add corner points to complete the loop
  // if(polyPts.size()>0)
  //   AddCorners(polyPts.front(), polyPts);
  //
  // // std::cout << "PolyPts 4:" << std::endl;
  // // for(auto& p : polyPts)
  // //   std::cout << "  " << p << std::endl;
  //
  // // std::cout << std::endl << std::endl;
  // if(polyPts.size()==0){
  //   return nullptr;
  // }else{
  //   return std::make_shared<Polygon_2>(polyPts.begin(), polyPts.end());
  // }

}


PolygonRasterizeIter& PolygonRasterizeIter::Increment()
{
  indices.first++;
  UpdateCellBounds();

  if(indices.first>xInd_end){

    // Increment the yInd and find the new x values
    unsigned int yInd = indices.second+1;
    UpdateEdges(yInd);

    double y = grid->yMin + grid->dy*yInd;

    if((edgesCW.size()==0)||(edgesCCW.size()==0)){
      isValid = false;
      return *this;
    }
    // If there's only one edge left and we're at the top of current edge
    if(edgesCW.size()==1){

      // And we're at the top of the current edge
      if(y>CGAL::to_double(edgesCW.front().second->y())-compTol){
        isValid = false;
        return *this;
      }
    }

    unsigned int xInd_begin;
    std::tie(xInd_begin, maxLeftIndBC, minRightIndBC, xInd_end) = GetXInds(y);

    if(xInd_end<xInd_begin){
      isValid=false;
      return *this;
    }

    indices = std::make_pair(xInd_begin, yInd);

    UpdateCellBounds();
  }

  if((indices.first<=maxLeftIndBC)||(indices.first>=minRightIndBC)){
    overlapPoly = BoundaryIntersection();
  }else{
    overlapPoly = nullptr;
  }

  return *this;
};


bool PolygonRasterizeIter::IsValid() const
{
  return isValid;
}

std::pair<int,int> const& PolygonRasterizeIter::Indices() const
{
  return indices;
}


std::pair<PolygonRasterizeIter::Polygon_2::Vertex_const_iterator,
          PolygonRasterizeIter::Polygon_2::Vertex_const_iterator> PolygonRasterizeIter::NextEdge(Polygon_2::Vertex_const_iterator const& currIter,
                                                                                                 Direction                               dir,
                                                                                                 bool                                    horizCheck)
{
 const double horizTol = 1e-12;

 auto srcVert = currIter;
 auto tgtVert = NextVertex(srcVert, dir);

 if(horizCheck) {
   // Iterate until we find a non-horizontal edge
   double ydiff = CGAL::to_double( tgtVert->y() - srcVert->y() );

   while( std::abs(ydiff) < horizTol ){
     srcVert = tgtVert;
     tgtVert = NextVertex(tgtVert,dir);
     ydiff = CGAL::to_double( tgtVert->y() - srcVert->y() );
   }
 }

 return std::make_pair(srcVert, tgtVert);
}


PolygonRasterizeIter::Polygon_2::Vertex_const_iterator PolygonRasterizeIter::NextVertex(PolygonRasterizeIter::Polygon_2::Vertex_const_iterator const& currIter,
                                                                                        Direction                               dir)
{
 if(dir==CCW){
   Polygon_2::Vertex_const_iterator result = currIter;
   result++;

   if(result==poly->vertices_end())
     return poly->vertices_begin();

   return result;

 }else{

   Polygon_2::Vertex_const_iterator result = currIter;

   if(currIter==poly->vertices_begin())
     result = poly->vertices_end();

   if(result != poly->vertices_begin())
     result--;

   return result;
 }
}
//
// std::pair<unsigned int, unsigned int> PolygonRasterizeIter::GetXInds(double xcw,  double xcw_next,
//                                                                      double xccw, double xccw_next) const
// {
//  unsigned int tempInd1 = grid->LeftNode(xcw);
//  unsigned int tempInd2 = grid->LeftNode(xcw_next);
//  unsigned int xInd_begin = std::min( tempInd1, tempInd2);
//
//  tempInd1 = grid->RightNode(xccw);
//  tempInd2 = grid->RightNode(xccw_next);
//  unsigned int xInd_end = std::max( tempInd1, tempInd2)-1;
//
//  return std::make_pair(xInd_begin, xInd_end);
// }
