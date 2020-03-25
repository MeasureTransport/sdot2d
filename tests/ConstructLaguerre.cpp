
#include <Eigen/Core>

#include <vector>
#include <memory>

#include <CGAL/Boolean_set_operations_2.h>

// standard includes
#include <iostream>

#include "SDOT/PolygonRasterize.h"
#include "SDOT/RegularGrid.h"
#include "SDOT/LaguerreDiagram.h"

using namespace sdot;



/** Works somewhat like an iterator over cells that intersect the boundary of a
    polygon.
*/
// class BoundaryPolygonGrid{
// public:
//
//   BoundaryPolygonGrid(std::shared_ptr<RegularGrid> const& gridIn,
//                       std::shared_ptr<Polygon_2>   const& polyIn);
//
//
//   BoundaryPolygonGrid& Increment();
//
//   bool IsValid() const{return isValid;};
//
//   std::pair<int,int> const& Indices() const{
//     return indices;
//   }
//
//   std::shared_ptr<Polygon_2> const& Intersection(){
//     // if(needsOverlapUpdate){
//     //   UpdateOverlap();
//     //   needsOverlapUpdate = false;
//     // }
//     return overlapPoly;
//   }
//
//
//
// private:
//   const double compTol = 1e-14;
//
//   //bool needsOverlapUpdate;
//
//   bool isValid;
//
//   Polygon_2::Edge_const_circulator startEdge;
//   Polygon_2::Edge_const_circulator currEdge;
//   Polygon_2::Edge_const_circulator nextEdge;
//
//   // x and y indices of the current grid cell in the regular grid
//   std::pair<int,int> indices;
//
//   // x and y indices where we started
//   std::pair<int,int> startIndices;
//
//   // Points defining the boundary within the current cell
//   //std::vector<Point_2> crossPts;
//
//   // Polygon describing the intersection of the laguerre cell and the grid cell
//   std::shared_ptr<Polygon_2> overlapPoly;
//
//   // The regular grid
//   std::shared_ptr<RegularGrid> grid;
//
//   // The Laguerre cell
//   std::shared_ptr<Polygon_2> lagCell;
//
//   Polygon_2 gridCell;
//
//   void UpdateOverlap();
//   void IncrementIndices();
//
//   /** Computes the points where a particular Laguerre cell edge crosses the
//       boundary of the current grid cell.  The first result is where the Laguerre
//       edge enters the cell.  The second result is where the Laguerre cells leaves
//       the cell.   If only one crossing exists (i.e., the  edge starts  or ends
//       in the grid cell), then a nullptr is used to indicate the missing crossing.
//       The `includeEdges` input gives control over what points are considered "inside"
//       the grid cell.  If includeEdges==true, then the grid cell edges are inclusive.  Otherwise
//       the grid cell edges are not inclusive.  For example, if includeEdges==true
//       and an Laguerre edge target point lies on the grid cell edge, then no
//       exit edge crossing will be returned because the Laguerre edge doesn't actually
//       leave the grid cell.  If includeEdges==false, then the Laguerre edge target
//       point would be considered outside the grid cell and it would be returned
//       as an exit point.
//   */
//   std::pair<std::shared_ptr<Point_2>, std::shared_ptr<Point_2>> EdgeCrossings(Polygon_2::Edge_const_circulator const& edge,
//                                                                               bool                                    includeEdges);
//
//   /** Adds corners of the current gridcell to a vector of points to complete the
//       intersection polygon.  This function assumes the first point in the polyPts
//       vector is the point where an edge enters the grid cell and the last  point
//       in the vector is the point where the edge leaves the grid cell.
//   */
//   void AddCellPts(std::vector<Point_2> &polyPts);
//
//
//   /** Checks if two points are equal up to the `compTol` tolerance.
//   */
//   bool Equal(Point_2 const& pt1, Point_2 const& pt2) const;
//
//   /** Checks if two points are equal up to the `compTol` tolerance.
//   */
//   bool Equal(double x1, double y1, double x2, double y2) const;
//
//   /** Checks whether a point is inside the current grid cell or not.  Boundaries
//       are inclusive in this function.
//   */
//   bool InGridCell(Point_2 const& pt) const;
// };
//
// BoundaryPolygonGrid::BoundaryPolygonGrid(std::shared_ptr<RegularGrid> const& gridIn,
//                                          std::shared_ptr<Polygon_2>   const& polyIn) : grid(gridIn),
//                                                                                        lagCell(polyIn)
// {
//   // Start by getting any edge on the boundary
//   startEdge = lagCell->edges_circulator();
//   currEdge = startEdge;
//   isValid = true;
//
//   // Figure out which grid cell contains the source vertex of the edge
//   double xsrc = CGAL::to_double( currEdge->source().x() );
//   double ysrc = CGAL::to_double( currEdge->source().y() );
//
//   unsigned int xInd = std::floor( (xsrc-grid->xMin)/grid->dx );
//   unsigned int yInd = std::floor( (ysrc-grid->yMin)/grid->dy );
//
//   indices = std::make_pair(xInd,yInd);
//   startIndices = indices;
//
//   // Find the edges entering and leaving this cell
//   std::vector<Point_2> polyPts;
//
//   std::shared_ptr<Point_2> enterPt, exitPt, dummy;
//   auto tempIter = currEdge;
//   std::tie(enterPt, exitPt) = EdgeCrossings(tempIter, true);
//
//   if(!enterPt){
//
//     // Find the edge that enters the current grid cell
//     tempIter--;
//
//     std::tie(enterPt, dummy) = EdgeCrossings(tempIter,true);
//     while((!enterPt)){
//       tempIter--;
//       std::tie(enterPt, dummy) = EdgeCrossings(tempIter,true);
//     }
//
//     // Add the entrance point and any internal points along the way
//     polyPts.push_back(*enterPt);
//     for(; tempIter!=startEdge; ++tempIter){
//       polyPts.push_back(tempIter->target());
//     }
//
//   }else{
//
//     // Check if the source and entrance points are the same
//     if( (std::abs(CGAL::to_double(currEdge->source().x()-enterPt->x()))>compTol) &&(std::abs(CGAL::to_double(currEdge->source().y()-enterPt->y()))>compTol) ){
//       polyPts.push_back(*enterPt);
//     }
//     polyPts.push_back(currEdge->source());
//   }
//
//   //  Find the exit point
//   tempIter = currEdge;
//   while(!exitPt){
//     polyPts.push_back(tempIter->target());
//     tempIter++;
//     std::tie(dummy, exitPt) = EdgeCrossings(tempIter,true);
//   }
//
//   nextEdge = tempIter;
//   polyPts.push_back(*exitPt);
//
//   AddCellPts(polyPts);
//   overlapPoly = std::make_shared<Polygon_2>(polyPts.begin(), polyPts.end());
//
//   //needsOverlapUpdate = true;
// }
//
// std::pair<std::shared_ptr<Point_2>, std::shared_ptr<Point_2>> BoundaryPolygonGrid::EdgeCrossings(Polygon_2::Edge_const_circulator const& edge,
//                                                                                                  bool includeEdges)
// {
//   // locations of grid cell boundaries
//   double xleft = indices.first * grid->dx + grid->xMin;
//   double xright = (indices.first+1) * grid->dx + grid->xMin;
//   double ybot = indices.second * grid->dy + grid->yMin;
//   double ytop = (indices.second + 1) * grid->dy + grid->yMin;
//
//   // locations of the source and target vertices of the edge
//   double tx = CGAL::to_double( edge->target().x() );
//   double ty = CGAL::to_double( edge->target().y() );
//   double sx = CGAL::to_double( edge->source().x() );
//   double sy = CGAL::to_double( edge->source().y() );
//
//
//   std::shared_ptr<Point_2> enterPt, exitPt;
//
//   double enterX, enterY;
//   double exitX, exitY;
//
//   // If it's vertical...
//   if(std::abs(tx-sx)<1e-12){
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
//   }else if(std::abs(ty-sy)<1e-12){
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
//   // If it's not vertical or horizontal...
//   }else{
//     double slope = (ty-sy)/(tx-sx);
//
//     // could enter from left and exit from right
//     if(tx>sx){
//       enterY = std::min(ytop, std::max(ybot, sy + slope*(xleft-sx)));
//       exitY = std::min(ytop, std::max(ybot, sy+slope*(xright-sx)));
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
//
// void BoundaryPolygonGrid::UpdateOverlap()
// {
//   unsigned int xInd = indices.first;
//   unsigned int yInd = indices.second;
//
//   // Compute the polygon intersection
//   gridCell.clear();
//   gridCell.push_back( Point_2(xInd*grid->dx + grid->xMin, yInd*grid->dy + grid->yMin) );
//   gridCell.push_back( Point_2((xInd+1)*grid->dx + grid->xMin, yInd*grid->dy + grid->yMin) );
//   gridCell.push_back( Point_2((xInd+1)*grid->dx + grid->xMin, (yInd+1)*grid->dy + grid->yMin) );
//   gridCell.push_back( Point_2(xInd*grid->dx + grid->xMin, (yInd+1)*grid->dy + grid->yMin) );
//
//   std::list<Polygon_with_holes_2> temp;
//   CGAL::intersection(*lagCell, gridCell, std::back_inserter(temp));
//   assert(temp.size()>0);
//
//   overlapPoly = std::make_shared<Polygon_2>(temp.begin()->outer_boundary());
// }
//
// void BoundaryPolygonGrid::IncrementIndices()
// {
//   double xleft = indices.first * grid->dx + grid->xMin;
//   double xright = (indices.first+1) * grid->dx + grid->xMin;
//   double ybot = indices.second * grid->dy + grid->yMin;
//   double ytop = (indices.second + 1) * grid->dy + grid->yMin;
//
//   // Now, we know that the source of the edge is inside the grid cell and the target is not
//   double tx = CGAL::to_double( currEdge->target().x() );
//   double ty = CGAL::to_double( currEdge->target().y() );
//   double sx = CGAL::to_double( currEdge->source().x() );
//   double sy = CGAL::to_double( currEdge->source().y() );
//
//   // Figure out which side of the grid cell the Laguerre cell edge exits
//   int xInd_next, yInd_next;
//
//   // If it's a vertical edge...
//   if(std::abs(tx - sx) < 1e-12) {
//
//     xInd_next = indices.first;
//
//     // Crosses out the top
//     if(ty > sy){
//       yInd_next = indices.second + 1;
//
//     // Crosses out the bottom
//     }else{
//       yInd_next = indices.second - 1;
//     }
//
//   }else{
//     double slope = (ty-sy)/(tx-sx);
//
//     xInd_next = indices.first;
//     yInd_next = indices.second;
//
//     // could exit out the right
//     if(tx > sx){
//       double rightY = sy + slope*(xright - sx);
//       if((rightY<=ytop)&&(rightY>=ybot)){
//         xInd_next++;
//       }
//
//     // could exit out the left
//     }else{
//       double leftY = sy + slope*(xleft - sx);
//       if((leftY<ytop+compTol)&&(leftY>ybot-compTol)){
//         xInd_next--;
//       }
//     }
//
//     // could exit out the top
//     if(ty > sy){
//       double topX = sx + (ytop-sy)/slope;
//       if((topX<xright+compTol)&&(topX>xleft-compTol)){
//         yInd_next++;
//       }
//
//     // could exit out the bottom
//     }else{
//       double botX = sx + (ybot-sy)/slope;
//       if((botX<xright+compTol)&&(botX>xleft-compTol)){
//         yInd_next--;
//       }
//     }
//   }
//
//   // Make sure we've moved
//   assert((xInd_next != indices.first) || (yInd_next != indices.second));
//
//   indices = std::make_pair(xInd_next,yInd_next);
// }
//
// bool BoundaryPolygonGrid::Equal(Point_2 const& pt1, Point_2 const& pt2) const
// {
//   return (std::abs(CGAL::to_double(pt1.x()-pt2.x()))<compTol)&&(std::abs(CGAL::to_double(pt1.y()-pt2.y()))<compTol);
// }
//
// bool BoundaryPolygonGrid::Equal(double x1, double y1, double x2, double y2) const
// {
//   return (std::abs(x1-x2)<compTol)&&(std::abs(y1-y2)<compTol);
// }
//
// bool BoundaryPolygonGrid::InGridCell(Point_2 const& pt) const
// {
//   double xleft = indices.first * grid->dx + grid->xMin;
//   double xright = (indices.first+1) * grid->dx + grid->xMin;
//   double ybot = indices.second * grid->dy + grid->yMin;
//   double ytop = (indices.second + 1) * grid->dy + grid->yMin;
//
//   double x = CGAL::to_double( pt.x() );
//   double y = CGAL::to_double( pt.y() );
//
//   bool inx = (x > xleft-compTol) && (x < xleft+compTol);
//   bool iny = (y > ybot-compTol) && (y < ytop+compTol);
//
//   return inx && iny;
// }
//
// void BoundaryPolygonGrid::AddCellPts(std::vector<Point_2> &polyPts)
// {
//   double xIn = CGAL::to_double( polyPts.at(0).x() );
//   double yIn = CGAL::to_double( polyPts.at(0).y() );
//   double xOut = CGAL::to_double( polyPts.at(polyPts.size()-1).x() );
//   double yOut = CGAL::to_double( polyPts.at(polyPts.size()-1).y() );
//
//   double xleft = indices.first * grid->dx + grid->xMin;
//   double xright = (indices.first+1) * grid->dx + grid->xMin;
//   double ybot = indices.second * grid->dy + grid->yMin;
//   double ytop = (indices.second + 1) * grid->dy + grid->yMin;
//
//   // Check to see if the edge lies completely on the boundary
//   if( Equal(xIn,yIn,xleft,ytop) && Equal(xOut,yOut,xleft,ybot)){
//     polyPts.push_back( Point_2(xright,ybot));
//     polyPts.push_back( Point_2(xright,ytop));
//     return;
//   }else if( Equal(xIn,yIn,xleft,ybot) && Equal(xOut,yOut,xright,ybot)){
//     polyPts.push_back( Point_2(xright,ytop));
//     polyPts.push_back( Point_2(xleft,ytop));
//     return;
//   }else if( Equal(xIn,yIn,xright,ybot) && Equal(xOut,yOut,xright,ytop)){
//     polyPts.push_back( Point_2(xleft,ytop));
//     polyPts.push_back( Point_2(xleft,ybot));
//     return;
//   }else  if( Equal(xIn,yIn,xright,ytop) && Equal(xOut,yOut,xleft,ytop)){
//     polyPts.push_back( Point_2(xleft,ybot));
//     polyPts.push_back( Point_2(xright,ybot));
//     return;
//   }
//
//   // Is the exit point along the bottom edge?
//   if(std::abs(yOut-ybot)<compTol){
//
//     // is the entrace point along the bottom edge?
//     if(std::abs(yIn-ybot)<compTol){
//       return;
//     }
//
//     // If not, and the exit is not the bottom right, add the bottom right point
//     if(std::abs(xOut-xright)>compTol){
//       polyPts.push_back( Point_2(xright,ybot));
//     }
//
//     // is the entrance point along the right edge?
//     if(std::abs(xIn-xright)<compTol){
//       return;
//     }
//
//     // if not, add the top right point
//     polyPts.push_back( Point_2(xright, ytop));
//
//     // is the entrance point along the top edge?
//     if(std::abs(yIn-ytop)<compTol){
//       return;
//     }
//     // if not, add the top left point
//     polyPts.push_back( Point_2(xleft, ytop));
//
//     // is the entrance point along the left edge?
//     if(std::abs(xIn-xleft)<compTol){
//       return;
//     }
//
//     // if not, and the exit point is not the bottom left corner, add the bottom left corner
//     if(std::abs(xOut-xleft)>compTol){
//       polyPts.push_back( Point_2(xleft, ybot));
//     }
//
//   // Is the exit point along the right boundary of the grid cell
//   }else if(std::abs(xOut-xright)<compTol){
//
//     // is the entrance point along the right boundary?
//     if(std::abs(xIn-xright)<compTol){
//       return;
//     }
//
//     // If not, and the exit is not the top right, add the top right point
//     if(std::abs(yOut-ytop)>compTol){
//       polyPts.push_back( Point_2(xright,ytop));
//     }
//
//     // is the entrance point along the top edge?
//     if(std::abs(yIn-ytop)<compTol){
//       return;
//     }
//
//     // if not, add the top left point
//     polyPts.push_back( Point_2(xleft, ytop));
//
//     // is the entrance point along the left edge?
//     if(std::abs(xIn-xleft)<compTol){
//       return;
//     }
//     // if not, add the bottom left point
//     polyPts.push_back( Point_2(xleft, ybot));
//
//     // is the entrance point along the bottom edge?
//     if(std::abs(yIn-ybot)<compTol){
//       return;
//     }
//
//     // if not, and the exit point is not the bottom right corner, add the bottom right corner
//     if(std::abs(yOut-ybot)>compTol){
//       polyPts.push_back( Point_2(xright, ybot));
//     }
//
//   // Is the exit point along the top boundary of the grid cell?
//   }else if(std::abs(yOut-ytop)<compTol){
//
//     // is the entrance point along the top boundary?
//     if(std::abs(yIn-ytop)<compTol){
//       return;
//     }
//
//     // If not, and the exit is not the top left, add the top left point
//     if(std::abs(xOut-xleft)>compTol){
//       polyPts.push_back( Point_2(xleft,ytop));
//     }
//
//     // is the entrance point along the left edge?
//     if(std::abs(xIn-xleft)<compTol){
//       return;
//     }
//
//     // if not, add the bottom left point
//     polyPts.push_back( Point_2(xleft, ybot));
//
//     // is the entrance point along the bottom edge?
//     if(std::abs(yIn-ybot)<compTol){
//       return;
//     }
//     // if not, add the bottom right point
//     polyPts.push_back( Point_2(xright, ybot));
//
//     // is the entrance point along the right edge?
//     if(std::abs(xIn-xright)<compTol){
//       return;
//     }
//
//     // if not, and the exit point is not the top right corner, add the top right corner
//     if(std::abs(xOut-xright)>compTol){
//       polyPts.push_back( Point_2(xright, ytop));
//     }
//
//   // Is the exit  point along the left boundary?
//   }else if(std::abs(xOut-xleft)<compTol){
//
//     // is the entrance point along the left boundary?
//     if(std::abs(xIn-xleft)<compTol){
//       return;
//     }
//
//     // If not, and the exit is not the bottom left, add the bottom left point
//     if(std::abs(yOut-ybot)>compTol){
//       polyPts.push_back( Point_2(xleft,ybot));
//     }
//
//     // is the entrance point along the bottom edge?
//     if(std::abs(yIn-ybot)<compTol){
//       return;
//     }
//
//     // if not, add the bottom right point
//     polyPts.push_back( Point_2(xright, ybot));
//
//     // is the entrance point along the right edge?
//     if(std::abs(xIn-xright)<compTol){
//       return;
//     }
//
//     // if not, add the top right point
//     polyPts.push_back( Point_2(xright, ytop));
//
//     // is the entrance point along the top edge?
//     if(std::abs(yIn-ytop)<compTol){
//       return;
//     }
//
//     // if not, and the exit point is not the top left corner, add the top left corner
//     if(std::abs(yOut-ytop)>compTol){
//       polyPts.push_back( Point_2(xleft, ytop));
//     }
//   }
// }
//
// BoundaryPolygonGrid& BoundaryPolygonGrid::Increment()
// {
//
//   // Increment the indices
//   currEdge = nextEdge;
//   IncrementIndices();
//
//   // Compute the overlaping polygon AND update the edge iterator if needed
//
//   // Points in the intersection polynomial
//   std::vector<Point_2> polyPts;
//
//   // Figure out where the current Laguerre edge enters and leaves the current grid cell
//   std::shared_ptr<Point_2> enterPt, exitPt;
//   std::tie(enterPt, exitPt) = EdgeCrossings(nextEdge,true);
//
//   // If things are working properly, the entrance point will only not exist if
//   // the source of the edge coincides with a corner.  So just add the source.
//   if(!enterPt){
//     polyPts.push_back(nextEdge->source());
//   }else{
//     polyPts.push_back(*enterPt);
//   }
//
//   /* To figure out what the next cell is, iterate over the Laguerre edges until
//      we find one that leaves the current grid cell
//   */
//   while(!exitPt){
//     polyPts.push_back(nextEdge->target());
//
//     nextEdge++;
//
//     // Check to see if we've made it all the way around the Laguerre cell
//     if(nextEdge==startEdge){
//       isValid = false;
//       break;
//     }
//
//     std::tie(enterPt, exitPt) = EdgeCrossings(nextEdge,true);
//   }
//
//   if(!isValid)
//     return *this;
//
//   polyPts.push_back(*exitPt);
//
//   // Add corners from the grid cell to complete the intersection polygon
//   AddCellPts(polyPts);
//   overlapPoly = std::make_shared<Polygon_2>(polyPts.begin(), polyPts.end());
//
//   // Check to see if we've made it all the way around the Laguerre cell
//   if(indices==startIndices){
//     isValid = false;
//     return *this;
//   }
//
//   return *this;
// }

/** Works somewhat like an iterator over grid cells contained entirely inside a
    polygon.  Works from the bottom of the polygon to the top, with the fast index
    over x.
*/


/** Class for solving semi-discrete optimal transport problems. */
class SemidiscreteOT
{

public:

  /**
    @param[in] gridIn Regular 2d grid defining cells of "continuous" distribution.
    @param[in] gridProbs Matrix of probabilities for each cell.  Rows correspond to X locations, columns to Y locations.
    @param[in] discrPts Matrix of discrete points.
    @param[in] discrProbs Probabilities assigned to each discrete point.
  */
  SemidiscreteOT(std::shared_ptr<RegularGrid> const& gridIn,
                 Eigen::MatrixXd              const& gridProbsIn,
                 Eigen::Matrix2Xd             const& discrPtsIn,
                 Eigen::VectorXd              const& discrProbsIn);

  double Objective(Eigen::VectorXd const& prices);

  Eigen::VectorXd Gradient(Eigen::VectorXd const& prices);

  Eigen::MatrixXd Hessian(Eigen::VectorXd const& prices);

private:

  std::shared_ptr<RegularGrid> grid;
  Eigen::MatrixXd  gridProbs;

  Eigen::Matrix2Xd discrPts;
  Eigen::VectorXd  discrProbs;

  Eigen::Matrix2Xd domain; // <- Points defining a polygon surrounding the domain of interest
};

SemidiscreteOT::SemidiscreteOT(std::shared_ptr<RegularGrid> const& gridIn,
                               Eigen::MatrixXd              const& gridProbsIn,
                               Eigen::Matrix2Xd             const& discrPtsIn,
                               Eigen::VectorXd              const& discrProbsIn) : grid(gridIn),
                                                                                   gridProbs(gridProbsIn),
                                                                                   discrPts(discrPtsIn),
                                                                                   discrProbs(discrProbsIn)
{
  domain.resize(2,4);
  domain << grid->xMin, grid->xMax, grid->xMax, grid->xMin,
            grid->yMin, grid->yMin, grid->yMax, grid->yMax;
}



double SemidiscreteOT::Objective(Eigen::VectorXd const& prices)
{
    // Notes:
    //   - The cost c(x,y) is the squared distance between x and y
    //   - See (17) of https://arxiv.org/pdf/1710.02634.pdf

    const int numCells = discrPts.cols();
    assert(numCells==prices.size());

    // Construct the Laguerre diagram
    LaguerreDiagram lagDiag(discrPts, prices, domain);

    // Holds the part of the objective for each cell in the Laguerre diagram
    Eigen::VectorXd cellParts;
    for(int cellInd=0; cellInd<numCells; ++cellInd)
      cellParts(cellInd) = prices(cellInd)*discrProbs(cellInd);

    // Loop over  all of the grid cells and add contributions to the relevant Laguerre cells
    // for(int  xInd=0; xInd<grid->NumCells(0); ++xInd){
    //   for(int  yInd=0; yInd<grid->NumCells(1); ++yInd){
    //
    //     //lagDiag.GetCells(gridCell);
    //
    //   }
    // }
    return 0.0;
}


int main(int argc, char* argv[])
{

  int numPts = 3;
  Eigen::VectorXd costs(numPts);
  costs << 1.0, 1.0, 1.0;

  Eigen::Matrix2Xd pts(2,numPts);
  pts << 0.1, 0.2, 0.1,
         0.1, 0.1, 0.2;

  Eigen::Matrix2Xd domain(2,4);
  domain << 0.0, 2.0, 2.0, 0.0,
            0.0, 0.0, 1.0, 1.0;

  LaguerreDiagram diag(pts, costs, domain);

  auto grid = std::make_shared<RegularGrid>(domain(0,0),domain(1,0), domain(0,2), domain(1,2), 200, 2000);

  double area = 0.0;
  Eigen::VectorXd localAreas = Eigen::VectorXd::Zero(3);
  Eigen::MatrixXd cellAreas = Eigen::MatrixXd::Zero(grid->NumCells(0), grid->NumCells(1));

  for(int polyInd=0; polyInd<3; ++polyInd){
    std::cout << "Polygon " << polyInd << std::endl;
    std::shared_ptr<PolygonRasterizeIter::Polygon_2> poly = diag.GetCell(polyInd);

    PolygonRasterizeIter gridIter(grid,poly);
    //
    // unsigned int oldYInd = gridIter.Indices().second;
    // std::cout << "yind = " << oldYInd << std::endl;
    // std::cout << "yval = " << grid->yMin + oldYInd*grid->dy << std::endl;
    // std::cout << "    ";

    while(gridIter.IsValid()){


      // if(gridIter.Indices().second != oldYInd){
      //   oldYInd = gridIter.Indices().second;
      //   std::cout << "\nyind = " << oldYInd << std::endl;
      //   std::cout << "yval = " << grid->yMin + oldYInd*grid->dy << std::endl;
      //   std::cout << "    ";
      // }

      double cellArea;

      // std::cout << gridIter.Indices().first << ",  " << std::flush;

      if(gridIter.IsBoundary()){
        cellArea = CGAL::to_double( gridIter.OverlapPoly()->area() );

        // unsigned int indX = gridIter.Indices().first;
        // unsigned int indY = gridIter.Indices().second;
        //
        // PolygonRasterizeIter::Polygon_2 tempPoly;
        // tempPoly.push_back(PolygonRasterizeIter::Point_2(grid->xMin + indX*grid->dx, grid->yMin+indY*grid->dy));
        // tempPoly.push_back(PolygonRasterizeIter::Point_2(grid->xMin + (indX+1)*grid->dx, grid->yMin+indY*grid->dy));
        // tempPoly.push_back(PolygonRasterizeIter::Point_2(grid->xMin + (indX+1)*grid->dx, grid->yMin+(indY+1)*grid->dy));
        // tempPoly.push_back(PolygonRasterizeIter::Point_2(grid->xMin + indX*grid->dx, grid->yMin+(indY+1)*grid->dy));
        //
        // // Compute the intersection of P and Q.
        // std::list<PolygonRasterizeIter::Polygon_with_holes_2> interList;
        // CGAL::intersection(*poly, tempPoly, std::back_inserter(interList));
        //
        // double trueArea = CGAL::to_double( interList.begin()->outer_boundary().area() );
        //
        // //std::cout << *gridIter.OverlapPoly() << std::endl;
        // double error = cellArea-trueArea;
        // if(std::abs(error)>std::numeric_limits<double>::epsilon()){
        //   // std::cout << "  Cell area = " << cellArea << " with error " << error << std::endl;
        // }
      }else{
        cellArea = grid->dx*grid->dy;
      }

      cellAreas(gridIter.Indices().first, gridIter.Indices().second) += cellArea;
      localAreas(polyInd) += std::abs(cellArea);
      area += std::abs(cellArea);

      gridIter.Increment();
    }
    // std::cout << "\n\n";
  }

  Eigen::VectorXd trueAreas(3);
  trueAreas << 0.15*0.15, (domain(0,2)-0.15)*0.15 + 0.5*(domain(1,2)-0.15)*(domain(0,2)-0.15), 0.15*(domain(1,2)-0.15) + 0.5*(domain(1,2)-0.15)*(domain(0,2)-0.15);

  // List all cells with errors
  for(unsigned int yInd=0; yInd<grid->NumCells(1); ++yInd){
    for(unsigned int xInd=0; xInd<grid->NumCells(0); ++xInd){
      if(std::abs(cellAreas(xInd,yInd) - grid->dx*grid->dy)>1e-15){
        std::cout << "Cell " << xInd << "," << yInd << " error = " << cellAreas(xInd,yInd) - grid->dx*grid->dy << std::endl;
      }
    }
  }

  std::cout << "Local areas = " << localAreas.transpose() << std::endl;

  //std::cout << "dx*dy = " << grid->dx * grid->dy << std::endl;
  std::cout << "True total area = " << grid->dx*grid->dy*grid->NumCells() << std::endl;
  std::cout << "Total Area = " << area << std::endl;
  std::cout << "Error = " << area - grid->dx*grid->dy*grid->NumCells() << std::endl;
  return 0;
}
