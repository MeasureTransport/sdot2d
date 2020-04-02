#include "SDOT/BoundingBox.h"

using namespace sdot;

BoundingBox::BoundingBox(double xmin_, double xmax_,
                         double ymin_, double ymax_) : xMin(xmin_), xMax(xmax_),
                                                       yMin(ymin_), yMax(ymax_)
{
  assert(xMin<xMax);
  assert(yMin<yMax);
}


bool BoundingBox::ClipSegment(Point_2 &srcPt, Point_2 &tgtPt) const
{
  double sx = CGAL::to_double(srcPt.x());
  double sy = CGAL::to_double(srcPt.y());
  double tx = CGAL::to_double(tgtPt.x());
  double ty = CGAL::to_double(tgtPt.y());

  // Check if it's vertical
  if(std::abs(sx-tx)<compTol){

    // Check to see if it crosses the domain at all
    if((sx>xMin-compTol)&&(sx<xMax+compTol)&&(std::max(sy,ty)>yMin-compTol)&&(std::min(sy,ty)<yMax+compTol)){
      srcPt = Point_2(sx, std::min(yMax, std::max(yMin, sy)));
      tgtPt = Point_2(tx, std::min(yMax, std::max(yMin, ty)));
      return true;
    }else{
      return false;
    }

  // Check if it's horizontal
  }else if(std::abs(sy-ty)<compTol){

    if((sy>yMin-compTol)&&(sy<yMax-compTol)&&(std::max(sx,tx)>xMin-compTol)&&(std::min(sx,tx)<xMax+compTol)){
      srcPt = Point_2(std::min(xMax, std::max(xMin, sx)), sy);
      tgtPt = Point_2(std::min(xMax, std::max(xMin, tx)), ty);
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

    qs(0) = (sx-xMin);///(tx-sx);
    qs(1) = (xMax-sx);//(tx-sx);
    qs(2) = (sy-yMin);///(ty-sy);
    qs(3) = (yMax-sy);///(ty-sy);


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

std::shared_ptr<BoundingBox::Polygon_2> BoundingBox::ClipPolygon(std::shared_ptr<BoundingBox::Polygon_2> const& poly) const
{

  std::vector<Point_2> newPts;

  // Get an edge circulator
  auto firstEdge = poly->edges_circulator();
  auto currEdge = firstEdge;

  // Loop over the edges
  do {

    bool srcInside = IsInside(currEdge->source());

    if(srcInside){
      newPts.push_back(currEdge->source());
    }

    // If the src is inside but the target is not, then we need to figure out where the edges cross
    if((!IsInside(currEdge->target())) && srcInside){

      Point_2 src = currEdge->source();
      Point_2 tgt = currEdge->target();

      // Clip the target to the bounding box
      ClipSegment(src,tgt);
      newPts.push_back(tgt);

      // Figure out the next edge that has a target inside
      auto nextEdge = currEdge;
      nextEdge++;
      while(!IsInside(nextEdge->target()))
        nextEdge++;

      // Clip this edge to the bounding box
      src = nextEdge->source();
      tgt = nextEdge->target();
      ClipSegment(src,tgt);

      // Add any necessary corners and the point where the polygon reenters
      AddCorners(src, newPts);
      newPts.push_back(src);
    }

  } while((++currEdge) != firstEdge);

  return std::make_shared<BoundingBox::Polygon_2>(newPts.begin(), newPts.end());
}

bool BoundingBox::IsInside(Point_2 const& pt) const{
  double x = CGAL::to_double(pt.x());
  double y = CGAL::to_double(pt.y());

  return (x>xMin-compTol)&&(x<xMax+compTol)&&(y>yMin-compTol)&&(y<yMax+compTol);
}

void BoundingBox::AddCorners(Point_2 const& nextPt, std::vector<Point_2> &polyPts) const
{
  double backX = CGAL::to_double( polyPts.back().x() );
  double backY = CGAL::to_double( polyPts.back().y() );

  double enterX = CGAL::to_double( nextPt.x() );
  double enterY = CGAL::to_double( nextPt.y() );

  // Check to see if we need to add add any corners.
  // If the xs and the ys are different OR if the points are far enough apart, then we need to add corner
  bool sameXdiffY = (std::abs(backX-enterX)<compTol)&&(std::abs(backY-enterY)>compTol);
  bool sameYdiffX = (std::abs(backY-enterY)<compTol)&&(std::abs(backX-enterX)>compTol);
  bool sameEdge = sameXdiffY || sameYdiffX;

  if(sameEdge)
    return;

  for(unsigned int i=0; i<4; ++i){
  //while(!sameEdge){
    // Add corner point

    // Currently on the right
    if((std::abs(backX-xMax)<compTol)&&(backY<yMax-compTol)){
      polyPts.push_back( Point_2(xMax,yMax) );

    // Currently on the top
    }else if((std::abs(backY-yMax)<compTol)&&(backX>xMin+compTol)){
      polyPts.push_back( Point_2(xMin, yMax) );

    // Currently on the left
    }else if((std::abs(backX-xMin)<compTol)&&(backY>yMin+compTol)){
      polyPts.push_back( Point_2(xMin, yMin) );

    // On the bottom
    }else if((std::abs(backY-yMin)<compTol)&&(backX<xMax-compTol)){
      polyPts.push_back( Point_2(xMax, yMin) );

    }else{
      std::cerr << "I shouldn't be here... " << std::endl;
      for(auto& pt : polyPts){
        std::cout << "  " << pt << std::endl;
      }
      assert(false);
    }

    backX = CGAL::to_double( polyPts.back().x() );
    backY = CGAL::to_double( polyPts.back().y() );

    sameXdiffY = (std::abs(backX-enterX)<compTol)&&(std::abs(backY-enterY)>compTol);
    sameYdiffX = (std::abs(backY-enterY)<compTol)&&(std::abs(backX-enterX)>compTol);
    sameEdge = sameXdiffY || sameYdiffX;

    if(sameEdge){
      break;
    }
  }
}
