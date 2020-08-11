#include "SDOT/BoundingBox.h"

#include <CGAL/Boolean_set_operations_2.h>


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

    if((sy>yMin-compTol)&&(sy<yMax+compTol)&&(std::max(sx,tx)>xMin-compTol)&&(std::min(sx,tx)<xMax+compTol)){
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

  Polygon_2 boxPoly;
  boxPoly.push_back(Point_2(xMin,yMin));
  boxPoly.push_back(Point_2(xMax,yMin));
  boxPoly.push_back(Point_2(xMax,yMax));
  boxPoly.push_back(Point_2(xMin,yMax));

  std::list<Polygon_with_holes_2> temp;
  CGAL::intersection(*poly, boxPoly, std::back_inserter(temp));
  if(temp.size()==0){
    return nullptr;
  }
  // // There should only be one polygon intersection because everything is simple and convex
  assert(std::distance(temp.begin(),temp.end())==1);

  auto outPoly = std::make_shared<Polygon_2>(temp.begin()->outer_boundary());

  return outPoly;

  ///////// TODO: FIX BUG BELOW HERE!

  std::vector<Point_2> newPts;

  // Get an edge circulator
  auto firstEdge = poly->edges_circulator();
  auto currEdge = firstEdge;

  // Loop over the edges, clipping them to the boundary and adding corners if necessary
  do {

    bool srcInside = IsInside(currEdge->source());
    bool tgtInside = IsInside(currEdge->target());

    if(srcInside){
      newPts.push_back(currEdge->source());

      // If the src is inside but the target is not, then we need to figure out where the edges cross
      if(!tgtInside){


        Point_2 src = currEdge->source();
        Point_2 tgt = currEdge->target();
        bool clipWorked = ClipSegment(src,tgt);

        if(clipWorked)
          newPts.push_back(tgt);

        // Figure out the next edge that intersects the bounding box
        auto nextEdge = currEdge;
        nextEdge++;
        src = nextEdge->source();
        tgt = nextEdge->target();
        while(!ClipSegment(src,tgt)){
          nextEdge++;
          src = nextEdge->source();
          tgt = nextEdge->target();
        }

        // Add any necessary corners and the point where the polygon reenters
        AddCorners(src, newPts);
        newPts.push_back(src);
        //currEdge = nextEdge;
      }

    // source and target are outside
    }else{

      Point_2 src = currEdge->source();
      Point_2 tgt = currEdge->target();

      if(ClipSegment(src,tgt)){

        newPts.push_back(src);

        // If the target is not inside, we might also have to add some corners
        if(!tgtInside){
          newPts.push_back(tgt);

          auto nextEdge = currEdge;
          nextEdge++;
          src = nextEdge->source();
          tgt = nextEdge->target();
          while(!ClipSegment(src,tgt)){
            nextEdge++;
            src = nextEdge->source();
            tgt = nextEdge->target();
          }

          AddCorners(src, newPts);
        }
      }
    }

  } while((++currEdge) != firstEdge);

  if(newPts.size()<3){
    return nullptr;
  }else{
    return std::make_shared<BoundingBox::Polygon_2>(newPts.begin(), newPts.end());
  }
}

bool BoundingBox::IsInside(Point_2 const& pt) const{
  double x = CGAL::to_double(pt.x());
  double y = CGAL::to_double(pt.y());

  return (x>xMin-compTol)&&(x<xMax+compTol)&&(y>yMin-compTol)&&(y<yMax+compTol);
}

bool BoundingBox::OnEdge(Point_2 const& pt) const{
  double x = CGAL::to_double(pt.x());
  double y = CGAL::to_double(pt.y());

  return (std::abs(x-xMin)<compTol)||(std::abs(x-xMax)<compTol)||(std::abs(y-yMin)<compTol)||(std::abs(y-yMax)<compTol);
}


void BoundingBox::AddCorners(Point_2 const& nextPt, std::vector<Point_2> &polyPts) const
{
  double backX = CGAL::to_double( polyPts.back().x() );
  double backY = CGAL::to_double( polyPts.back().y() );


  if((backX>xMax+compTol)||(backX<xMin-compTol)||(backY<yMin-compTol)||(backY>yMax+compTol)){
    std::stringstream error;
    error << "ERROR: In BoundingBox::AddCorners, the last polyPt is not inside the grid cell." << std::endl;
    error << "  last polyPt = " << backX << ", " << backY << std::endl;
    error << "  xmin, xmax = " << xMin << ", " << xMax << std::endl;
    error << "  ymin, ymax = " << yMin << ", " << yMax << std::endl;
    throw std::runtime_error(error.str());
  }

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
      std::cerr << "In BoundingBox::AddCorners, the point doesn't seem to be on an edge... " << std::endl;
      for(auto& pt : polyPts){
        std::cerr << "  " << pt << std::endl;
      }
      std::cerr << "Edges:" << std::endl;
      std::cerr << "  xmin, xmax = " << xMin << ", " << xMax << std::endl;
      std::cerr << "  ymin, ymax = " << yMin << ", " << yMax << std::endl;

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
