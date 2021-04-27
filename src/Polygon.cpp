#include "SDOT/Polygon.h"
#include "SDOT/Assert.h"

#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>

#include <CGAL/ch_graham_andrew.h>

using namespace sdot;

Polygon::Polygon(Eigen::Matrix2Xd const& pts)
{
  xmin = std::numeric_limits<double>::infinity();
  ymin = std::numeric_limits<double>::infinity();
  xmax = -std::numeric_limits<double>::infinity();
  ymax = -std::numeric_limits<double>::infinity();

  cgalPoly = std::make_shared<Polygon_2>();
  for(int colInd=0; colInd<pts.cols()-1; ++colInd){
    cgalPoly->push_back(Point_2(pts(0,colInd),pts(1,colInd)));

    xmin = std::min(xmin, pts(0,colInd));
    ymin = std::min(ymin, pts(1,colInd));

    xmax = std::max(xmax, pts(0,colInd));
    ymax = std::max(ymax, pts(1,colInd));
  }

  // Only add the last point if the first and last points are different
  if((pts.col(0)-pts.col(pts.cols()-1)).cwiseAbs().maxCoeff()>1e-10){
    cgalPoly->push_back(Point_2(pts(0,pts.cols()-1),pts(1,pts.cols()-1)));
  }

  if(cgalPoly->orientation()==CGAL::CLOCKWISE){
    cgalPoly->reverse_orientation();
  }
  SDOT_ASSERT(cgalPoly->orientation()==CGAL::COUNTERCLOCKWISE);
  SDOT_ASSERT(cgalPoly->is_simple());

  // Compute the centroid and area of the polygon if it's convex
  if(cgalPoly->is_convex()){
    _ComputeCentroidArea();
    _ComputeRadius();
  }

}

Polygon::Polygon(std::shared_ptr<Polygon_2> const& cgalPolyIn) : cgalPoly(cgalPolyIn)
{
  if(cgalPoly->size()>0){
    if(cgalPoly->orientation()==CGAL::CLOCKWISE){
      cgalPoly->reverse_orientation();
    }
    assert(cgalPoly->orientation()==CGAL::COUNTERCLOCKWISE);
    assert(cgalPoly->is_simple());
  }

  xmin = CGAL::to_double(  cgalPoly->left_vertex()->x() );
  xmax = CGAL::to_double( cgalPoly->right_vertex()->x()  );
  ymin = CGAL::to_double( cgalPoly->bottom_vertex()->y() );
  ymax = CGAL::to_double( cgalPoly->top_vertex()->y() );

  if(cgalPoly->is_convex()){
    _ComputeCentroidArea();
    _ComputeRadius();
  }
}


Eigen::Matrix2Xd Polygon::GetVertices() const
{
  Eigen::Matrix2Xd points(2,cgalPoly->size());

  // Copy the vertices from the polygon to the Eigen::Matrix
  unsigned int col = 0;
  for(auto vertIter = cgalPoly->vertices_begin(); vertIter!=cgalPoly->vertices_end(); ++vertIter, ++col){
    points(0,col) = CGAL::to_double( vertIter->x() );
    points(1,col) = CGAL::to_double( vertIter->y() );
  }

  return points;
}

std::shared_ptr<Polygon> Polygon::ConvexHull() const
{
  if(IsConvex()){
    return std::make_shared<Polygon>(cgalPoly);
  }else{

    std::shared_ptr<Polygon_2> outPoly = std::make_shared<Polygon_2>();
    CGAL::ch_graham_andrew( cgalPoly->vertices_begin(),
                            cgalPoly->vertices_end(),
                            std::back_inserter(*outPoly));

    return std::make_shared<Polygon>(outPoly);
  }
}


std::vector<std::shared_ptr<Polygon>> Polygon::ConvexPartition() const
{
  typedef CGAL::Partition_traits_2<Polygon::K>                Traits;
  typedef Traits::Point_2                                     Part_Point_2;
  typedef Traits::Polygon_2                                   Part_Polygon_2;
  typedef std::list<Part_Polygon_2>                           Part_Polygon_list;

  if(IsConvex()){
    return std::vector<std::shared_ptr<Polygon>>(1,std::make_shared<Polygon>(cgalPoly));

  }else{

    Part_Polygon_list partition_cgal;
    CGAL::optimal_convex_partition_2(cgalPoly->vertices_begin(),
                                    cgalPoly->vertices_end(),
                                    std::back_inserter(partition_cgal));

   assert(CGAL::convex_partition_is_valid_2(cgalPoly->vertices_begin(),
                                            cgalPoly->vertices_end(),
                                            partition_cgal.begin(),
                                            partition_cgal.end()));

    std::vector<std::shared_ptr<Polygon>> output;

    for(const Part_Polygon_2& poly : partition_cgal){
      auto polyPtr = std::make_shared<Polygon_2>();
      for(auto& p : poly.container()){
         polyPtr->push_back( Point_2(p.x(), p.y()) );
       }
      assert(polyPtr->is_convex());
      output.push_back(std::make_shared<Polygon>(polyPtr));
    }
    return output;
  }
}


void Polygon::_ComputeRadius(){

  maxRadius = 0.0;

  for(auto vertIter = cgalPoly->vertices_begin(); vertIter!=cgalPoly->vertices_end(); ++vertIter){
    maxRadius = std::max(maxRadius, std::pow(CGAL::to_double(vertIter->x())-centroid(0),2.0)+std::pow(CGAL::to_double(vertIter->y())-centroid(1),2.0));
  }

  maxRadius = std::sqrt(maxRadius);
}

void Polygon::_ComputeCentroidArea(){

  // Loop over triangles making up the convex polygon
  auto iter1 = cgalPoly->vertices_circulator();
  auto iter2 = iter1;
  iter2++;
  auto iter3 = iter2;
  iter3++;

  Eigen::Vector2d triCenter(2);
  triCenter << (1.0/3.0)*(CGAL::to_double(iter1->x())+CGAL::to_double(iter2->x())+CGAL::to_double(iter3->x())),
               (1.0/3.0)*(CGAL::to_double(iter1->y())+CGAL::to_double(iter2->y())+CGAL::to_double(iter3->y()));

  double triArea = std::abs(0.5*CGAL::to_double(iter1->x()*(iter2->y()-iter3->y()) + iter2->x()*(iter3->y()-iter1->y()) + iter3->x()*(iter1->y()-iter2->y())));
  area = triArea;
  centroid = triArea*triCenter;

  // Loop over the remaining triangles
  iter2++;
  iter3++;
  for(; iter3!=iter1; iter2++, iter3++){
    triCenter << (1.0/3.0)*CGAL::to_double(iter1->x()+iter2->x()+iter3->x()),
                 (1.0/3.0)*CGAL::to_double(iter1->y()+iter2->y()+iter3->y());

    triArea = std::abs(0.5*CGAL::to_double(iter1->x()*(iter2->y()-iter3->y()) + iter2->x()*(iter3->y()-iter1->y()) + iter3->x()*(iter1->y()-iter2->y())));
    centroid = triArea*triCenter;
  }

}

/** Returns the intersection polygon of this polygon with another. */
Polygon Polygon::Intersection(Polygon const& otherPoly) const
{
  std::list<Polygon_with_holes_2> intR;
  CGAL::intersection(*cgalPoly, *otherPoly.cgalPoly, std::back_inserter(intR));

  return Polygon(std::make_shared<Polygon_2>(intR.begin()->outer_boundary()));
}

/** Returns true if this polygon intersects the other polygon and false otherwise. */
bool Polygon::Intersects(Polygon const& otherPoly) const
{
  // Check the bounding boxes
  if((xmax>=otherPoly.xmin)&&(otherPoly.xmax>=xmin)&&(ymax>=otherPoly.ymin)&&(otherPoly.ymax>=ymin)){
    return CGAL::do_intersect(*cgalPoly, *otherPoly.cgalPoly);
  }else{
    return false;
  }
}



/** Translate all vertices in the polygon. */
Polygon& Polygon::Translate(Eigen::Vector2d const& step)
{
  xmin += step(0);
  xmax += step(0);
  ymin += step(1);
  ymax += step(1);

  centroid += step;

  for(auto vertIter = cgalPoly->vertices_begin(); vertIter!=cgalPoly->vertices_end(); ++vertIter){
    cgalPoly->set(vertIter, Point_2(vertIter->x()+step(0), vertIter->y()+step(1)));
  }

  return *this;
}

/** Rotate the polygon around its centroid. */
Polygon& Polygon::Rotate(double angle)
{
  Eigen::Vector2d delta(2), newPt(2);
  Eigen::Matrix2d rotMat(2,2);
  rotMat << cos(angle), -sin(angle),
            sin(angle), cos(angle);

  xmin = std::numeric_limits<double>::infinity();
  xmax = -std::numeric_limits<double>::infinity();
  ymin = std::numeric_limits<double>::infinity();
  ymax = -std::numeric_limits<double>::infinity();

  for(auto vertIter = cgalPoly->vertices_begin(); vertIter!=cgalPoly->vertices_end(); ++vertIter){
    delta(0) = CGAL::to_double(vertIter->x()) - centroid(0);
    delta(1) = CGAL::to_double(vertIter->y()) - centroid(1);
    newPt = centroid + rotMat*delta;

    cgalPoly->set(vertIter, Point_2(newPt(0), newPt(1)));

    xmin = std::min(xmin, newPt(0));
    xmax = std::max(xmax, newPt(0));
    ymin = std::min(ymin, newPt(1));
    ymax = std::min(ymax, newPt(1));
  }

  return *this;
}


double Polygon::SecondMoment() const
{
  auto vert1 = cgalPoly->vertices_begin();
  auto vert2 = vert1;
  vert2++;
  double Jx = 0;
  double Jy = 0;

  for(; vert2!=cgalPoly->vertices_end(); ++vert1, ++vert2){
    Jx += CGAL::to_double((vert1->x()*vert2->y()-vert2->x()*vert1->y())*(vert1->y()*vert1->y() + vert1->y()*vert2->y()+vert2->y()*vert2->y()));
    Jy += CGAL::to_double((vert1->x()*vert2->y()-vert2->x()*vert1->y())*(vert1->x()*vert1->x() + vert1->x()*vert2->x()+vert2->x()*vert2->x()));
  }

  double Jz = (Jx+Jy)/12.0;

  // apply parallel axis theorem to get MOI around centroid
  return Jz - Area()*centroid.squaredNorm();
}
