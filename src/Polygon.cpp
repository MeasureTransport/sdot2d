#include "SDOT/Polygon.h"

#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>

#include <CGAL/ch_graham_andrew.h>

using namespace sdot;

Polygon::Polygon(Eigen::Matrix2Xd const& pts)
{
  cgalPoly = std::make_shared<Polygon_2>();
  for(int colInd=0; colInd<pts.cols()-1; ++colInd){
    cgalPoly->push_back(Point_2(pts(0,colInd),pts(1,colInd)));
  }

  // Only add the last point if the first and last points are different
  if((pts.col(0)-pts.col(pts.cols()-1)).cwiseAbs().maxCoeff()>1e-10){
    cgalPoly->push_back(Point_2(pts(0,pts.cols()-1),pts(1,pts.cols()-1)));
  }

  if(cgalPoly->orientation()==CGAL::CLOCKWISE){
    cgalPoly->reverse_orientation();
  }
  assert(cgalPoly->orientation()==CGAL::COUNTERCLOCKWISE);
  assert(cgalPoly->is_simple());
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
